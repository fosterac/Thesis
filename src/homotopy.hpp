#include <algorithm>

#include "HomotopyTypes.h"

//Problem API
#include "Problems.h"

//Include source files
#include "evaluator.hpp"
#include "Scalarization.hpp"
#include "Constraint.hpp"
#include "mesh.hpp"
#include "MeshSet.hpp"
#include "CommInterface.hpp"
#include "JobQueue.hpp"

//Optimizer interface
#include "nlopt.hpp"
#include "NloptAdapt.hpp"
#include "optimizer.hpp"

using namespace Homotopy;

namespace Pareto {

	class homotopy {
	private:
		Problem::Interface *Prob;

        typedef Communication::Interface comm_t;
        comm_t& Comm;

        //The evaluation chain is agnostic to the
        //local vs. remote evaluation strategy
        typedef JobQueue< comm_t > queue_t;
        queue_t Queue;

        //typedef Evaluator< EvaluationStrategy::Local< functionSet_t > > eval_t;
        typedef Evaluator< EvaluationStrategy::Cached< EvaluationStrategy::ReferenceTo < queue_t > > > eval_t;

        DynamicScalarization< eval_t > Scal;
		optimizer * Opt;

		double tolerance;
        double fd_step;
        FiniteDifferences::FD_TYPE fd_type;

        //Default values
        static const double FDstep = 1e-6;
        static const FiniteDifferences::FD_TYPE FDtype = FiniteDifferences::CENTRAL ;

		//For holding the mesh corners
		std::vector< std::vector< double > > Design;
		std::vector< std::vector< double > > Objective;
		std::vector< std::vector< double > > Lambda;

		static std::vector< double > GetF(const std::vector< Problem::FUNCTION > &f, const std::vector< double > &x){
			std::vector< double > result;
			int i;
			for(i=0; i<f.size(); i++) { result.push_back((f[i])(x)); }
			return result;
		}

        void GetCorners(){
            //Find the individual optimae using a fixed scalarization
            FixedScalarization< eval_t > S(Prob, Queue);
            //FixedScalarization< Evaluator< EvaluationStrategy::Local< functionSet_t > > > S(Prob, Prob->Objectives);

            //Establish finite difference parameters
            FiniteDifferences::Params_t FDpar = { this->fd_step, this->fd_type };

			int i;
			for( i=0; i<Prob->Objectives.size(); i++){

                optimizer * op = new OptNlopt( &S, tolerance, FDpar);

				//Optimize the objective
				std::vector<double> lam(Prob->Objectives.size(), 0.0);
				lam[i] = 1.0;
				S.SetWeights( &lam );
				std::vector<double> start(Prob->dimDesign, 0.0);
                int j;
                for(j=0; j<Prob->dimDesign; j++){
                    start[j] = (Prob->upperBounds[j] - Prob->lowerBounds[j])/2.0 + Prob->lowerBounds[j];
                }

                std::vector<double> x(start);
                this->Queue.NewGroup( i );
                optimizer::EXIT_COND flag = (optimizer::EXIT_COND) op->RunFrom( x );
                while( flag == optimizer::RERUN ) {
                    this->Queue.Poll();
                    this->Queue.NewGroup( i );

                    x = start;
                    flag = (optimizer::EXIT_COND) op->RunFrom( x );
                }

                //Make sure we have a successful optimization
                //(Should probably print the results after they're available)
                if( flag != OptNlopt::SUCCESS ) throw std::runtime_error( "Unable to optimize objective." );

				std::vector<double> f( S.e.eval( x ) );

				//Cache the results
				this->Design.push_back( x );
				this->Objective.push_back( f );
				this->Lambda.push_back( lam );
                    
                //TODO: Extra jobs are somehow lingering in the list...
                this->Queue.Clear();
                delete op;
			}
        }

	public:
        homotopy( Problem::Interface *P, double tolerance, double fd_step, Communication::Interface & c) : Prob(P), Comm( c ), Queue( Comm ), Scal( Prob, Queue ), tolerance(tolerance), fd_step(FDstep), fd_type(FDtype), Opt( NULL ){
			this->GetCorners();
		}

        typedef FEqDistanceConstraint< boost::function<objVars_t (const designVars_t&)> > FunctionSpaceEqDistConstr;
		void GetFront(int NumPoints, int Iterations, int id, int worldSize){
			
            //Get the list of ID's this instance is responsible for
            int Subsets = worldSize;
            std::vector< ind_t > IDs;
            IDs.push_back( id );
            //int id;
            //for(id=0; id<Mesh::Simplex::eta( this->Prob->dimObj - 1, Subsets ); id++){ IDs.push_back( id ); }

            //Instantiate mesh
			//Mesh::Simplex mesh( this->Design, this->Objective, this->Lambda, NumPoints);
            MeshSet< Mesh::SimplexSubset, comm_t > mesh( this->Design, this->Objective, this->Lambda, NumPoints, IDs , Subsets, Comm);
            mesh.Generate();

			/*
			TODO: This point needs further investigation as starting from
			the Obj Space guess location seems to be a winning prospect for
			quickly converging on convex fronts.  However, in the real FON
			example, leads to roundoff error for large #'s of points.
			Starting, however, from the real locations of the interpolated
			Pareto Set (not front) allows much larger mesh sizes without
			roundoff errors and contains points rooted in some kind of truth.
			
			//Project the design space mesh onto the objective space
			int m;
			for(m=0;m<mesh.Points.size();m++){
				//mesh.Points[ m ].ObjectiveCoords = GetF( Prob->Objectives, mesh.Points[ m ].DesignCoords ) ;
			}
            */

			//Constraints

			//Sum constraints on the Lambda parameters
			BoundSumConstraint LT( BoundSumConstraint::LESS_THAN, 1.0, Prob->dimDesign, Scal.dimDesign);
			BoundSumConstraint GT( BoundSumConstraint::GREATER_THAN, 0.0, Prob->dimDesign, Scal.dimDesign);
			Scal.InequalityConstraints.push_back( LT.function );
			Scal.InequalityConstraints.push_back( GT.function );

			//Incorporate constraint
            boost::function<objVars_t (const designVars_t&)> f = boost::bind( &eval_t::eval, &(Scal.e), _1);
			std::vector< FunctionSpaceEqDistConstr* > NeighborConstraints( mesh.MeshDim );
			int n;
			for(n=0;n<mesh.MeshDim;n++) {
				NeighborConstraints[n] = new FunctionSpaceEqDistConstr (f, this->Prob->dimDesign, NULL, NULL) ;
            }

            //Establish finite difference parameters
            FiniteDifferences::Params_t FDpar = { this->fd_step, this->fd_type };

			//Construct optimizer
			this->Opt = new OptNlopt(&Scal, tolerance, FDpar);

			//Run a set of updates
			int j;
			for(j=0; j<Iterations; j++){
                //Send ghost nodes
                //Receive ghost nodes
                mesh.Refresh();

                for(id=0;id<mesh.GetSize();id++){
                    //Clear out the backlog of evaluations
                    Queue.Clear();

                    Mesh::MeshBase *m = mesh.Get(id);

                    //Status variables
                    std::vector< bool > flags(m->Points.size(), false);
                    optimizer::EXIT_COND ec;

                    //Initial run of the optimizer for all the points
                    int i;
				    for(i=0; i<m->Points.size(); i++){

                        //Specify a group of evaluations
                        Queue.NewGroup( i );

                        //Refine the point
                        ec = RefinePoint(i, this->Opt, *m, NeighborConstraints);
                        if( ec != optimizer::RERUN ) flags[i] = true;
                    }

                    //Do we have everything?
                    bool alldone = ( std::find( flags.begin(), flags.end(), false ) == flags.end() );

                    //If not, run the poll loop
                    while( !alldone ){

                        //Run poll_loop() to get i
                        i = Queue.Poll();
                        Queue.NewGroup( i );

                        //Run the update
                        ec = RefinePoint(i, this->Opt, *m, NeighborConstraints);

                        if( ec != optimizer::RERUN ) flags[i] = true;

                        //Again, do we have everything?
				        alldone = ( std::find( flags.begin(), flags.end(), false ) == flags.end() );
				    }
                }
			}
			delete this->Opt;

			mesh.Print();
			mesh.WriteOut( "front.txt" );
		}

        optimizer::EXIT_COND RefinePoint( int i, optimizer * opt, Mesh::MeshBase &mesh, std::vector< FunctionSpaceEqDistConstr * > &NeighborConstraints ) {
            if( !mesh.Points[i].Neighbors.empty() ){

				//NOTE: Pass only the required number of constraints to the
				//optimizer and then establish their neighbors.
                std::vector< std::vector< double >* > n (mesh.GetNeighborLocOf( i ));
				int iter;
				for(iter=0; iter<n.size()/2; iter++) {
                    std::vector< double > *left = n[2*iter];
                    std::vector< double > *right = n[2*iter+1];

					NeighborConstraints[iter]->UpdateFrom( left, right );
					Scal.EqualityConstraints.push_back( NeighborConstraints[iter]->function );
				}

                opt->RefreshConstraints();

				//Get Design points
				std::vector< double > x( mesh.Points[i].DesignCoords );

				//Add the lambda values
				int k;
				for(k=0;k<mesh.Points[i].LambdaCoords.size()-1;k++) {
					x.push_back( mesh.Points[i].LambdaCoords[k] );
				}

                optimizer::EXIT_COND flag = (optimizer::EXIT_COND) opt->RunFrom( x );

				if( flag == OptNlopt::SUCCESS ) {
					//Update design points
                    std::vector< double > d( x.begin(), x.begin() + Prob->dimDesign );

					//Update obj points
					std::vector< double > f( this->Scal.e.eval( mesh.Points[i].DesignCoords ) );

					//Update lam points
					std::vector< double > l ( x.begin() + Prob->dimDesign, x.end() );
					l.push_back( 1.0 - std::accumulate( l.begin(), l.end(), 0.0 ) );

                    //Run the update
                    mesh.UpdatePoint( i, d, f, l );
				}

				//Pop the EQDist Constraints
				Scal.EqualityConstraints.erase( Scal.EqualityConstraints.end() - iter, Scal.EqualityConstraints.end());

                return flag;
			}

            //Neighborless points always succeed
            return optimizer::SUCCESS;
        }
	};
}
