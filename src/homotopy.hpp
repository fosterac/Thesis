#include <algorithm>

#include "HomotopyTypes.h"

//Problem API
#include "Problems.h"

//Include source files
#include "evaluator.hpp"
#include "Scalarization.hpp"
#include "Constraint.hpp"
#include "mesh.hpp"
#include "CommInterface.hpp"
#include "JobQueue.hpp"

//Optimizer interface
#include "nlopt.hpp"
#include "NloptAdapt.hpp"
#include "optimizer.hpp"

namespace Pareto {

	class Homotopy {
	private:
		Problem::Interface *Prob;
        
        typedef JobQueue< Communication::SimulatedRemote< functionSet_t > > queue_t;

        //typedef Communication::AdHoc< typename Communication::CommImpl::Iface > comm_t;
        //typedef JobQueue< comm_t > queue_t;

        queue_t Queue;
        
        //typedef Evaluator< EvaluationStrategy::Local< functionSet_t > > eval_t;
        typedef Evaluator< EvaluationStrategy::Cached< EvaluationStrategy::ReferenceTo < queue_t > > > eval_t;
        
        DynamicScalarization< eval_t > Scal;
		Optimizer * Opt;

		double tolerance;
        static const double FDstep = 1e-7;
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

	public:
        Homotopy( Problem::Interface *P, double tolerance) : Prob(P), Queue( Prob->Objectives ), Scal( Prob, Queue ), tolerance(tolerance), 
        //Homotopy( Problem::Interface *P, double tolerance) : Prob(P), Scal( Prob, Queue ), tolerance(tolerance), 
			Opt( NULL ) {
				
				//Find the individual optimae using a fixed scalarization
                FixedScalarization< eval_t > S(Prob, Queue);

                //Establish finite difference parameters
                FiniteDifferences::Params_t FDpar = { FDstep, FDtype };
				Optimizer * op = new OptNlopt( &S, tolerance, FDpar);

				int i;
				for( i=0; i<Prob->Objectives.size(); i++){
					//Optimize the objective
					std::vector<double> lam(Prob->Objectives.size(), 0.0);
					lam[i] = 1.0;
					S.SetWeights( &lam );
					std::vector<double> x(Prob->dimDesign, 0.0);
                    int j;
                    for(j=0; j<Prob->dimDesign; j++){
                        x[j] = (Prob->upperBounds[j] - Prob->lowerBounds[j])/2.0 + Prob->lowerBounds[j];
                    }

                    Queue.NewGroup( 0 );
                    OptNlopt::EXIT_COND flag = (OptNlopt::EXIT_COND) op->RunFrom( x );
                    while( flag == OptNlopt::RERUN ) {
                        Queue.Poll();
                        Queue.NewGroup( 0 );
                        flag = (OptNlopt::EXIT_COND) op->RunFrom( x );
                    }

					std::vector<double> f( GetF( Prob->Objectives, x) );

					//Cache the results
					Design.push_back( x );
					Objective.push_back( f );
					Lambda.push_back( lam );
				}
				delete op;
		}

        //TODO: Hack-ish, but only for now
        //void SetDispatcherAndHandler( typename comm_t::dispatcher_t d, typename comm_t::handler_t h ){
            //this->Queue.RemoteEvaluator.Dispatcher = d;
            //this->Queue.RemoteEvaluator.Handler = h;
        //}

        typedef FEqDistanceConstraint< boost::function<objVars_t (const designVars_t&)> > FunctionSpaceEqDistConstr;

		void GetFront(int NumPoints, int Iterations){
			//Instantiate mesh
			Mesh::Simplex mesh( this->Design, this->Objective, this->Lambda, NumPoints);

			/*
			TODO: This point needs further investigation as starting from
			the Obj Space guess location seems to be a winning prospect for 
			quickly converging on convex fronts.  However, in the real FON 
			example, leads to roundoff error for large #'s of points.
			Starting, however, from the real locations of the interpolated 
			Pareto Set (not front) allows much larger mesh sizes without 
			roundoff errors and contains points rooted in some kind of truth.
			*/
			//Project the design space mesh onto the objective space
			int m;
			for(m=0;m<mesh.Points.size();m++){
				//mesh.Points[ m ].ObjectiveCoords = GetF( Prob->Objectives, mesh.Points[ m ].DesignCoords ) ;
			}

			//Constraints

			//Sum constraints on the Lambda parameters
			BoundSumConstraint LT( BoundSumConstraint::LESS_THAN, 1.0, Prob->dimDesign, Scal.dimDesign);
			BoundSumConstraint GT( BoundSumConstraint::GREATER_THAN, 0.0, Prob->dimDesign, Scal.dimDesign);
			Scal.InequalityConstraints.push_back( LT.function );
			Scal.InequalityConstraints.push_back( GT.function );

			//Incorporate constraint
            
            bool placeholder = false;
            boost::function<objVars_t (const designVars_t&)> f = boost::bind( &eval_t::eval, &(Scal.e), _1, placeholder);
			std::vector< FunctionSpaceEqDistConstr* > NeighborConstraints( mesh.MeshDim );
			int n;
			for(n=0;n<mesh.MeshDim;n++) {
				NeighborConstraints[n] = new FunctionSpaceEqDistConstr (f, NULL, NULL) ;
            }

            //Establish finite difference parameters
            FiniteDifferences::Params_t FDpar = { FDstep, FDtype };

			//Construct optimizer
			//this->Opt = new OptNlopt(&Scal, tolerance, FDpar);

			//Run a set of updates
			int j;
			for(j=0; j<Iterations; j++){
                //Set of optimizers
                std::vector< Optimizer* > opts( mesh.Points.size(), NULL ) ;

                //Status variables
                std::vector< bool > flags(mesh.Points.size(), false);
                Optimizer::EXIT_COND ec;

                //Initial run of the optimizer for all the points
                int i;
                int k;
                for(k=0;k<2;k++){
				for(i=k%2; i<NumPoints; i+=2){
				//for(i=0; i<mesh.Points.size(); i++){
                    //Define a local optimizer
                    opts[i] = new OptNlopt(&Scal, tolerance, FDpar);
                    
                    //Specify a group of evaluations
                    Queue.NewGroup( i );
                    
                    //Refine the point
                    //ec = RefinePoint(i, this->Opt, mesh, NeighborConstraints);
					ec = RefinePoint(i, opts[i], mesh, NeighborConstraints);
                    if( ec != Optimizer::RERUN ) flags[i] = true;
                }
                

                //Do we have everything?
                bool alldone = ( std::find( flags.begin(), flags.end(), false ) == flags.end() );

                //If not, run the poll loop
                while( !alldone ){
                    
                    //Run poll_loop() to get i
                    //i = Comm.PollLoop();
                    i = Queue.Poll();
                    Queue.NewGroup( i );

                    //Run the update
                    //ec = RefinePoint(i, this->Opt, mesh, NeighborConstraints);
                    ec = RefinePoint(i, opts[i], mesh, NeighborConstraints);
                    
                    if( ec != Optimizer::RERUN ) flags[i] = true;


                    //Again, do we have everything?
				    //alldone = ( std::find( flags.begin(), flags.end(), false ) == flags.end() );
                    alldone = ( std::count( flags.begin(), flags.end(), false ) <= NumPoints/2 );
				}
                }

                std::vector< Optimizer* >::iterator iter;
                for(iter=opts.begin(); iter!=opts.end(); iter++) delete *iter;
			}
			delete this->Opt; 
			
			mesh.Print();			
			mesh.WriteOut( "front.txt" );	
		}

        Optimizer::EXIT_COND RefinePoint( int i, Optimizer * opt, Mesh::MeshBase &mesh, std::vector< FunctionSpaceEqDistConstr * > &NeighborConstraints ) {
            if( !mesh.Points[i].Neighbors.empty() ){
						
				//NOTE: Pass only the required number of constraints to the 
				//optimizer and then establish their neighbors.
				int iter;
				for(iter=0; iter<mesh.Points[i].Neighbors.size()/2; iter++) {
					std::vector< double > *left = &mesh.Points[ mesh.Points[i].Neighbors[2*iter] ].ObjectiveCoords;
					std::vector< double > *right = &mesh.Points[ mesh.Points[i].Neighbors[2*iter+1] ].ObjectiveCoords;

					NeighborConstraints[iter]->UpdateFrom( left, right );
					Scal.EqualityConstraints.push_back( NeighborConstraints[iter]->function );
				}

				//this->Opt->RefreshConstraints();
                opt->RefreshConstraints();

				//Get Design points
				std::vector< double > x( mesh.Points[i].DesignCoords );

				//Add the lambda values
				int k;
				for(k=0;k<mesh.Points[i].LambdaCoords.size()-1;k++) { 
					x.push_back( mesh.Points[i].LambdaCoords[k] ); 
				}

                //OptNlopt::EXIT_COND flag = OptNlopt::RERUN;
                //while( flag == OptNlopt::RERUN ) flag = (OptNlopt::EXIT_COND)this->Opt->RunFrom( x );
                //Optimizer::EXIT_COND flag = (Optimizer::EXIT_COND) this->Opt->RunFrom( x );
                Optimizer::EXIT_COND flag = (Optimizer::EXIT_COND) opt->RunFrom( x );

				if( flag == OptNlopt::SUCCESS ) {
					//Update design points
					mesh.Points[i].DesignCoords.assign( x.begin(), x.begin() + Prob->dimDesign);

					//Update obj points
					std::vector< double > f( GetF( Prob->Objectives, mesh.Points[i].DesignCoords ) );
					mesh.Points[i].ObjectiveCoords.assign( f.begin(), f.end() );

					//Update lam points
					std::vector< double > l ( x.begin() + Prob->dimDesign, x.end() );
					l.push_back( 1.0 - std::accumulate( l.begin(), l.end(), 0.0 ) );
					mesh.Points[i].LambdaCoords.assign( l.begin(), l.end() );
				}

				//Pop the EQDist Constraints
				Scal.EqualityConstraints.erase( Scal.EqualityConstraints.end() - iter, Scal.EqualityConstraints.end());

                return flag;
			}

            //Neighborless points always succeed
            return Optimizer::SUCCESS;
        }
	};
}