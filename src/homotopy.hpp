namespace Pareto {

	void PrintVec( std::vector< double > &vec ){
		int i;
		for(i=0;i<vec.size();i++){ 
			printf("%lf ", vec[i] );
		}
		printf("\n");
	}

	class Homotopy {
	private:
		Problem::Interface *Prob;
		//DynamicScalarization< typename Problem::FUNCTION > Scal;
        DynamicScalarization< Evaluator< EvaluationStrategy::Cached< EvaluationStrategy::Local< functionSet_t > > > > Scal;
		Optimizer * Opt;
		//OptNlopt * Opt;

		double tolerance;

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
		Homotopy( Problem::Interface *P, double tolerance) : Prob(P), Scal( Prob ), tolerance(tolerance), 
			//Opt( OptNlopt(Scal.f, &Scal, tolerance) ) {
			Opt( NULL ) {
				
				//Find the individual optimae using a fixed scalarization
				//FixedScalarization< typename Problem::FUNCTION > S(Prob);
                FixedScalarization< Evaluator<EvaluationStrategy::Local< functionSet_t > > > S(Prob);
                FiniteDifferences::Params_t fd_par;
                fd_par.step = 1e-6;
                fd_par.type = FiniteDifferences::CENTRAL;
				Optimizer * op = new OptNlopt(S.f, &S, 1e-4, fd_par);

				int i;
				for( i=0; i<Prob->Objectives.size(); i++){
					//Optimize the objective
					std::vector<double> lam(Prob->Objectives.size(), 0.0);
					lam[i] = 1.0;
					S.SetWeights( &lam );
					std::vector<double> x(Prob->dimDesign, (double)(i+1) / (Prob->Objectives.size()+2));
					op->RunFrom( x );
					std::vector<double> f( GetF( Prob->Objectives, x) );

					//Cache the results
					Design.push_back( x );
					Objective.push_back( f );
					Lambda.push_back( lam );
				}
				delete op;
		}

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
				mesh.Points[ m ].ObjectiveCoords = GetF( Prob->Objectives, mesh.Points[ m ].DesignCoords ) ;
			}

			//Constraints

			//Sum constraints on the Lambda parameters
			BoundSumConstraint< typename Problem::FUNCTION > LT( BoundSumConstraint< typename Problem::FUNCTION >::LESS_THAN, 1.0, Prob->dimDesign, Scal.dimDesign);
			BoundSumConstraint< typename Problem::FUNCTION > GT( BoundSumConstraint< typename Problem::FUNCTION >::GREATER_THAN, 0.0, Prob->dimDesign, Scal.dimDesign);
			Scal.InequalityConstraints.push_back( LT.function );
			Scal.InequalityConstraints.push_back( GT.function );

			//Incorporate constraint
			std::vector< FEqDistanceConstraint< typename Problem::FUNCTION >* > NeighborConstraints( mesh.MeshDim );
			int n;
			for(n=0;n<mesh.MeshDim;n++) {
				NeighborConstraints[n] = new FEqDistanceConstraint< typename Problem::FUNCTION > (Prob->Objectives, NULL, NULL) ;
				//Scal.EqualityConstraints.push_back( NeighborConstraints[n]->function );
			}

            //Establish finite difference parameters
            FiniteDifferences::Params_t FDpar = { 1e-6, FiniteDifferences::CENTRAL };

			//Construct optimizer
			this->Opt = new OptNlopt(Scal.f, &Scal, tolerance, fd_par);

			//Run a set of updates
			int j;
			for(j=0; j<Iterations; j++){

				//for each point
				int i;
				//for(i=j%2; i<NumPoints; i+=2){
				for(i=0; i<mesh.Points.size(); i++){
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

						this->Opt->RefreshConstraints();

						//Get Design points
						std::vector< double > x( mesh.Points[i].DesignCoords );

						//Add the lambda values
						int k;
						for(k=0;k<mesh.Points[i].LambdaCoords.size()-1;k++) { 
							x.push_back( mesh.Points[i].LambdaCoords[k] ); 
						}

						if( this->Opt->RunFrom( x ) ) {
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
					}
				}
			}
			delete this->Opt; 
			
			mesh.Print();			
			mesh.WriteOut( "front.txt" );	
		}
	};
}