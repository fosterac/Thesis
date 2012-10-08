namespace Pareto {
	class Homotopy {
	private:
		Problem::Interface *Prob;
		DynamicScalarization< typename Problem::FUNCTION > Scal;
		Optimizer * Opt;

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
				FixedScalarization< typename Problem::FUNCTION > S(Prob);
				Optimizer * op = new OptNlopt(S.f, &S, 1e-4);

				int i;
				for( i=0; i<Prob->Objectives.size(); i++){
					//Optimize the objective
					std::vector<double> lam(Prob->Objectives.size(), 0.0);
					lam[i] = 1.0;
					S.SetWeights( &lam );
					std::vector<double> x(Prob->dimDesign, 0.5);
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

			//Incorporate constraint
			FEqDistanceConstraint< typename Problem::FUNCTION > C(Prob->Objectives, NULL, NULL);
			Scal.EqualityConstraints.push_back( C.function );

			//Construct optimizer
			this->Opt = new OptNlopt(Scal.f, &Scal, tolerance);

			//Run a set of updates
			int j;
			for(j=0; j<Iterations; j++){

				//for each point
				int i;
				for(i=j%2; i<NumPoints; i+=2){
				//for(i=0; i<NumPoints; i++){
					if( !mesh.Points[i].Neighbors.empty() ){
						
						//Get the neighbor locations
						std::vector< double > *left		= &mesh.Points[ mesh.Points[i].Neighbors[0] ].ObjectiveCoords;
						std::vector< double > *right	= &mesh.Points[ mesh.Points[i].Neighbors[1] ].ObjectiveCoords;

						//Update equidistant constraint
						C.UpdateFrom( left, right );

						//Get Design points
						std::vector< double > x( mesh.Points[i].DesignCoords );
						//Add the lambda values
						int k;
						for(k=0;k<mesh.Points[i].LambdaCoords.size()-1;k++) { 
							x.push_back( mesh.Points[i].LambdaCoords[k] ); 
						}

						if( Opt->RunFrom( x ) ) {
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
					}
				}
			}
			delete this->Opt; 
			mesh.Print();
		}
	};
}