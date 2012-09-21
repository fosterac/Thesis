#include <math.h>
#include <numeric>
#include <algorithm>
#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <cassert>
#include <stdio.h>

namespace Mesh {

	struct Interface {

	};

	class MeshBase : public Interface {
	public:
		class MeshPoint {
		public:
			std::vector< double > DesignCoords;
			std::vector< double > ObjectiveCoords;

			std::vector< MeshPoint* > Neighbors;

			std::vector< double > LambdaCoords;

			MeshPoint(	std::vector< double > &Design, 
						std::vector< double > &Objective, 
						std::vector< double > &Lambdas ) :
						DesignCoords( Design ), 
						ObjectiveCoords( Objective ), 
						LambdaCoords( Lambdas ) { }
			void Print() {
				int i;
				printf("Design: (");
				for(i=0;i<DesignCoords.size();i++) printf("%lf ", DesignCoords[i]);
				printf(") Obj: (");
				for(i=0;i<ObjectiveCoords.size();i++) printf("%lf ", ObjectiveCoords[i]);
				printf(") Lam: (");
				for(i=0;i<LambdaCoords.size();i++) printf("%lf ", LambdaCoords[i]);
				printf(")\n");
			}
		};

		std::vector< MeshPoint > Corners;
		std::vector< MeshPoint > InternalPoints;

		int DesignDim;
		int ObjectiveDim;

		//Verbose Constructor
		MeshBase(	std::vector< std::vector< double > > DesignSpace,
					std::vector< std::vector< double > > ObjectiveSpace,
					std::vector< std::vector< double > > Lambdas ) : 
					DesignDim( DesignSpace[0].size() ), 
					ObjectiveDim( ObjectiveSpace[0].size() )
		{
			
			//verify that we have the right number of points
			assert( DesignSpace.size() == ObjectiveSpace.size() );
			assert( ObjectiveSpace.size() == Lambdas.size() );

			int i;
			for(i=0; i<DesignSpace.size(); i++) {

				//verify the dimensionality of each point
				assert( DesignDim == DesignSpace[i].size() );
				assert( ObjectiveDim == ObjectiveSpace[i].size() );
				assert( ObjectiveDim == Lambdas[i].size() );

				MeshPoint p( DesignSpace[i], ObjectiveSpace[i], Lambdas[i] );
				this->Corners.push_back( p );
			}
		}
		void Print() {
			int i;
			printf("Corners: \n");
			for(i=0;i<Corners.size();i++) Corners[i].Print();
			printf("Points: \n");
			for(i=0;i<InternalPoints.size();i++) InternalPoints[i].Print();
		}
	};

	class Simplex : public MeshBase {
	private:
	public:
		static double eta(int dim, int n) {
			return eta_num(dim, n) / eta_den(dim);
		}
		static double eta_num(int dim, int n) {
			if (dim == 1) return n;
			else return (n+dim -1) * eta_num(dim-1, n);
		}
		static double eta_den(int dim) {
			if (dim == 1) return 1;
			else return dim * eta_den(dim - 1);
		}
		static double coord_to_ind(std::vector< int > &coords, int n){
			return coord_to_ind_aux(coords, n, 0);
		}
		static double coord_to_ind_aux(std::vector< int > &coords, int n, int ind) {
			if ( ind == coords.size()-1 ) return coords[ind];
			else {
				int this_coord = coords[ind];
				int dim = coords.size() - ind;
				return ( eta(dim, n) - eta(dim, n - this_coord) ) + coord_to_ind_aux(coords, n-this_coord, ind+1);
			}
		}

		typedef boost::function<int ( int )> func;
		static int bisection_search(func & f, int a, int b){
			assert ( (f)(a) * (f)(b) <= 0 );
			if( (f)(a) == 0 ) return a;
			if( (f)(b) == 0 ) return b;
			while ( fabs(b - a) > 1.0 ) {
				int sum = a+b;
				int c = ((sum)/2) + ( (sum % 2) ? 1 : 0 ) ;
				if( (f)(c) == 0 ) return c;
				if( (f)(a) * (f)(c) <= 0 ) b = c;
				else a = c;
			}
			return a;
		}
		static int eta_helper(int dim, int ind, int x){ return eta(dim, x) - ind; }
		static std::vector< int > ind_to_coord(int ind, int dim, int n){
			//flip the coordinate system
			ind = eta(dim, n) - 1 - ind;

			std::vector< int > coords(dim);
			int i;
			for(i=0; i<dim; i++){
				func f ( boost::bind( &eta_helper, dim-i, ind, _1 ) );
				int c = bisection_search(f, 0, n);
				coords[i] = c;
				ind -= eta(dim - i, c);
			}

			std::vector< int >::iterator it;
			for(it = coords.begin(); it != coords.end(); it++){
				*it = n - std::accumulate(coords.begin(), it, 0) - *it -1;
			}
			return coords;
		}

		typedef boost::function<double ( int )> baryTransFunc;
		static double barycentric_trans_helper(double n, int x){ return x*n; }
		static double barycentric_trans_helper_2(double n, double x, double y){ 
			return y + n*x; 
			//return y;
		}
		//double barycentric_trans_helper(int x){ return double(x)/this->PointsPerSide; }
	//public:
		int PointsPerSide;

		Simplex(	std::vector< std::vector< double > > DesignSpace,
					std::vector< std::vector< double > > ObjectiveSpace,
					std::vector< std::vector< double > > Lambdas, 
					int NumberOfPoints ) : MeshBase( DesignSpace, ObjectiveSpace, Lambdas),
					PointsPerSide( NumberOfPoints ) {
			
			//Generate points
			baryTransFunc baryTrans( boost::bind( &barycentric_trans_helper, 1.0/this->PointsPerSide, _1 ) );
			int i;
			for(i=0; i<this->PointsPerSide; i++){

				//Get the point n-index
				std::vector< int > v( Mesh::Simplex::ind_to_coord( i, this->Corners.size(), this->PointsPerSide) ) ;

				//Generate the barycentric coordinates of this point
				v.push_back( PointsPerSide - std::accumulate(v.begin(), v.end(), 0) );
				std::vector< double > barycentric( v.size() );
				std::transform( v.begin(), v.end(), barycentric.begin(), baryTrans );

				//Generate point location in various spaces from
				//barycentric coordinates
				std::vector< double > design( this->DesignDim, 0.0 );
				std::vector< double > objective( this->ObjectiveDim, 0.0 );
				std::vector< double > lambda( this->ObjectiveDim, 0.0 );

				int j;
				for(j=0; j<barycentric.size(); j++){
					printf("%lf\n", barycentric[j]);
					boost::function<double ( double, double )> comb ( 
						boost::bind( &barycentric_trans_helper_2, 
						barycentric[j], _1, _2) );
					std::transform( Corners[j].DesignCoords.begin(), Corners[j].DesignCoords.end(), design.begin(), design.begin(), comb);
					std::transform( Corners[j].ObjectiveCoords.begin(), Corners[j].ObjectiveCoords.end(), objective.begin(), objective.begin(), comb);
					std::transform( Corners[j].LambdaCoords.begin(), Corners[j].LambdaCoords.end(), lambda.begin(), lambda.begin(), comb);
				}
				MeshPoint m( design, objective, lambda );
				m.Print();
				this->InternalPoints.push_back( m );
			}
		}
	};

}