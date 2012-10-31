#include <math.h>
#include <numeric>
#include <algorithm>
#include "boost/bind.hpp"
#include "boost/function.hpp"

//for mesh writeout
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <cassert>
#include <stdio.h>

namespace Mesh {

	//Placeholder for a standard interface
	struct Interface {	};

	//Basic function for outputting vector to a stream
	template< typename T>
	std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec){
		os << "[ " ;
		int i;
		for(i=0;i<vec.size()-1;i++){ os << vec[i] << ", " ;	}
		return os << vec[i] << " ]" ;
	}

	class MeshBase : public Interface {
	public:
		class MeshPoint {
		public:
			//The three basic pieces of each meshpoint
			std::vector< double > DesignCoords;
			std::vector< double > ObjectiveCoords;
			std::vector< double > LambdaCoords;

			//Maintain neighbor references
			//std::vector< MeshPoint* > Neighbors;
			//By index in the global array
			std::vector< int > Neighbors;

			MeshPoint(	std::vector< double > &Design, 
						std::vector< double > &Objective, 
						std::vector< double > &Lambdas ) :
						DesignCoords( Design ), 
						ObjectiveCoords( Objective ), 
						LambdaCoords( Lambdas ) { }

			//Simple print function that dumps all values
			void Print() {
				int i;
				printf("Design: ( ");
				for(i=0;i<DesignCoords.size();i++)		printf("%lf ", DesignCoords[i]);
				printf(") Obj: ( ");
				for(i=0;i<ObjectiveCoords.size();i++)	printf("%lf ", ObjectiveCoords[i]);
				printf(") Lam: ( ");
				for(i=0;i<LambdaCoords.size();i++)		printf("%lf ", LambdaCoords[i]);
				printf(") Links: ( ");
				for(i=0;i<Neighbors.size();i++)			printf("%d ", Neighbors[i]);
				printf(")\n");
			}

			//For outputting a meshpoint to the front file
			friend std::ostream& operator<<(std::ostream& os, const MeshPoint& mp){
				return os << mp.ObjectiveCoords;
			}
		};

		std::vector< MeshPoint > Corners;
		std::vector< MeshPoint > Points;

		int DesignDim;
		int ObjectiveDim;
		int MeshDim;

		//Verbose Constructor
		MeshBase(	std::vector< std::vector< double > > DesignSpace,
					std::vector< std::vector< double > > ObjectiveSpace,
					std::vector< std::vector< double > > Lambdas ) : 
					DesignDim( DesignSpace[0].size() ), 
					ObjectiveDim( ObjectiveSpace[0].size() ),
					MeshDim( ObjectiveDim - 1 )
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
		virtual void Print() {
			int i;
			printf("Corners: \n");
			for(i=0;i<Corners.size();i++) Corners[i].Print();
			printf("Points: \n");
			for(i=0;i<Points.size();i++) Points[i].Print();
		}
		void WriteOut(const char* filename) {
			std::ofstream os( filename );
			os << "{" << std::endl;

			std::vector<std::string> cols;
			int i;
			for(i=0;i<this->ObjectiveDim;i++) {
				std::stringstream s;
				s << "\"" << i << "\"" ;
				cols.push_back( std::string( s.str() ) );
			}
			os << "\"columns\" : " << cols << "," << std::endl;

			os << "\"data\" : ";
			os << this->Points;

			os << "}" << std::endl;
			os.close();
		}
	};

	class Simplex : public MeshBase {
	private:
	public:
		static int eta(int dim, int n) {
			return eta_num(dim, n) / eta_den(dim);
		}
		static int eta_num(int dim, int n) {
			if (dim == 1) return n;
			else return (n+dim -1) * eta_num(dim-1, n);
		}
		static int eta_den(int dim) {
			if (dim == 1) return 1;
			else return dim * eta_den(dim - 1);
		}
		static int coord_to_ind(std::vector< int > &coords, int n){
			return coord_to_ind_aux(coords, n, 0);
		}
		static int coord_to_ind_aux(std::vector< int > &coords, int n, int ind) {
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
		}
		static bool isBoundaryPoint( std::vector< int > &coords, int n) {
			return ( std::accumulate( coords.begin(), coords.end(), 0 ) == (n - 1) );
		}
		static std::vector< std::vector< int > > getNeighborsAux( std::vector< int > coords, int n){
			std::vector< std::vector< int > > result;

			if( ! isBoundaryPoint( coords, n ) ) {
				int i;
				for(i=0;i<coords.size();i++) {
					if( coords[i] == 0 ) continue;

					coords[i] += 1;
					result.push_back( coords );					
					coords[i] -= 2;
					result.push_back( coords );
					coords[i] += 1;
				}
			}
			else{
				coords.pop_back();

				if( !coords.empty() ) {
					result = getNeighborsAux( coords, n );
					int i;
					for(i=0;i<result.size();i++){
						int sum = std::accumulate( result[i].begin(), result[i].end(), 0);
						result[i].push_back( n - 1 - sum );
					}
				}
			}
			return result;
		}
		static std::vector< int > getNeighbors( std::vector< int > coords, int n ){
			std::vector< int > neighbors;
			std::vector< std::vector< int > > NCoords(getNeighborsAux(coords, n) );
			int i;
			for(i=0;i<NCoords.size();i++){ neighbors.push_back( coord_to_ind(NCoords[i], n) ); }
			
			return neighbors;
		}

	//public:
        std::set< int > GetNeighborIDs( int ID ){
            //Get a vector of neighbors
            std::vector< int > vec (getNeighbors( ind_to_coord( ID, this->MeshDim, this->PointsPerSide ), this->PointsPerSide ));
            //Copy to a set
            std::set< int > retval(vec.begin(), vec.end()); 
            //std::vector< int >::iterator i;
            //for(i=vec.begin(); i!=vec.end(); i++) retval.insert( *i );
            return retval;
        }

		int PointsPerSide;

		Simplex(	std::vector< std::vector< double > > DesignSpace,
					std::vector< std::vector< double > > ObjectiveSpace,
					std::vector< std::vector< double > > Lambdas, 
					int NumberOfPoints ) : MeshBase( DesignSpace, ObjectiveSpace, Lambdas),
					PointsPerSide( NumberOfPoints ) {
			
			//Generate points
			baryTransFunc baryTrans( boost::bind( &barycentric_trans_helper, 1.0/(this->PointsPerSide-1), _1 ) );
			int i;
			for(i=0; i<eta(this->MeshDim, this->PointsPerSide); i++){

				//Get the point n-index
				std::vector< int > v( Mesh::Simplex::ind_to_coord( i, this->MeshDim, this->PointsPerSide) ) ;

				//Get the neighbor indicies
				std::vector< int > neighbors ( Mesh::Simplex::getNeighbors( v, this->PointsPerSide ) );

				//Generate the barycentric coordinates of this point
				v.push_back( PointsPerSide - 1 - std::accumulate(v.begin(), v.end(), 0) );
				std::vector< double > barycentric( v.size() );
				std::transform( v.begin(), v.end(), barycentric.begin(), baryTrans );

				//Generate point location in various spaces from
				//barycentric coordinates
				std::vector< double > design( this->DesignDim, 0.0 );
				std::vector< double > objective( this->ObjectiveDim, 0.0 );
				std::vector< double > lambda( this->ObjectiveDim, 0.0 );
				int j;
				for(j=0; j<barycentric.size(); j++){
					//printf("%lf\t", barycentric[j]);
					boost::function<double ( double, double )> comb ( 
						boost::bind( &barycentric_trans_helper_2, 
						barycentric[j], _1, _2) );
					std::transform( Corners[j].DesignCoords.begin(), Corners[j].DesignCoords.end(), design.begin(), design.begin(), comb);
					std::transform( Corners[j].ObjectiveCoords.begin(), Corners[j].ObjectiveCoords.end(), objective.begin(), objective.begin(), comb);
					std::transform( Corners[j].LambdaCoords.begin(), Corners[j].LambdaCoords.end(), lambda.begin(), lambda.begin(), comb);
				}
				//printf("\n");

				MeshPoint m( design, objective, lambda );
				m.Neighbors = neighbors;

				this->Points.push_back( m );
			}
		}
	};
}