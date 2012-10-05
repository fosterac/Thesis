#include <gtest/gtest.h>

#include <vector>
#include "mesh.hpp"

#include <stdio.h>

namespace {

	using namespace std;

	template< typename T >
		void PrintV(std::vector< T > v){
			int i;
			for(i=0; i<v.size(); i++) { cout << v[i] << "\t"; }
			cout << std::endl;
		}

	class SimplexTest : public ::testing::Test {
	protected:
		//vector< vector<double> > Design;
		virtual void SetUp(){
		}		
	};

	//Verify that the coordinate transforms are
	//stable and reversible
	TEST(SimplexTest, CoordinateTransforms){
		int dim = 8;
		int n = 8;

		int i;
		for(i = 0; i < Mesh::Simplex::eta(dim, n); i++){
			std::vector< int > v( Mesh::Simplex::ind_to_coord( i, dim, n) ) ;
			EXPECT_EQ( i,  Mesh::Simplex::coord_to_ind( v , n ) );
		}
	}
	//Build a 1D mesh (line)
	TEST(SimplexTest, 1D){

		std::vector< std::vector < double > > v, L;
		double vals[3][3] = { {1, 0}, {0, 1} };
		double *a;

		a = &vals[0][0];
		std::vector< double > p;
		p.assign(a, a+2);
		v.push_back(p);

		a = &vals[1][0];
		p.assign(a, a+2);
		v.push_back(p);

		Mesh::Simplex mesh( v, v, v, 5);
		mesh.Print();
	}

	//Build a 2D mesh (trianglar plane)
	TEST(SimplexTest, 2D){

		std::vector< std::vector < double > > D, O, L;
		double Ovals[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 0} };
		double Lvals[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
		double *a;
		double *b;
		int dim = 3;
		
		std::vector< double > p;
		int i;
		for(i=0;i<dim;i++){
			a = &Ovals[i][0];
			p.assign(a, a+dim);
			D.push_back(p);

			p.assign(a, a+dim);
			O.push_back(p);

			a = &Lvals[i][0];
			p.assign(a, a+dim);
			L.push_back(p);
		}

		Mesh::Simplex mesh( D, O, L, 3);
		mesh.Print();
	}

	//Build a 3D mesh (quarter pyramid)
	TEST(SimplexTest, 3D){

		std::vector< std::vector < double > > D, O, L;
		double Ovals[4][4] = { {0, 0, 0, 0}, {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0} };
		double Lvals[4][4] = { {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} };
		double *a;
		double *b;
		int dim = 4;
		
		std::vector< double > p;
		int i;
		for(i=0;i<dim;i++){
			a = &Ovals[i][0];
			p.assign(a, a+dim);
			D.push_back(p);

			p.assign(a, a+dim);
			O.push_back(p);

			a = &Lvals[i][0];
			p.assign(a, a+dim);
			L.push_back(p);
		}

		Mesh::Simplex mesh( D, O, L, 3);
		mesh.Print();
	}
}