//Basic marching homotopy test to validate additional
//problem constraints

#include <gtest/gtest.h>

#include <vector>
#include <boost/function.hpp>

#include "Problems.h"
#include "Scalarization.hpp"
#include "Constraint.hpp"
#include <nlopt.hpp>
#include "NloptAdapt.hpp"
#include "optimizer.h"
#include "mesh.hpp"

#include <stdio.h>

namespace {

	void Print(const std::vector< double > &x){
		int i;
		for(i=0; i<x.size(); i++) {printf("%lf ", x[i]);}
	}
	std::vector< double > GetF(const std::vector< Problem::FUNCTION > &f, const std::vector< double > &x){
		std::vector< double > result;
		int i;
		for(i=0; i<f.size(); i++) { result.push_back((f[i])(x)); }
		return result;
	}
	double FDist(const std::vector< Problem::FUNCTION > &f, double step, const std::vector< double > &y, const std::vector< double > &x){
		std::vector< double > fy, fx;
		fy = GetF(f, y);
		fx = GetF(f, x);
		int i;
		double dist = 0.0;
		for(i=0; i<f.size(); i++) {
			dist += (fy[i] - fx[i]) * (fy[i] - fx[i]);
		}
		return fabs( sqrt(dist) - step );
	}
	double Dist(double step, const std::vector< double > &y, const std::vector< double > &x){
		int i;
		double dist = 0.0;
		for(i=0; i<y.size(); i++) {
			dist += (y[i] - x[i]) * (y[i] - x[i]);
		}
		return fabs( sqrt(dist) - step );
	}
	void PrintF(std::vector< double > &x, std::vector< Problem::FUNCTION > &f){ 
		Print(x);
		printf(" -> ");
		Print(GetF(f, x));	
		printf("\n");
	}

	TEST(MARCHINGHOMOTOPY, 2D){

		//Set up the problem
		Problem::Interface * P = Problem::Factory("FON", 2, 3);
		
		//Get the trade-off endpoints
		FixedScalarization< typename Problem::FUNCTION > S(P);
		Optimizer * op = new OptNlopt(S.f, &S, 1e-4);

		//Objective one
		std::vector<double> w1(P->Objectives.size());
		w1[0] = 1.0; w1[1] = 0.0;
		S.SetWeights(&w1);
		std::vector<double> x1(P->dimDesign, 0.5);
		op->RunFrom(x1);
		std::vector<double> f1( GetF( P->Objectives, x1) );
		//Print(x1);

		//Objective two
		std::vector<double> w2(P->Objectives.size());
		w2[0] = 0.0; w2[1] = 1.0;
		S.SetWeights(&w2);
		std::vector<double> x2(P->dimDesign, 0.5);
		op->RunFrom(x2);
		std::vector<double> f2( GetF( P->Objectives, x2) );
		//Print(x2);

		//Setup the mesh
		int NumPoints = 146;
		std::vector< std::vector< double > > Des;
		Des.push_back( x1 ); Des.push_back( x2 );
		std::vector< std::vector< double > > Obj;
		Obj.push_back( f1 ); Obj.push_back( f2 );
		std::vector< std::vector< double > > Lam;
		Lam.push_back( w1 ); Lam.push_back( w2 );
		
		Mesh::Simplex mesh( Des, Obj, Lam, NumPoints);

		//Optimize the meshpoints
		DynamicScalarization< typename Problem::FUNCTION > D(P);
		FEqDistanceConstraint< typename Problem::FUNCTION > C(P->Objectives, NULL, NULL);
		D.EqualityConstraints.push_back( C.function );

		op = new OptNlopt(D.f, &D, 1e-6);

		//Run a set of updates
		int j;
		for(j=0; j<5; j++){
			//for each point
			int i;
			for(i=0; i<NumPoints; i++){
				if( !mesh.Points[i].Neighbors.empty() ){
					
					//Get the neighbor locations

					/*
					TODO: This point needs further investigation as starting from
					the Obj Space guess location seems to be a winning prospect for 
					quickly converging on convex fronts.  However, in the real FON 
					example, leads to roundoff error for large #'s or points.
					Starting, however, from the real locations of the interpolated 
					Pareto Set (not front) allows much larger mesh sizes without 
					roundoff errors and contains points rooted in some kind of truth.
					*/

					//std::vector< double > *left		= &mesh.Points[ mesh.Points[i].Neighbors[0] ].ObjectiveCoords;
					//std::vector< double > *right	= &mesh.Points[ mesh.Points[i].Neighbors[1] ].ObjectiveCoords;
					std::vector< double > *left		= new std::vector<double> ( GetF( P->Objectives, mesh.Points[ mesh.Points[i].Neighbors[0] ].DesignCoords ) );
					std::vector< double > *right	= new std::vector<double> ( GetF( P->Objectives, mesh.Points[ mesh.Points[i].Neighbors[1] ].DesignCoords ) );

					//Update equidistant constraint
					C.UpdateFrom( left, right );

					std::vector< double > x( mesh.Points[i].DesignCoords );
					x.push_back( mesh.Points[i].LambdaCoords[0] );

					op->RunFrom( x );

					delete left;
					delete right;

					//Update design points
					mesh.Points[i].DesignCoords.assign( x.begin(), x.end() - 1);

					//Update obj points
					std::vector< double > f( GetF( P->Objectives, mesh.Points[i].DesignCoords ) );
					mesh.Points[i].ObjectiveCoords.assign( f.begin(), f.end() );

					//Update lam points
					mesh.Points[i].LambdaCoords[0] = x.back();
					mesh.Points[i].LambdaCoords[1] = 1 - x.back();
				}
			}
		}

		mesh.Print();
	}
}