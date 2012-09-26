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
	void PrintF(std::vector< double > &x, std::vector< Problem::FUNCTION > &f){ 
		Print(x);
		printf(" -> ");
		Print(GetF(f, x));	
		printf("\n");
	}

	TEST(MARCHINGHOMOTOPY, 2D){

		int Points = 146;
		double step = 0.01;

		//Set up the objects
		Problem::Interface * P = Problem::Factory("FON", 2, 3);

		//Before we can get this working, we need to find a way
		//(decorator?) to reliably strip the weights from the 
		//solution vector.  
		//Problem::Interface * P = Problem::Factory("WFG2", 2, 4);

		FixedScalarization< typename Problem::FUNCTION > S(P);
		
		//Get the starting point
		Optimizer * op = new OptNlopt(P->Objectives[0], &S, 1e-4);
		std::vector<double> x1(P->dimDesign, 0.3);
		op->RunFrom(x1);
		PrintF(x1, P->Objectives);

		DynamicScalarization< typename Problem::FUNCTION > D(P);
		//StepConstraint< typename Problem::FUNCTION > C(NULL, step);
		FStepConstraint< typename Problem::FUNCTION > C(P->Objectives, NULL, step);
		D.EqualityConstraints.push_back( C.function );

		op = new OptNlopt(D.f, &D, 1e-4);
		x1.push_back(1.0);
		std::vector<double> x(x1);

		//Get the rest of the points
		int i;
		for(i=0; i<Points; i++){
			std::vector< double > last_x(x.begin(), x.begin() + P->dimDesign);
			std::vector< double > last_f( GetF(P->Objectives, last_x) );

			//Build the additional constraint
			C.UpdateFrom( &last_f );
			//C.UpdateFrom( &last_x );

			op->RunFrom(x);
			PrintF(x, P->Objectives);

			//printf("DistErr: %lf\n", FDist(P->Objectives, 0, last_x, x));
		}
	}
}