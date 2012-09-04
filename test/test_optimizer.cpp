//Simple test to verify Optimizer class functionality

#include <gtest/gtest.h>

#include <stdio.h>

#include <vector>
#include <boost/function.hpp>
#include "Problems.h"

#include "Scalarization.hpp"

#include <nlopt.hpp>
#include "NloptAdapt.hpp"


#include "optimizer.h"

namespace {

// Tests that the optimizer works
TEST(OPTIMIZERTest, Alive) {
	
	int DesignVars = 3;

	//Problem::Interface * P = Problem::Factory("BASIN", 1, DesignVars);
	Problem::Interface * P = Problem::Factory("FON", 2, DesignVars);
	
	Optimizer * op = new OptNlopt(P, 1e-4);
	
	std::vector<double> x(DesignVars);
	int i;
	for(i=0; i<x.size(); i++) { x[i] = 0.3; }
	
	printf("starting at (%lf,%lf,%lf) \n", x[0], x[1], x[2]);

	std::vector<double> w(DesignVars);
	for(i=0; i<w.size(); i++) { w[i] = 1.0; }

	op->SetWeights(w);
	double result = op->RunFrom(x);
	
	printf("minimum at (%lf,%lf,%lf) ", x[0], x[1], x[2]);
	printf(" =  %lf\n", result);
}

//Sample a simple Pareto front
TEST(OPTIMIZERTest, Shotgun) {
	
	//Choose a Problem
	int DesignVars = 3;
	Problem::Interface * P = Problem::Factory("FON", 2, DesignVars);
	//Problem::Interface * P = Problem::Factory("WFG2", 2, DesignVars);
	
	//Instantiate an optimizer
	Optimizer * op = new OptNlopt(P, 1e-4);

	//Set up an initial point
	std::vector<double> x(DesignVars);
	int i;
	for(i=0; i<DesignVars; i++) { x[i] = 0.5; }

	//Set up the weight vector
	std::vector<double> w(P->Objectives.size());

	//Sample the curve
	int NumPts = 10;
	int j;
	for(j=0; j<NumPts; j++){

		//Set the weights
		w[0] =  j* 1.0/(double)(NumPts-1);
		w[1] = 1.0 - w[0];
		op->SetWeights(w);
		
		//Re-initialize the starting points
		for(i=0; i<DesignVars; i++) { x[i] = 0.5; }

		//Run the Optimization
		double result = op->RunFrom(x);

		//Print the results
		printf("Point %d: \t %lf, %lf -> (%lf, %lf) = %lf\n", j, w[0], w[1], (P->Objectives[0])(x), (P->Objectives[1])(x), result);
	}
}

}  //namespace
