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
	//Problem::Interface * P = Problem::Factory("FON", 2, DesignVars);
	Problem::Interface * P = Problem::Factory("CONST_TEST", 1, 2);

	//DynamicScalarization< typename Problem::FUNCTION > S(P);
	FixedScalarization< typename Problem::FUNCTION > S(P);

	Optimizer * op = new OptNlopt(S.f, &S, 1e-4);
	
	printf("starting at: ");
	//std::vector<double> x(S.dimDesign);
	std::vector<double> x(2);
	int i;
	for(i=0; i<x.size(); i++) { 
		x[i] = 0.3;
		printf("%lf ", x[i]);
	}
	printf("\n");

	std::vector<double> w(P->Objectives.size());
	for(i=0; i<w.size(); i++) { w[i] = 0.5; }

	//S.SetWeights(&w);
	double result = op->RunFrom(x);
	
	printf("minimum at: ");
	for(i=0; i<x.size(); i++) { printf("%lf ", x[i]); }

	printf(" =  %lf\n", result);
}

}  //namespace
