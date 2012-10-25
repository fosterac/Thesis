//Simple test to verify Optimizer class functionality

#include "gtest/gtest.h"

#include <stdio.h>

#include <vector>
#include "boost/function.hpp"

#include "Problems.h"

#include "evaluator.hpp"
#include "Scalarization.hpp"
using namespace Homotopy;

#include "nlopt.hpp"
#include "NloptAdapt.hpp"


#include "optimizer.hpp"

namespace {

// Tests that the optimizer works
TEST(OPTIMIZERTest, Alive) {
	
	int DesignVars = 3;

	//Problem::Interface * P = Problem::Factory("BASIN", 1, DesignVars);
	Problem::Interface * P = Problem::Factory("FON", 2, DesignVars);
	//Problem::Interface * P = Problem::Factory("CONST_TEST", 1, 2);

	//DynamicScalarization< typename Problem::FUNCTION > S(P);
	//FixedScalarization< typename Problem::FUNCTION > S(P);
    FixedScalarization< Evaluator<EvaluationStrategy::Local< functionSet_t > > > S(P);

    FiniteDifferences::Params_t fd_par;
    fd_par.step = 1e-6;
    fd_par.type = FiniteDifferences::CENTRAL;

	Optimizer * op = new OptNlopt(S.f, &S, 1e-4, fd_par);
	
	printf("starting at: ");
	std::vector<double> x(S.dimDesign);
	//std::vector<double> x(2);
	int i;
	for(i=0; i<x.size(); i++) { 
		x[i] = 0.3;
		printf("%lf ", x[i]);
	}
	printf("\n");

	std::vector<double> w(P->Objectives.size());
	for(i=0; i<w.size(); i++) { w[i] = 0.5; }

	//S.SetWeights(&w);
	op->RunFrom(x);
    double result = S(x);
	
	printf("minimum at: ");
	for(i=0; i<x.size(); i++) { printf("%lf ", x[i]); }

	printf(" =  %lf\n", result);
}

}  //namespace
