//Simple test to verify Problem class functionality

#include <gtest/gtest.h>

#include <math.h>
#include <cstdio>

#include <nlopt.hpp>

#include <boost/function.hpp>
#include "Problems.h"

//namespace {

// Tests that the problem interface works
#include "NloptAdapt.hpp"

TEST(PROBLEMSTest, SingleObjective) {

	int DesignVars = 3;

	Problem::Interface * P = Problem::Factory("BASIN", 1, DesignVars);
	//Problem::Interface * P = Problem::Factory("WFG2", 2, DesignVars);

	nlopt::opt opt(nlopt::LD_SLSQP, DesignVars);

	opt.set_lower_bounds(P->lowerBounds);
	opt.set_upper_bounds(P->upperBounds);

	NloptAdapt<Problem::FUNCTION> NA(P->Objectives[0], .000001);

	opt.set_min_objective(&NloptAdapt<typename Problem::FUNCTION>::iface, (void*)&NA);

	opt.set_xtol_rel(1e-4);

	std::vector<double> x(DesignVars);
	int i;
	for(i=0; i<x.size(); i++) { x[i] = 0.3; }
	double minf;
	nlopt::result result = opt.optimize(x, minf);
	printf("found minimum at f(%lf,%lf) = %lf\n", x[0], x[1], minf);
	}

#include "Scalarization.hpp"

TEST(PROBLEMSTest, MultiObjective) {

	int DesignVars = 3;

	//Problem::Interface * P = Problem::Factory("BASIN", 1, DesignVars);
	Problem::Interface * P = Problem::Factory("FON", 2, DesignVars);

	nlopt::opt opt(nlopt::LD_SLSQP, DesignVars);
	opt.set_lower_bounds(P->lowerBounds);
	opt.set_upper_bounds(P->upperBounds);

	Scalarization<Problem::FUNCTION> s( P->Objectives );
	std::vector<double> weights;
	weights.push_back(1.0);
	weights.push_back(1.0);
	s.SetWeights(&weights);

	NloptAdapt< Scalarization<Problem::FUNCTION> > NA(s, .000001);

	opt.set_min_objective(&NloptAdapt< Scalarization<Problem::FUNCTION> >::iface, (void*)&NA);

	opt.set_xtol_rel(1e-4);

	std::vector<double> x(DesignVars);
	int i;
	for(i=0; i<x.size(); i++) { x[i] = 0.3; }
	printf("starting at (%lf,%lf,%lf) \n", x[0], x[1], x[2]);
	double minf=0;
	nlopt::result result = opt.optimize(x, minf);
	printf("found minimum at f(%lf,%lf,%lf) = %lf\n", x[0], x[1], x[2], minf);

	
	}
//}  //namespace
