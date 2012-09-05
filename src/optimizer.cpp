#include <vector>

#include <boost/function.hpp>

#include "Problems.h"
#include "Scalarization.hpp"

#include <nlopt.hpp>
#include "NloptAdapt.hpp"

#include <stdio.h>

#include "optimizer.h"

//Nlopt-based optimizer
OptNlopt::OptNlopt(Problem::Interface *p, double tolerance) : 
				P(p), opt(nlopt::LD_SLSQP, P->dimDesign), 
				s( P->Objectives ), NA(s, .000001)	{
	//Set the state of the optimizer:
	//require upper & lower bounds
	opt.set_lower_bounds(P->lowerBounds);
	opt.set_upper_bounds(P->upperBounds);
	opt.set_xtol_rel(tolerance);

	//Pass a scalarized function through the 
	//Nlopt Adapter to the Nlopt object
	opt.set_min_objective(&NloptAdapt< Scalarization<Problem::FUNCTION> >::iface, (void*)(&(this->NA)));
}

double OptNlopt::RunFrom(std::vector< double > &x) {
	double minf;
	nlopt::result result = opt.optimize(x, minf);
	return minf;
}

void OptNlopt::SetWeights(std::vector< double > &x) { s.SetWeights(&x); }