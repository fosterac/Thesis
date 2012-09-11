#include <vector>

#include <boost/function.hpp>

#include "Problems.h"
#include "Scalarization.hpp"

#include <nlopt.hpp>
#include "NloptAdapt.hpp"

#include <stdio.h>

#include "optimizer.h"

//Nlopt-based optimizer
OptNlopt::OptNlopt(Scalarization * s, double tolerance) : 
				S(s), NA(*S, .000001), opt(nlopt::LD_SLSQP, S->dimDesign) {
	//Set the state of the optimizer:
	//require upper & lower bounds
	opt.set_lower_bounds(S->lowerBounds);
	opt.set_upper_bounds(S->upperBounds);
	opt.set_xtol_rel(tolerance);

	//Pass a scalarized function through the 
	//Nlopt Adapter to the Nlopt object
	opt.set_min_objective(&NloptAdapt< Scalarization >::iface, (void*)(&(this->NA)));
}

double OptNlopt::RunFrom(std::vector< double > &x) {
	double minf;
	nlopt::result result = opt.optimize(x, minf);
	return minf;
}
