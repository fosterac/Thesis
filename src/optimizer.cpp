#include <vector>

#include <boost/function.hpp>

#include "Problems.h"
#include "Scalarization.hpp"
#include <nlopt.hpp>
#include "NloptAdapt.hpp"
#include "optimizer.h"

#include <stdio.h>

//Nlopt-based optimizer
OptNlopt::OptNlopt(Problem::FUNCTION &Obj, Scalarization< typename Problem::FUNCTION > *s, double tolerance) : 
				S(s), NA(Obj, &S->EqualityConstraints, &S->InequalityConstraints, 1e-6), 
				opt(nlopt::LD_SLSQP, S->dimDesign), 
				EqTolerances(S->EqualityConstraints.size(), tolerance),
				InEqTolerances(S->InequalityConstraints.size(), tolerance){

	//Set the state of the optimizer:
	
	//Boundary values
	if (!S->lowerBounds.empty()) opt.set_lower_bounds(S->lowerBounds);
	if (!S->upperBounds.empty()) opt.set_upper_bounds(S->upperBounds);
	
	//Set the stop conditions
	//this requires some attention
	opt.set_xtol_rel(tolerance);
	opt.set_ftol_abs(tolerance);

	//Pass a scalarized function through the 
	//Nlopt Adapter to the Nlopt object
	opt.set_min_objective(&NloptAdapt< typename Problem::FUNCTION >::ObjIface, (void*)(&(this->NA)));

	//Likewise pass the constraints to the optimizer
	if ( !S->EqualityConstraints.empty() ) {
		opt.add_equality_mconstraint(&NloptAdapt< typename Problem::FUNCTION >::EqConstrIface, (void*)(&(this->NA)), this->EqTolerances);
	}
	if ( !S->InequalityConstraints.empty() ) {
		opt.add_inequality_mconstraint(&NloptAdapt< typename Problem::FUNCTION >::InEqConstrIface, (void*)(&(this->NA)), this->InEqTolerances);
	}
}

/*double OptNlopt::RunFrom(std::vector< double > &x) {
	double minf;
	nlopt::result result = opt.optimize(x, minf);
	printf("NLOPT status: %d\n", result);
	return minf;
}*/

int OptNlopt::RunFrom(std::vector< double > &x) {
	int result = 0;
	double minf;
	try {
		//nlopt::result result = opt.optimize(x, minf);
		result = opt.optimize(x, minf);
	}
	//catch (const std::exception& ex) {
	catch (const std::runtime_error& ex) {
		result = 0;
	}

	return result;
}
