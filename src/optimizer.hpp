#ifndef optimizer_hpp
#define optimizer_hpp

#include "HomotopyTypes.h"

using namespace Homotopy;

class Optimizer {
public:
	//virtual double RunFrom(std::vector< double > &) = 0;
	virtual int RunFrom(std::vector< double > &) = 0;
	virtual void RefreshConstraints() = 0;
};

//NloptBased optimizer
class OptNlopt : public Optimizer {
private:
	ScalarizationInterface *S;
	NloptAdapt< typename Problem::FUNCTION > NA;
	nlopt::opt opt;

	double tolerance;
	std::vector< double > EqTolerances;
	std::vector< double > InEqTolerances;

public:
    OptNlopt(function_t &Obj, ScalarizationInterface *s, double tolerance, FiniteDifferences::Params_t fd_par) : 
				    Optimizer(), S(s), NA(Obj, &S->EqualityConstraints, &S->InequalityConstraints, fd_par), 
				    opt(nlopt::LD_SLSQP, S->dimDesign), tolerance(tolerance), 
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

	int RunFrom(std::vector< double > &x) {
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

    void RefreshConstraints(){
	    //printf("Called right func\n");
	    opt.remove_equality_constraints();
	    opt.remove_inequality_constraints();
	    //printf("Called right func\n");
	    std::vector< double > E( this->S->EqualityConstraints.size(), this->tolerance) ;
	    this->EqTolerances.assign(E.begin(), E.end());
	    std::vector< double > I( this->S->InequalityConstraints.size(), this->tolerance);
	    this->InEqTolerances.assign(I.begin(), I.end() );
	    //printf("Called right func\n");
	    if ( !S->EqualityConstraints.empty() ) {
		    this->opt.add_equality_mconstraint(&NloptAdapt< typename Problem::FUNCTION >::EqConstrIface, (void*)(&(this->NA)), this->EqTolerances);
	    }
	    if ( !S->InequalityConstraints.empty() ) {
		    this->opt.add_inequality_mconstraint(&NloptAdapt< typename Problem::FUNCTION >::InEqConstrIface, (void*)(&(this->NA)), this->InEqTolerances);
	    }
	    //printf("Called right func\n");
    }
};

#endif