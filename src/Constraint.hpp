/*
General constraint interface
*/

#include <boost/bind.hpp>
#include <cassert>

template< typename T >
class GeneralConstraint {
protected:
	virtual double eval(const std::vector< double > &x) =0;
public:
	T function;
	GeneralConstraint() :
		function( boost::bind( &GeneralConstraint::eval, this, _1 ) ) {}
};

#include <numeric>

template< typename T >
class BoundSumConstraint : public GeneralConstraint< T > {
protected:
	double value;
	int sign;
	int from;
	int to;

	double eval(const std::vector< double > &x) {
		return this->sign * ( this->value - std::accumulate(x.begin() + this->from, x.begin() + this->to, 0.0) );
	}

public:
	enum SIGN { LESS_THAN = -1, GREATER_THAN = 1 };
	BoundSumConstraint( SIGN s, double value, int from, int to ) : 
	GeneralConstraint< T > (), sign(s), value(value), from(from), to(to) {}
};

#include <math.h>

template< typename T >
class StepConstraint : public GeneralConstraint< T >  {
protected:
	const std::vector< double > * x1;
	double step;
	static double L2Dist(const std::vector< double > &x1, const std::vector< double > &x2 ) {
		int i;
		double sum = 0.0;
		for(i=0;i<x1.size();i++){
			sum += pow( x1[i] - x2[i], 2.0 );
		}
		return sqrt(sum);
	}
	virtual double eval(const std::vector< double > &x){
		return L2Dist(*this->x1, x) - this->step;
	}
public:
	StepConstraint( const std::vector< double >  * x1, double step ) :
	  GeneralConstraint< T > (), x1(x1), step(step) { }
	void UpdateFrom( const std::vector< double > * x1) { this->x1 = x1; }
};
template< typename T >
class FStepConstraint : public StepConstraint< T >  {
protected:
	const std::vector< T > & F_trans;
	std::vector< double > f;
	void Feval(const std::vector< double > &x){
		int i;
		for(i=0; i<this->F_trans.size(); i++) { this->f[i] = ( (F_trans[i])(x) ); }
	}
	virtual double eval(const std::vector< double > &x){
		this->Feval(x);
		return L2Dist(*this->x1, this->f) - this->step ;
	}
public:
	FStepConstraint( const std::vector< T > &F_trans, const std::vector< double >  * f1, double step) :
	  F_trans(F_trans), f(F_trans.size(), 0.0), StepConstraint< T >( f1, step) {}
};

template< typename T >
class FEqDistanceConstraint : public FStepConstraint< T >  {
protected:
	std::vector< double > const * x2;

	double eval(const std::vector< double > &x){
		this->Feval(x);
		return L2Dist(*this->x1, this->f) - L2Dist(*this->x2, this->f) ;
	}
public:
	FEqDistanceConstraint( const std::vector< T > &F_trans, const std::vector< double >  * f1, const std::vector< double >  * f2) :
	  x2(f2), FStepConstraint< T >(F_trans, f1, 0.0) {}
	void UpdateFrom( std::vector< double > * const f1, std::vector< double > * const f2) { 
		this->x1 = f1; 
		this->x2 = f2; 
	}
};

