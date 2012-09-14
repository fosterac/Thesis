/*
Scalarizing adapter allowing dynamic weight changes
*/
#include <boost/bind.hpp>
#include <cassert>
#include <stdio.h>

//The public interface is quite similar to the Problem::Interface
//and could probably be refactored.
template< typename T >
class Scalarization {
protected:
	virtual double eval(const std::vector<double> &x) =0;
	Problem::Interface * P;
public:
	T f;
	std::vector< T > EqualityConstraints;
	std::vector< T > InequalityConstraints;

	int dimDesign;
	int dimObj;

	std::vector< double > lowerBounds;
	std::vector< double > upperBounds;

	Scalarization(Problem::Interface * p) : P(p), f(boost::bind( &Scalarization<T>::eval, this, _1) ),
											lowerBounds(P->lowerBounds), upperBounds(P->upperBounds),
											EqualityConstraints(P->EqualityConstraints),
											InequalityConstraints(P->InequalityConstraints){
		this->dimObj = 1;
		this->dimDesign = P->dimDesign;
	}

	//operator () overload
	double operator() (const std::vector<double> &x){
		return (this->eval(x));
	}
};

template< typename T >
class FixedScalarization : public Scalarization< T > {
protected:
	std::vector<double> * weights;

	//Could probably get rid of this when 
	//the interface is stabilized
	std::vector<double> default_weights;

	double eval(const std::vector<double> &at){	
		double result = 0.0;
		int i;
		for(i=0;i<this->P->Objectives.size();i++){
			result += (*(this->weights))[i] * (this->P->Objectives[i])(at);
		}
		return result;
	}

public:
	FixedScalarization(Problem::Interface * p) : Scalarization<T>(p), weights(NULL) {
		int i;
		for(i=0; i<this->P->Objectives.size(); i++) { this->default_weights.push_back(1.0); }
		this->weights = &(this->default_weights);
	}
	
	//Careful: No validation of pointer
	void SetWeights(std::vector<double> * w) { 
		assert ( w->size() == this->P->Objectives.size() );
		this->weights = w; 
	}

};

template< typename T >
class DynamicScalarization : public Scalarization < T > {
protected:
	double eval(const std::vector<double> &at){	
		//Allow for the weights to be part of the optimization
		//by assuming theyre tacked on to the end of the at vector

		//NOTE: the last weight is the one implied by the 1-sum(others)

		double result = 0.0;
		double weight_sum = 0.0;

		//Determine where the weights start
		int offset = at.size() - this->P->Objectives.size() + 1;

		//Separate the set of "Problem" variables
		const std::vector< double > vars(at.begin(), at.begin() + offset);

		//Evaluate the explicitly weighted components
		int i;
		for(i=0;i<(this->P->Objectives.size() - 1);i++){
			result += at[i + offset] * (this->P->Objectives[i])(vars);
			weight_sum += at[i + offset];
		}
		result += (1.0 - weight_sum) * (this->P->Objectives[i])(vars);
		return result;
	}
public:
	DynamicScalarization(Problem::Interface * p) : Scalarization<T>(p) {
		this->dimDesign = this->P->dimDesign + this->P->Objectives.size() - 1;
		int i;
		for(i=0; i<this->P->Objectives.size() - 1; i++) {
			this->lowerBounds.push_back(0.0);
			this->upperBounds.push_back(1.0);
		}
	}
};
