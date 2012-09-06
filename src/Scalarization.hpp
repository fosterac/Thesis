/*
Scalarizing adapter allowing dynamic weight changes
*/

#include <cassert>

template<typename T>
class Scalarization{
private:
	//Templated function object:  could be boost::function, lambda, functor, etc.
	std::vector<T> &obj;
	std::vector<double> * weights;

	//Could probably get rid of this when 
	//the interface is stabilized
	std::vector<double> default_weights;

	double eval(const std::vector<double> &at){	
		double result = 0.0;
		int i;
		for(i=0;i<obj.size();i++){
			result += ((*weights)[i]) * (obj[i])(at);
		}
		return result;
	}
	double eval_dyn(const std::vector<double> &at){	
		//Allow for the weights to be part of the optimization
		//by assuming theyre tacked on to the end of the at vector

		//NOTE: the last weight is the one implied by the 1-sum(others)

		double result = 0.0;
		double weight_sum = 0.0;

		int offset = at.size() - obj.size() + 1;

		int i;
		for(i=0;i<(obj.size() - 1);i++){
			result += at[i + offset] * (obj[i])(at);
		}
		result += (1.0 - weight_sum) * (obj[i])(at);
		return result;
	}

public:
	Scalarization(std::vector<T> &Obj) : obj(Obj), weights(NULL) {
		int i;
		for(i=0; i<this->obj.size(); i++) { this->default_weights.push_back(1.0); }
		this->weights = &(this->default_weights);
	}

	//Careful: No validation of pointer
	void SetWeights(std::vector<double> * w) { 
		assert ( w->size() == obj.size() );
		this->weights = w; 
	}

	//operator () overload
	double operator() (const std::vector<double> &x){
		//return (this->eval(x));
		return (this->eval_dyn(x));
	}
};