/*
Adapter wrapping a function object in the Nlopt library interface
*/

template<typename T>
class NloptAdapt{
private:
	//Templated function object:  could be boost::function, lambda, functor, etc.
	T &obj;
	double eval(const std::vector<double> &at){	return (this->obj)(at);	}

	//Elements to do manually perform finite differences
	//Should probably offload to another class...
	double fd_step;
	//Central difference quotient
	double finiteDifference(int dim, std::vector<double> &at){
		at[dim] -= this->fd_step; 
		double left = this->eval(at);
		at[dim] += 2*this->fd_step; 
		double right = this->eval(at);
		at[dim] -= this->fd_step; 

		return (right - left) / (2*this->fd_step);
	}
	void gradEval(std::vector<double> &grad, const std::vector<double> &at){
		std::vector<double> local_at(at);
		int i;
		for(i=0; i<grad.size(); i++){	grad[i] = this->finiteDifference(i, local_at); }
	}

public:
	NloptAdapt(T &Obj, double FD_Step) : obj(Obj), fd_step(FD_Step) {}

	//C-style callback interface expected by Nlopt
	static double iface(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data){
		NloptAdapt *nl = NULL;
		if(my_func_data) nl = reinterpret_cast<NloptAdapt *>(my_func_data);
		if(!grad.empty()){
			nl->gradEval(grad, x);
		}
		return nl->eval(x);
	}
};