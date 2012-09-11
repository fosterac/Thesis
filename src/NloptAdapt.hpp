/*
Adapter wrapping a function object in the Nlopt library interface
*/

template<typename T>
class NloptAdapt{
private:
	//Templated function object:  could be boost::function, lambda, functor, etc.
	T & obj;
	std::vector< T > * constr;

	//Elements to manually perform finite differences
	//Should probably offload to another class...
	double fd_step;
	//Central difference quotient
	static double finiteDifference(const int dim, std::vector<double> &at, const T& func, const double fd_step){
		at[dim] -= fd_step; 
		double left = eval(func, at);
		at[dim] += 2*fd_step; 
		double right = eval(func, at);
		at[dim] -= fd_step; 

		return (right - left) / (2*fd_step);
	}

	static double eval(const T& func, const std::vector<double> &at){	return (func)(at);	}

	static void gradEval(std::vector<double> &grad, const std::vector<double> &at, const T &func, const double fd_step){
		std::vector<double> local_at(at);
		int i;
		for(i=0; i<grad.size(); i++){	grad[i] = finiteDifference(i, local_at, func, fd_step); }
	}

public:
	NloptAdapt(T &Obj, double FD_Step) : obj(Obj), constr(NULL), fd_step(FD_Step) {}
	NloptAdapt(T &Obj, std::vector<T> &Constr, double FD_Step) : obj(Obj), constr(&Constr), fd_step(FD_Step) {}

	//C-style callback interfaces expected by Nlopt
	static double iface(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data){
		NloptAdapt *nl = NULL;
		if(my_func_data) nl = reinterpret_cast<NloptAdapt *>(my_func_data);
		if(!grad.empty()){
			gradEval(grad, x, nl->obj, nl->fd_step);
		}
		return eval(nl->obj, x);
	}
	/*static double ConstrIface(const std::vector<double> &x, std::vector<double> &grad, void *data){
		NloptAdapt *nl = NULL;
		if(my_func_data) nl = reinterpret_cast<NloptAdapt *>(my_func_data);
		if(!grad.empty()){
			nl->gradEval(grad, x);
		}
		return nl->eval(x);
	}*/
};