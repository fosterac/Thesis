/*
Adapter wrapping a function object in the Nlopt library interface
*/

#include "FiniteDifferences.hpp"

template<typename T>
class NloptAdapt{
private:
	//Templated function object:  could be boost::function, lambda, functor, etc.
	T & obj;
	std::vector< T > * EqConstr;
	std::vector< T > * InEqConstr;

	//Elements to manually perform finite differences
	//Should probably offload to another class...
	double fd_step;
	//Central difference quotient
	/*static double finiteDifference(const int dim, std::vector<double> &at, T& func, const double fd_step){
		at[dim] -= fd_step; 
		double left = eval(func, at);
		at[dim] += 2*fd_step; 
		double right = eval(func, at);
		at[dim] -= fd_step; 

		return (right - left) / (2*fd_step);
	}*/

	static double eval( T& func, const std::vector<double> &at){	return (func)(at);	}

	static void gradEval(std::vector<double> &grad, const std::vector<double> &at, T &func, const double fd_step){
		std::vector<double> local_at(at);
		int i;
		//for(i=0; i<grad.size(); i++){	grad[i] = finiteDifference(i, local_at, func, fd_step); }
		for(i=0; i<at.size(); i++){	grad[i] = FiniteDifferences::Forward<T>(i, local_at, func, -fd_step); }
	}
	static void gradEval(double* grad, const std::vector< double > &at, T & func, const double fd_step){
		std::vector<double> local_at(at);
		int i;
		//for(i=0; i<at.size(); i++){	grad[i] = finiteDifference(i, local_at, func, fd_step); }
		for(i=0; i<at.size(); i++){	grad[i] = FiniteDifferences::Forward<T>(i, local_at, func, -fd_step); }
	}

	static void ConstrIface(unsigned m, double* result, unsigned n, const double* x, double* grad, std::vector< T > &Constr, NloptAdapt *nl){
		//Copy array data to vector
		std::vector< double > at(n);
		int i;
		for(i=0; i<n; i++){ at[i] = x[i]; }
		
		for(i=0; i<m; i++){
			result[i] = eval(Constr[i], at);
			if(grad) gradEval(grad+i*n, at, Constr[i], nl->fd_step);
		}
	}

public:
	NloptAdapt(T &Obj, double FD_Step) : obj(Obj), EqConstr(NULL), InEqConstr(NULL), fd_step(FD_Step) {}
	NloptAdapt(T &Obj, std::vector<T> *EqConstr, std::vector<T> *InEqConstr, double FD_Step) : obj(Obj), EqConstr(EqConstr), InEqConstr(InEqConstr), fd_step(FD_Step) {}

	//C-style callback interfaces expected by Nlopt
	static double ObjIface(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data){
		assert ( my_func_data != NULL );
		NloptAdapt *nl = reinterpret_cast<NloptAdapt *>(my_func_data);
		if(!grad.empty()){
			gradEval(grad, x, nl->obj, nl->fd_step);
		}
		return eval(nl->obj, x);
	}

	//C-style interface for single constraints
	/*static double ConstrIface(const std::vector<double> &x, std::vector<double> &grad, void *data){

	}*/

	//C-style interface for multiple constriaints
	static void EqConstrIface(unsigned m, double* result, unsigned n, const double* x, double* grad, void* data){
		assert ( data != NULL );
		NloptAdapt *nl = reinterpret_cast<NloptAdapt *>(data);
		ConstrIface( m, result, n, x, grad, *nl->EqConstr, nl);
	}
	static void InEqConstrIface(unsigned m, double* result, unsigned n, const double* x, double* grad, void* data){
		assert ( data != NULL );
		NloptAdapt *nl = reinterpret_cast<NloptAdapt *>(data);
		ConstrIface( m, result, n, x, grad, *nl->InEqConstr, nl);
	}
};