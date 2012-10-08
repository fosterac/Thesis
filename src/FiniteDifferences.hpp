namespace FiniteDifferences {
	struct CENTRAL {};
	template< typename T >
	static double Central(const int dim, std::vector<double> &at, T& func, const double fd_step){
		double tmp = at[dim];

		at[dim] -= fd_step; 
		double left = (func)(at);

		at[dim] = tmp + fd_step; 
		double right = (func)(at);

		at[dim] = tmp;

		return (right - left) / (2*fd_step);
	}
	struct FORWARD {};
	template< typename T >
	static double Forward(const int dim, std::vector<double> &at, T& func, const double fd_step){
		double tmp = at[dim];

		double center = (func)(at);

		at[dim] += fd_step; 
		double forward = (func)(at);

		at[dim] = tmp;

		return (forward - center) / (fd_step);
	}
	template< typename T >
	static double ForwardOpt(const int dim, std::vector<double> &at, double center, T& func, const double fd_step){
		double tmp = at[dim];

		at[dim] += fd_step; 
		double forward = (func)(at);

		at[dim] = tmp;

		return (forward - center) / (fd_step);
	}
	
	template<typename Grad, typename T>
	void GradEval(Grad grad, T& func, const std::vector< double > & at, double fd_step, CENTRAL c){
		std::vector<double> local_at(at);
		int i;
		for(i=0; i<at.size(); i++){	grad[i] = Central<T>(i, local_at, func, fd_step); }
	}
	template<typename Grad, typename T>
	void GradEval(Grad grad, T& func, const std::vector< double > & at, double fd_step, FORWARD f){
		std::vector<double> local_at(at);
		double f_at = (func) (local_at);
		int i;
		for(i=0; i<at.size(); i++){	grad[i] = ForwardOpt<T>(i, local_at, f_at, func, fd_step); }
	}
}