namespace FiniteDifferences {

	template< typename T >
	static double Central(const int dim, std::vector<double> &at, T& func, const double fd_step){
		at[dim] -= fd_step; 
		double left = (func)(at);
		at[dim] += 2*fd_step; 
		double right = (func)(at);
		at[dim] -= fd_step; 

		return (right - left) / (2*fd_step);
	}

	template< typename T >
	static double Forward(const int dim, std::vector<double> &at, T& func, const double fd_step){
		double center = (func)(at);

		at[dim] += fd_step; 
		double forward = (func)(at);

		at[dim] -= fd_step; 

		return (forward - center) / (fd_step);
	}

}