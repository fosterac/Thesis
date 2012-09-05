namespace Problem{

	//Use variant here...
	typedef boost::function<double (const std::vector<double> &)> FUNCTION;

	struct Interface {

		std::vector< FUNCTION > Objectives;
		std::vector< FUNCTION > Constraints;

		int dimDesign;
		int dimObj;

		std::vector< double > lowerBounds;
		std::vector< double > upperBounds;

		//Func eval odometer?
	};

	Interface * Factory( std::string, int, int );
}