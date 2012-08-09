class Problem {
public:

	//typedef double (*FUNCTION)(const std::vector<double> &);
	typedef boost::function<double (const std::vector<double> &)> FUNCTION;

	std::vector< FUNCTION > Objectives;
	std::vector< FUNCTION > Constraints;

	int dimDesign;
	int dimObj;

	std::vector< double > lowerBounds;
	std::vector< double > upperBounds;

	//Func eval odometer?
};

/*
class Factory{
public:
	static 
};*/

class Basin : public Problem{
private:
	double obj(const std::vector< double > &);
public:
	Basin(int DimObj, int DimDesign);
};


/*
class Basin : public Problem {
public:
	Basin(int DimObj);
	~Basin();
private:
	BasinImpl * impl;
}*/