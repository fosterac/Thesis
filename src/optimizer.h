class Optimizer {
public:
	virtual double RunFrom(std::vector< double > &) = 0;
	//AdjustConstraints
};

class OptNlopt : public Optimizer {
private:
	Scalarization *S;
	NloptAdapt< Scalarization > NA;
	nlopt::opt opt;

public:
	OptNlopt(Scalarization *s, double tolerance);
	double RunFrom(std::vector< double > &);
};