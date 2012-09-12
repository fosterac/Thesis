class Optimizer {
public:
	virtual double RunFrom(std::vector< double > &) = 0;
	//AdjustConstraints
};

class OptNlopt : public Optimizer {
private:
	Scalarization *S;
	NloptAdapt< typename Problem::FUNCTION > NA;
	nlopt::opt opt;

	std::vector< double > EqTolerances;
	std::vector< double > InEqTolerances;

public:
	OptNlopt(Scalarization *s, double tolerance);
	double RunFrom(std::vector< double > &);
};