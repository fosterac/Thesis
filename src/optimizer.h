class Optimizer {
public:
	virtual double RunFrom(std::vector< double > &) = 0;
	virtual void SetWeights(std::vector< double > &) = 0;
	//AdjustConstraints
};

class OptNlopt : public Optimizer {
private:
	Problem::Interface * P;
	nlopt::opt opt;
	Scalarization<Problem::FUNCTION> s;
	NloptAdapt< Scalarization<Problem::FUNCTION> > NA;

public:
	OptNlopt(Problem::Interface *p, double tolerance);
	double RunFrom(std::vector< double > &);
	void SetWeights(std::vector< double > &);
};