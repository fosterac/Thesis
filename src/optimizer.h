class Optimizer {
public:
	virtual double RunFrom(std::vector< double > &) = 0;
};

class OptNlopt : public Optimizer {
private:
	Scalarization< typename Problem::FUNCTION > *S;
	NloptAdapt< typename Problem::FUNCTION > NA;
	nlopt::opt opt;

	std::vector< double > EqTolerances;
	std::vector< double > InEqTolerances;

public:
	OptNlopt(Problem::FUNCTION &Obj, Scalarization< typename Problem::FUNCTION > *s, double tolerance);
	double RunFrom( std::vector< double > & );
};