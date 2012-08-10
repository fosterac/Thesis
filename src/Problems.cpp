#include <vector>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <string>

#include "Problems.h"

#include <numeric>
namespace helpers {
	double square( double l, double r )	{ return l + (r) * (r); }
}

//Simple quadratic problem
class Basin : public Problem::Interface{
private:
	double obj(const std::vector< double > &x){
		return std::accumulate(x.begin(), x.end(), 0.0, helpers::square);
	}
public:
	Basin(int DimObj, int DimDesign)	{

		//Problem::FUNCTION f = obj;
		boost::function<double (const std::vector<double> &)> f( boost::bind(&Basin::obj, this, _1) );
		this->Objectives.push_back(	f );

		this->dimObj = DimObj;
		this->dimDesign = DimDesign;

		int i;
		for(i=0; i<dimDesign; i++){
			this->lowerBounds.push_back( 0.0 );
			this->upperBounds.push_back( 1.0 );
		}
	}
};

Problem::Interface * Problem::Factory( std::string s, int DimObj, int DimDesign){
	Interface * toReturn;
	//Instantiate a simple Basin problem
	if ( s.compare(std::string("BASIN")) == 0 ){
		toReturn = new Basin(DimObj, DimDesign);
	}
	else toReturn = NULL;

	return toReturn;
}
