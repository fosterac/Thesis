#include <vector>

#include <boost/bind.hpp>
#include <boost/function.hpp>


#include "Problems.h"

#include <numeric>
namespace helpers {
	double square( double l, double r )	{ return l + (r) * (r); }
}

//using namespace boost::bind;

//class Basin : public Problem {

//private:
	double Basin::obj(const std::vector< double > &x){
		return std::accumulate(x.begin(), x.end(), 0.0, helpers::square);
	}

//public:
	Basin::Basin(int DimObj, int DimDesign)	{

		//Problem::FUNCTION f = obj;
		boost::function<double (const std::vector<double> &)> FUNC( boost::bind(&Basin::obj, this, _1) );
		//Problem::FUNCTION f = *FUNC.target<Problem::FUNCTION>();
		this->Objectives.push_back(	FUNC	);
		//this->Objectives.push_back(	f	);

		this->dimObj = DimObj;
		this->dimDesign = DimDesign;

		int i;
		for(i=0; i<dimDesign; i++){
			this->lowerBounds.push_back( 0.0 );
			this->upperBounds.push_back( 1.0 );
		}
	}

//};
/*
Basin::Basin(int DimObj){
	this->impl = new BasinImpl(DimObj);
}
Basin::~Basin(){
	delete this->impl;
}*/