#include <vector>

#include <boost/bind.hpp>
#include <boost/function.hpp>

#include <string>
#include <cassert>
#include <stdio.h>

#include <iostream>

#include "Problems.h"
#include <numeric>
namespace helpers {
	double addsquare( double l, double r )	{ return l + (r) * (r); }
	double square( double x )	{ return (x) * (x); }
}

//Simple quadratic problem
class Basin : public Problem::Interface{
private:
	double obj(const std::vector< double > &x){
		return std::accumulate(x.begin(), x.end(), 0.0, helpers::addsquare);
	}
public:
	Basin(int DimObj, int DimDesign)	{

		typename Problem::FUNCTION f( boost::bind(&Basin::obj, this, _1) );
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

#include <math.h>

//Problem FON from Pereyra 2009
class FON : public Problem::Interface{
private:
	double obj1(const std::vector< double > &x){
		return 1.0 - exp( -1.0 * ( helpers::square(x[0] - 1.0/sqrt(3.0)) + helpers::square(x[1] - 1.0/sqrt(3.0)) + helpers::square(x[2] - 1.0/sqrt(3.0)) ));
	}
	double obj2(const std::vector< double > &x){
		return 1.0 - exp( -1.0 * ( helpers::square(x[0] + 1.0/sqrt(3.0)) + helpers::square(x[1] + 1.0/sqrt(3.0)) + helpers::square(x[2] + 1.0/sqrt(3.0)) ));
	}
public:
	FON(int DimObj, int DimDesign)	{

		typename Problem::FUNCTION f( boost::bind(&FON::obj1, this, _1) );
		this->Objectives.push_back(	f );
		typename Problem::FUNCTION g( boost::bind(&FON::obj2, this, _1) );
		this->Objectives.push_back(	g );

		this->dimObj = DimObj;
		this->dimDesign = DimDesign;

		int i;
		for(i=0; i<dimDesign; i++){
			this->lowerBounds.push_back( -5.0 );
			this->upperBounds.push_back( 5.0 );
		}
	}
};



//// Toolkit includes. //////////////////////////////////////////////////////

#include "wfgProblems/Toolkit/ExampleProblems.h"
#include "wfgProblems/Toolkit/TransFunctions.h"


//// Used namespaces. ///////////////////////////////////////////////////////

using namespace WFG::Toolkit;
using namespace WFG::Toolkit::Examples;

//WFG2
class WFG2 : public Problem::Interface{
private:
	double obj(const std::vector< double > &x, const int i, const int k, const int M){
		//We're seeing an assertion in the WFGProblems code fail 
		//as the input vectors are "out of bounds" (<0), which is
		//a result of the finite differences scheme.

		//int j;
		//for(j=0;j<x.size();j++){ printf("%lf ", x[j]);}
		//printf("\n");
		return Problems::WFG2(x, k, M).at(i);
	}
public:
	WFG2(int DimObj, int DimDesign)	{
		this->dimObj = DimObj;
		this->dimDesign = DimDesign;

		int k = 1*(DimObj - 1);

		int i;
		for(i=0; i<dimObj; i++){
			typename Problem::FUNCTION f( boost::bind(&WFG2::obj, this, _1, i, k, DimObj) );
			this->Objectives.push_back(	f );
		}
		for(i=0; i<dimDesign; i++){
			this->lowerBounds.push_back( 0.0 );
			this->upperBounds.push_back( 1.0 );
		}
	}
};

Problem::Interface * Problem::Factory( std::string s, int DimObj, int DimDesign){
	Interface * toReturn;

	toReturn = NULL;
	
	//Instantiate a simple Basin problem
	if ( s.compare(std::string("BASIN")) == 0 ){
		//This only allows one objective
		assert (DimObj == 1);
		toReturn = new Basin(DimObj, DimDesign);
	}

	//Instantiate a FON problem
	if ( s.compare(std::string("FON")) == 0 ){
		//This only allows one objective
		assert (DimDesign == 3);
		toReturn = new FON(DimObj, DimDesign);
	}
	
	//Instantiate a WFG2 Problem
	if ( s.compare(std::string("WFG2")) == 0 ){
		assert (DimDesign % 2 != 0 );
		toReturn = new WFG2(DimObj, DimDesign);
	}

	return toReturn;
}
