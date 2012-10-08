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

//Trivial quadratic problem
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

//Trivial constrainted problem (from nlopt tutorial)
class ConstrainedSR : public Problem::Interface{
private:
	int a1, a2, b1, b2;

	double obj(const std::vector< double > &x){
		return sqrt(x[1]);
	}
	double constr1(const std::vector< double > &x){
		return -1 * x[1];
	}
	double constr2(const std::vector< double > &x){
		return pow(a1*x[0] + b1, 3) - x[1];
	}
	double constr3(const std::vector< double > &x){
		return pow(a2*x[0] + b2, 3) - x[1];
	}
public:
	ConstrainedSR(int DimObj, int DimDesign){

		a1 = 2;
		a2 = -1;

		b1 = 0;
		b2 = 1;

		this->dimObj = DimObj;
		this->dimDesign = DimDesign;

		typename Problem::FUNCTION f( boost::bind(&ConstrainedSR::obj, this, _1) );
		this->Objectives.push_back(	f );

		typename Problem::FUNCTION c1( boost::bind(&ConstrainedSR::constr1, this, _1) );
		this->InequalityConstraints.push_back(	c1 );

		typename Problem::FUNCTION c2( boost::bind(&ConstrainedSR::constr2, this, _1) );
		this->InequalityConstraints.push_back(	c2 );

		typename Problem::FUNCTION c3( boost::bind(&ConstrainedSR::constr3, this, _1) );
		this->InequalityConstraints.push_back(	c3 );
	}
};

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
			this->lowerBounds.push_back( -1.0 );
			this->upperBounds.push_back( 1.0 );
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
class WFG5 : public Problem::Interface{
private:
	double obj(const std::vector< double > &x, const int i, const int k, const int M){
		//We're seeing an assertion in the WFGProblems code fail 
		//as the input vectors are "out of bounds" (<0), which is
		//a result of the finite differences scheme.
		
		//return (Problems::WFG2(x, k, M))[i];

		//This should be done with a decorator anyway
		std::vector< double > var_x(x.begin(), x.begin() + this->dimDesign);
		return (Problems::WFG5(var_x, k, M))[i];
	}
public:
	WFG5(int DimObj, int DimDesign)	{
		this->dimObj = DimObj;
		this->dimDesign = DimDesign;

		int k = 1*(DimObj - 1);

		assert ( k == 1);

		int i;
		for(i=0; i<dimObj; i++){
			typename Problem::FUNCTION f( boost::bind(&WFG5::obj, this, _1, i, k, DimObj) );
			this->Objectives.push_back(	f );
		}
		for(i=0; i<dimDesign; i++){
			//The zero bound is enforced in the WFG code
			//so we perturb it to allow the central difference
			//scheme to work.
			//this->lowerBounds.push_back( 0.0 );
			this->lowerBounds.push_back( 1e-4 );
			this->upperBounds.push_back( 1.0 );
		}
	}
};

Problem::Interface * Problem::Factory( std::string s, int DimObj, int DimDesign){
	
	Interface * toReturn = NULL;
	
	//Instantiate a simple Basin problem
	if ( s.compare(std::string("BASIN")) == 0 ){
		//This only allows one objective
		assert (DimObj == 1);
		toReturn = new Basin(DimObj, DimDesign);
	}

	//Instantiate a problem to test constraints
	if ( s.compare(std::string("CONST_TEST")) == 0 ){
		//This only allows one objective
		assert (DimObj == 1);
		assert (DimDesign == 2);
		toReturn = new ConstrainedSR(DimObj, DimDesign);
	}

	//Instantiate a FON problem
	if ( s.compare(std::string("FON")) == 0 ){
		assert (DimDesign == 3);
		toReturn = new FON(DimObj, DimDesign);
	}
	
	//Instantiate a WFG2 Problem
	if ( s.compare(std::string("WFG5")) == 0 ){
		//assert (DimDesign % 2 != 0 );
		toReturn = new WFG5(DimObj, DimDesign);
	}

	return toReturn;
}
