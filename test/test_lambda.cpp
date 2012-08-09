#include <gtest/gtest.h>

#include <boost/lambda/lambda.hpp>
#include <boost/function.hpp>

#include <iostream>
#include <vector>
#include <algorithm>

#include "Problems.h"

using namespace std;
using namespace boost::lambda;

class LambdaTest : public ::testing::Test {
protected:
	vector<double> v;

	virtual void SetUp(){
		v.push_back( 1.0 );
		v.push_back( 2.0 );
		v.push_back( 3.0 );
	}
};


// Tests boost lambda functionality
TEST_F(LambdaTest, Alive) {
	for_each(v.begin(), v.end(), cout << _1 << '\n');
	std::for_each(v.begin(), v.end(), _1 = _1 * _1);
	for_each(v.begin(), v.end(), cout << _1 << '\n');
}

//Test that the Problems wrapping technique works
TEST_F(LambdaTest, WrapProblem) {
	for_each(v.begin(), v.end(), cout << _1 << '\n');

	Problem * P = new Basin(1, 3);
	cout << (P->Objectives[0]) (v) << endl;

}
