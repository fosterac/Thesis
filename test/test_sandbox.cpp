#include "gtest/gtest.h"

#include "boost/lambda/lambda.hpp"
#include "boost/function.hpp"

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

	//Problem * P = new Basin(1, 3);
	Problem::Interface * P = Problem::Factory("BASIN", 1, 3);
	cout << (P->Objectives[0]) (v) << endl;

}

class VerticalIterTest : public ::testing::Test {
public:
	template< typename T >
	class vertical_iterator : public T::const_iterator {
		int column;
	public:
		vertical_iterator( int c ) : column( c ) {}
		typename T::value_type::value_type operator*() const { return (this->T::const_iterator::operator*())[column]; }
		vertical_iterator& operator= (const typename T::const_iterator & rhs) { this->T::const_iterator::operator=(rhs); }
	};
	typedef std::vector< std::vector< double > > Data;
};

#include <algorithm>

// Tests vertical iterator concept
TEST_F(VerticalIterTest, Alive) {
	const std::vector< std::vector< double > > data;
	std::vector< double > v(3);
	v[0] = 0.0; v[1] = 1.0; v[2] = 2.0;
	//data.push_back(v);
	v[0]++ ; v[1]++ ; v[2]++ ;
	//data.push_back(v);

	vertical_iterator< Data > vi(1);
	for(vi = data.begin(); vi != data.end(); vi++){
		printf("%lf\n", *vi);
	}

	vertical_iterator< Data > vi1(1);
	vertical_iterator< Data > vi2(1);
	vi1 = data.begin();
	vi2 = data.end();
	printf("%lf\n", *std::min_element(vi1, vi2));
}