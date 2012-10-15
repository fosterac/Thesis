//Simple test to verify interpolator functionality

#include "gtest/gtest.h"

#include <vector>

#include "interpolator.hpp"

namespace {
	TEST(InterpolatorTest, Alive) {
		std::vector< std::vector< double > > data;
		std::vector< double > v(3);
		v[0] = -1.0; v[1] = 0.0; v[2] = 2.0;
		data.push_back( v );
		v[0] = 1.0; v[1] = 0.0; v[2] = 3.0;
		data.push_back( v );

		Interpolation::RBF rbf(data);
		std::vector< double > point(2, 0.0);
		double result = rbf.evaluate( point );
		printf("%lf\n", result);
		EXPECT_NEAR( 2.5, double(result), 1e-3 );
	}
}