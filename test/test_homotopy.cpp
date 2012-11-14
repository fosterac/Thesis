//Basic marching homotopy test to validate additional
//problem constraints

#include "gtest/gtest.h"

#include "homotopy.hpp"

namespace {

	TEST(HOMOTOPY, 2D){
		//Set up the problem
		Problem::Interface * P = Problem::Factory("FON", 2, 3);
		//Problem::Interface * P = Problem::Factory("WFG5", 3, 30);
		//Problem::Interface * P = Problem::Factory("DTLZ2", 3, 30);

		//Instantiate the homotopy
		Pareto::Homotopy h( P, 1e-4 );

		//Homotopically deform the ansatz
		h.GetFront(9, 10);
	}
}