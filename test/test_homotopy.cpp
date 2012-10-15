//Basic marching homotopy test to validate additional
//problem constraints

#include <gtest/gtest.h>

#include <vector>
#include "boost/function.hpp"

#include "Problems.h"
#include "Scalarization.hpp"
#include "Constraint.hpp"
#include <nlopt.hpp>
#include "NloptAdapt.hpp"
#include "optimizer.h"
#include "mesh.hpp"

#include <stdio.h>

#include "homotopy.hpp"

namespace {

	TEST(HOMOTOPY, 2D){
		//Set up the problem
		Problem::Interface * P = Problem::Factory("FON", 2, 3);
		//Problem::Interface * P = Problem::Factory("WFG5", 3, 5);
		//Problem::Interface * P = Problem::Factory("DTLZ2", 3, 10);

		//Instantiate the homotopy
		Pareto::Homotopy h( P, 1e-4 );

		//Homotopically deform the ansatz
		h.GetFront(10, 10);
	}
}