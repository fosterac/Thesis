//Basic marching homotopy test to validate additional
//problem constraints

#include <gtest/gtest.h>

#include <vector>
#include <boost/function.hpp>

#include "Problems.h"
#include "Scalarization.hpp"
#include "Constraint.hpp"
#include <nlopt.hpp>
#include "NloptAdapt.hpp"
#include "optimizer.h"
#include "mesh.hpp"

#include "homotopy.hpp"

#include <stdio.h>

namespace {

	TEST(HOMOTOPY, 2D){
		//Set up the problem
		Problem::Interface * P = Problem::Factory("FON", 2, 3);
		//Problem::Interface * P = Problem::Factory("WFG5", 2, 5);

		//Instantiate the homotopy
		Pareto::Homotopy h( P, 1e-4 );

		//Homotopically deform the ansatz
		h.GetFront(10, 5);
	}
}