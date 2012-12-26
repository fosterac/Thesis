//Basic marching homotopy test to validate additional
//problem constraints

#include "gtest/gtest.h"
#include "homotopy.hpp"

#include <algorithm>

namespace {

    //Problem::Interface * P = Problem::Factory("FON", 2, 3);
	//Problem::Interface * P = Problem::Factory("WFG5", 3, 30);
    //Problem::Interface * P = Problem::Factory("SURROGATE", 3, 2);

	TEST(HOMOTOPY, DISABLED_MOTTA1){
        Problem::Interface * P = Problem::Factory("MOTTA1", 3, 3);

        //Set up the communication framework
        Homotopy::Communication::SimulatedRemote< functionSet_t > Comm( P->Objectives );

		//Instantiate the homotopy
		//Pareto::homotopy h( P, 1e-3, 1e-6, Comm, Design, Objective, Lambda );
        Pareto::homotopy h( P, 1e-3, 1e-6, Comm );

		//Homotopically deform the ansatz
		h.GetFront(15, 5, 0, 1);
	}

    TEST(HOMOTOPY, DTLZ2){
        Problem::Interface * P = Problem::Factory("DTLZ2", 3, 3);

        std::vector< std::vector< double > > Design;
        std::vector< std::vector< double > > Objective;
        std::vector< std::vector< double > > Lambda;

        std::vector< double > d(3, 0.0);
        d[0] = 1.0; d[1] = 1.0 ; d[2] = 0.5;
        Design.push_back( d );
        d[0] = 0.0; d[1] = 1.0 ; d[2] = 0.5;
        Design.push_back( d );
        d[0] = 0.0; d[1] = 0.0 ; d[2] = 0.5;
        Design.push_back( d );

        d[0] = 0.0; d[1] = 0.0 ; d[2] = 1.0;
        Objective.push_back( d );
        d[0] = 0.0; d[1] = 1.0 ; d[2] = 0.0;
        Objective.push_back( d );
        d[0] = 1.0; d[1] = 0.0 ; d[2] = 0.0;
        Objective.push_back( d );

        d[0] = 1.0; d[1] = 1.0 ; d[2] = 0.0;
        Lambda.push_back( d );
        d[0] = 1.0; d[1] = 0.0 ; d[2] = 1.0;
        Lambda.push_back( d );
        d[0] = 0.0; d[1] = 1.0 ; d[2] = 1.0;
        Lambda.push_back( d );

        //Set up the communication framework
        Homotopy::Communication::SimulatedRemote< functionSet_t > Comm( P->Objectives );

		//Instantiate the homotopy
		Pareto::homotopy h( P, 1e-3, 1e-6, Comm, Design, Objective, Lambda );
        //Pareto::homotopy h( P, 1e-3, 1e-6, Comm );

		//Homotopically deform the ansatz
		h.GetFront(15, 5, 0, 1);
	}
}