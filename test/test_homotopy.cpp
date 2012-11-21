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
        //Problem::Interface * P = Problem::Factory("SURROGATE", 2, 2);

        //Set up the communication framework
        Homotopy::Communication::SimulatedRemote< functionSet_t > Comm( P->Objectives );
        
        /*Homotopy::Communication::AdHoc< Homotopy::Communication::CommImpl::Simulator< functionSet_t > > Comm;
        Comm.comm_.P = P->Objectives;
        Comm.Dispatcher = boost::bind( &Homotopy::Communication::CommImpl::Simulator< functionSet_t >::dispatcher,
            &Comm.comm_, _1, _2 );
        Comm.Handler = boost::bind( &Homotopy::Communication::CommImpl::Simulator< functionSet_t >::handler,
            &Comm.comm_, _1, _2, _3 );*/

		//Instantiate the homotopy
		Pareto::homotopy h( P, 1e-3, Comm );

		//Homotopically deform the ansatz
		h.GetFront(5, 1);
	}
}