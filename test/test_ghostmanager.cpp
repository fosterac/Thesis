//Validate the evaluator

#include "gtest/gtest.h"

#include "boost/bind.hpp"
#include "boost/function.hpp"

#include "Problems.h"
#include "HomotopyTypes.h"

#include "GhostManager.hpp"

#include <stdio.h>

using namespace Homotopy;

namespace {

    void Print(const std::vector< double > & x ){
        int i;
	    for(i=0; i<x.size(); i++) { 
		    printf("%lf ", x[i]);
	    }
	    printf("\n");
    }

    class TestExchanger{
        ind_t rec_ind;
        objVars_t obj;

    public:
        TestExchanger( ind_t ind ) : rec_ind( ind ) {
            this->obj = objVars_t (1, 0.0);
        }
        void exchanger(  std::queue< NeighborMessage_t > & toSend, std::queue< nodeEnvelope_t > & received ){
            //Mock send
            while( !toSend.empty() ){
                NeighborMessage_t& m = toSend.front();
                printf("To %d: \n", m.first );
                std::vector< nodeEnvelope_t >::iterator i;
                for( i=m.second.begin(); i!=m.second.end(); i++){
                    Print( i->second );
                }
                toSend.pop();
            }

            //Mock receive
            nodeEnvelope_t e ( this->rec_ind, this->obj );
            received.push( e );
        }
    };
    

	TEST(GHOSTMANAGER, Alive){
        TestExchanger te( 1 );

        std::vector< double > v(3, 1.0);
        Point_t point0( 0, v, v, v);
        Point_t point1( 1, v, v, v);

        std::vector< Point_t* > pointvec;
        pointvec.push_back( &point0 );

        std::map< ind_t, std::vector< Point_t* > > ns;
        ns[0] = pointvec;

        std::map< ind_t, Point_t* > gg;
        gg[1] = &point1;

        GhostManager< TestExchanger > gm( te, gg, ns );

        gm.Exchange();

        Print(point1.ObjectiveCoords);
	}
}