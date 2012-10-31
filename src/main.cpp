/*

Test program showing the general state of functionality

*/

#include <stdio.h>
#include <vector>

#include "boost/function.hpp"

#include "Problems.h"
#include "Scalarization.hpp"
#include "nlopt.hpp"
#include "NloptAdapt.hpp"

#include "optimizer.h"

int main(int argc, char** argv){

	//Choose a Problem
	int DesignVars = 3;
	Problem::Interface * P = Problem::Factory("FON", 2, DesignVars);
	//Problem::Interface * P = Problem::Factory("WFG2", 2, DesignVars);
	
	//Scalarize the problem
	FixedScalarization< typename Problem::FUNCTION > * S = new FixedScalarization< typename Problem::FUNCTION >(P);

	//Instantiate an optimizer
	Optimizer * op = new OptNlopt(S->f, S, 1e-4);

	//Set up an initial point
	std::vector<double> x(S->dimDesign);
	int i;
	for(i=0; i<DesignVars; i++) { x[i] = 0.5; }

	//Set up the weight vector
	std::vector<double> w(P->Objectives.size());

	//Sample the curve
	int NumPts = 10;
	int j;
	for(j=0; j<NumPts; j++){

		//Set the weights
		w[0] =  j* 1.0/(double)(NumPts-1);
		w[1] = 1.0 - w[0];
		S->SetWeights(&w);
		
		//Re-initialize the starting points
		//for(i=0; i<DesignVars; i++) { x[i] = 0.5; }

		//Run the Optimization
		double result = op->RunFrom(x);

		//Print the results
		//printf("Point %d: \t %lf, %lf -> (%lf, %lf) = %lf\n", j, w[0], w[1], (P->Objectives[0])(x), (P->Objectives[1])(x), result);
		printf("%lf, %lf\n", (P->Objectives[0])(x), (P->Objectives[1])(x));
	}

	return 0;
}