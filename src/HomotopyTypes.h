#ifndef HomotopyTypes_h
#define HomotopyTypes_h

//Standard includes
#include <vector>

//Boost includes
#include "boost/bind.hpp"
#include "boost/function.hpp"

//Switch for local/remote evaluation mechanisms
#ifdef HAS_MPI
#define REMOTE_EVAL
#else
#define LOCAL_EVAL
#endif

namespace Homotopy {

    //Point variables
    typedef std::vector< double > designVars_t;
    typedef std::vector< double > objVars_t;
    typedef std::vector< double > lamVars_t;

    //Function types
    typedef boost::function<double (const designVars_t &)> function_t;
    typedef std::vector< function_t > functionSet_t;
    
}

#endif