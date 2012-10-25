#ifndef HomotopyTypes_h
#define HomotopyTypes_h

#include <vector>

#include "Problems.h"

namespace Homotopy {

    typedef std::vector< double > designVars_t;
    typedef std::vector< double > objVars_t;
    typedef std::vector< double > lamVars_t;

    typedef Problem::FUNCTION function_t;
    typedef std::vector< function_t > functionSet_t;

}

#endif