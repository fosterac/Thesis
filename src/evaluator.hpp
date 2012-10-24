#ifndef evaluator_hpp
#define evaluator_hpp

#include <boost/bind.hpp>
#include <boost/thread.hpp>

#include "HomotopyTypes.h"

#include <map>

#include <stdio.h>

namespace Homotopy {

    //Evaluator interface
    template< typename T >
    class Evaluator {
        T impl_;
    public:
        Evaluator( Problem::Interface * p ) : impl_(p) {}
        objVars_t eval( designVars_t &x ) {
            return impl_.eval( x );
        }
    };

    namespace EvaluationStrategy {
        //Simple local, sequential evaluation
        class Local {
        protected:
            Problem::Interface * P;
        public:
            Local ( Problem::Interface * p ) : P(p) {}
            objVars_t eval( designVars_t &x ) {
                objVars_t results( P->dimObj );
                int i;
                for(i=0; i<P->dimObj; i++){
                    results[i] = (P->Objectives[i])( x );
                }
                return results;
            }
        };

        //Map-based result caching (unbound) 
        template< typename T >
        class Cached {
            T impl_;
            //Establish a map-based cache
            typedef std::map< designVars_t, objVars_t > CacheType;
            CacheType cache_;
        public:
            Cached ( Problem::Interface * p ) : impl_(p) {}
            objVars_t eval( designVars_t &x ) {
                //search the cache for the designVars
                CacheType::iterator i = this->cache_.find( x );
                //spit out the answer if we already have it
                if( i != this->cache_.end() ) {
                    return i->second;
                }
                else {
                    //evaluate the designVars
                    objVars_t result ( impl_.eval( x ) ) ;
                    //cache the result
                    this->cache_[ x ] = result;
                    return result;
                }
            }
        };

        /*
        Decided not to support this until I can think of a better way to
        distribute thread-capable code from source.  
        */
        /*
        //Threaded, local evaluation
        class AsyncLocal {
        protected:
            Problem::Interface * P;
            static void exeAndStore( Problem::FUNCTION &f, designVars_t &x, double* result ) {
                (*result) = (f)( x );
            }
        public:
            AsyncLocal( Problem::Interface * p ) : P(p) {}
            objVars_t eval( designVars_t &x ) {
                //Thread container
                std::vector< boost::thread* > threads;
                //Results store
                objVars_t results( P->dimObj );

                int i;
                for(i=0; i<P->dimObj; i++){
                    //Launch threads
                    threads.push_back( new boost::thread( boost::bind( &AsyncLocal::exeAndStore,  P->Objectives[i], x, &results[i] ) ) );
                }
                for(i=0; i<P->dimObj; i++){
                    //Join thread group
                    threads[i]->join();
                }
                return results;
            }
        };
        */
    }
}
#endif