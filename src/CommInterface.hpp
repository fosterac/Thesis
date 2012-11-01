#ifndef CommInterface_hpp
#define CommInterface_hpp

#include "HomotopyTypes.h"

#include <queue>

namespace Homotopy {
    namespace Communication {

        class MockCommunicator {
            int i;
            int size;
        public:
            MockCommunicator(int size) : i(0), size( size ) {}
            int PollLoop() { return (i++)%size; }
        };

        class Interface{
        public:
            virtual void dispatch( const designVars_t &, int ) =0;
            virtual void collect(std::queue< std::pair< int, objVars_t> >& ) =0;
        };

        template< typename T>
        class SimulatedRemote : public Interface {
            struct req {
                designVars_t x;
                int id;
            };
            const int EVALS_AT_ONCE;
            T& P;
            std::vector< req > requests;
        public:
            typedef T BaseType;
            SimulatedRemote( T & p ) : P(p), EVALS_AT_ONCE( 100 ) {
                srand ( unsigned ( time (NULL) ) );
            }
            void dispatch( const designVars_t &x, int id ) {
                req r = { x, id };
                requests.push_back( r );
            }
            void collect(std::queue< std::pair< int, objVars_t> >& results ){
                int evals = ( (requests.size() >= EVALS_AT_ONCE) ? (EVALS_AT_ONCE) : requests.size() );

                std::vector< req > v( requests.begin(), requests.begin() + evals);
                requests.erase( requests.begin(), requests.begin() + evals );

                //std::random_shuffle( v.begin(), v.end() );

                typename std::vector< req >::iterator it;
                for(it=v.begin(); it!=v.end(); it++){
                
                    objVars_t r( this->P.size() );
                    int i;
                    for(i=0; i<this->P.size(); i++){
                        r[i] = (this->P[i])( it->x );
                    }

                    std::pair< int, objVars_t > p( it->id, r );
                    results.push( p );
                }
            }
        };

        //To get this to compile on MPI-less machines
        //typedef MPI_Status Status_t;
        typedef int Status_t;

        template< typename T >
        class AdHoc : public Interface {
            //CommImpl::MPI comm_;
            T comm_;

            bool HandleMessage( std::queue< std::pair< int, objVars_t> >& q ) {
                return (this->Handler )( comm_.GetStatus(), comm_.GetValue(), q );
            }

            void initialize() {
                comm_.Init();
                comm_.AsyncRecv();
            }
        public:
            typedef int BaseType;

            AdHoc() { 
                this->initialize();
            }

            typedef boost::function<void (size_t, const designVars_t &)> dispatcher_t;
            typedef boost::function<bool (Status_t status, size_t length, std::queue< std::pair< int, objVars_t> >& q)> handler_t;

            dispatcher_t Dispatcher;
            handler_t Handler;

            void dispatch( const designVars_t &x , int id ) { 
                (this->Dispatcher) (id, x);
            }
            void collect(std::queue< std::pair< int, objVars_t> >& q ) {

                //Wait for a message
                while( ! comm_.HasMessage() ) { }

                //Receive all incomming messages
                while( comm_.HasMessage() ) {
                    if( comm_.ShouldStop() ) {
                        comm_.Shutdown();
                        //TODO: Figure out shutdown protocol
                        return;
                    } 
                    else {
                        if( this->HandleMessage( q ) )
                            comm_.AsyncRecv();
                        else
                            break;
                    }
                }
            }
        };

        
        //Communication layer implementations
        //that should allow writing simulated interfaces for testing
        namespace CommImpl {

            //Interface to implement
            class Iface {
            public:
                virtual void Init() {} 
                virtual void AsyncRecv() {} 
                virtual bool HasMessage() { return false;} 
                virtual bool ShouldStop() { return false; } 
                virtual void Shutdown() {} 
                virtual Status_t GetStatus(){ return 0; }
                virtual size_t GetValue(){ return 0; }
            };

            /*
            class MPI : public Iface {
            private:
                MPI_Status status;
                MPI_Request req;
                size_t recv_value;
                int flag;
                bool is_running_;

            public:
                MPI() : is_running( false ), recv_value( 0 ), flag( 0 ) {}
                void Init() {
                    is_running_ = true;
                }
                void AsyncRecv() {
                    MPI_Irecv(&recv_value, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, comm_m, &req);
                }
                bool HasMessage() {
                    if (req != MPI_REQUEST_NULL){
                        MPI_Test(&req, &flag, &status);
                        if(flag) return true;
                    }
                    return false;

                }
                bool ShouldStop() {
                    return status.MPI_TAG == MPI_STOP_TAG;
                }
                void Shutdown() {

                    is_running = false;
                }

                MPI_Status GetStatus(){ return req; }
                size_t GetValue(){ return recv_value; }
            };
            */
        }

        /*
        template< typename T >
        class Interface {
            T CommImpl_;

            void pushEvaluation(designVars_t &x){
            }
            void pushNeighbor(designVars_t &x, objVars_t &f, lamVars_t &l){
            }

            void runPollLoop() {

                setupPoll();

                //Init
                CommImpl_.Init()
                //Get message
                CommImpl_.AsyncRecv();
                while(true) {

                    prePoll();

                    if( CommImpl_.HasMessage() ) {
                        if( CommImpl_.ShouldStop() ) {
                            CommImpl_.Shutdown();
                            return;
                        } 
                        else {
                            if( CommImpl_.HandleMessage() )
                                //Get message
                                CommImpl_.AsyncRecv();
                            else
                                break;
                        }
                    }

                    postPoll();
                }
            }
        };
        */
    }
}
#endif
