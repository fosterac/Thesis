#ifndef CommInterface_hpp
#define CommInterface_hpp

#include "HomotopyTypes.h"

namespace Homotopy {
    namespace Communication {

        class SimpleMockCommunicator {
            int i;
            int size;
        public:
            MockCommunicator(int size) : i(0), size( size ) {}
            int PollLoop() { return (i++)%size; }
        };


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

        //Communication layer implementations
        //that should allow writing simulated interfaces for testing
        namespace CommImpl {

            //Interface to implement
            class Iface {
            public:
                virtual void Init() {} =0;
                virtual void AsyncRecv() {} =0;
                virtual bool HasMessage() {} =0;
                virtual bool HandleMessage() {} =0;
                virtual bool ShouldStop() {} =0;
                virtual void ShutDown() {} =0;
            };

            class MPI : public Iface {
            private:
                MPI_Status status;
                MPI_Request req;
                size_t recv_value;
                int flag;
                bool is_running_;
            public:
                MPI() : is_running( false ), recv_value( 0 ), flag( 0 ) {}
                virtual void Init() {
                    is_running_ = true;
                }
                virtual void AsyncRecv() {
                    MPI_Irecv(&recv_value, 1, MPI_LONG, MPI_ANY_SOURCE, MPI_ANY_TAG, comm_m, &req);
                }
                virtual bool HasMessage() {
                    if (req != MPI_REQUEST_NULL){
                        MPI_Test(&req, &flag, &status);
                        if(flag) return true;
                    }
                    return false;

                }
                virtual bool HandleMessage() {} =0;
                virtual bool ShouldStop() {
                    return status.MPI_TAG == MPI_STOP_TAG;
                }
                virtual void Shutdown() {

                    is_running = false;
                }
            };
            class Simulator : public Iface {};
        }
        */
    }
}
#endif
