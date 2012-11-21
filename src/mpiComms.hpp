#ifdef mpiComms_hpp
#define mpiComms_hpp

#include "CommInterface.hpp"

namespace CommImpl {

    typedef MPI_Status Status_t;

    class MPI : public Iface {
    private:
        MPI_Status status;
        MPI_Request req;
        size_t recv_value;
        int flag;
        bool is_running_;

    public:
        MPI_Comm comm_m;

        MPI() : is_running_( false ), recv_value( 0 ), flag( 0 ) {}
        void Init() {
            is_running_ = true;
        }
        void AsyncRecv() {
            MPI_Irecv(&recv_value, 1, MPI_LONG, MPI_ANY_SOURCE,
                        MPI_ANY_TAG, comm_m, &req);
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

            is_running_ = false;
        }

        MPI_Status GetStatus(){ return status; }
        size_t GetValue(){ return recv_value; }
    };
}

#endif