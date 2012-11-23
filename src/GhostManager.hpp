#ifndef JobQueue_hpp
#define JobQueue_hpp

#include "HomotopyTypes.h"
#include <map>
#include <queue>
#include "mesh.hpp"

namespace Homotopy {

    typedef Mesh::ind_t ind_t;
    typedef Mesh::MeshBase::MeshPoint Point_t;
    typedef std::pair< ind_t, objVars_t > nodeEnvelope_t;
    typedef std::pair< ind_t, std::vector< nodeEnvelope_t > > NeighborMessage_t;
    
    template< typename T>
    class GhostManager {

        T& Exchanger;
        std::map< ind_t, Point_t* > &GlobalToGhost;
        std::map< ind_t, std::vector< Point_t* > > &NeighborSubsetsToLocals;
 
    public:
        GhostManager( T& exchanger, std::map< ind_t, Point_t* > &globaltoghost, std::map< ind_t, std::vector< Point_t* > > &neighborsubsetstolocals  ) : 
          Exchanger( exchanger ), GlobalToGhost( globaltoghost ), NeighborSubsetsToLocals( neighborsubsetstolocals ) {}

        void Exchange() {
            //Container of send requests format: pair( To, vector Content )
            //Complicated, but should be automatically serializable...
            std::queue< NeighborMessage_t > ToSend;

            //iterate over the neighbors
            std::map< ind_t, std::vector< Point_t* > >::iterator n;
            for(n=NeighborSubsetsToLocals.begin(); n!=NeighborSubsetsToLocals.end(); n++){
                std::vector< nodeEnvelope_t > toNode;

                //gather the points to send to these neighbors
                std::vector< Point_t* >::iterator v;
                for(v=n->second.begin(); v!=n->second.end(); v++){
                    nodeEnvelope_t node( (*v)->ID, (*v)->ObjectiveCoords );
                    toNode.push_back( node );
                }
                ToSend.push( std::pair< ind_t, std::vector< nodeEnvelope_t > > ( n->first, toNode ) );
            }

            //Container for results
            std::queue< nodeEnvelope_t > newNeighbors;
            Exchanger.exchanger( ToSend, newNeighbors );

            //Process ghost updates
            while( !newNeighbors.empty() ){
                nodeEnvelope_t &neighbor = newNeighbors.front();

                //Do we have this registered as a ghost
                std::map< ind_t, Point_t* >::iterator ghost_pos = GlobalToGhost.find( neighbor.first );
                if( ghost_pos != GlobalToGhost.end() ) {
                    (*ghost_pos).second->ObjectiveCoords = neighbor.second;
                }
                else {
                    //Something is wrong...
                    printf("GhostManager received erroneous node!\n");
                }
                newNeighbors.pop();
            }
        }
    };

}

#endif
