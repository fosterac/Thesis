#ifndef JobQueue_hpp
#define JobQueue_hpp

#include "HomotopyTypes.h"
#include <map>
#include <queue>
#include <algorithm>
#include "mesh.hpp"

namespace Homotopy {

    typedef Mesh::ind_t ind_t;
    typedef Mesh::MeshBase::MeshPoint Point_t;
    typedef std::pair< ind_t, objVars_t > nodeEnvelope_t;
    typedef std::pair< ind_t, std::vector< nodeEnvelope_t > > NeighborMessage_t;

    template< typename T>
    class GhostManager {

        T& Exchanger;
        std::vector< ind_t > IDs;
        std::map< ind_t, std::vector< Point_t* > > GlobalToGhosts;
        std::map< ind_t, std::vector< Point_t* > > SubsToLocals;
 
    public:
        GhostManager( T& exchanger ) : Exchanger( exchanger ) {}

        void Add( ind_t ID, std::map< ind_t, Point_t* > &globaltoghost, std::map< ind_t, std::vector< Point_t* > > &neighborsubsetstolocals ) {
            this->IDs.push_back( ID );

            std::map< ind_t, Point_t* >::iterator gg;
            for(gg=globaltoghost.begin(); gg!=globaltoghost.end(); gg++){
                if( this->GlobalToGhosts.find( gg->first ) == this->GlobalToGhosts.end() ) {
                    this->GlobalToGhosts[ gg->first ] = std::vector< Point_t* > ();
                }
                this->GlobalToGhosts[ gg->first ].push_back( gg->second );
            }
            std::map< ind_t, std::vector< Point_t* > >::iterator nstl;
            for(nstl=neighborsubsetstolocals.begin(); nstl!=neighborsubsetstolocals.end(); nstl++){
                if( this->SubsToLocals.find( nstl->first ) == this->SubsToLocals.end() ) {
                    this->SubsToLocals[ nstl->first ] = std::vector< Point_t* > ();
                }
                this->SubsToLocals[ nstl->first ].insert( this->SubsToLocals[ nstl->first ].end(), nstl->second.begin(), nstl->second.end() );
            }
        }          

        void Exchange() {
            //Container of send requests format: pair( To, vector Content )
            //Complicated, but should be automatically serializable...
            std::queue< NeighborMessage_t > ToSend;

            //iterate over the neighbors
            std::map< ind_t, std::vector< Point_t* > >::iterator n;
            for(n=SubsToLocals.begin(); n!=SubsToLocals.end(); n++){
                
                //If the node is not managed by this instance,
                //then prepare a message to send
                if( std::find( this->IDs.begin(), this->IDs.end(), n->first ) == this->IDs.end() ){
                    std::vector< nodeEnvelope_t > toNode;
                    //gather the points to send to these neighbors
                    std::vector< Point_t* >::iterator v;
                    for(v=n->second.begin(); v!=n->second.end(); v++){
                        nodeEnvelope_t node( (*v)->ID, (*v)->ObjectiveCoords );
                        toNode.push_back( node );
                    }
                    ToSend.push( std::pair< ind_t, std::vector< nodeEnvelope_t > > ( n->first, toNode ) );
                }
                //Process locally managed updates
                else {
                    std::vector< Point_t* >::iterator p;
                    for(p=n->second.begin(); p!=n->second.end(); p++){
                        //Maybe we should check this?
                        std::vector< Point_t* >& v = this->GlobalToGhosts[ (*p)->ID ];

                        std::vector< Point_t* >::iterator it;
                        for( it=v.begin(); it!=v.end(); it++ ){
                            (*it)->ObjectiveCoords = (*p)->ObjectiveCoords;
                        }
                    }
                }
            }

            //Container for results
            std::queue< nodeEnvelope_t > newNeighbors;
            Exchanger.exchanger( ToSend, newNeighbors );

            //Process received node updates
            while( !newNeighbors.empty() ){
                nodeEnvelope_t &neighbor = newNeighbors.front();

                //Do we have this registered as a ghost
                std::map< ind_t, std::vector< Point_t* > >::iterator ghost_pos = GlobalToGhosts.find( neighbor.first );
                if( ghost_pos != GlobalToGhosts.end() ) {
                    std::vector< Point_t* >::iterator it;
                    for( it=ghost_pos->second.begin(); it!=ghost_pos->second.end(); it++ ){
                        (*it)->ObjectiveCoords = neighbor.second;
                    }
                }
                else {
                    //Something is wrong...
                    printf("GhostManager received erroneous node!\n");
                }
                newNeighbors.pop();
            }
        }
    };


    template< typename T, typename C >
    class MeshSet : Mesh::MeshBase {
    private:
        std::vector< T* > meshes;
        GhostManager< C > gman;
    public:
        MeshSet(    std::vector< Mesh::point_t > DesignSpace,
                    std::vector< Mesh::point_t > ObjectiveSpace,
					std::vector< Mesh::point_t > Lambdas, 
					int NumberOfPoints, 
                    std::vector<Mesh::ind_t> IDs, 
                    Mesh::ind_t SubsetsPerSide) : 
                                            MeshBase( DesignSpace, ObjectiveSpace, Lambdas ) {
            std::vector< ind_t >::iterator i;
            for(i=IDs.begin(); i!=IDs.end(); i++){
                meshes.push_back( new T(    this->DesignSpace, this->ObjectiveSpace, this->Lambdas, 
                                            NumberOfPoints, *i, SubsetsPerSide ) );
                gman.Add( *i, meshes.back().GlobalToGhost, meshes.back().NeighborSubsetsToLocals );
            }
        }
        virtual void Generate() {
            typename std::vector< T* >::iterator i;
            for(i=meshes.begin(); i!=meshes.end(); i++){
                i->Generate();
            }
        }
        virtual void Refresh() {
            this->gman.Exchange();
            typename std::vector< T* >::iterator i;
            for(i=meshes.begin(); i!=meshes.end(); i++){
                i->Refresh();
            }
        }
        virtual void Print() {
            typename std::vector< T* >::iterator i;
            for(i=meshes.begin(); i!=meshes.end(); i++){
                i->Print();
            }
        }
        virtual void WriteOut( const char* f ) {
            //Don't do anything for now, since the outfile
            //will need to persist.
        }
    };
}

#endif
