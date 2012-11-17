cmake_minimum_required(VERSION 2.8)

find_package ( MPI )
if ( MPI_FOUND )
        add_definitions( -DHAS_MPI )
	include_directories(${MPI_INCLUDE_PATH})
endif( MPI_FOUND )

add_subdirectory( dependencies )
add_subdirectory( src )
add_subdirectory( test )
