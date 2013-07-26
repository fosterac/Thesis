add_subdirectory( nlopt-2.3 )
if(EXISTS 
"${CMAKE_CURRENT_SOURCE_DIR}/lib/CMakeLists.txt")
   add_subdirectory( lib )
endif()

add_subdirectory( alglib-3.6.0 )

add_subdirectory( wfgProblems )
add_subdirectory( wfgHypervolume )

