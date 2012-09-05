cmake_minimum_required(VERSION 2.8)

SET(COMMON_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
SET(COMMON_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/libs")

add_subdirectory( dependencies )
add_subdirectory( src )
add_subdirectory( test )