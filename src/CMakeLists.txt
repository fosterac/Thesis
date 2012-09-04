include_directories(${CMAKE_SOURCE_DIR}/dependencies)
include_directories(${CMAKE_SOURCE_DIR}/dependencies/nlopt-2.3/api)

add_library (Problems STATIC Problems.cpp)
target_link_libraries(Problems wfgProblems)

add_library (Optimizer STATIC optimizer.cpp)
target_link_libraries(Optimizer nlopt)