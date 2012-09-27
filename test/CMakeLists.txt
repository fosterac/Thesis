include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/dependencies)
include_directories(${CMAKE_SOURCE_DIR}/dependencies/gtest-1.6.0/include)
include_directories(${CMAKE_SOURCE_DIR}/dependencies/nlopt-2.3/api)

enable_testing()

add_executable(testSandbox
    test_sandbox.cpp
)
target_link_libraries(testSandbox gtest gtest_main nlopt Problems )
add_test(
    NAME testSandbox
    COMMAND testSandbox
)

add_executable(testNlopt
    test_nlopt.cpp
)
target_link_libraries(testNlopt gtest gtest_main nlopt )
add_test(
    NAME testNlopt
    COMMAND testNlopt
)

file( COPY data/ DESTINATION data/)
add_executable(testHypervolume
    test_hypervolume.cpp
)
target_link_libraries(testHypervolume gtest gtest_main nlopt wfgHypervolume)
add_test(
    NAME testHypervolume
    COMMAND testHypervolume
)

add_executable(testWFGPROBLEMS
    test_wfgproblems.cpp
)
target_link_libraries(testWFGPROBLEMS gtest gtest_main nlopt wfgProblems)
add_test(
    NAME testWFGPROBLEMS
    COMMAND testWFGPROBLEMS
)

add_executable(testProblems
    test_problems.cpp
)
target_link_libraries(testProblems gtest gtest_main nlopt Problems Optimizer)
add_test(
    NAME testProblems
    COMMAND testProblems
)

add_executable(testMesh
    test_mesh.cpp
)
target_link_libraries(testMesh gtest gtest_main )
add_test(
    NAME testMesh
    COMMAND testMesh
)

add_executable(testOptimizer
    test_optimizer.cpp
)
target_link_libraries(testOptimizer gtest gtest_main nlopt Problems Optimizer )
add_test(
    NAME testOptimizer
    COMMAND testOptimizer
)

add_executable(testMarching
    test_marching.cpp
)
target_link_libraries(testMarching gtest gtest_main nlopt Problems Optimizer )
add_test(
    NAME testMarching
    COMMAND testMarching
)