include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/dependencies)
include_directories(${CMAKE_SOURCE_DIR}/dependencies/gtest-1.6.0/include)
include_directories(${CMAKE_SOURCE_DIR}/dependencies/nlopt-2.3/api)

enable_testing()

add_executable(testNlopt
    test_opt.cpp
)
target_link_libraries(testNlopt gtest gtest_main nlopt Problems)
add_test(
    NAME testNlopt
    COMMAND testNlopt
)


add_executable(testLambda
    test_lambda.cpp
)
target_link_libraries(testLambda gtest gtest_main nlopt Problems )
add_test(
    NAME testLambda
    COMMAND testLambda
)


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
