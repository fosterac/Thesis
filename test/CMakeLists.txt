include_directories(${CMAKE_SOURCE_DIR}/src)

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

