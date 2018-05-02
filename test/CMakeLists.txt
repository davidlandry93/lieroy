cmake_minimum_required(VERSION 3.0)

find_package(GTest REQUIRED)

enable_testing()

set(SHORT_TESTS
  test_se3.cpp
  test_so3.cpp
  test_se3_gaussian_distribution.cpp
  test_se3_uniform_distribution.cpp)
add_executable(short_tests test_main.cpp ${SHORT_TESTS})
target_link_libraries(short_tests GTest::GTest GTest::Main pthread lieroy)

set(MEDIUM_TESTS
  test_se3_medium.cpp
  test_so3_medium.cpp
  test_se3_gaussian_distribution_medium.cpp)
add_executable(medium_tests test_main.cpp ${MEDIUM_TESTS})
target_link_libraries(medium_tests GTest::GTest GTest::Main pthread lieroy)
