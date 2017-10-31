project(pylie)
cmake_minimum_required(VERSION 3.1)

find_package(Eigen3 REQUIRED)

option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
option(BUILD_TESTS "Build tests" ON)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${EIGEN3_INCLUDE_DIR}
  )

if(BUILD_PYTHON_BINDINGS)
  add_subdirectory(python_bindings)
endif(BUILD_PYTHON_BINDINGS)

if(BUILD_TESTS)
  add_subdirectory(test)
endif(BUILD_TESTS)
