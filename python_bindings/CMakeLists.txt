cmake_minimum_required(VERSION 3.7)

option(USE_PYTHON3 "Wether to use Python 3 instead of Python 2" ON)

if(${USE_PYTHON3})
  set(Python_ADDITIONAL_VERSIONS 3)
  find_package(Boost 1.63 COMPONENTS python3 numpy3)
else(${USE_PYTHON3})
  find_package(Boost 1.63 COMPONENTS python numpy)
endif(${USE_PYTHON3})

find_package(Eigen3 REQUIRED)
find_package(PythonLibs)

message(STATUS "Python Version Found: ${PYTHONLIBS_VERSION_STRING}")

include_directories(
  ${Boost_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/src)

add_library(pylie_bindings SHARED pylie.cpp)
target_link_libraries(pylie_bindings pylie ${Boost_LIBRARIES})
set_target_properties(pylie_bindings PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  OUTPUT_NAME "pylie"
  PREFIX "")
