cmake_minimum_required(VERSION 3.5)

project(lieroy VERSION 0.1)

find_package(Eigen3 REQUIRED)
find_package(OpenMP)

option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
option(BUILD_TESTS "Build tests" ON)


add_library(lieroy INTERFACE)
target_include_directories(lieroy BEFORE INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
  ${EIGEN3_INCLUDE_DIR})

# Even though it is header only, lieroy still requires the consumer to use a c++11 compiler.
target_compile_features(lieroy INTERFACE cxx_auto_type cxx_range_for cxx_rvalue_references)

if(OpenMP_FOUND)
  target_compile_options(lieroy INTERFACE ${OpenMP_CXX_FLAGS})
  target_link_libraries(lieroy INTERFACE ${OpenMP_CXX_FLAGS})
  message(STATUS Using OpenMP)
endif(OpenMP_FOUND)


add_executable(distribution_of_se3 src/distribution_of_se3.cpp)
target_link_libraries(distribution_of_se3 lieroy)

if(BUILD_PYTHON_BINDINGS)
  add_subdirectory(python)
endif(BUILD_PYTHON_BINDINGS)

if(BUILD_TESTS)
  add_subdirectory(test)
endif(BUILD_TESTS)


# === INSTALLATION ===
include(CMakePackageConfigHelpers)
set(LIEROY_INSTALL_DIR lib/lieroy)

# Generate the lieroy-config file which makes lieroy available to find_package commands.
configure_package_config_file(lieroy-config.cmake.in ${CMAKE_BINARY_DIR}/lieroy-config.cmake
  INSTALL_DESTINATION ${LIEROY_INSTALL_DIR}/cmake)

# Install the executables.

# Install the header files of the library.
install(DIRECTORY include/lieroy
  DESTINATION include
  FILES_MATCHING REGEX ".h$|.hpp$")

# Define lieroy_export as the export for target lieroy.
install(TARGETS lieroy distribution_of_se3
  EXPORT lieroy-targets
  RUNTIME DESTINATION bin)

# Install an export file which allows consumers to use the target lieroy::lieroy in their project.
install(EXPORT lieroy-targets
  NAMESPACE lieroy::
  DESTINATION ${LIEROY_INSTALL_DIR}/cmake)

# Installe the auto generated files to make the targets available through find_package.
install(FILES ${CMAKE_BINARY_DIR}/lieroy-config.cmake
  DESTINATION ${LIEROY_INSTALL_DIR}/cmake)
