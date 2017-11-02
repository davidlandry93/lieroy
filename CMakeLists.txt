cmake_minimum_required(VERSION 3.8)

project(pylie VERSION 0.1)

find_package(Eigen3 REQUIRED)

option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
option(BUILD_TESTS "Build tests" ON)


add_library(pylie INTERFACE)
target_include_directories(pylie INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>
  ${EIGEN3_INCLUDE_DIR})
# Even though it is header only, pylie still requires the consumer to use a c++11 compiler.
target_compile_features(pylie INTERFACE cxx_std_11)


add_executable(distribution_of_se3 src/distribution_of_se3.cpp)
target_link_libraries(distribution_of_se3 pylie)

if(BUILD_PYTHON_BINDINGS)
  add_subdirectory(python)
endif(BUILD_PYTHON_BINDINGS)

if(BUILD_TESTS)
  add_subdirectory(test)
endif(BUILD_TESTS)


# === INSTALLATION ===
include(CMakePackageConfigHelpers)
set(PYLIE_INSTALL_DIR lib/pylie)

# Generate the pylie-config file which makes pylie available to find_package commands.
configure_package_config_file(pylie-config.cmake.in ${CMAKE_BINARY_DIR}/pylie-config.cmake
  INSTALL_DESTINATION ${PYLIE_INSTALL_DIR}/cmake)

# Install the executables.

# Install the header files of the library.
install(DIRECTORY include/pylie
  DESTINATION include/pylie
  FILES_MATCHING REGEX ".h$|.hpp$")

# Define pylie_export as the export for target pylie.
install(TARGETS pylie distribution_of_se3
  EXPORT pylie-targets
  RUNTIME DESTINATION bin)

# Install an export file which allows consumers to use the target pylie::pylie in their project.
install(EXPORT pylie-targets
  NAMESPACE pylie::
  DESTINATION ${PYLIE_INSTALL_DIR}/cmake)

# Installe the auto generated files to make the targets available through find_package.
install(FILES ${CMAKE_BINARY_DIR}/pylie-config.cmake
  DESTINATION ${PYLIE_INSTALL_DIR}/cmake)
