cmake_minimum_required(VERSION 3.7)

option(USE_PYTHON3 "Wether to use Python 3 instead of Python 2" ON)

if(USE_PYTHON3)
  set(PYTHON_COMPONENTS_SUFFIX 3)
else(USE_PYTHON3)
  set(PYTHON_COMPONENTS_SUFFIX "")
endif(USE_PYTHON3)

find_package(Boost 1.63 COMPONENTS python${PYTHON_COMPONENTS_SUFFIX} numpy${PYTHON_COMPONENTS_SUFFIX})
find_package(Eigen3 REQUIRED)
find_package(PythonLibs)

set(PYLIE_BINDINGS_OUTPUT_NAME "pylie")
add_library(pylie_bindings SHARED pylie/pylie.cpp)
target_include_directories(pylie_bindings PRIVATE
  ${Boost_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_DIRS})
target_link_libraries(pylie_bindings pylie ${Boost_LIBRARIES})
set_target_properties(pylie_bindings PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/pylie
  OUTPUT_NAME ${PYLIE_BINDINGS_OUTPUT_NAME}
  PREFIX ""
  SUFFIX ".so")

# === INSTALLATION ===

find_program(PYTHON "python")

if(PYTHON)
  # Placeholder file to remember every time we built the python package.
  set(PACKAGE_WAS_BUILT_FILE ${CMAKE_BINARY_DIR}/timestamp)
  set(SETUP_FILE ${CMAKE_CURRENT_SOURCE_DIR}/setup.py)

  # Generate a setup.py file from template.
  get_target_property(PYLIE_INCLUDE_DIRS pylie INTERFACE_INCLUDE_DIRECTORIES)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in ${SETUP_FILE})

  # Add a custom command that builds the package.
  add_custom_command(OUTPUT ${PACKAGE_WAS_BUILT_FILE}
    COMMAND ${PYTHON} ${SETUP_FILE} build
    COMMAND ${CMAKE_COMMAND} -E touch ${PACKAGE_WAS_BUILT_FILE})

  install(CODE "execute_process(COMMAND ${PYTHON} ${SETUP_FILE} install)")
endif(PYTHON)
