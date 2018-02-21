cmake_minimum_required(VERSION 3.7)

option(USE_PYTHON3 "Wether to use Python 3 instead of Python 2" ON)
option(AUTO_INSTALL_BINDINGS "Install the python package of the bindings at install time" OFF)

if(USE_PYTHON3)
  set(PYTHON_COMPONENTS_SUFFIX 3)
else(USE_PYTHON3)
  set(PYTHON_COMPONENTS_SUFFIX "")
endif(USE_PYTHON3)

find_package(Boost 1.63 COMPONENTS python${PYTHON_COMPONENTS_SUFFIX} numpy${PYTHON_COMPONENTS_SUFFIX})
find_package(Eigen3 REQUIRED)
find_package(PythonLibs)

add_library(lieroy_bindings SHARED lieroy/lieroy.cpp)
target_include_directories(lieroy_bindings PRIVATE
  ${Boost_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_DIRS})
target_link_libraries(lieroy_bindings lieroy ${Boost_LIBRARIES})
set_target_properties(lieroy_bindings PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lieroy
  OUTPUT_NAME "lieroy"
  PREFIX ""
  SUFFIX ".so")

add_custom_command(TARGET lieroy_bindings PRE_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
  ${CMAKE_CURRENT_SOURCE_DIR}/lieroy/__init__.py
  ${CMAKE_CURRENT_BINARY_DIR}/lieroy/__init__.py)

set(SETUP_FILE ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in ${SETUP_FILE})

# === INSTALLATION ===

find_program(PIP "pip")

if(PIP)
  set(PYTHON_BINDINGS_INSTALLATION_COMMAND  ${PIP} install -I .)

  # Install the package at install time.
  install(CODE "execute_process(COMMAND ${PYTHON_BINDINGS_INSTALLATION_COMMAND} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})")

  if(AUTO_INSTALL_BINDINGS)
    message(STATUS "Auto install bindings is ON")

    # Add a custom command that installs the package every time the bindings a rebuilt.
    add_custom_command(TARGET lieroy_bindings POST_BUILD
      COMMAND ${PYTHON_BINDINGS_INSTALLATION_COMMAND}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  endif(AUTO_INSTALL_BINDINGS)
endif(PIP)

