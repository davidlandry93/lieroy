cmake_minimum_required(VERSION 3.7)

option(USE_PYTHON3 "Wether to use Python 3 instead of Python 2" ON)
option(AUTO_INSTALL_BINDINGS "Install the python package of the bindings at install time" OFF)

if(USE_PYTHON3)
  set(PYTHON_COMPONENTS_SUFFIX 3)
else(USE_PYTHON3)
  set(PYTHON_COMPONENTS_SUFFIX "")
endif(USE_PYTHON3)

find_package(Boost 1.63 COMPONENTS python${PYTHON_COMPONENTS_SUFFIX} numpy${PYTHON_COMPONENTS_SUFFIX} REQUIRED)
find_package(Eigen3 REQUIRED)

set(Python_ADDITIONAL_VERSIONS 3.5)
find_package(PythonLibs 3)

message(STATUS ${Boost_INCLUDE_DIRS})

add_library(lieroy_bindings SHARED lieroy/lieroy_core.cpp)
target_include_directories(lieroy_bindings PRIVATE
  ${PYTHON_INCLUDE_OVERRIDE}
  ${Boost_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  ${PYTHON_INCLUDE_DIRS})
target_link_libraries(lieroy_bindings lieroy ${Boost_LIBRARIES})
set_target_properties(lieroy_bindings PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lieroy
  OUTPUT_NAME "lieroy_core"
  PREFIX ""
  SUFFIX ".so")


# Copy the python files to the build tree where the actual python package is installed from.
file(GLOB PY_SOURCE_FILES
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/lieroy
  "lieroy/*.py")

foreach(PY_FILE ${PY_SOURCE_FILES})
  add_custom_command(TARGET lieroy_bindings POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/lieroy/${PY_FILE}
    ${CMAKE_CURRENT_BINARY_DIR}/lieroy/${PY_FILE})
endforeach(PY_FILE)


set(SETUP_FILE ${CMAKE_CURRENT_BINARY_DIR}/setup.py)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in ${SETUP_FILE})

# === INSTALLATION ===

find_program(PIP "pip")

if(PIP)
  set(PYTHON_BINDINGS_INSTALLATION_COMMAND  ${PIP} install --ignore-installed .)

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

