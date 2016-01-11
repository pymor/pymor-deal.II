cmake_minimum_required(VERSION 2.8.8)
FIND_PACKAGE(deal.II 8.1 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()
DEAL_II_INITIALIZE_CACHED_VARIABLES()

project(DEALII_example)

# Add a CMake parameter for choosing a desired Python version
set(REQUIRED_PYTHON_VERSION "2.7" CACHE STRING "Python version to use for compiling the example library")

# Set a default build configuration if none is specified. 'MinSizeRel' produces the smallest binaries
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'MinSizeRel' as none was specified.")
  set(CMAKE_BUILD_TYPE MinSizeRel CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# Try to autodetect Python (can be overridden manually if needed)
find_package(PythonLibs ${REQUIRED_PYTHON_VERSION} REQUIRED)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
include_directories(${PYTHON_INCLUDE_DIR})
# pybind11 path
include_directories(/home/r_milk01/fast_projekte/uni/paper/15_pymor_paper/local/include)

foreach(_target dealii_bindings)
    # Create the binding library
    add_library(${_target} SHARED ${_target}.cc )
    # Don't add a 'lib' prefix to the shared library
    set_target_properties(${_target} PROPERTIES PREFIX "")
    DEAL_II_SETUP_TARGET(${_target})
    target_link_libraries(${_target} ${PYTHON_LIBRARY})
endforeach(_target example elliptic)
