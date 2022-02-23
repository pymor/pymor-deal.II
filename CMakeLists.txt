cmake_minimum_required(VERSION 3.2)

# this should only be necessary while the deal setup macro is copied
set(CMAKE_CXX_STANDARD 17 CACHE STRING "")

project(pymor_dealii)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/lib)

FIND_PACKAGE(deal.II 9.3
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

# Set a default build configuration if none is specified. 'MinSizeRel' produces the smallest binaries
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
endif()
# Add a CMake parameter for choosing a desired Python version
set(PYBIND11_PYTHON_VERSION "3" CACHE STRING "Python version to use for compiling the example library")

add_subdirectory(lib)
add_subdirectory(src)
