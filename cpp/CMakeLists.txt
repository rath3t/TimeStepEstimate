cmake_minimum_required(VERSION 3.16)
project(timeStepEstimate CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
  ${dune-common_MODULE_PATH})

#include the dune macros
include(DuneMacros)
add_definitions(-DDUNE_LOCALFEFUNCTIONS_USE_EIGEN=1)
add_definitions(-DEIGEN_DEFAULT_TO_ROW_MAJOR=1)
add_definitions(-DBUILD_SHARED_LIBS=1)
add_definitions(-DOMPI_SKIP_MPICXX=0)
# start a dune project with information from dune.module
dune_project()

dune_enable_all_packages()

add_subdirectory(src)
add_subdirectory(dune)
add_subdirectory(doc)
add_subdirectory(cmake/modules)
add_subdirectory(python)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project()
