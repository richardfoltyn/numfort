# src/CMakeLists.txt
# Top-level cmake configuration for numfort

cmake_minimum_required(VERSION 3.0)

# include C in project languages, needed for LAPACK detection
project(numfort C Fortran)
set(LIBRARY_NAME ${CMAKE_PROJECT_NAME} CACHE STRING "File name of main library")

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

################################################################################
# User-definable options
if (WIN32)
    option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
else ()
    option(BUILD_SHARED_LIBS "Build shared libraries" ON)
endif ()
option(USE_MKL "Build against Intel MKL" ON)
option(BUILD_TESTS "Build (unit) test executables" OFF)
option(BUILD_EXAMPLES "Build examples/demo code" OFF)
option(BUILD_BENCHMARKS "Build benchmarks" OFF)

################################################################################
# Additional CMAKE modules
# include FindMKL module taken from ARMA project

list(APPEND CMAKE_MODULE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/submodules/cmake/Modules")

if (USE_MKL)
    find_package(MKL REQUIRED)
    set(${PROJECT_NAME}_BLA_LIBRARIES ${MKL_LIBRARIES})
else ()
    find_package(LAPACK REQUIRED)
    set(${PROJECT_NAME}_BLA_LIBRARIES ${LAPACK_LIBRARIES})
    set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
    set(CMAKE_SHARED_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
endif ()

################################################################################
# Testing

if(BUILD_TESTS)
    enable_testing()
    find_package(FCore)
endif()

###############################################################################
# Path suffixes beneath CMAKE_INSTALL_PREFIX

# if INSTALL_LIBRARY_PREFIX or INSTALL_INCLUDE_PREFIX were specified on the
# command line, use those as prefixes beneath CMAKE_INSTALL_PREFIX where
# library and include files should be placed.

if (NOT INSTALL_LIBRARY_PREFIX)
    set(INSTALL_LIBRARY_PREFIX lib)
endif()

if (NOT INSTALL_INCLUDE_PREFIX)
    set(INSTALL_INCLUDE_PREFIX include/numfort)
endif()

################################################################################
# Set Fortran compiler FLAGS

set(_FFLAGS)
# FFLAGS for external libraries (these are most likely F77 code, avoid
# enforcing any Fortran standards, etc.)
set(_FFLAGS_EXTERNAL)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # explicitly add -cpp to run C preprocessor if required
    set(GNU_WARN "-Wall -Wextra -Wimplicit-interface -Wimplicit-procedure \
      -Warray-temporaries -Wrealloc-lhs -pedantic")
    set(_FFLAGS "${GNU_WARN} -std=f2008")
    # enable optimization that does not break debugging
    # set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -Og")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (WIN32)
        set (_FFLAGS "/warn:all /standard-semantics /stand:f08")
    else ()
        set (_FFLAGS "-warn all -standard-semantics -std08")
    endif (WIN32)
endif ()

set (${PROJECT_NAME}_FFLAGS "${_FFLAGS}")
set (${PROJECT_NAME}_FFLAGS_EXTERNAL "${_FFLAGS_EXTERNAL}")

# append additional FLAGS to default cmake FLAGS
# set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_FFLAGS}")

###############################################################################
# Report options and settings
if (BUILD_SHARED_LIBS)
    message(STATUS "Building SHARED libraries")
else(BUILD_SHARED_LIBS)
    message(STATUS "Building STATIC libraries")
endif(BUILD_SHARED_LIBS)

if (USE_MKL)
	message(STATUS "Building against Intel MKL")
endif (USE_MKL)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Default compiler flags: ${CMAKE_Fortran_FLAGS}")
message(STATUS "Add. compiler flags: ${${PROJECT_NAME}_FFLAGS}")
message(STATUS "Add. compiler flags for external libs: ${${PROJECT_NAME}_FFLAGS_EXTERNAL}")
# Make sure build type is upper case, otherwise the wrong flags will be shown
string(TOUPPER "${CMAKE_BUILD_TYPE}" _BUILD_TYPE)
message(STATUS "Add. compiler flags for ${CMAKE_BUILD_TYPE}: ${CMAKE_Fortran_FLAGS_${_BUILD_TYPE}}")

message(STATUS "Library install directory: ${CMAKE_INSTALL_PREFIX}/${INSTALL_LIBRARY_PREFIX}")
message(STATUS "MOD file install directory: ${CMAKE_INSTALL_PREFIX}/${INSTALL_INCLUDE_PREFIX}")

################################################################################
# Output directories

set(${PROJECT_NAME}_MODDIR "${CMAKE_BINARY_DIR}/include")
set(CMAKE_Fortran_MODULE_DIRECTORY "${${PROJECT_NAME}_MODDIR}")

# store test executables in this root directory
set(${PROJECT_NAME}_TESTS_DIRECTORY "${CMAKE_BINARY_DIR}/bin/tests")
set(${PROJECT_NAME}_EXAMPLES_DIRECTORY "${CMAKE_BINARY_DIR}/bin/examples")
set(${PROJECT_NAME}_BENCHMARKS_DIRECTORY "${CMAKE_BINARY_DIR}/bin/benchmarks")

add_subdirectory(numfort)

if (BUILD_TESTS)
    add_subdirectory(tests)
endif()

if (BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

if (BUILD_EXAMPLES)
    add_subdirectory(examples)
endif ()
