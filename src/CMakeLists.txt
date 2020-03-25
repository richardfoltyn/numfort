# src/CMakeLists.txt
# Top-level cmake configuration for numfort

cmake_minimum_required(VERSION 3.0)

project(numfort C Fortran)
set(NUMFORT_LIBRARY_NAME ${CMAKE_PROJECT_NAME})

# number of version components
set(NUMFORT_VERSION_COUNT 3)
set(NUMFORT_VERSION_MAJOR 0)
set(NUMFORT_VERSION_MINOR 1)
set(NUMFORT_VERSION_PATCH 0)
set(NUMFORT_VERSION
    "${NUMFORT_VERSION_MAJOR}.${NUMFORT_VERSION_MINOR}.${NUMFORT_VERSION_PATCH}"
)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE
    )
endif (NOT CMAKE_BUILD_TYPE)

################################################################################
# User-definable options
if (WIN32)
    option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
else ()
    option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
endif ()
option(USE_MKL "Build against Intel MKL" ON)
option(BUILD_BLAS95 "Build Fortran 95 wrappers for BLAS" OFF)
option(BUILD_TESTS "Build (unit) tests" OFF)
option(BUILD_EXAMPLES "Build examples/demo code" OFF)
option(OPT_REPORT "Generate compiler optimization reports" OFF)
option(USE_OPENMP "Use some OpenMP features" OFF)

# If enabled, pass -assume std_mod_proc_name to ifort such that module
# names will use _MP_ as module name / procedure name separator.
# Off by default, this causes major pain since most third-party libraries
# are compiled without this option.
option(IFORT_STD_MOD_PROC_NAME
    "Assume standard-compliant name mangling for module procedures" OFF
)

################################################################################
# Additional CMAKE modules

list(APPEND CMAKE_MODULE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/../submodules/cmake/Modules")

################################################################################
# Detect features supported by Fortran compiler

include(CheckFortranKinds)
check_Fortran_kinds()

add_definitions(
    -D__NUMFORT_DEFAULT_INT_KIND__=${Fortran_DEFAULT_INT_KIND}
    -D__NUMFORT_DEFAULT_REAL_KIND__=${Fortran_DEFAULT_REAL_KIND}
    -D__NUMFORT_DEFAULT_INT_SIZE__=${Fortran_DEFAULT_INT_SIZE}
)

################################################################################
# Testing

if(BUILD_TESTS)
    enable_testing()
endif()

################################################################################
# BLAS / LAPACK libraries

# 1. Determine if we need a BLAS/LAPACK library
# Note that to compile a static NUMFORT library without any unit tests,
# examples, etc. does not require BLAS/LAPACK.

set(NEED_BLAS_LAPACK OFF)

if (BUILD_TESTS OR BUILD_EXAMPLES)
    set(NEED_BLAS_LAPACK ON)
endif()

# Dynamically-linked libraries need to be linked
if (BUILD_SHARED_LIBS)
    set(NEED_BLAS_LAPACK ON)
endif()

# 2. Try to find a BLAS / LAPACK library
if (NEED_BLAS_LAPACK)
    if (USE_MKL)
        find_package(MKL OPTIONAL_COMPONENTS SEQUENTIAL BLAS95 LAPACK95)

        if (NOT MKL_FOUND)
            message(WARNING "Intel MKL requested but NOT FOUND.")
        else()
            set(BLA_LIBRARIES ${MKL_LIBRARIES})
        endif()

        if (MKL_BLAS95_FOUND)
            list(APPEND BLA_INCLUDE_DIRS ${MKL_BLAS95_INCLUDE_DIR})
        endif()

        if (MKL_LAPACK95_FOUND)
            list(APPEND BLA_INCLUDE_DIRS ${MKL_LAPACK95_INCLUDE_DIR})
        endif()

        if (NOT MKL_BLAS95_FOUND)
            message(WARNING "No BLAS Fortran 95 interfaces found, building NUMFORT wrapper")
            set(BUILD_BLAS95 ON)
        endif()
    endif ()

    # Fallback: if MKL was not found or no MKL requested, try to find whatever
    # BLAS/LAPACK library is present on the system.
    if (NOT USE_MKL OR NOT MKL_FOUND)
        find_package(LAPACK REQUIRED)
        set(BLA_LIBRARIES ${LAPACK_LIBRARIES})
        set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
        set(CMAKE_SHARED_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
        # Examples / unit tests use blas95 interfaces, so build these
        # in case we are not using MKL Fortran95 interfaces.
        set(BUILD_BLAS95 ON)
    endif ()
endif()


###############################################################################
# Path suffixes beneath CMAKE_INSTALL_PREFIX

# if CMAKE_INSTALL_LIBDIR or CMAKE_INSTALL_INCLUDEDIR were specified on the
# command line, use those as prefixes beneath CMAKE_INSTALL_PREFIX where
# library and include files should be placed.

# store values in tempoprary variables, will be overwritten by GNUInstallDirs
set(_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}")
set(_INSTALL_INCLUDEDDIR "${CMAKE_INSTALL_INCLUDEDIR}")

include(GNUInstallDirs)

# fill in any missing values
if (_INSTALL_LIBDIR)
    set(CMAKE_INSTALL_LIBDIR "${_INSTALL_LIBDIR}")
elseif (NOT CMAKE_INSTALL_LIBDIR)
    set(CMAKE_INSTALL_LIBDIR lib)
endif()

if (_INSTALL_INCLUDEDDIR)
    set(CMAKE_INSTALL_INCLUDEDIR "${_INSTALL_INCLUDEDDIR}")
elseif (NOT CMAKE_INSTALL_INCLUDEDIR)
    set(CMAKE_INSTALL_INCLUDEDIR include)
endif()

################################################################################
# Installation path


# Allow for installing ifort and gfortran version of libraries side-by-side by
# appending compiler suffix
string(TOLOWER ${CMAKE_Fortran_COMPILER_ID} Fortran_COMPILER_ID)

# Find major compiler version
string(REGEX REPLACE "(^[0-9]+).*" "\\1" Fortran_COMPILER_VERSION_MAJOR
    "${CMAKE_Fortran_COMPILER_VERSION}")

set(Fortran_COMPILER_SUFFIX "${Fortran_COMPILER_ID}-${Fortran_COMPILER_VERSION_MAJOR}")

# Installation directory within <prefix>/lib, <prefix>/include, etc.
# If ont defined by -DNUMFORT_INSTALL_DIR by the user, create
# version- and compiler-specific directory name.
if (NOT DEFINED NUMFORT_INSTALL_DIR)
    set(_version ${NUMFORT_VERSION_MAJOR}.${NUMFORT_VERSION_MINOR})
    set(NUMFORT_INSTALL_DIR
        "${NUMFORT_LIBRARY_NAME}-${_version}-${Fortran_COMPILER_SUFFIX}"
        CACHE STRING "Install directory name"
        FORCE
    )
endif()


################################################################################
# Output directories

set(NUMFORT_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/include/numfort")
set(CMAKE_Fortran_MODULE_DIRECTORY "${NUMFORT_MODULE_DIRECTORY}")

# include directory containing CPP header files
set(NUMFORT_HEADER_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/numfort/include")

# set include directory to directory where MOD files are stored
list(APPEND NUMFORT_INCLUDE_DIRS
    "${NUMFORT_MODULE_DIRECTORY}"
    "${CMAKE_BINARY_DIR}/include/numfort/external"
    "${NUMFORT_HEADER_INCLUDE_DIR}"
)

# store test executables in this root directory
set(NUMFORT_TESTS_DIRECTORY "${CMAKE_BINARY_DIR}/bin/tests")
set(NUMFORT_EXAMPLES_DIRECTORY "${CMAKE_BINARY_DIR}/bin/examples")
set(NUMFORT_BENCHMARKS_DIRECTORY "${CMAKE_BINARY_DIR}/bin/benchmarks")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")


###############################################################################
# Report options and settings
if (BUILD_SHARED_LIBS)
    message(STATUS "Building SHARED libraries")
else(BUILD_SHARED_LIBS)
    message(STATUS "Building STATIC libraries")
endif(BUILD_SHARED_LIBS)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Library install directory: ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
message(STATUS "MOD file install directory: ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")

################################################################################
# Compiler FLAGS

# Notes on the floating point model used by ifort:
# https://software.intel.com/en-us/articles/consistency-of-floating-point-results-using-the-intel-compiler/

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR MINGW)
    set(_FFLAGS "-frecursive")
    set(_FFLAGS_DEBUG "-fcheck=all -O0")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (WIN32)
        set(_FFLAGS "/warn:all,nounused /auto")
        set(_FFLAGS_DEBUG "/check:all")

        if (NOT IFORT_STD_MOD_PROC_NAME)
            set(_FFLAGS "${_FFLAGS} /assume:nostd_mod_proc_name")
        else ()
            set(_FFLAGS "${_FFLAGS} /assume:std_mod_proc_name")
        endif ()
    else ()
        set(_FFLAGS " -warn all,nounused -auto")
        set(_FFLAGS_DEBUG "-check all")

        if (NOT IFORT_STD_MOD_PROC_NAME)
            set(_FFLAGS "${_FFLAGS} -assume nostd_mod_proc_name")
        else ()
            set(_FFLAGS "${_FFLAGS} -assume std_mod_proc_name")
        endif ()
    endif (WIN32)
endif ()

#set_property(DIRECTORY PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_FFLAGS}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${_FFLAGS_DEBUG}")


################################################################################
# Required targets

# External third-party projects
add_subdirectory(external)

# Stand-alone interface definitions for F77 routines and F95 wrappers
add_subdirectory(interfaces)
add_subdirectory(wrappers)

# Numfort library code
add_subdirectory(numfort)


################################################################################
# Optional targets: tests, examples, etc.

if (BUILD_TESTS)
    add_subdirectory(tests)
endif()

if (BUILD_EXAMPLES)
    add_subdirectory(examples)
endif ()
