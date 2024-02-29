
cmake_minimum_required(VERSION 3.12)

project(numfort VERSION 0.1.0 LANGUAGES C Fortran)
set(NUMFORT_LIBRARY_NAME ${CMAKE_PROJECT_NAME})

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE
    )
endif ()

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
option(USE_OPENMP "Use some OpenMP features" ON)
option(USE_LTO "Use link-time optimization" OFF)


################################################################################
# Additional CMAKE modules

list(APPEND CMAKE_MODULE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/submodules/cmake-modules/Modules")

################################################################################
# Testing

if(BUILD_TESTS)
    enable_testing()
endif()

################################################################################
# LTO / IPO

if (USE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT SUPPORTS_LTO)
    if (NOT SUPPORTS_LTO)
        set(USE_LTO OFF)
        message(WARNING "LTO not supported")
    endif()
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
# Report options and settings
if (BUILD_SHARED_LIBS)
    message(STATUS "Building SHARED libraries")
else(BUILD_SHARED_LIBS)
    message(STATUS "Building STATIC libraries")
endif(BUILD_SHARED_LIBS)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
if (USE_LTO)
    message(STATUS "Link-time optimization ENABLED")
endif()

################################################################################
# Output directories

set(NUMFORT_TESTS_DIRECTORY "${CMAKE_BINARY_DIR}/bin/tests")
set(NUMFORT_EXAMPLES_DIRECTORY "${CMAKE_BINARY_DIR}/bin/examples")
set(NUMFORT_BENCHMARKS_DIRECTORY "${CMAKE_BINARY_DIR}/bin/benchmarks")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/include/numfort")
set(EXTERNAL_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/include/external")

# Output directory for additional static library used for linking tests
# and examples irrespective of whether main build is shared/static.
# This is not meant to be installed anywhere, so keep it separate.
set(NUMFORT_STATIC_DIRECTORY "${CMAKE_BINARY_DIR}/lib-static")


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

# Allow for installing ifort and gfortran version of libraries side-by-side by
# appending compiler suffix
string(TOLOWER ${CMAKE_Fortran_COMPILER_ID} Fortran_COMPILER_ID)

# Find major compiler version
string(REGEX REPLACE "(^[0-9]+).*" "\\1" Fortran_COMPILER_VERSION_MAJOR
    "${CMAKE_Fortran_COMPILER_VERSION}")

set(Fortran_COMPILER_SUFFIX "${Fortran_COMPILER_ID}-${Fortran_COMPILER_VERSION_MAJOR}")

# Installation directory within <prefix>/lib, <prefix>/include, etc.
# If ont defined by -DINSTALL_DIR_NAME by the user, create
# version- and compiler-specific directory name.
if (NOT DEFINED INSTALL_DIR_NAME)
    set(_version ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
    set(INSTALL_DIR_NAME
        "${NUMFORT_LIBRARY_NAME}-${_version}-${Fortran_COMPILER_SUFFIX}"
        CACHE STRING "Install directory name"
        FORCE
    )
endif()


################################################################################
# Universal compiler FLAGS for all targets

unset(_FFLAGS)
unset(_FFLAGS_DEBUG)

# Notes on the floating point model used by ifort:
# https://software.intel.com/en-us/articles/consistency-of-floating-point-results-using-the-intel-compiler/

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR MINGW)
    set(_FFLAGS "-frecursive")
    set(_FFLAGS_DEBUG "-fcheck=all -O0")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel(LLVM)?")
    if (WIN32)
        set(_FFLAGS "/auto /assume:recursion")
        set(_FFLAGS_DEBUG "/check:all")
        # Disable seriously annoying warning that ifort is deprecated
        if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
            set(_FFLAGS "${_FFLAGS} /Qdiag-disable:10448")
        endif ()
    else ()
        set(_FFLAGS "-auto -assume recursion")
        set(_FFLAGS_DEBUG "-check all")

        # Disable seriously annoying warning that ifort is deprecated
        if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
            set(_FFLAGS "${_FFLAGS} -diag-disable=10448")
        endif ()
    endif (WIN32)
endif ()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_FFLAGS}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${_FFLAGS_DEBUG}")


################################################################################
# Process subdirectories

add_subdirectory(src)

if (BUILD_TESTS)
    add_subdirectory(tests)
endif()

if (BUILD_EXAMPLES)
    add_subdirectory(examples)
endif ()


################################################################################
# Installation

# Path there cmake config files will be installed
set(NUMFORT_INSTALL_CONFIG_PATH
    ${CMAKE_INSTALL_LIBDIR}/${INSTALL_DIR_NAME}/cmake
)

set(NUMFORT_INSTALL_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/install")
set(NUMFORT_INSTALL_NAMESPACE ${NUMFORT_LIBRARY_NAME})

set(NUMFORT_CONFIG_FILE ${NUMFORT_LIBRARY_NAME}-config.cmake)
set(NUMFORT_CONFIG_VERSION_FILE ${NUMFORT_LIBRARY_NAME}-config-version.cmake)
set(NUMFORT_TARGETS_FILE ${NUMFORT_LIBRARY_NAME}Targets.cmake)

configure_file("${NUMFORT_INSTALL_SOURCE_DIR}/config-version.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${NUMFORT_CONFIG_VERSION_FILE}"
    @ONLY
)

include(CMakePackageConfigHelpers)

configure_package_config_file("${NUMFORT_INSTALL_SOURCE_DIR}/config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${NUMFORT_CONFIG_FILE}"
    INSTALL_DESTINATION "${NUMFORT_INSTALL_CONFIG_PATH}"
    PATH_VARS NUMFORT_INSTALL_CONFIG_PATH
)

# Note: INCLUDES DESTINATION is equivalent to adding
# INTERFACE_INCLUDE_DIRECTORIES property on every exported target
install(TARGETS ${NUMFORT_LIBRARY_NAME}
    EXPORT ${NUMFORT_LIBRARY_NAME}_targets
    RUNTIME DESTINATION "${CMAKE_INSTALL_LIBDIR}/${INSTALL_DIR_NAME}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}/${INSTALL_DIR_NAME}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}/${INSTALL_DIR_NAME}"
    INCLUDES DESTINATION
        "${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_DIR_NAME}"
)

export(EXPORT ${NUMFORT_LIBRARY_NAME}_targets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/${NUMFORT_TARGETS_FILE}"
    NAMESPACE ${NUMFORT_INSTALL_NAMESPACE}
)

install(EXPORT ${NUMFORT_LIBRARY_NAME}_targets
    FILE "${NUMFORT_TARGETS_FILE}"
    NAMESPACE ${NUMFORT_INSTALL_NAMESPACE}::
    DESTINATION "${NUMFORT_INSTALL_CONFIG_PATH}"
)


# Install CMake config files that were generated above
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${NUMFORT_CONFIG_FILE}"
    "${CMAKE_CURRENT_BINARY_DIR}/${NUMFORT_CONFIG_VERSION_FILE}"
    "${NUMFORT_INSTALL_SOURCE_DIR}/compiler_test.f90"
    DESTINATION "${NUMFORT_INSTALL_CONFIG_PATH}"
    COMPONENT devel
)

# Install only numfort's MOD files, ignore any modules potentially created
# by external projects
# TODO: Do not install MOD files that are used only internally.
install(DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${INSTALL_DIR_NAME}"
    FILES_MATCHING PATTERN "numfort*.mod"
    PATTERN "external" EXCLUDE
)