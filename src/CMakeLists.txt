# src/CMakeLists.txt
# Top-level cmake configuration for numfort

cmake_minimum_required(VERSION 3.0)

# include C in project languages, needed for LAPACK detection
project(numfort C Fortran)
set(${PROJECT_NAME}_LIBRARY_NAME ${CMAKE_PROJECT_NAME})
set(NUMFORT_LIBRARY_NAME ${CMAKE_PROJECT_NAME})

# number of version components
set(NUMFORT_VERSION_COUNT 3)
set(NUMFORT_VERSION_MAJOR 0)
set(NUMFORT_VERSION_MINOR 1)
set(NUMFORT_VERSION_PATCH 0)
set(NUMFORT_VERSION
    "${NUMFORT_VERSION_MAJOR}.${NUMFORT_VERSION_MINOR}.${NUMFORT_VERSION_PATCH}"
)

# list of libraries for targets that link against the final numfort library
# (for examples, benchmarks, etc.)
list(APPEND ${PROJECT_NAME}_LIBRARIES ${${PROJECT_NAME}_LIBRARY_NAME}_static)

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
option(BUILD_TESTS "Build (unit) tests" OFF)
option(BUILD_EXAMPLES "Build examples/demo code" OFF)
option(BUILD_BENCHMARKS "Build benchmarks" OFF)

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

include(CheckFortranPDT)

# check_Fortran_PDT()
# if (Fortran_SUPPORTS_PDT_KIND)
#   add_definitions(-D__NUMFORT_SUPPORTS_PDT_KIND__)
# endif()
#
# if (Fortran_SUPPORTS_PDT_LEN)
#   add_definitions(-D__NUMFORT_SUPPORTS_PDT_LEN__)
# endif()

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
# MKL
if (USE_MKL)
    # find_package(MKL REQUIRED)
    find_package(MKL REQUIRED SEQUENTIAL LAPACK95 BLAS95)
    set(${PROJECT_NAME}_BLA_LIBRARIES ${MKL_LIBRARIES})
else ()
    find_package(LAPACK REQUIRED)
    set(${PROJECT_NAME}_BLA_LIBRARIES ${LAPACK_LIBRARIES})
    set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
    set(CMAKE_SHARED_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
endif ()

list(APPEND ${PROJECT_NAME}_LIBRARIES ${${PROJECT_NAME}_BLA_LIBRARIES})

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
# Output directories

set(${PROJECT_NAME}_MODDIR "${CMAKE_BINARY_DIR}/include")
set(CMAKE_Fortran_MODULE_DIRECTORY "${${PROJECT_NAME}_MODDIR}")

# set include directory to directory where MOD files are stored
list(APPEND ${PROJECT_NAME}_INCLUDE_DIRS
    "${${PROJECT_NAME}_MODDIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/numfort/include"
)

# store test executables in this root directory
set(${PROJECT_NAME}_TESTS_DIRECTORY "${CMAKE_BINARY_DIR}/bin/tests")
set(${PROJECT_NAME}_EXAMPLES_DIRECTORY "${CMAKE_BINARY_DIR}/bin/examples")
set(${PROJECT_NAME}_BENCHMARKS_DIRECTORY "${CMAKE_BINARY_DIR}/bin/benchmarks")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")



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

if (BUILD_TESTS)
    message(STATUS "Building unit tests")
endif ()

if (BUILD_BENCHMARKS)
    message(STATUS "Building benchmarks")
endif ()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Library install directory: ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
message(STATUS "MOD file install directory: ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")

################################################################################
# Compiler FLAGS

# Notes on the floating point model used by ifort:
# https://software.intel.com/en-us/articles/consistency-of-floating-point-results-using-the-intel-compiler/

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR MINGW)
    # explicitly add -cpp to run C preprocessor if required
    # By default set -Wall, but disable warnings when comparing reals
    # and -Wsurpsiring, as this generates numerous warnings about finalizers.
    set(_FFLAGS "${_FFLAGS} -Wall -Wno-compare-reals -Wno-surprising")
    set(_FFLAGS_DEBUG "-fcheck=all -Og")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (WIN32)
#         set(_FFLAGS_FP "/fp:precise /fp:source")
        set(_FFLAGS "/warn:all /auto ${_FFLAGS_FP}")
        set(_FFLAGS_DEBUG "/check:all")

        if (NOT IFORT_STD_MOD_PROC_NAME)
            set(_FFLAGS "${_FFLAGS} /assume:nostd_mod_proc_name")
        else ()
            set(_FFLAGS "${_FFLAGS} /assume:std_mod_proc_name")
        endif ()
    else ()
#         set(_FFLAGS_FP "-fp-model precise -fp-model source")
        set(_FFLAGS " -warn all -auto ${_FFLAGS_FP}")
        set(_FFLAGS_DEBUG "-check all")

        if (NOT IFORT_STD_MOD_PROC_NAME)
            set(_FFLAGS "${_FFLAGS} -assume nostd_mod_proc_name")
        else ()
            set(_FFLAGS "${_FFLAGS} -assume std_mod_proc_name")
        endif ()
    endif (WIN32)
endif ()

set_property(DIRECTORY PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_FFLAGS}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${_FFLAGS_DEBUG}")

################################################################################
# External projects

# Add these subdirectories before setting any compiler numfort-specifc
# compiler flags.

add_subdirectory(external)
add_subdirectory(interfaces)

################################################################################
# Targets
# Create single numfort library by merging all object libraries

include_directories(${${PROJECT_NAME}_INCLUDE_DIRS})

add_subdirectory(numfort)

# Update fortran version module with current version information
configure_file(numfort/numfort.f90.in
    ${CMAKE_BINARY_DIR}/numfort.f90
)

# Add a few general source files, makes it easier in case generator does
# not support targets w/o source files.
set(SOURCE_FILES
    ${CMAKE_BINARY_DIR}/numfort.f90
)

# get list of targets the ultimate library depends on. This list is populated
# in the sub-components' directories
get_property(_TARGETS GLOBAL PROPERTY ${CMAKE_PROJECT_NAME}_COMPONENT_TARGETS)

set(COMPONENT_OBJECTS)
foreach (_comp IN LISTS _TARGETS)
    list(APPEND COMPONENT_OBJECTS $<TARGET_OBJECTS:${_comp}>)
endforeach ()

add_library(${NUMFORT_LIBRARY_NAME} ${SOURCE_FILES} ${COMPONENT_OBJECTS})
if (BUILD_SHARED_LIBS)
    target_link_libraries(${NUMFORT_LIBRARY_NAME} ${${PROJECT_NAME}_BLA_LIBRARIES})
endif ()

set_target_properties(${NUMFORT_LIBRARY_NAME} PROPERTIES
    Fortran_MODULE_DIRECTORY "${${PROJECT_NAME}_MODDIR}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

set_target_properties(${NUMFORT_LIBRARY_NAME} PROPERTIES
    VERSION ${NUMFORT_VERSION}
)

# for benchmarks, examples, etc
add_library(${NUMFORT_LIBRARY_NAME}_static STATIC ${SOURCE_FILES} ${COMPONENT_OBJECTS})



################################################################################
# INSTALLATION

include(CMakePackageConfigHelpers)

# Allow for installing ifort and gfortran version of libraries side-by-side by
# appending compiler suffix
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(Fortran_COMPILER_SUFFIX "gf")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(Fortran_COMPILER_SUFFIX "ifort")
else ()
    message(FATAL_ERROR "Unsupported Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif ()

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

# Path there cmake config files will be installed
set(NUMFORT_INSTALL_CONFIG_PATH
    ${CMAKE_INSTALL_LIBDIR}/${NUMFORT_INSTALL_DIR}/cmake
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

configure_package_config_file("${NUMFORT_INSTALL_SOURCE_DIR}/config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${NUMFORT_CONFIG_FILE}"
    INSTALL_DESTINATION "${NUMFORT_INSTALL_CONFIG_PATH}"
    PATH_VARS NUMFORT_INSTALL_CONFIG_PATH
)

# Note: INCLUDES DESTINATION is equivalent to adding
# INTERFACE_INCLUDE_DIRECTORIES property on every exported target
install(TARGETS ${NUMFORT_LIBRARY_NAME}
    EXPORT ${NUMFORT_LIBRARY_NAME}_targets
    RUNTIME DESTINATION "${CMAKE_INSTALL_LIBDIR}/${NUMFORT_INSTALL_DIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}/${NUMFORT_INSTALL_DIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}/${NUMFORT_INSTALL_DIR}"
    INCLUDES DESTINATION
        "${CMAKE_INSTALL_INCLUDEDIR}/${NUMFORT_INSTALL_DIR}"
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
install(DIRECTORY "${${PROJECT_NAME}_MODDIR}/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${NUMFORT_INSTALL_DIR}"
    FILES_MATCHING PATTERN "numfort*.mod"
)

install(FILES
        "${${PROJECT_NAME}_MODDIR}/blas_interfaces.mod"
        "${${PROJECT_NAME}_MODDIR}/lapack_interfaces.mod"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${NUMFORT_INSTALL_DIR}"
)




################################################################################
# Tests, examples, etc.

if (BUILD_TESTS)
    add_subdirectory(tests)
endif()

if (BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

if (BUILD_EXAMPLES)
    add_subdirectory(examples)
endif ()
