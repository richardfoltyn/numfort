# src/CMakeLists.txt
# Top-level cmake configuration for numfort

cmake_minimum_required(VERSION 3.0)

# include C in project languages, needed for LAPACK detection
project(numfort C Fortran)
set(${PROJECT_NAME}_LIBRARY_NAME ${CMAKE_PROJECT_NAME})

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
# include FindMKL module taken from ARMA project

list(APPEND CMAKE_MODULE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/submodules/cmake/Modules")

################################################################################
# Testing

if(BUILD_TESTS)
    enable_testing()
endif()

################################################################################
# MKL
if (USE_MKL)
    find_package(MKL REQUIRED)
    set(${PROJECT_NAME}_BLA_LIBRARIES ${MKL_LIBRARIES})
else ()
    find_package(LAPACK REQUIRED)
    set(${PROJECT_NAME}_BLA_LIBRARIES ${LAPACK_LIBRARIES})
    set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
    set(CMAKE_SHARED_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
endif ()

list(APPEND ${PROJECT_NAME}_LIBRARIES ${${PROJECT_NAME}_BLA_LIBRARIES})

################################################################################
# Compiler FLAGS

unset(_FFLAGS)
unset(_FFLAGS_DEBUG)
# global _FFLAGS will be added to CMAKE_Fortran_FLAGS and thus effect
# compilation of F77 external code
unset(_FFLAGS_GLOBAL)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR MINGW)
    # explicitly add -cpp to run C preprocessor if required
    set(GNU_WARN "-Wall -Wextra -Wimplicit-interface -Wimplicit-procedure \
      -Warray-temporaries -Wrealloc-lhs -pedantic")
    set(_FFLAGS "${GNU_WARN} -std=f2008")
    set(_FFLAGS_DEBUG "-fcheck=all -Og")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (WIN32)
        set (_FFLAGS "/warn:all /standard-semantics /stand:f08")
		set(_FFLAGS_DEBUG "/check:all")

		if (NOT IFORT_STD_MOD_PROC_NAME)
            set(_FFLAGS_GLOBAL "/assume:nostd_mod_proc_name")
        endif ()
    else ()
        set (_FFLAGS "-warn all -standard-semantics -std08")
        set(_FFLAGS_DEBUG "-check all")

		if (NOT IFORT_STD_MOD_PROC_NAME)
            set(_FFLAGS_GLOBAL "-assume nostd_mod_proc_name")
        endif ()
    endif (WIN32)
endif ()

# add standard semantics to compiler flags even for F77 files, in case
# these reference any submodules
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_FFLAGS_GLOBAL}")

message(STATUS "Default compiler flags: ${CMAKE_Fortran_FLAGS}")
# Make sure build type is upper case, otherwise the wrong flags will be shown
string(TOUPPER "${CMAKE_BUILD_TYPE}" _CONFIG)
if (CMAKE_Fortran_FLAGS_${_CONFIG})
    message(STATUS "Default compiler flags for ${CMAKE_BUILD_TYPE}: ${CMAKE_Fortran_FLAGS_${_CONFIG}}")
endif ()
if (_FFLAGS)
    message(STATUS "Add. compiler flags: ${_FFLAGS}")
endif ()
if (_FFLAGS_${_CONFIG})
    message(STATUS "Add. compiler flags for ${CMAKE_BUILD_TYPE}: ${_FFLAGS_${_CONFIG}}")
endif ()

set(${CMAKE_PROJECT_NAME}_FFLAGS "${_FFLAGS}")
set(${CMAKE_PROJECT_NAME}_FFLAGS_DEBUG "${_FFLAGS_DEBUG}")

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
message(STATUS "Library install directory: ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
message(STATUS "MOD file install directory: ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")

################################################################################
# Output directories

set(${PROJECT_NAME}_MODDIR "${CMAKE_BINARY_DIR}/include")
set(CMAKE_Fortran_MODULE_DIRECTORY "${${PROJECT_NAME}_MODDIR}")

# set include directory to directory where MOD files are stored
list(APPEND ${PROJECT_NAME}_INCLUDE_DIRS "${${PROJECT_NAME}_MODDIR}")

# store test executables in this root directory
set(${PROJECT_NAME}_TESTS_DIRECTORY "${CMAKE_BINARY_DIR}/bin/tests")
set(${PROJECT_NAME}_EXAMPLES_DIRECTORY "${CMAKE_BINARY_DIR}/bin/examples")
set(${PROJECT_NAME}_BENCHMARKS_DIRECTORY "${CMAKE_BINARY_DIR}/bin/benchmarks")

################################################################################
# External projects

add_subdirectory(external)
add_subdirectory(interfaces)

################################################################################
# Detect features supported by Fortran compiler
# Do this only after adding external projects to avoid any possible
# re-definition

include(CheckFortranPDT)
include(CheckFortranKinds)

check_Fortran_PDT()
check_Fortran_kinds()

if (Fortran_SUPPORTS_PDT_KIND)
    add_definitions(-D__SUPPORTS_PDT_KIND__)
endif()

if (Fortran_SUPPORTS_PDT_LEN)
    add_definitions(-D__SUPPORTS_PDT_LEN__)
endif()

add_definitions(-D__DEFAULT_INT_KIND__=${Fortran_DEFAULT_INT_KIND}
	-D__DEFAULT_INT_SIZE__=${Fortran_DEFAULT_INT_SIZE}
	-D__DEFAULT_REAL_KIND__=${Fortran_DEFAULT_REAL_KIND}
	-D__DEFAULT_REAL_SIZE__=${Fortran_DEFAULT_REAL_SIZE}
)

################################################################################
# Targets
# Create single numfort library by merging all object libraries

add_subdirectory(numfort)

# Update fortran version module with current version information
configure_file(numfort/version.f90.in
	${CMAKE_BINARY_DIR}/version.f90
)

# Add a few general source files, makes it easier in case generator does
# not support targets w/o source files.
set(SOURCE_FILES
	${CMAKE_BINARY_DIR}/version.f90
	numfort/numfort.f90
)

include_directories("${${PROJECT_NAME}_MODDIR}")

# get list of targets the ultimate library depends on. This list is populated
# in the sub-components' directories
get_property(_TARGETS GLOBAL PROPERTY ${CMAKE_PROJECT_NAME}_COMPONENT_TARGETS)

set(COMPONENT_OBJECTS)
foreach (_comp IN LISTS _TARGETS)
    list(APPEND COMPONENT_OBJECTS $<TARGET_OBJECTS:${_comp}>)
endforeach ()

set(LIBRARY_NAME ${${PROJECT_NAME}_LIBRARY_NAME})
add_library(${LIBRARY_NAME} ${SOURCE_FILES} ${COMPONENT_OBJECTS})
target_link_libraries(${LIBRARY_NAME} ${${PROJECT_NAME}_BLA_LIBRARIES})
set_target_properties(${LIBRARY_NAME} PROPERTIES
	Fortran_MODULE_DIRECTORY "${${PROJECT_NAME}_MODDIR}"
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
	LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
	ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

# for benchmarks, examples, etc
add_library(${LIBRARY_NAME}_static STATIC ${SOURCE_FILES} ${COMPONENT_OBJECTS})
add_dependencies(${LIBRARY_NAME}_static ${_TARGETS})

install(TARGETS ${LIBRARY_NAME}
    RUNTIME DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

# Install only numfort's MOD files, ignore any modules potentially created
# by external projects
install(DIRECTORY ${${PROJECT_NAME}_MODDIR}/
	DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/numfort"
	FILES_MATCHING PATTERN "numfort*.mod"
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
