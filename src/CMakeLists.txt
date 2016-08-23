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

################################################################################
# Additional CMAKE modules
# include FindMKL module taken from ARMA project

list(APPEND CMAKE_MODULE_PATH
  "${CMAKE_CURRENT_SOURCE_DIR}/submodules/cmake/Modules")

################################################################################
# Testing

if(BUILD_TESTS)
    enable_testing()
    find_package(FCore)
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
message(STATUS "Library install directory: ${CMAKE_INSTALL_PREFIX}/${INSTALL_LIBRARY_PREFIX}")
message(STATUS "MOD file install directory: ${CMAKE_INSTALL_PREFIX}/${INSTALL_INCLUDE_PREFIX}")

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

add_subdirectory(external)
add_subdirectory(numfort)

################################################################################
# Targets
# Create single numfort library by merging all object libraries

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
    RUNTIME DESTINATION "${INSTALL_LIBRARY_PREFIX}"
    LIBRARY DESTINATION "${INSTALL_LIBRARY_PREFIX}"
    ARCHIVE DESTINATION "${INSTALL_LIBRARY_PREFIX}"
)

# Install only numfort's MOD files, ignore any modules potentially created
# by external projects
install(DIRECTORY ${${PROJECT_NAME}_MODDIR}/
	DESTINATION "${INSTALL_INCLUDE_PREFIX}"
	FILES_MATCHING PATTERN "numfort*.mod"
)
