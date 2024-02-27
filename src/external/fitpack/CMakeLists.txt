# src/external/fitpack/CMakeLists.txt
# CMake configuration for FITPACK library

cmake_minimum_required(VERSION 3.12)

project(fitpack Fortran)
set(PROJECT_NAME_PRETTY FITPACK)
set(LIBRARY_NAME ${PROJECT_NAME})

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif ()

################################################################################
# User-definable options
if (WIN32)
    option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
else ()
    option(BUILD_SHARED_LIBS "Build shared libraries" ON)
endif ()
option(BUILD_EXAMPLES "Build examples/demo code" OFF)

################################################################################
# Default standalone build

set(BUILD_STANDALONE ON)
set(EXAMPLES_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/include")

################################################################################
# Numfort subproject setup

if (${CMAKE_PROJECT_NAME} STREQUAL "numfort")
    set(BUILD_STANDALONE OFF)
	set(EXAMPLES_DIRECTORY
		"${NUMFORT_EXAMPLES_DIRECTORY}/external/${PROJECT_NAME}"
	)
    # add to list of object libraries that will be used to create ultimate
    # numfort library
    set_property(GLOBAL APPEND
        PROPERTY ${CMAKE_PROJECT_NAME}_COMPONENT_TARGETS
        ${LIBRARY_NAME}
    )

    # Store module files in global module directory
    set(MODULE_DIRECTORY "${EXTERNAL_MODULE_DIRECTORY}/${LIBRARY_NAME}")

    message(STATUS "Building ${PROJECT_NAME_PRETTY} as sub-project")
else ()
	message(STATUS "Building ${PROJECT_NAME_PRETTY} as standalone project")
endif ()

################################################################################
# Additional compiler flags for Fortran 2008 source files

unset(F18_FLAGS)
unset(F18_FLAGS_Debug)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR MINGW)
    set(GNU_WARN -Wall -Wextra -Warray-temporaries -Wrealloc-lhs -pedantic -Wuninitialized)
    set(F18_FLAGS ${GNU_WARN} -std=f2018)

elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel(LLVM)?")
    if (WIN32)
        set(F18_FLAGS /warn:all /stand:f18)
    else ()
        set(F18_FLAGS -warn all -std18)
    endif (WIN32)
endif ()

################################################################################
# Target definitions

set(SOURCE_FILES
	fitpack_real64.f
)

# additional files replacing original FITPACK implementation, written in F90
# These FLAGS are appended to any other previously defined FLAGS.
list(APPEND F18_SOURCES constraints_tree.f90)
list(APPEND SOURCE_FILES ${F18_SOURCES})

set(_name ${LIBRARY_NAME})

# build shared/static library from object library only if standalone build.
# Otherwise, parent project will use object library directly.
if (BUILD_STANDALONE)
    add_library(${_name} ${SOURCE_FILES})
    if (BUILD_SHARED_LIBS)
        set_target_properties(${_name} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    endif ()
else ()
    add_library(${_name} OBJECT ${SOURCE_FILES})
endif ()

set_target_properties(${_name} PROPERTIES
    Fortran_MODULE_DIRECTORY ${MODULE_DIRECTORY}
)

target_compile_options(${_name} PRIVATE
    ${EXTERNAL_FFLAGS}
    $<$<CONFIG:Debug>:${EXTERNAL_FFLAGS_Debug}>
)

target_include_directories(${_name} PUBLIC
    "${MODULE_DIRECTORY}"
)

# These do not support generator expressions
set_property(SOURCE ${F18_SOURCES} APPEND PROPERTY COMPILE_OPTIONS ${F18_FLAGS})
set_property(SOURCE ${F18_SOURCES} APPEND PROPERTY COMPILE_OPTIONS 
    ${F18_FLAGS_${CMAKE_BUILD_TYPE}}
)


# Building original example files shipped with FITPACK does not work ATM
#if (BUILD_EXAMPLES)
#    add_subdirectory(examples)
#endif ()
