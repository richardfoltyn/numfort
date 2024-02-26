# src/external/fitpack/CMakeLists.txt
# CMake configuration for FITPACK library

cmake_minimum_required(VERSION 3.10)

project(fitpack Fortran)
set(PROJECT_NAME_PRETTY FITPACK)
set(LIBRARY_NAME ${PROJECT_NAME})

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
    set(MODULE_DIRECTORY ${NUMFORT_MODULE_DIRECTORY}/external)

    message(STATUS "Building ${PROJECT_NAME_PRETTY} as sub-project")
else ()
	message(STATUS "Building ${PROJECT_NAME_PRETTY} as standalone project")
endif ()

################################################################################
# Additional compiler flags for Fortran 2008 source files

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR MINGW)
    # explicitly add -cpp to run C preprocessor if required
    set(GNU_WARN "-Wall -Wextra -Warray-temporaries -Wrealloc-lhs -pedantic -Wuninitialized")
    set(F18_FLAGS "${GNU_WARN} -std=f2018")

    # Disable "surprising" warnings for older GFORTRAN version as this
    # produces endless warnings about array finalizers
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 7)
        set(F18_FLAGS "${F18_FLAGS} -Wno-surprising")
    endif()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel(LLVM)?")
    if (WIN32)
        set(F18_FLAGS "/warn:all /standard-semantics /stand:f18")
        # NOTE: When setting -standard-semantics, resets any previously
        # given -assume std_mod_proc_name, so we need to append it again!
        if (NOT IFORT_STD_MOD_PROC_NAME)
            set(F18_FLAGS "${F18_FLAGS} /assume:nostd_mod_proc_name")
        else ()
            set(F18_FLAGS "${F18_FLAGS} /assume:std_mod_proc_name")
        endif ()
    else ()
        set(F18_FLAGS "-warn all -standard-semantics -std18")
        # NOTE: When setting -standard-semantics, resets any previously
        # given -assume std_mod_proc_name, so we need to append it again!
        if (NOT IFORT_STD_MOD_PROC_NAME)
            set(F18_FLAGS "${F18_FLAGS} -assume nostd_mod_proc_name")
        else ()
            set(F18_FLAGS "${F18_FLAGS} -assume std_mod_proc_name")
        endif ()
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

if (F18_FLAGS)
    set_source_files_properties(${F18_SOURCES} PROPERTIES
        COMPILE_FLAGS "${F18_FLAGS}"
    )
endif ()

string(TOUPPER ${CMAKE_BUILD_TYPE} _build)
if (F18_FLAGS_${_build})
    set_source_files_properties(${F18_SOURCES} PROPERTIES
        COMPILE_FLAGS "${F18_FLAGS_${_build}}"
    )
endif ()


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

# Building original example files shipped with FITPACK does not work ATM
#if (BUILD_EXAMPLES)
#    add_subdirectory(examples)
#endif ()
