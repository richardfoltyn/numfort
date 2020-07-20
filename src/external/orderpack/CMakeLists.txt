# src/external/orderpack/CMakeLists.txt
# Build configuration for ORDERPACK library

cmake_minimum_required(VERSION 3.10)

project(orderpack Fortran)
set(PROJECT_NAME_PRETTY ORDERPACK)
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

################################################################################
# Default standalone build

set(BUILD_STANDALONE ON)
set(EXAMPLES_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/include")

################################################################################
# Numfort subproject setup

if (${CMAKE_PROJECT_NAME} STREQUAL "numfort")
    set(BUILD_STANDALONE OFF)
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
# Target definitions

set(SOURCE_FILES
    ctrper.f90
    fndnth.f90
    indmed.f90
    indnth.f90
    inspar.f90
    inssor.f90
    median.f90
    mrgref.f90
    mrgrnk.F90
    mulcnt.f90
    rapknr.f90
    refpar.f90
    refsor.f90
    rinpar.f90
    rnkpar.f90
    uniinv.f90
    unipar.f90
    unirnk.f90
    unista.f90
    valmed.f90
    valnth.f90
)

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

target_include_directories(${_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
if (USE_LTO)
    set_target_properties(${_name} PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
