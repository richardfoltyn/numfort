# orderpack/CMakeLists.txt
# Build configuration for ORDERPACK library

cmake_minimum_required(VERSION 3.0)

project(orderpack Fortran)
set(PROJECT_NAME_PRETTY ORDERPACK)
set(${PROJECT_NAME}_LIBRARY_NAME ${PROJECT_NAME})

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release." FORCE
    )
endif (NOT CMAKE_BUILD_TYPE)

################################################################################
# User-definable options
if (WIN32)
    option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
else ()
    option(BUILD_SHARED_LIBS "Build shared libraries" ON)
endif ()

################################################################################
# Numfort subproject setup

if (${CMAKE_PROJECT_NAME} STREQUAL "numfort")
	set(BUILD_OBJECT_LIBS ON)
	set(${PROJECT_NAME}_EXAMPLES_DIRECTORY
		"${${CMAKE_PROJECT_NAME}_EXAMPLES_DIRECTORY}/external/${PROJECT_NAME}"
	)
    # add to list of object libraries that will be used to create ultimate
    # numfort library
    set_property(GLOBAL APPEND
        PROPERTY ${CMAKE_PROJECT_NAME}_COMPONENT_TARGETS
        ${${PROJECT_NAME}_LIBRARY_NAME}
    )
    message(STATUS "Building ${PROJECT_NAME_PRETTY} as sub-project")
else ()
	set(BUILD_OBJECT_LIBS OFF)
	set(${PROJECT_NAME}_EXAMPLES_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
	message(STATUS "Building ${PROJECT_NAME_PRETTY} as standalone project")
endif ()

# Link demo executables against these libraries
# set(${PROJECT_NAME}_LIBRARIES ${${PROJECT_NAME}_LIBRARY_NAME}_static)

################################################################################
# Set Fortran compiler FLAGS

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # enable optimization that does not break debugging
    set(CMAKE_Fortran_FFLAGS_DEBUG "-Og")
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
    mrgrnk.f90
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

set(_name ${${PROJECT_NAME}_LIBRARY_NAME})
add_library(${_name} OBJECT ${SOURCE_FILES})
if (BUILD_SHARED_LIBS)
    set_target_properties(${_name} PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif ()

# build shared/static library from object library only if standalone build.
# Otherwise, parent project will use object library directly.
if (${PROJECT_NAME}_STANDALONE)
    add_library(${_name} $<TARGET_OBJECTS:${_name}>)
endif ()