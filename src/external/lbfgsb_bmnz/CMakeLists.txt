# lbfgsb_bmnz/CMakeLists.txt
# Build configuration for L-BFGS-B library
cmake_minimum_required(VERSION 3.0)

project(lbfgsb_bmnz Fortran)
set(PROJECT_NAME_PRETTY L-BFGS-B)
set(${PROJECT_NAME}_LIBRARY_NAME ${PROJECT_NAME})

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
option(USE_MKL "Build against Intel MKL" ON)


################################################################################
# Numfort subproject setup

if (${CMAKE_PROJECT_NAME} STREQUAL "numfort")
    set(${PROJECT_NAME}_STANDALONE FALSE)
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
	set(${PROJECT_NAME}_STANDALONE TRUE)
	set(${PROJECT_NAME}_EXAMPLES_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
	message(STATUS "Building ${PROJECT_NAME_PRETTY} as standalone project")
endif ()

################################################################################
# BLAS libraries

if (${PROJECT_NAME}_STANDALONE)
    if (USE_MKL)
        find_package(MKL REQUIRED)
        list(APPEND ${PROJECT_NAME}_LIBRARIES ${MKL_LIBRARIES})
    else ()
        find_package(BLAS REQUIRED)
        list(APPEND ${PROJECT_NAME}_LIBRARIES ${LAPACK_LIBRARIES})
        set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
        set(CMAKE_SHARED_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
    endif ()
endif ()

################################################################################
# Target definitions

set(SOURCE_FILES
    lbfgsb_real64.f
    # lbfgsb.f
    # linpack.f
    # timer.f
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

# if (BUILD_EXAMPLES)
# 	add_subdirectory(examples)
# endif ()
