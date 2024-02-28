# src/external/CDFLIB90/CMakeLists.txt
# Build configuration for CDFLIB90 library

cmake_minimum_required(VERSION 3.12)

project(cdflib90 Fortran)
set(PROJECT_NAME_PRETTY CDFLIB90)
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

################################################################################
# Default standalone build

set(BUILD_STANDALONE ON)
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
    set(MODULE_DIRECTORY "${EXTERNAL_MODULE_DIRECTORY}/${LIBRARY_NAME}")

    message(STATUS "Building ${PROJECT_NAME_PRETTY} as sub-project")
else ()
	message(STATUS "Building ${PROJECT_NAME_PRETTY} as standalone project")
endif ()

################################################################################
# Target definitions

set(SOURCE_FILES
    biomath_constants_mod.f90
    biomath_mathlib_mod.f90
    biomath_strings_mod.f90
    cdf_beta_mod.f90
    cdf_chisq_mod.f90
    cdf_gamma_mod.f90
    cdf_nc_f_mod.f90
    cdf_neg_binomial_mod.f90
    cdf_poisson_mod.f90
    zero_finder.f90
    biomath_interface_mod.f90
    biomath_sort_mod.f90
    cdf_aux_mod.f90
    cdf_binomial_mod.f90
    cdf_f_mod.f90
    cdf_nc_chisq_mod.f90
    cdf_nc_t_mod.f90
    cdf_normal_mod.f90
    cdf_t_mod.f90
)


# build shared/static library from object library only if standalone build.
# Otherwise, parent project will use object library directly.
if (BUILD_STANDALONE)
    add_library(${LIBRARY_NAME} ${SOURCE_FILES})
    if (BUILD_SHARED_LIBS)
        set_target_properties(${LIBRARY_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    endif ()
else ()
    add_library(${LIBRARY_NAME} OBJECT ${SOURCE_FILES})
endif ()

set_target_properties(${LIBRARY_NAME} PROPERTIES
    Fortran_MODULE_DIRECTORY "${MODULE_DIRECTORY}"
)

target_include_directories(${LIBRARY_NAME} PUBLIC
    "${MODULE_DIRECTORY}"
)

target_compile_options(${LIBRARY_NAME} PRIVATE
    ${EXTERNAL_FFLAGS}
    $<$<CONFIG:Debug>:${EXTERNAL_FFLAGS_DEBUG}>
)

if (USE_LTO)
    set_target_properties(${LIBRARY_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

