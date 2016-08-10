# src/CMakeLists.txt
# Top-level cmake configuration for numfort

cmake_minimum_required(VERSION 3.0)

# include C in project languages, needed for LAPACK detection
project(numfort C Fortran)
set(LIBRARY_NAME ${CMAKE_PROJECT_NAME} CACHE STRING "File name of main library")

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set(CMAKE_VERBOSE_MAKEFILE ON)

################################################################################
# User-definable options
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(USE_MKL "Build against Intel MKL" ON)
option(BUILD_TESTS "Build unit tests" ON)
option(BUILD_BENCHMARKS "Build benchmarks" OFF)

################################################################################
# Additional CMAKE modules
# include FindMKL module taken from ARMA project

list(APPEND CMAKE_MODULE_PATH
"${CMAKE_CURRENT_SOURCE_DIR}/submodules/cmake-modules/Modules")
# set(FC gfortran)
include(FindLAPACK)
include(FindMKL)

# if (LAPACK_FOUND)
	# message("LAPACK_LINKER_FLAGS: ${LAPACK_LINKER_FLAGS}")
	# message("LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")
# endif (LAPACK_FOUND)

# if (MKL_FOUND)
	# message("MLK_LIBRARIES: ${MKL_LIBRARIES}")
# endif (MKL_FOUND)

if (USE_MKL AND (NOT MKL_FOUND))
	message(FATAL_ERROR "MKL build requested but MKL not found")
endif()

if ((NOT LAPACK_FOUND) AND (NOT MKL_FOUND))
	message(FATAL_ERROR "Neither MKL nor other LAPACK library found.")
else()
	if (MKL_FOUND)
		set(NUMFORT_LAPACK_LIBRARIES ${MKL_LIBRARIES})
	elseif (LAPACK_FOUND)
		set(NUMFORT_LAPACK_LIBRARIES ${MKL_LIBRARIES})
		set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
		set(CMAKE_SHARED_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS})
	endif()
endif()

################################################################################
# Testing

if(BUILD_TESTS)
    enable_testing()
    find_package(fcore)
endif()

################################################################################
# Detect target architecture
include(FindTargetArch)
find_target_arch()

# set installation directories (relative to PREFIX)
if (NOT INSTALL_INCLUDE_ROOT OR NOT INSTALL_LIBRARY_ROOT)
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        set(INCLUDE_SUBDIR "gfortran")
    elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
        if (TARGET_ARCH_BITS STREQUAL "64")
            set(INCLUDE_SUBDIR "intel64")
        else()
            set(INCLUDE_SUBDIR "intel")
        endif()
    endif ()
    set(INSTALL_INCLUDE_ROOT "include/${INCLUDE_SUBDIR}"
        CACHE STRING "Relative installation directory for MOD files" FORCE
    )
    set(INSTALL_LIBRARY_ROOT "lib/${INCLUDE_SUBDIR}"
        CACHE STRING "Relative installation directory for libraries" FORCE
    )
endif()

################################################################################
# Set Fortran compiler FLAGS

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set (${PROJECT_NAME}_Fortran_FLAGS "-Wall -Wextra -Wimplicit-interface -Wimplicit-procedure -std=f2008")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (WIN32)
        set (${PROJECT_NAME}_Fortran_FLAGS "/warn:all /standard-semantics /stand:f08")
		if (USE_MKL)
			set (${PROJECT_NAME}_Fortran_FLAGS "${${PROJECT_NAME}_Fortran_FLAGS} /Qmkl")
		endif (USE_MKL)
    elseif (WIN32)
        set (${PROJECT_NAME}_Fortran_FLAGS "-warn all -standard-semantics -std08")
    endif (WIN32)
else ()
    message(FATAL_ERROR "Unsupported Fortran compiler: ${Fortran_COMPILER_NAME}")
endif ()

# append additional FLAGS to default cmake FLAGS
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${${PROJECT_NAME}_Fortran_FLAGS}")
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${${PROJECT_NAME}_Fortran_FLAGS_RELEASE}")
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${${PROJECT_NAME}_Fortran_FLAGS_DEBUG}")

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
message(STATUS "Fortran compiler flags: ${CMAKE_Fortran_FLAGS}")
message(STATUS "Additional Fortran compiler flags for ${CMAKE_BUILD_TYPE}: ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}")

message(STATUS "Library install directory: ${CMAKE_INSTALL_PREFIX}/${INSTALL_LIBRARY_ROOT}")
message(STATUS "MOD file install directory: ${CMAKE_INSTALL_PREFIX}/${INSTALL_INCLUDE_ROOT}")

################################################################################
# Output directories

set(NUMFORT_MODDIR "${CMAKE_BINARY_DIR}/include")

# Set output directory for static and shared libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_Fortran_MODULE_DIRECTORY "${NUMFORT_MODDIR}")

# store test executables in this root directory
set(${PROJECT_NAME}_TESTS_DIRECTORY "${CMAKE_BINARY_DIR}/bin/tests")

add_subdirectory(numfort)

if (BUILD_TESTS)
    add_subdirectory(tests)
endif()

if (BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()
