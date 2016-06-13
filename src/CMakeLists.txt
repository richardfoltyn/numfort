# Top-level cmake configuration
# src/CMakeLists.txt

cmake_minimum_required(VERSION 2.8)

# include C in project languages, needed for LAPACK detection
project(numfort C Fortran)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set(CMAKE_VERBOSE_MAKEFILE ON)

# include FindMKL module taken from ARMA project
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")
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

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
	set (CMAKE_Fortran_FLAGS "-Wall -Wextra -Wimplicit-interface -Wimplicit-procedure -std=f2008")
	# append any FFLAGS (in particular for machine-specific optimization)
	set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} $ENV{FFLAGS}")
	# set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
	#set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
	if (WIN32)
		set (CMAKE_Fortran_FLAGS "/nologo /standard-semantics /stand:f08 /warn:all $ENV{FFLAGS} ")
		if (MKL_FOUND)
			# with ifort on Windows just use the appropriate option to link against MKL
			set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /Qmkl ")
		endif (MKL_FOUND)
		# set (CMAKE_Fortran_FLAGS_RELEASE "/O2")
		# set (CMAKE_Fortran_FLAGS_DEBUG " /debug:full /Od /check:all")
		# set (USE_MKL ON CACHE TYPE BOOL "Link against Intel MKL")
	elseif (WIN32)
		set (CMAKE_Fortran_FLAGS "-standard-semantics -std08 $ENV{FFLAGS}")
		# set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
		# set (CMAKE_Fortran_FLAGS_DEBUG "-g -O0")
	endif (WIN32)
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

set(NUMFORT_MODDIR "${CMAKE_BINARY_DIR}/include")

# Set output directory for static and shared libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

set(CMAKE_Fortran_MODULE_DIRECTORY "${NUMFORT_MODDIR}")

add_subdirectory(numfort)
add_subdirectory(tests)
add_subdirectory(benchmarks)
