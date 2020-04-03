
################################################################################
# Add compiler flags specific to external projects

unset(_FFLAGS)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR MINGW)
    # there shouldn't be many issues with LHS reallocation or temporary arrays
    # in old F77 code, but still add this.
    set(_FFLAGS "-Wall -Wno-uninitialized -Wno-conversion -Wno-unused -Wno-tabs -Wno-unused-dummy-argument")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (WIN32)
        set(_FFLAGS "/warn:all,nounused,nousage")
        #set(_FLAGS_RELEASE "/Qipo")
    else ()
        set(_FFLAGS "-warn all,nounused,nousage ${OpenMP_Fortran_FLAGS}")
        #set(_FFLAGS_RELEASE "-ipo")
    endif (WIN32)
endif ()

# Store a copy of "original" FFLAGS without any additional warning options
# applied to external libraries. These will be used for external examples
# which otherwise generate a lot of warnings that we want to ignore.
set(CMAKE_Fortran_FLAGS_DEBUG_ORIG "${CMAKE_Fortran_FLAGS_DEBUG}")
set(CMAKE_Fortran_FLAGS_ORIG "${CMAKE_Fortran_FLAGS}")

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_FFLAGS}")

string(TOUPPER "${CMAKE_BUILD_TYPE}" _CONFIG)
message(STATUS "Flags for EXTERNAL code: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${_CONFIG}}")

################################################################################

# add_subdirectory(ARPACK)
add_subdirectory(fitpack)
add_subdirectory(lbfgsb_bmnz)
add_subdirectory(simplex_csiro)
add_subdirectory(minpack)
add_subdirectory(newuoa2)
add_subdirectory(orderpack)
add_subdirectory(quadpack)
add_subdirectory(CDFLIB90)
add_subdirectory(random)
add_subdirectory(sobol)
add_subdirectory(slsqp)
