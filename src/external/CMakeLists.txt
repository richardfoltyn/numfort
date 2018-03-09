
################################################################################
# Add compiler flags specific to external projects

unset(_FFLAGS_DEBUG)
unset(_FFLAGS)
unset(_FFLAGS_RELEASE)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR MINGW)
    # there shouldn't be many issues with LHS reallocation or temporary arrays
    # in old F77 code, but still add this.
    set(_FFLAGS "-Warray-temporaries -Wrealloc-lhs -Wno-unused -Wno-tabs -Wno-unused-dummy-argument")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (WIN32)
        set(_FFLAGS "/warn:all,nounused")
        #set(_FLAGS_RELEASE "/Qipo")
    else ()
        set(_FFLAGS "-warn all,nounused")
        #set(_FFLAGS_RELEASE "-ipo")
    endif (WIN32)
endif ()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_FFLAGS}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${_FFLAGS_DEBUG}")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${_FFLAGS_RELEASE}")

string(TOUPPER "${CMAKE_BUILD_TYPE}" _CONFIG)
message(STATUS "Flags for EXTERNAL code: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${_CONFIG}}")

################################################################################

add_subdirectory(ARPACK)
add_subdirectory(fitpack)
add_subdirectory(lbfgsb_bmnz)
add_subdirectory(simplex_csiro)
add_subdirectory(minpack)
add_subdirectory(orderpack)
add_subdirectory(quadpack)
add_subdirectory(CDFLIB90)
add_subdirectory(random)
add_subdirectory(sobol)
add_subdirectory(slsqp)
