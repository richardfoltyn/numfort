
################################################################################
# Add compiler flags specific to external projects

unset(EXTERNAL_FFLAGS)
unset(EXTERNAL_FFLAGS_DEBUG)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR MINGW)
    set(GNU_WARN -Wall -Wno-uninitialized -Wno-conversion -Wno-unused -Wno-tabs 
        -Wno-unused-dummy-argument)
    set(EXTERNAL_FFLAGS ${GNU_WARN} -std=legacy)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel(LLVM)?")
    if (WIN32)
        set(EXTERNAL_FFLAGS /warn:all,nounused,nousage)
    else ()
        set(EXTERNAL_FFLAGS -warn all,nounused,nousage ${OpenMP_Fortran_FLAGS})
    endif (WIN32)
endif ()


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
