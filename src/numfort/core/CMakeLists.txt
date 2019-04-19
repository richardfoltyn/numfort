# src/numfort/numfort/interp/CMakeLists.txt

set(COMPONENT_NAME core)

set(SOURCE_FILES
    core.F90
    constants.f90
    combinatorics.F90
    cephes_stats.F90
    libcwrapper.F90
    logexp.F90
    misc.F90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
