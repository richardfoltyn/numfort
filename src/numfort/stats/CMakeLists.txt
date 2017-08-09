
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(SOURCE_FILES
    random.F90
    stats_core.F90
    dcont.F90
    ddisc.F90
    dnorm.F90
    duniform.F90
    drandint.F90
    dgenpareto.F90
    stats.f90
    stats_lm.F90
    stats_qrng.F90
)

include_directories(
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    ${${CMAKE_PROJECT_NAME}_INCLUDE_DIRS}
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common core)
# external libraries
add_dependencies(${COMPONENT_NAME} cdflib90 random sobol_direction_num)
