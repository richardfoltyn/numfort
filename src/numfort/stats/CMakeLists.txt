
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(SOURCE_FILES
    common.f90
    random.F90
    stats_core.F90
    dlogistic.F90
    dnorm.F90
    dmvnorm.F90
    duniform.F90
    dfrechet.F90
    dgenpareto.F90
    dlognorm.F90
    markov.F90
    stats.f90
    stats_lm.F90
    stats_lm_common.f90
    stats_lm_real32.F90
    stats_lm_real64.F90
    stats_qrng.F90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} arrays common core linalg interpolate)
# external libraries
add_dependencies(${COMPONENT_NAME} cdflib90 random sobol_direction_num)
