
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(SOURCE_FILES
    common.f90
    data_helpers.f90
    data_helpers_logical.F90
    data_helpers_real32.F90
    data_helpers_real64.F90
    random.F90
    stats_core.f90
    stats_core_common.f90
    stats_core_real32.F90
    stats_core_real64.F90
    dlogistic.F90
    dnorm.F90
    dmvnorm.F90
    duniform.F90
    dfrechet.F90
    dgenpareto.F90
    dlognorm.F90
    glmnet.f90
    glmnet_real32.F90
    glmnet_real64.F90
    markov.f90
    markov_real64.F90
    markov_simulate_real64.F90
    markov_simulate_real64_int8.F90
    stats.f90
    stats_ineq_real32.F90
    stats_ineq_real64.F90
    stats_lm.F90
    stats_lm_common.f90
    stats_lm_real32.F90
    stats_lm_real64.F90
    stats_np.f90
    stats_np_real32.F90
    stats_np_real64.F90
    stats_qrng.F90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} arrays common core linalg interpolate)

# External dependencies
target_link_libraries(${COMPONENT_NAME} PRIVATE
    cdflib90 random sobol_direction_num
)
