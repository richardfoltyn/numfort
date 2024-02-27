get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(${COMPONENT_NAME}_FILES
    glmnet.f90
    stats_pcr.f90
    sobol_seq.F90
    PARENT_SCOPE
)
