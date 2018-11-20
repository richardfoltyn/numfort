
get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(${COMPONENT_NAME}_FILES
    random_integer.f90
    descriptive.f90
    dfrechet.f90
    dlogistic.f90
    dmvnorm.f90
    dnorm.f90
    markov.f90
    PARENT_SCOPE
)
