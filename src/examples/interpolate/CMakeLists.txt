get_filename_component(COMPONENT_NAME "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

set(${COMPONENT_NAME}_FILES
    concon.f90
    curfit.f90
    splder.f90
    PARENT_SCOPE
)

set(${COMPONENT_NAME}_FILES_OMP
    curfit_omp.f90
    concon_omp.f90
    PARENT_SCOPE
)
