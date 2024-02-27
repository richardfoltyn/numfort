set(COMPONENT_NAME interfaces)

set(SOURCE_FILES
    blas.f90
    lapack.f90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})

# append to list of object library targets
set_property(GLOBAL APPEND
    PROPERTY ${CMAKE_PROJECT_NAME}_COMPONENT_TARGETS ${COMPONENT_NAME}
)

# Manually add MOD files which don't have a numfort_ prefix
install(FILES
    "${CMAKE_Fortran_MODULE_DIRECTORY}/blas_interfaces.mod"
    "${CMAKE_Fortran_MODULE_DIRECTORY}/lapack_interfaces.mod"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${NUMFORT_INSTALL_DIR}"
)
