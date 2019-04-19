set(COMPONENT_NAME interfaces)

set(SOURCE_FILES
    blas.f90
    lapack.f90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})

install(FILES
    "${NUMFORT_MODULE_DIRECTORY}/blas_interfaces.mod"
    "${NUMFORT_MODULE_DIRECTORY}/lapack_interfaces.mod"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${NUMFORT_INSTALL_DIR}"
)
