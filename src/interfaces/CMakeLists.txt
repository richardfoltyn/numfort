set(COMPONENT_NAME interfaces)

set(SOURCE_FILES
    blas.f90
    lapack.f90
)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
