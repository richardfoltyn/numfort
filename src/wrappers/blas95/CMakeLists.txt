

set(COMPONENT_NAME "wrappers_blas95")


set(SOURCE_FILES
    blas95.F90
)

message(STATUS "Building BLAS95 wrappers")

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} interfaces)

target_include_directories(${COMPONENT_NAME}
    PRIVATE
    "${CMAKE_SOURCE_DIR}/src/numfort/include"
)

# append to list of object library targets
set_property(GLOBAL APPEND
    PROPERTY ${CMAKE_PROJECT_NAME}_COMPONENT_TARGETS ${COMPONENT_NAME}
)

if (BUILD_BLAS95)
    # Install blas95 wrappers mod file if not built against MKL; otherwise
    # client code is expected to use the interfaces provided by MKL.
    install(FILES
        "${CMAKE_Fortran_MODULE_DIRECTORY}/blas95.mod"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${NUMFORT_INSTALL_DIR}"
    )
endif()
