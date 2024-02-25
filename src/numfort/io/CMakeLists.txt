# src/numfort/numfort/interp/CMakeLists.txt

set(COMPONENT_NAME io)

set(SOURCE_FILES
    io.f90
    io_common.F90
    io_fixed.F90
)

unset(CPP_DEFINITIONS)
if(UNIX)
    list(APPEND CPP_DEFINITIONS _LINUX)
elseif (WIN32)
    list(APPEND CPP_DEFINITIONS _WIN32)
endif ()

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
target_compile_definitions(${COMPONENT_NAME} PRIVATE ${CPP_DEFINITIONS})

add_dependencies(${COMPONENT_NAME} common)
