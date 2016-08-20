
set(COMPONENT_NAME optimize)

set(SOURCE_FILES)

# subdirectories add files to SOURCE_FILES
add_subdirectory(lbfgsb_bmnz)

set(F08_SOURCE_FILES
    optimize.f90
    common.f90
    optim_result.f90
    minimize_lbfgsb.f90
)

set(SOURCE_FILES
    ${SOURCE_FILES}
    ${F08_SOURCE_FILES}
)

set_source_files_properties(${F08_SOURCE_FILES}
    PROPERTIES COMPILE_FLAGS "${${PROJECT_NAME}_FFLAGS}"
)

include_directories(${${PROJECT_NAME}_MODDIR})

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
add_dependencies(${COMPONENT_NAME} common)
# set_target_properties(${COMPONENT_NAME}
#	PROPERTIES COMPILE_FLAGS "${${PROJECT_NAME}_FFLAGS}"
# )
