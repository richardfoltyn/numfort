
set(COMPONENT_NAME optimize)

set(SOURCE_FILES)

# subdirectories add files to SOURCE_FILES
add_subdirectory(lbfgsb_bmnz)

set(F08_SOURCE_FILES
    minimize_lbfgsb.f90
)

set(SOURCE_FILES
    ${SOURCE_FILES}
    ${F08_SOURCE_FILES}
)

set_source_files_properties(${F08_SOURCE_FILES}
    PROPERTIES COMPILE_FLAGS "${${PROJECT_NAME}_FFLAGS}"
)

include_directories(${${PROJECT_NAME}_MODDIR}/common)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})

# set_target_properties(${COMPONENT_NAME}
#	PROPERTIES COMPILE_FLAGS "${${PROJECT_NAME}_FFLAGS}"
# )
