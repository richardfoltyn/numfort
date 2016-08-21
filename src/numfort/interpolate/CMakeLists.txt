# src/numfort/numfort/interp/CMakeLists.txt

set(COMPONENT_NAME interpolate)

set(SOURCE_FILES
	interpolate.f90
	common.f90
	interp_result.f90
	interp_linear.f90
	fitpack_interfaces.f90
	fitpack08.f90
)

set_source_files_properties(${SOURCE_FILES}
	PROPERTIES COMPILE_FLAGS "${${PROJECT_NAME}_FFLAGS}"
)

# Add directory with F77 FITPACK code; Adds F77 files to SOURCE_FILES
add_subdirectory(fitpack)

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
