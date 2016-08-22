# src/fitpack/fitpack/CMakeLists.txt
# CMake configuration for numfort library

set(SOURCE_FILES
	mnbisp.f
	mncloc.f
	mncoco.f
	mnconc.f
	mncosp.f
	mncual.f
	mncuev.f
	mncurf.f
	mndbin.f
	mnevpo.f
	mnfour.f
	mninst.f
	mnpade.f
	mnparc.f
	mnpasu.f
	mnperc.f
	mnpogr.f
	mnpola.f
	mnprof.f
	mnregr.f
	mnspal.f
	mnspde.f
	mnspev.f
	mnspgr.f
	mnsphe.f
	mnspin.f
	mnspro.f
	mnsuev.f
	mnsurf.f)

# set output directory to test directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${${PROJECT_NAME}_EXAMPLES_DIRECTORY}")

# Add all tests as separate executables
foreach (_file IN LISTS SOURCE_FILES)
	get_filename_component(exe_name "${_file}" NAME_WE)
	set(_target dierckx_${exe_name})

	add_executable(${_target} "${_file}")
	target_link_libraries(${_target} ${${PROJECT_NAME}_LIBRARIES})
	set_target_properties(${_target} PROPERTIES OUTPUT_NAME ${exe_name})

endforeach (_file)

# copy sample data from samples/ subdirectory to test directory
file(GLOB FITPACK_SAMPLE_DATA "da*.txt")
file(COPY ${FITPACK_SAMPLE_DATA} DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
