# fitpack/CMakeLists.txt
# CMake configuration for FITPACK library

cmake_minimum_required(VERSION 3.0)

project(fitpack Fortran)
set(PROJECT_NAME_PRETTY FITPACK)
set(${PROJECT_NAME}_LIBRARY_NAME ${PROJECT_NAME})

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

################################################################################
# User-definable options
if (WIN32)
    option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
else ()
    option(BUILD_SHARED_LIBS "Build shared libraries" ON)
endif ()
option(BUILD_EXAMPLES "Build examples/demo code" OFF)

################################################################################
# Numfort subproject setup

if (${CMAKE_PROJECT_NAME} STREQUAL "numfort")
	set(BUILD_OBJECT_LIBS ON)
	set(${PROJECT_NAME}_EXAMPLES_DIRECTORY
		"${${CMAKE_PROJECT_NAME}_EXAMPLES_DIRECTORY}/external/${PROJECT_NAME}"
	)
    # add to list of object libraries that will be used to create ultimate
    # numfort library
    set_property(GLOBAL APPEND
        PROPERTY ${CMAKE_PROJECT_NAME}_COMPONENT_TARGETS
        ${${PROJECT_NAME}_LIBRARY_NAME}
    )
    message(STATUS "Building ${PROJECT_NAME_PRETTY} as sub-project")
else ()
	set(BUILD_OBJECT_LIBS OFF)
	set(${PROJECT_NAME}_EXAMPLES_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
	message(STATUS "Building ${PROJECT_NAME_PRETTY} as standalone project")
endif ()

# Link demo executables against these libraries
set(${PROJECT_NAME}_LIBRARIES ${${PROJECT_NAME}_LIBRARY_NAME}_static)

################################################################################
# Set Fortran compiler FLAGS

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # enable optimization that does not break debugging
    set(CMAKE_Fortran_FFLAGS_DEBUG "-Og")
endif ()

################################################################################
# Target definitions

set(SOURCE_FILES
	bispeu.f
	bispev.f
	clocur.f
	cocosp.f
	concon.f
	concur.f
	cualde.f
	curev.f
	curfit.f
	dblint.f
	evapol.f
	fourco.f
	fpader.f
	fpadno.f
	fpadpo.f
	fpback.f
	fpbacp.f
	fpbfout.f
	fpbisp.f
	fpbspl.f
	fpchec.f
	fpched.f
	fpchep.f
	fpclos.f
	fpcoco.f
	fpcons.f
	fpcosp.f
	fpcsin.f
	fpcurf.f
	fpcuro.f
	fpcyt1.f
	fpcyt2.f
	fpdeno.f
	fpdisc.f
	fpfrno.f
	fpgivs.f
	fpgrdi.f
	fpgrpa.f
	fpgrre.f
	fpgrsp.f
	fpinst.f
	fpintb.f
	fpknot.f
	fpopdi.f
	fpopsp.f
	fporde.f
	fppara.f
	fppasu.f
	fpperi.f
	fppocu.f
	fppogr.f
	fppola.f
	fprank.f
	fprati.f
	fpregr.f
	fprota.f
	fprppo.f
	fprpsp.f
	fpseno.f
	fpspgr.f
	fpsphe.f
	fpsuev.f
	fpsurf.f
	fpsysy.f
	fptrnp.f
	fptrpe.f
	insert.f
	parcur.f
	parder.f
	pardeu.f
	parsur.f
	percur.f
	pogrid.f
	polar.f
	profil.f
	regrid.f
	spalde.f
	spgrid.f
	sphere.f
	splder.f
	splev.f
	splint.f
	sproot.f
	surev.f
	surfit.f
)

set(_name ${${PROJECT_NAME}_LIBRARY_NAME})
add_library(${_name} OBJECT ${SOURCE_FILES})
if (BUILD_SHARED_LIBS)
    set_target_properties(${_name} PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif ()

# build shared/static library from object library only if standalone build.
# Otherwise, parent project will use object library directly.
if (${PROJECT_NAME}_STANDALONE)
    add_library(${_name} $<TARGET_OBJECTS:${_name}>)
endif ()

if (BUILD_EXAMPLES)
    # create static library which examples can link against
    add_library(${_name}_static $<TARGET_OBJECTS:${_name}>)
	add_subdirectory(examples)
endif ()