
################################################################################
# Additional compiler flags for core numfort library (can be more restrictive
# than flags for external code, examples, unit tests, etc.

unset(_FFLAGS)
unset(_FFLAGS_Debug)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR MINGW)
    set(GNU_WARN -Wall -Wextra -Wno-compare-reals -Warray-temporaries 
        -Wrealloc-lhs -pedantic -Wuninitialized
    )

    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 9)
        # Skip warning about out-of-bounds array indices in do loops since
        # compiler cannot detect conditional execution.
        set(GNU_WARN ${GNU_WARN} -Wno-do-subscript)
    endif()
    set(_FFLAGS ${GNU_WARN} --std=f2018)
    set(_FFLAGS ${_FFLAGS} -finline-functions -fipa-cp-clone -fgraphite-identity 
        -floop-nest-optimize -ftree-loop-distribution -ftree-vectorize
    )
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel(LLVM)?")
    if (WIN32)
        set(_FFLAGS /warn:all /stand:f18 /auto)
    else ()
        set(_FFLAGS -warn all -std18 -auto)
        set(_FFLAGS_REPORT -qopt-report=2 -qopt-report-phase=ipo,loop,par,vec)
        if (OPT_REPORT)
            set(_FFLAGS ${_FFLAGS} ${_FFLAGS_REPORT})
        endif()
    endif (WIN32)
endif ()

if (USE_OPENMP)
    find_package(OpenMP)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()


################################################################################
# Target definitions

# add individual numfort components as object libraries
set(COMPONENTS
    core
    common
    arrays
    interpolate
    integrate
    io
    linalg
    optimize
    polynomial
    stats
)

foreach (_comp IN LISTS COMPONENTS)
    add_subdirectory(${_comp})

    # Add component-specific include directory, if applicable
    if (IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${_comp}/include")
        target_include_directories(${_comp} PRIVATE
            "${CMAKE_CURRENT_SOURCE_DIR}/${_comp}/include"
        )
    endif()

    target_include_directories(${_comp} 
        PRIVATE
            "${CMAKE_CURRENT_SOURCE_DIR}/include"
        PUBLIC
            "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${NUMFORT_INSTALL_DIR}>"
            "$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>"
    )

    if (USE_LTO)
        set_target_properties(${_comp} PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()

    # we need to make sure that -fPIC etc. is passed to compiler when building
    # object libraries that will be used to create a shared lib
    if (BUILD_SHARED_LIBS)
        set_target_properties(${_comp} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    endif ()

    target_compile_options(${_comp} BEFORE PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:${_FFLAGS}>
        $<$<COMPILE_LANGUAGE:Fortran>:$<$<CONFIG:Debug>:${_FFLAGS_Debug}>>
    )

    # append to list of object library targets
    set_property(GLOBAL APPEND
        PROPERTY ${CMAKE_PROJECT_NAME}_COMPONENT_TARGETS ${_comp}
    )
endforeach()


# Update fortran version module with current version information
configure_file(numfort.f90.in
    ${CMAKE_BINARY_DIR}/numfort.f90
)

# Additional sources that are not part of any object library
set(SOURCE_FILES
    ${CMAKE_BINARY_DIR}/numfort.f90 
)

# get list of targets the ultimate library depends on. This list is populated
# in the sub-components' directories
get_property(_TARGETS GLOBAL PROPERTY ${CMAKE_PROJECT_NAME}_COMPONENT_TARGETS)

set(COMPONENT_OBJECTS)
foreach (_comp IN LISTS _TARGETS)
    list(APPEND COMPONENT_OBJECTS $<TARGET_OBJECTS:${_comp}>)
endforeach ()

add_library(${NUMFORT_LIBRARY_NAME} 
    ${SOURCE_FILES}
    ${COMPONENT_OBJECTS}
)

set_target_properties(${NUMFORT_LIBRARY_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

set_target_properties(${NUMFORT_LIBRARY_NAME} PROPERTIES
    VERSION ${PROJECT_VERSION}
)

if (BUILD_TESTS OR BUILD_EXAMPLES)
    # build static library so tests & examples can link against it
    add_library(${NUMFORT_LIBRARY_NAME}_static STATIC 
        ${SOURCE_FILES}
        ${COMPONENT_OBJECTS}
    )

    set_target_properties(${NUMFORT_LIBRARY_NAME}_static PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${NUMFORT_STATIC_DIRECTORY}"
        LIBRARY_OUTPUT_DIRECTORY "${NUMFORT_STATIC_DIRECTORY}"
        ARCHIVE_OUTPUT_DIRECTORY "${NUMFORT_STATIC_DIRECTORY}"
        OUTPUT_NAME ${NUMFORT_LIBRARY_NAME}
    )

    # Add public include directories, this will propagage to tests and 
    # examples.
    target_include_directories(${NUMFORT_LIBRARY_NAME}_static
        PUBLIC
            "$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>"
            "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/numfort/include>"
    )

endif ()