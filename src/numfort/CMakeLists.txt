# src/numfort/CMakeLists.txt
# CMake configuration for numfort library

################################################################################
# Additional compiler flags for core numfort library (can be more restrictive
# than flags for external code, examples, unit tests, etc.

unset(_FFLAGS)
unset(_FFLAGS_DEBUG)
unset(_FFLAGS_RELEASE)


if (USE_OPENMP)
    find_package(OpenMP)
endif()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" OR MINGW)
    # explicitly add -cpp to run C preprocessor if required
    set(GNU_WARN "-Wall -Wextra -Wno-compare-reals -Warray-temporaries -Wrealloc-lhs -pedantic -Wuninitialized")

    # Disable "surprising" warnings for older GFORTRAN version as this
    # produces endless warnings about array finalizers
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 7)
        set(GNU_WARN "${GNU_WARN} -Wno-surprising")
    endif()
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 6)
        # breaks compilation with gfortran 5.x, possibly with 6.x as well.
        set(GNU_WARN "${GNU_WARN} -Wno-unused-dummy-argument")
    endif()
    set(_FFLAGS "${GNU_WARN} --std=f2008 ${OpenMP_Fortran_FLAGS}")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (WIN32)
        set(_FFLAGS "/warn:all /standard-semantics /stand:f08 ${OpenMP_Fortran_FLAGS}")
        # NOTE: When setting -standard-semantics, resets any previously
        # given -assume std_mod_proc_name, so we need to append it again!
        if (NOT IFORT_STD_MOD_PROC_NAME)
            set(_FFLAGS "${_FFLAGS} /assume:nostd_mod_proc_name")
        else ()
            set(_FFLAGS "${_FFLAGS} /assume:std_mod_proc_name")
        endif ()
    else ()
        set(_FFLAGS "-warn all -standard-semantics -std08 ${OpenMP_Fortran_FLAGS}")
        set(_FFLAGS_REPORT "-qopt-report=2 -qopt-report-phase=ipo,loop,par,vec")
        if (OPT_REPORT)
            set(_FFLAGS "${_FFLAGS} ${_FFLAGS_REPORT}")
        endif()
        # NOTE: When setting -standard-semantics, resets any previously
        # given -assume std_mod_proc_name, so we need to append it again!
        if (NOT IFORT_STD_MOD_PROC_NAME)
            set(_FFLAGS "${_FFLAGS} -assume nostd_mod_proc_name")
        else ()
            set(_FFLAGS "${_FFLAGS} -assume std_mod_proc_name")
        endif ()
    endif (WIN32)
endif ()

# Set FFLAGS for this and all sub-directories which are part of the core
# NUMFORT library
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${_FFLAGS}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${_FFLAGS_DEBUG}")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${_FFLAGS_RELEASE}")

#set_property(DIRECTORY PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)

string(TOUPPER "${CMAKE_BUILD_TYPE}" _CONFIG)
message(STATUS "Flags for NUMFORT library: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${_CONFIG}}")

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
    # set mod-file output directory
    set_target_properties(${_comp} PROPERTIES
        Fortran_MODULE_DIRECTORY "${NUMFORT_MODULE_DIRECTORY}"
    )

    # Add component-specific include directory, if applicable
    if (IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/${_comp}/include")
        target_include_directories(${_comp} PRIVATE
            "${CMAKE_CURRENT_SOURCE_DIR}/${_comp}/include"
        )
    endif()

    target_include_directories(${_comp} PRIVATE
        ${NUMFORT_INCLUDE_DIRS}
    )

    # we need to make sure that -fPIC etc. is passed to compiler when building
    # object libraries that will be used to create a shared lib
    if (BUILD_SHARED_LIBS)
        set_target_properties(${_comp} PROPERTIES POSITION_INDEPENDENT_CODE ON)
    endif ()

    # append to list of object library targets
    set_property(GLOBAL APPEND
        PROPERTY ${CMAKE_PROJECT_NAME}_COMPONENT_TARGETS ${_comp}
    )
endforeach()


# Update fortran version module with current version information
configure_file(numfort.f90.in
    ${CMAKE_BINARY_DIR}/numfort.f90
)

# Add a few general source files, makes it easier in case generator does
# not support targets w/o source files.
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

add_library(${NUMFORT_LIBRARY_NAME} ${SOURCE_FILES} ${COMPONENT_OBJECTS})

set_target_properties(${NUMFORT_LIBRARY_NAME} PROPERTIES
    Fortran_MODULE_DIRECTORY "${NUMFORT_MODULE_DIRECTORY}"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)

set_target_properties(${NUMFORT_LIBRARY_NAME} PROPERTIES
    VERSION ${NUMFORT_VERSION}
)


################################################################################
# Installation

include(CMakePackageConfigHelpers)

# Path there cmake config files will be installed
set(NUMFORT_INSTALL_CONFIG_PATH
    ${CMAKE_INSTALL_LIBDIR}/${NUMFORT_INSTALL_DIR}/cmake
)

set(NUMFORT_INSTALL_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../install")
set(NUMFORT_INSTALL_NAMESPACE ${NUMFORT_LIBRARY_NAME})

set(NUMFORT_CONFIG_FILE ${NUMFORT_LIBRARY_NAME}-config.cmake)
set(NUMFORT_CONFIG_VERSION_FILE ${NUMFORT_LIBRARY_NAME}-config-version.cmake)
set(NUMFORT_TARGETS_FILE ${NUMFORT_LIBRARY_NAME}Targets.cmake)

configure_file("${NUMFORT_INSTALL_SOURCE_DIR}/config-version.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${NUMFORT_CONFIG_VERSION_FILE}"
    @ONLY
)

configure_package_config_file("${NUMFORT_INSTALL_SOURCE_DIR}/config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${NUMFORT_CONFIG_FILE}"
    INSTALL_DESTINATION "${NUMFORT_INSTALL_CONFIG_PATH}"
    PATH_VARS NUMFORT_INSTALL_CONFIG_PATH
)

# Note: INCLUDES DESTINATION is equivalent to adding
# INTERFACE_INCLUDE_DIRECTORIES property on every exported target
install(TARGETS ${NUMFORT_LIBRARY_NAME}
    EXPORT ${NUMFORT_LIBRARY_NAME}_targets
    RUNTIME DESTINATION "${CMAKE_INSTALL_LIBDIR}/${NUMFORT_INSTALL_DIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}/${NUMFORT_INSTALL_DIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}/${NUMFORT_INSTALL_DIR}"
    INCLUDES DESTINATION
        "${CMAKE_INSTALL_INCLUDEDIR}/${NUMFORT_INSTALL_DIR}"
)

export(EXPORT ${NUMFORT_LIBRARY_NAME}_targets
    FILE "${CMAKE_CURRENT_BINARY_DIR}/${NUMFORT_TARGETS_FILE}"
    NAMESPACE ${NUMFORT_INSTALL_NAMESPACE}
)

install(EXPORT ${NUMFORT_LIBRARY_NAME}_targets
    FILE "${NUMFORT_TARGETS_FILE}"
    NAMESPACE ${NUMFORT_INSTALL_NAMESPACE}::
    DESTINATION "${NUMFORT_INSTALL_CONFIG_PATH}"
)


# Install CMake config files that were generated above
install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/${NUMFORT_CONFIG_FILE}"
    "${CMAKE_CURRENT_BINARY_DIR}/${NUMFORT_CONFIG_VERSION_FILE}"
    "${NUMFORT_INSTALL_SOURCE_DIR}/compiler_test.f90"
    DESTINATION "${NUMFORT_INSTALL_CONFIG_PATH}"
    COMPONENT devel
)

# Install only numfort's MOD files, ignore any modules potentially created
# by external projects
# TODO: Do not install MOD files that are used only internally.
install(DIRECTORY "${NUMFORT_MODULE_DIRECTORY}/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${NUMFORT_INSTALL_DIR}"
    FILES_MATCHING PATTERN "numfort*.mod"
    PATTERN "external" EXCLUDE
)

