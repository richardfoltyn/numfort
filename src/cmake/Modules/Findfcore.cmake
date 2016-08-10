# Findfcore.cmake
# Attempts to find local installation of Fortran corelib (fcore) library
# Usage:
#   find_package(fcore)
#
# Upon successful completion defines the following cached variables:
#   fcore_FOUND : TRUE if fcore fcore
#   fcore_INCLUDE_DIRS : include directories for compiled module files
#   fcore_LIBRARIES : fcore library paths

if (WIN32)
    set(HOME $ENV{USERPROFILE})
else ()
    set(HOME $ENV{HOME})
endif ()

include(FindTargetArch)
find_target_arch()

# additional subdirectories to look through when searching for include files
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(FCORE_ARCH_SUBDIR "gfortran")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (TARGET_ARCH_BITS STREQUAL "64")
        set(FCORE_ARCH_SUBDIR "intel64")
    else()
        set(FCORE_ARCH_SUBDIR "intel")
    endif()
endif ()

set(FCORE_NAME fcore)
set(FCORE_MODFILES assertion_mod.mod
    collection_mod.mod
    foounit_mod.mod
    linked_list_mod.mod
    corelib_string.mod
    test_case_mod.mod)

find_library(${FCORE_NAME}_LIBRARY
    NAMES ${FCORE_NAME}
    HINTS
    ${HOME}/lib/${FCORE_ARCH_SUBDIR}
    ${HOME}/local/lib/${FCORE_ARCH_SUBDIR}
    ${HOME}/.local/lib/${FCORE_ARCH_SUBDIR}    
    PATHS
    ${CMAKE_SYSTEM_LIBRARY_PATH}
    ${HOME}/lib
    ${HOME}/local/lib
    ${HOME}/.local/lib
)

find_path(${FCORE_NAME}_INCLUDE_DIR NAMES ${FCORE_MODFILES}
    HINTS
    ${HOME}/include/${FCORE_ARCH_SUBDIR}
    ${HOME}/local/include/${FCORE_ARCH_SUBDIR}
    ${HOME}/.local/include/${FCORE_ARCH_SUBDIR}
    PATHS
    ${HOME}/include
    ${HOME}/local/include
    ${HOME}/.local/include
    PATH_SUFFIXES ${FCORE_NAME})

find_package_handle_standard_args(${FCORE_NAME} DEFAULT_MSG
    ${FCORE_NAME}_LIBRARY
    ${FCORE_NAME}_INCLUDE_DIR)

if(${FCORE_NAME}_FOUND)
    set(${FCORE_NAME}_LIBRARIES ${${FCORE_NAME}_LIBRARY}
        CACHE STRING "FCore libraries")
    set(${FCORE_NAME}_INCLUDE_DIRS ${${FCORE_NAME}_INCLUDE_DIR}
        CACHE PATH "FCore include directories")

    if(NOT ${FCORE_NAME}_FIND_QUIETLY)
        message(STATUS "Found ${FCORE_NAME} include directories: ${${FCORE_NAME}_INCLUDE_DIRS}")
    endif()
else()
    if(${FCORE_NAME}_FIND_REQUIRED)
        message(FATAL_ERROR "Could not find ${FCORE_NAME}")
    endif()
endif()
