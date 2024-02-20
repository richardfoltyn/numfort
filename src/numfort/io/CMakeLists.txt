# src/numfort/numfort/interp/CMakeLists.txt

option(USE_HDF5 "Build with HDF5 wrapper to load/save files" ON)
option(HDF5_STATIC "Use static linking against HDF5" ON)

################################################################################
# HDF5 library
#
# We attempt to find HDF5 installation using the following steps:
# 1.    If HDF5_ROOT is specified, try to locate hdf5-config.cmake in the
#       tree beneath the root folder.
# 2.    If no HDF5_ROOT is specified, search for hdf5-config.cmake in
#       default CMAKE search paths.
# 3.    Assume that no hdf5-config.cmake exists and use FindHDF5 that was
#       shipped with CMake.

if (USE_HDF5)

    # HDF5 package names used in variables created by FindHDF5
    string(TOUPPER hdf5 HDF5_PACKAGE_NAME)

    # Link against shared or static library
    set(lib_type "static")
    if (NOT HDF5_STATIC)
        set(lib_type "shared")
    endif ()

    # components to look for (does not include lib type, ie. static of shared)
    set(HDF5_FIND_COMPONENTS Fortran C)

    # store either -DHDF5_ROOT from command-line or env. variable HDF5_ROOT;
    # command-line spec. takes precedence!
    set(_hdf5_root "${HDF5_ROOT}")
    if (NOT _hdf5_root AND ENV{HDF5_ROOT})
        set(_hdf5_root "$ENV{HDF5_ROOT}")
    endif()

    # Invoke FIND_PACKAGE in config mode, ie directly including the
    # file hdf5-config.cmake that should have been created by the installer.
    if (_hdf5_root)
        find_package(HDF5 QUIET NO_MODULE
            COMPONENTS ${HDF5_FIND_COMPONENTS} ${lib_type}
            HINTS ${_hdf5_root}
            PATH_SUFFIXES cmake share/cmake
            NO_DEFAULT_PATH
        )
    else()
        find_package(HDF5 QUIET NO_MODULE
            COMPONENTS ${HDF5_FIND_COMPONENTS} ${lib_type}
        )
    endif()

    if (HDF5_VERSION_STRING)
        # this is usually set in the FindHDF5 file, not in hdf5-config.cmake
        set(HDF5_FOUND TRUE)

        unset(HDF5_LIBRARIES)
        foreach(comp IN LISTS HDF5_FIND_COMPONENTS)
            string(TOUPPER "${HDF5_PACKAGE_NAME}_${comp}_${lib_type}_LIBRARY" comp_lib)
            list(APPEND HDF5_LIBRARIES ${${comp_lib}})
        endforeach()

        list(APPEND HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIR_FORTRAN})

        # store version string for further processing below
        set(_hdf_version ${HDF5_VERSION_STRING})
    else()
        # Look for HDF5 using whatever came with CMake, most likely some FindHDF5
        # file.
        # Need to set environment variable to communicate desired
        # HDF5 root to FindHDF5 scripts
        set(_defined_hdf5_env FALSE)
        if (NOT ENV{HDF5_ROOT} AND HDF5_ROOT)
            set(_defined_hdf5_env TRUE)
            set(ENV{HDF5_ROOT} ${HDF5_ROOT})
        endif()

        # communicate to FindHDF5 that we want static libraries
        # Note: This defines the variables HDF5_LIBRARIES, HDF5_INCLUDE_DIRS
        # HDF5_Fortran_INCLUDE_DIRS
        set(HDF5_USE_STATIC_LIBRARIES ${HDF5_STATIC})
        find_package(HDF5 REQUIRED COMPONENTS ${HDF5_FIND_COMPONENTS})

        set(HDF5_LIBRARIES ${HDF5_LIBRARIES} ${HDF5_Fortran_LIBRARIES})
        set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS} ${HDF5_Fortran_INCLUDE_DIRS})

        # Append shared/ or static/ to HDF5 include path as this is present
        # on Windows installation but not correctly returned by FindHDF5
        if (WIN32)
            foreach(_lib IN LISTS HDF5_INCLUDE_DIRS)
                list(APPEND HDF5_INCLUDE_DIRS_SUFFIX "${_lib}/${lib_type}")
            endforeach()
            list(APPEND HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS_SUFFIX})
        endif ()

        # unset HDF5_ROOT if we defined it above to restore environment to
        # its previous state
        if (_defined_hdf5_env)
            unset(ENV{HDF5_ROOT})
        endif()

        set(_hdf_version ${HDF5_VERSION})
    endif()

    if (NOT WIN32)
        find_package(ZLIB REQUIRED)
        list(APPEND HDF5_LIBRARIES ${ZLIB_LIBRARIES})
    endif()
endif ()

if (USE_HDF5)
    message(STATUS "Building against HDF5 ${_hdf_version}")
    message(STATUS "HDF5 include path: ${HDF5_INCLUDE_DIRS}")
    message(STATUS "HDF5 link path: ${HDF5_LIBRARY_DIRS}")
    message(STATUS "HDF5 libraries: ${HDF5_LIBRARIES}")
endif ()


################################################################################

set(COMPONENT_NAME io)

set(SOURCE_FILES
    io.f90
    io_common.F90
    io_fixed.F90
    hdf5.F90
)

if (USE_HDF5)
    list(APPEND SOURCE_FILES
        hdf5.F90
    )
endif ()

unset(CPP_DEFINITIONS)
if(UNIX)
    list(APPEND CPP_DEFINITIONS _LINUX)
elseif (WIN32)
    list(APPEND CPP_DEFINITIONS _WIN32)
endif ()

add_library(${COMPONENT_NAME} OBJECT ${SOURCE_FILES})
target_compile_definitions(${COMPONENT_NAME} PRIVATE ${CPP_DEFINITIONS})

add_dependencies(${COMPONENT_NAME} common)

if (USE_HDF5)
    target_include_directories(${COMPONENT_NAME} PUBLIC
        ${HDF5_INCLUDE_DIRS}
    )

    target_link_libraries(${COMPONENT_NAME} PUBLIC
        ${HDF5_LIBRARIES}
    )

    # on Windows, add bin/ directory in HDF5_ROOT that contains debug symbols
    if (WIN32 AND EXISTS "${HDF5_ROOT}/bin")
        link_directories(${HDF5_ROOT}/bin)
    endif()
endif ()
