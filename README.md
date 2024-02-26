# Numfort #

_numfort_ implements some of the functionality of 
[NumPy](https://numpy.org/) and [Scipy](https://scipy.org/) in
Fortran, albeit on a much less ambitios scale. The main goal
is to provide routines that make solving quantitative macroeconomic
models in Fortran easier.

## Obtaining the code ##

1. To clone the git repository, run 
    
    ```bash
    git clone https://github.org/richardfoltyn/numfort.git --recursive
    ```
        
1.  _numfort_ contains a git submodule with CMake scripts required to build the 
    library. It needs to be initalized using
    ```bash
    git submodule init
    ```
    Submodules can be updated to the latest version by executing
    ```bash
    git submodule update --remote
    ```        
    Note that a simple `git pull` will not update the submodule.
        
## Build instructions ##

### Linux ###

_numfort_ is built using [CMake](https://cmake.org/).
The library itself does not have any compile-time dependencies, but 
requires BLAS and LAPACK libraries to be present when linking
any client application using _numfort_.

To compile the library, adapt the following to your environment:

```bash
# GCC compiler version
GCC_VERSION=12

# Define source directory
SRC_DIR=$HOME/repos/numfort/src

# Build directory
BUILD_DIR=$HOME/build/gnu/${GCC_VERSION}/numfort

# Installation prefix
INSTALL_PREFIX="$HOME/.local"

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

FC=gfortran-${GCC_VERSION} CC=gcc-${GCC_VERSION} \
cmake -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}" "${SRC_DIR}"
```

```bash
cmake --build .
cmake --install .
```

You may want to specify an alternative Fortran compiler or compile flags (`FFLAGS`).
For example, to build using Intel's Fortran compiler `ifx` or the now
deprecated `ifort` and optimize the
code for the host machine architecture, you could run
```bash
FC=ifort FFLAGS="-xHost" \
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/ /path/to/numfort/src
```
_numfort_ was tested to compile with the following compilers:

-   GNU `gfortran` 11.x, 12.x, and 13.x
-   Intel `ifort` 2021, `ifx` 2024


### Advanced build instructions ###

To build the unit tests or example files, several libraries are needed at 
compile time:
1.  The _fcore_ library available [here](https://github.com/richardfoltyn/fortran-corelib).
    The latter is build as a CMake project, so we use the 
    `CMAKE_PREFIX_PATH` variable to specify where CMake should look for it.

1.  BLAS and LAPACK libraries: _numfort_ works either with Intel MKL
    or the some other implementation such as OpenBLAS.
    1.  For Intel MKL, the variable `MKL_ROOT` should point to the 
        desired MKL installation directory.
        Optionally, `MKL_FORTRAN95_ROOT` can be specified if the Fortran 95
        wrappers for BLAS and LAPACK are available (this is only required
        for gfortran, MKL itself ships the required libraries for `ifort` and `ifx`).
        
        #### Linux ####
        
        The corresponding CMake invocation on Linux could be something like
        
            cmake -DCMAKE_PREFIX_PATH=/path/to/fcore \
                -DMKL_ROOT=/opt/intel/compilers_and_libraries_20xx/linux/mkl \
                -DMKL_FORTRAN95_ROOT=/path/to/mkl/fortran95 \
                -DBUILD_TESTS=ON \
                -DBUILD_EXAMPLES=ON \
                <NUMFORT_REPOSITORY>/src
                
    2.  If another BLAS/LAPACK implementation should be used, set
        `USE_MKL=OFF` and _numfort_ will use whichever library is found 
         by CMake.
         
         #### Linux ####
         
         On Linux, the CMake invocation could be something like
         
            cmake -DCMAKE_PREFIX_PATH=/path/to/fcore \
                -DUSE_MKL=OFF \
                -DBUILD_TESTS=ON \
                -DBUILD_EXAMPLES=ON \
                <NUMFORT_REPOSITORY>/src


## Usage

Add the following to a CMake-based project to use the _numfort_ library:

```CMake
find_package(numfort REQUIRED)

target_link_libraries(<target> PRIVATE numfort::numfort)
```

## License

This program is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software 
Foundation, either version 3 of the License, or (at your option) any later 
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE. See the GNU General Public License for more details.

_numfort_ includes several bundled projects in the `src/external` directory
which were distributed under various different licenses by their original authors.
Any changes and additions made to the code in `src/external` is licensed
under a project's respective original license.
See `LICENSES_bundled.txt` and the sub-directories in `src/external` for 
details.

## Authors

With the exception of the files in `src/external`, _numfort_ was written
by Richard Foltyn. See the `src/external` directories for information
about the authors of these bundled components.