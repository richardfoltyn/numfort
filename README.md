# Numfort #

## About ##

NUMFORT replicates some of the functionality of Numpy and Scipy in
Fortran, albeit on a much less ambitios scale. The main goal
is to provide routines that make solving quantitative macroeconomic
models in Fortran easier.

## Obtaining the code ##

1. To clone the git repository, run 
    
        git clone https://bitbucket.org/richardfoltyn/numfort.git --recursive
        
1.  NUMFORT contains a git submodule with CMake scripts required to build the 
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

The library itself does not have any compile-time dependencies, but 
requires BLAS and LAPACK libraries to be present when linking
any client application using NUMFORT.

To compile the library, create a build directory and run

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
make -j 16
make install
```

You may want to specify an alternative Fortran compiler or compile flags (`FFLAGS`).
For example, to build using Intel's Fortran compiler `ifort`and optimize the
code for the host machine architecture, you could run
```bash
FC=ifort FFLAGS="-xHost" \
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/ <NUMFORT_REPOSITORY>/src
```
NUMFORT was tested with the following compilers:

-   Intel ifort 2018 and 2019
-   GNU gfortran 7 and 8



### Advanced build instructions ###

To build the unit tests or example files, several libraries are needed at 
compile time:
1.  The FCORE library available [here](https://bitbucket.org/richardfoltyn/fortran-corelib).
    The latter is build as a CMake project, so we use the 
    CMAKE_PREFIX_PATH variable to specify where CMake should look for it.

1.  BLAS and LAPACK libraries: NUMFORT works either with Intel MKL
    or the some other implementation such as OpenBLAS.
    1.  For Intel MKL, the variable MKL_ROOT should point to the 
        desired MKL installation directory.
        Optionally, MKL_FORTRAN95_ROOT can be specified if the Fortran 95
        wrappers for BLAS and LAPACK are available (this is only required
        for gfortran, MKL itself ships the required libraries for ifort).
        
        #### Linux ####
        
        The corresponding CMake invocation on Linux could be something like
        
            cmake -DCMAKE_PREFIX_PATH=/path/to/fcore \
                -DMKL_ROOT=/opt/intel/compilers_and_libraries_20xx/linux/mkl \
                -DMKL_FORTRAN95_ROOT=/path/to/mkl/fortran95 \
                -DBUILD_TESTS=ON \
                -DBUILD_EXAMPLES=ON \
                <NUMFORT_REPOSITORY>/src
                
    2.  If another BLAS/LAPACK implementation should be used, set
        `USE_MKL=OFF` and NUMFORT will use whichever library is found 
         by CMake.
         
         #### Linux ####
         
         On Linux, the CMake invocation could be something like
         
            cmake -DCMAKE_PREFIX_PATH=/path/to/fcore \
                -DUSE_MKL=OFF \
                -DBUILD_TESTS=ON \
                -DBUILD_EXAMPLES=ON \
                <NUMFORT_REPOSITORY>/src