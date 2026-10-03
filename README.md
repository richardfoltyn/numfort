# Numfort

_numfort_ implements some of the functionality of 
[NumPy](https://numpy.org/) and [SciPy](https://scipy.org/) in
Fortran, albeit on a much less ambitious scale. The main goal
is to provide routines that make solving quantitative macroeconomic
models in Fortran easier.

## Table of Contents

- [Obtaining the code](#obtaining-the-code)
- [Build instructions](#build-instructions)
- [Usage](#usage)
- [License](#license)
- [Authors](#authors)

***
## Obtaining the code

1. To clone the git repository, run 
    
    ```bash
    git clone https://github.com/richardfoltyn/numfort.git --recursive
    ```
        
1.  _numfort_ contains a git submodule with CMake scripts required to build the 
    library. It needs to be initialized using
    ```bash
    git submodule update --init --recursive
    ```
    This checks out the version recorded by the repository. After a `git pull`,
    run the same command to update the submodule to the recorded revision.

***
## Build instructions

### Linux

_numfort_ is built using [CMake](https://cmake.org/).
The library itself does not have any compile-time dependencies, but 
requires BLAS and LAPACK libraries to be present when linking
any client application using _numfort_.

Set the GCC major version in your shell. On Fedora this sets the build path;
the default compiler must be the matching version.

Bash:
```bash
GCC_VERSION=16
```

Fish:
```fish
set GCC_VERSION 16
```

The remaining commands work in both shells and build outside the source tree.

#### Debian/Ubuntu ####

On Debian/Ubuntu, the compiler executables have version suffixes:

```sh
cmake -S "$HOME/repos/numfort" -B "$HOME/build/gnu/$GCC_VERSION/numfort" \
    -DCMAKE_Fortran_COMPILER=gfortran-$GCC_VERSION \
    -DCMAKE_C_COMPILER=gcc-$GCC_VERSION \
    -DCMAKE_INSTALL_PREFIX="$HOME/.local"
```

#### Fedora ####

On Fedora, the default GCC compiler does not have a version suffix:

```sh
cmake -S "$HOME/repos/numfort" -B "$HOME/build/gnu/$GCC_VERSION/numfort" \
    -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_INSTALL_PREFIX="$HOME/.local"
```

#### Building and installing ####

```sh
cmake --build "$HOME/build/gnu/$GCC_VERSION/numfort" --parallel 8
cmake --install "$HOME/build/gnu/$GCC_VERSION/numfort"
```

For GCC 16 on Fedora, the CMake package is installed under
`~/.local/lib64/numfort-0.1-gnu-16/cmake/` and modules under
`~/.local/include/numfort-0.1-gnu-16/`.

#### Using a specific compiler ####

You may want to specify an alternative Fortran compiler or compile flags (`FFLAGS`).
For example, to build using Intel's Fortran compiler `ifx` or the now-deprecated `ifort` and optimize the
code for the host machine architecture, you could run:

```sh
cmake -S "$HOME/repos/numfort" -B "$HOME/build/intel/numfort" \
    -DCMAKE_Fortran_COMPILER=ifort -DCMAKE_Fortran_FLAGS=-xHost \
    -DCMAKE_INSTALL_PREFIX=/path/to/install
```
_numfort_ was tested to compile with the following compilers:

-   GNU `gfortran` 11.x - 16.x
-   Intel `ifort` 2021 and `ifx` 2024

***
### Advanced build instructions

To build the unit tests or example files, several libraries are needed at
compile time:
1.  The _fcore_ library available [here](https://github.com/richardfoltyn/fortran-corelib).
    The latter is built as a CMake project, so we use the 
    `CMAKE_PREFIX_PATH` variable to specify where CMake should look for it.

1.  BLAS and LAPACK libraries: _numfort_ works either with Intel MKL
    or some other implementation such as OpenBLAS.
    1.  For Intel MKL, the variable `MKL_ROOT` should point to the
        desired MKL installation directory.
        Optionally, `MKL_FORTRAN95_ROOT` can be specified if the Fortran 95
        wrappers for BLAS and LAPACK are available (this is only required
        for gfortran, MKL itself ships the required libraries for `ifort` and `ifx`).
        
        On Fedora, after setting `GCC_VERSION` as above and installing fcore
        and MKL95, set the MKL version in your shell:

        Bash:
        ```bash
        MKL_VERSION=2026.1
        ```

        Fish:
        ```fish
        set MKL_VERSION 2026.1
        ```

        Then configure with either shell:

        ```sh
        cmake -S "$HOME/repos/numfort" -B "$HOME/build/gnu/$GCC_VERSION/numfort" \
            -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc \
            -DCMAKE_INSTALL_PREFIX="$HOME/.local" \
            -DCMAKE_PREFIX_PATH="$HOME/.local" \
            -DMKL_ROOT="/opt/intel/oneapi/mkl/$MKL_VERSION" \
            -DMKL_FORTRAN95_ROOT="$HOME/.local/share/mkl/$MKL_VERSION/gnu/$GCC_VERSION" \
            -DBUILD_TESTS=ON -DBUILD_EXAMPLES=ON
        ```
                
    2.  If another BLAS/LAPACK implementation should be used, set
        `USE_MKL=OFF` and _numfort_ will use whichever library is found
        by CMake.
         
        On Fedora, configure with:

        ```sh
        cmake -S "$HOME/repos/numfort" -B "$HOME/build/gnu/$GCC_VERSION/numfort" \
            -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_C_COMPILER=gcc \
            -DCMAKE_INSTALL_PREFIX="$HOME/.local" \
            -DCMAKE_PREFIX_PATH="$HOME/.local" \
            -DUSE_MKL=OFF \
            -DBUILD_TESTS=ON -DBUILD_EXAMPLES=ON
        ```

***
## Usage

Add the following to a CMake-based project to use the _numfort_ library:

```cmake
find_package(numfort REQUIRED)

target_link_libraries(<target> PRIVATE numfort::numfort)
```

***
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
Any changes and additions made to the code in `src/external` are licensed
under a project's respective original license.
See `LICENSES_bundled.txt` and the subdirectories in `src/external` for details.

## Authors

With the exception of the files in `src/external`, _numfort_ was written
by Richard Foltyn. See the `src/external` directories for information
about the authors of these bundled components.
