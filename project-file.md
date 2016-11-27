project: Numfort
summary:    Library providing a unified and convenient interface to
            existing and new numerical Fortran routines.
author: Richard Foltyn
project_bitbucket: https://bitbucket.org/richardfoltyn/numfort
src_dir: ./src/numfort
output_dir: ./docs
docmark: !
docmark_alt: *
predocmark: <
predocmark_alt: >
include:    src/numfort/core/include
            src/numfort/arrays/include
            src/numfort/common/include
            src/numfort/interpolate/include
            src/numfort/linalg/include
            src/numfort/optimize/include
exclude_dir:    src/numfort/core/include
                src/numfort/arrays/include
                src/numfort/common/include
                src/numfort/interpolate/include
                src/numfort/linalg/include
                src/numfort/optimize/include
exclude: common.F90
display: public
         protected
