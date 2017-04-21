# Exports direction numbers required to create Sobol sequences to Fortran.
# See http://web.maths.unsw.edu.au/%7Efkuo/sobol/

import os
import numpy as np
import math
from io import StringIO

IDX_D = 0
IDX_S = 1
IDX_A = 2
IDX_M = 3


FORTRAN_MODULE_HEADER = r'''
! Warning: Automatically generated file, do not edit manually!
module sobol_direction_num
    use, intrinsic :: iso_fortran_env
    implicit none
'''

FORTRAN_MODULE_FOOTER = r'''
end module
'''

FORTRAN_LINEWIDTH = 78

SOBOL_MAX_DIM = 1024


def read_data(fn, nrow, ncol):

    d = np.empty(nrow, dtype=np.int32)
    s = np.empty_like(d)
    a = np.empty_like(d)

    m = np.zeros((nrow, ncol), dtype=np.int32)

    with open(fn, 'r') as f:
        # skip header line
        f.readline()
        for i, line in enumerate(f):
            dat = line.split()
            d[i] = dat[IDX_D]
            s[i] = dat[IDX_S]
            a[i] = dat[IDX_A]
            m[i,:s[i]] = dat[IDX_M:]

    return d, s, a, m


def write_fortran(fn, d, s, a, m):

    with open(fn, 'wt', encoding='ascii') as f:
        print(FORTRAN_MODULE_HEADER, file=f)

        fmt = 'integer, parameter :: SOBOL_DEFAULT_MAX_DIM = {:d}'
        print(fmt.format(SOBOL_MAX_DIM), file=f)

        fmt = 'integer, parameter :: SOBOL_DEFAULT_MAX_M = {:d}'
        max_nm = np.amax(s[:SOBOL_MAX_DIM])
        print(fmt.format(max_nm), file=f)

        # write_array(f, 'SOBOL_D', d[:SOBOL_MAX_DIM])
        write_array(f, 'SOBOL_DEFAULT_S', s[:SOBOL_MAX_DIM])
        write_array(f, 'SOBOL_DEFAULT_A', a[:SOBOL_MAX_DIM])

        mflat = np.empty(np.sum(s[:SOBOL_MAX_DIM]), dtype=np.int32)
        istart = 0
        for i in range(SOBOL_MAX_DIM):
            iend = istart + s[i]
            mflat[istart:iend] = m[i,:s[i]]
            istart = iend

        write_array(f, 'SOBOL_DEFAULT_M', mflat)

        print(FORTRAN_MODULE_FOOTER, file=f)


def write_array(f, name, data):
    prefix = ' ' * 4

    print('{:s}integer, parameter :: {:s}(*) = &'.format(prefix, name), file=f)

    width = math.ceil(math.log10(np.amax(data)))
    n = (FORTRAN_LINEWIDTH - len(prefix)) // (width + 2)
    nlines = math.ceil(len(data) / n)
    fmt = '{{:{:d}d}}'.format(width)
    for i in range(nlines):
        print(prefix, file=f, end='')
        istart = i*n
        iend = istart + n
        s = ', '.join(fmt.format(x) for x in data[istart:iend])
        print(s, file=f, end='')
        if i < nlines - 1:
            print(' &', file=f)
        else:
            print(' ]\n\n', file=f)


def main():
    fn = 'new-joe-kuo-6.21201.txt'
    fn_fort = 'sobol_direction_num_mod.f90'

    with open(fn, 'r+') as f:
        # skip first line, contains header
        lines = -1
        for line in f:
            lines += 1

    # determine number of columns from last line
    fields = line.split()
    # First three fields containd D, S and A
    ncol = len(fields) - 3

    d, s, a, m = read_data(fn, lines, ncol)
    write_fortran(fn_fort, d, s, a, m)


if __name__ == '__main__':
    main()
