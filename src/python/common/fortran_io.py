"""
Module contains various functions used to convert Python data into Fortran
source code and write these to files.

Author: Richard Foltyn
"""

import numpy as np


FORTRAN_ARRAY_1D = '''
    real (PREC), parameter :: {name:s}(*) = [ real (PREC) :: &
        {data:s} ]
'''

FORTRAN_ARRAY_RESHAPE = '''
    real (PREC), parameter :: {name:s}({target_shape:s}) = &
        reshape([ real (PREC) :: &
        {data:s} ], &
        shape=[{target_shape:s}])
'''


def write_fortran_array(x, name, f, real_fmt='{:.12f}', indent=8, maxline=80):
    """
    Function convert numpy arrays to Fortran source code creating
    compile-time constant arrays.

    Parameters
    ----------
    x : array_like
    name : str
        Identifier or the Fortran array
    f : object
        File handle into which Fortran source code is printed
    real_fmt : str
        Python format descriptor used to convert floating-point numbers
        into text strings.
    indent : int
        Number of spaces used to indent Fortran array data
    maxline : int
        Max. number of characters per line, including any indentation white
        space, etc.
    """

    x = np.atleast_1d(x)
    shp = x.shape
    shp_F = np.array(shp)[::-1]

    # Store whether array was originally 1d before flattening
    is_1d = (x.ndim == 1)

    x = x.flatten()

    # Format array such that each row contains 4
    real_fmt = real_fmt+'_PREC'
    lines = []

    # Length if additional prefix and suffix characters for each line
    line_suffix = ', &'
    pad_len = len(line_suffix) + indent
    # within-line separator between numeric literals
    sep = ', '
    # current line length
    llen = pad_len
    # tokens on current line
    tokens = []
    for i, xi in enumerate(x):
        txt = real_fmt.format(xi)

        if len(tokens) == 0 or llen + len(txt) + len(sep) <= maxline:
            tokens.append(txt)
            llen += len(txt) + len(sep)
        else:
            line_txt = sep.join(tokens)
            lines.append(line_txt)

            # reset line tokens, append current token as first element
            llen = pad_len + len(txt) + len(sep)
            tokens = [txt]

    # Append any remaining tokens as last line
    if len(tokens) > 0:
        line_txt = sep.join(tokens)
        lines.append(line_txt)

    # merge individual rows using an indentation of 8 spaces
    sep = line_suffix + '\n' + ' ' * indent
    data = sep.join(lines)

    # Write constant array definition
    if is_1d:
        txt = FORTRAN_ARRAY_1D.format(name=name, data=data)
    else:
        target_shape_txt = ','.join('{:d}'.format(i) for i in shp_F)
        txt = FORTRAN_ARRAY_RESHAPE.format(
            name=name,
            data=data,
            target_shape=target_shape_txt)

    print(txt, file=f)
