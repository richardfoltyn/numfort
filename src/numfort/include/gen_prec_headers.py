from __future__ import print_function

import os
from os.path import join


BLAS_ROUTINES = ['GEMM']
LAPACK_ROUTINES = ['GESVD', 'GELSD', 'GESDD']

FILE_HEADER = '''
/*
    Include file for {{:s}} precision. 
    THIS FILE IS GENERATED AUTOMATICALLY BY {:s}, DO NOT MODIFY.
*/
'''.format(__file__)

ADD_DEFINITION_TEMPLATE = '''
#ifdef {name:s}
#undef {name:s}
#endif
#define {name:s} {value:s}
'''


def add_definition(name, value, file):
    s = ADD_DEFINITION_TEMPLATE.format(name=name, value=value)
    print(s, file=file)


def gen_prec(prec, prefix, filename):

    with open(filename, 'wt', encoding='ascii') as f:
        s = FILE_HEADER.format(prec)
        print(s, file=f)

        # Add precision definition
        add_definition('__PREC', prec, f)

        prefix = prefix.upper()

        s = '/* BLAS routines */'
        print(s, file=f)

        for name in BLAS_ROUTINES:
            defname = '__{:s}'.format(name.upper())
            value = '{:s}{:s}'.format(prefix, name)
            add_definition(defname, value, f)

        s = '/* LAPACK routines */'
        print(s, file=f)
        for name in LAPACK_ROUTINES:
            defname = '__{:s}'.format(name.upper())
            value = '{:s}{:s}'.format(prefix, name)
            add_definition(defname, value, f)


def main():

    prec = ['real32', 'real64']
    prefix = ['S', 'D']
    # Assume that script is run in include/ directory
    basedir = os.getcwd()

    for p, px in zip(prec, prefix):
        fn = 'numfort_{:s}.h'.format(p)
        fn = join(basedir, fn)

        gen_prec(p, px, fn)


if __name__ == '__main__':
    main()
