"""
Python code to create the CSV data used as input in Fortran unit tests
that compare the Gauss-Hermite quadrature sample points and weights
to the ones produces by Numpy's hermegauss() function.

"""

import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from os.path import join, expanduser


# list of number of quadrature nodes (= Hermite polynomial degree) for
# which quadrature data should be created.
degrees = np.hstack((np.arange(1, 11), 20, 31, 49, 50, 100, 101))

# Output directory where quadrature data should be stored (in fixed-format CSV)
outdir = join(expanduser('~'), 'tmp')

data = np.zeros((0, 3))
for n in degrees:
    x, w = hermegauss(n)
    xx = np.tile(n, reps=(n, 1))
    xx = np.hstack((xx, x[:, np.newaxis], w[:, np.newaxis]))
    data = np.vstack((data, xx))

# Store list of polynomial degrees as otherwise Fortran code does not know
# array size of actual data.
fn = join(outdir, 'hermegauss_n.csv')
fmt = '%3.0f'
np.savetxt(fn, degrees, fmt=fmt, delimiter='')

# Store quadrature nodes and weights
fn = join(outdir, 'hermegauss_data.csv')
fmt = ['%3.0f', '%50.30f', '%50.30f']
np.savetxt(fn, data, fmt=fmt, delimiter='')

