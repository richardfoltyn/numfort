"""
Create diagnostic data for PDF, CDF and PPF values for the log-normal
distribution evaluated and selected location/scale parameters, at various
points in the function domain.

This data is used to create Fortran source code that is included in Fortran
unit test code to check the correctness of the Fortran implementation.

Author: Richard Foltyn
"""

from scipy.stats import lognorm

import numpy as np
import os
from os.path import join

from common.env import env_setup
from common.fortran_io import write_fortran_array


def main(env):

    sigma_all = np.array([1.0, 0.1, 10.0])
    mu_all = np.array([0.0, -1.0, 2.0])

    xx = np.linspace(1.0e-10, 20, 20)
    # Do not compute quantile for CDF=1.0 as this will return inf, which
    # cannot be represented in a compile-time fortran array.
    qq = np.linspace(0.0, 1.0-1.0e-6, 20)

    nx = len(xx)
    pdf = np.zeros((len(sigma_all), len(mu_all), nx))
    cdf = np.zeros_like(pdf)
    ppf = np.zeros_like(pdf)

    for i, sigma in enumerate(sigma_all):
        s = sigma
        for j, mu in enumerate(mu_all):
            scale = np.exp(mu)

            fx = lognorm.pdf(xx, s=s, scale=scale)
            pdf[i, j] = fx

            fx = lognorm.cdf(xx, s=s, scale=scale)
            cdf[i, j] = fx

            fx = lognorm.ppf(qq, s=s, scale=scale)
            ppf[i, j] = fx

    destdir = join(env.fortran_root, 'tests', 'stats', 'include')
    if not os.path.exists(destdir):
        os.mkdir(destdir)

    fn = join(destdir, 'lognorm_data.F90')
    with open(fn, 'wt') as f:
        name = 'SCALE_PARAMS'
        write_fortran_array(sigma_all, name, f)

        name = 'LOC_PARAMS'
        write_fortran_array(mu_all, name, f)

        name = 'XDATA'
        write_fortran_array(xx, name, f)

        name = 'QUANTILE_DATA'
        write_fortran_array(qq, name, f)

        name = 'SCIPY_PDF_DATA'
        write_fortran_array(pdf, name, f)

        name = 'SCIPY_CDF_DATA'
        write_fortran_array(cdf, name, f)

        name = 'SCIPY_PPF_DATA'
        write_fortran_array(ppf, name, f)


if __name__ == '__main__':
    x = env_setup()
    main(x)


