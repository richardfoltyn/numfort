"""
SLSQP example using the Scipy wrapper around the Fortran library.
Used to check Fortran wrapper results vs. Scipy wrapper.

Uses Rosenbrock function with some nonlinear inequality constraints,
see: http://se.mathworks.com/help/optim/ug/example-nonlinear-constrained-minimization.html
"""

import numpy as np
from scipy.optimize import minimize


def fobj(x):
    x1, x2 = x
    fx = 100.0 * (x2-x1**2.0)**2.0 + (1.0-x1)**2.0
    return fx


def fobj_grad(x):
    x1, x2 = x
    dfdx1 = - 200.0 * (x2-x1**2.0) * 2.0 * x1 - 2.0 * (1.0-x1)
    dfdx2 = 200.0 * (x2-x1**2.0)
    return np.array([dfdx1, dfdx2])


def f_ieqcons(x):
    """
    Returns evaluated constraint such that C(x) >= 0 holds.
    """
    x1, x2 = x
    cx = -x1**2.0 - x2**2.0 + 1
    return cx


def f_ieqcons_grad(x):
    x1, x2 = x
    dx1 = -2.0 * x1
    dx2 = -2.0 * x2
    dx = np.array([dx1, dx2]).reshape((1, 2))
    return dx


def f_eqcons(x):
    x1, x2 = x
    cx = -x1**2.0 - x2**3.0 + 1.0
    return cx


def f_eqcons_grad(x):
    x1, x2 = x
    return np.array([-2.0*x1, - 3.0*x2**2.0]).reshape((1, 2))


def main():

    # Example 1 with inequality constraint
    x0 = np.array([0.1, 0.1])
    bounds = ((-1.0, 1.0), (-1.0, 1.0))
    constr = {'type': 'ineq', 'fun': f_ieqcons, 'jac': f_ieqcons_grad}
    res = minimize(fobj, x0, method='SLSQP', jac=fobj_grad, bounds=bounds,
                   constraints=constr, tol=1e-8)

    s = 'Min: {}; f(xmin)={:.8f}'.format(res.x, res.fun)
    print(s)

    s = 'Constr. value at xmin: {}'.format(f_ieqcons(res.x))
    print(s)

    # Example 2 with equality constraint
    x0 = np.array([0.1, 0.1])
    bounds = ((-1.0, 1.0), (-1.0, 1.0))
    constr = {'type': 'eq', 'fun': f_eqcons, 'jac': f_eqcons_grad}
    res = minimize(fobj, x0, method='SLSQP', jac=fobj_grad, bounds=bounds,
                   constraints=constr, tol=1e-8)

    s = 'Min: {}; f(xmin)={:.8f}'.format(res.x, res.fun)
    print(s)

    s = 'Constr. value at xmin: {}'.format(f_eqcons(res.x))
    print(s)


if __name__ == '__main__':
    main()
