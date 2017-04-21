Direction numbers for Sobol sequences
=====================================

We use the direction numbers suggested by Joe and Kuo:
    1.  S. Joe and F. Y. Kuo, Remark on Algorithm 659: Implementing Sobol's
        quasirandom sequence generator, ACM Trans. Math. Softw. 29, 49-57 (2003)
    2.  S. Joe and F. Y. Kuo, Constructing Sobol sequences with better
        two-dimensional projections, SIAM J. Sci. Comput. 30, 2635-2654 (2008)

These are provided on the authors'
[website](http://web.maths.unsw.edu.au/%7Efkuo/sobol/) and are stored
in the file `new-joe-kuo-6.21201.txt` for up to 21201 dimensions.

The Python program `process_direction.py` uses this input file to create
a Fortran file that makes these numbers available via the module
`sobol_direction_num`.
