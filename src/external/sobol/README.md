# Sobol sequence generator

This page contains the primitive polynomials and various sets of initial 
direction numbers for generating Sobol sequences.

This is a joint project between Stephen Joe and Frances Kuo. More details can 
be found in the following papers:

- S. Joe and F. Y. Kuo, Remark on Algorithm 659: Implementing Sobol's quasirandom sequence generator, ACM Trans. Math. Softw. 29, 49-57 (2003). 
  [Link to paper](http://doi.acm.org/10.1145/641876.641879).
- S. Joe and F. Y. Kuo, Constructing Sobol sequences with better two-dimensional projections, SIAM J. Sci. Comput. 30, 2635-2654 (2008). 
  [Link to paper](http://dx.doi.org/10.1137/070709359).

Here is a [https://web.maths.unsw.edu.au/%7Efkuo/sobol/joe-kuo-notes.pdf](3-page) 
notes on generating Sobol sequences. 

# Additions by Richard Foltyn

The Python program `process_direction.py` uses this input file to create
a Fortran file that makes these numbers available via the module
`sobol_direction_num`.

# License

This program and the accompanying direction numbers above are covered by this 
BSD-style license (see LICENSE.txt).