#!/bin/sed -rf

# Convert FORTRAN 77 code to FORTRAN 90.
# Author: Richard Foltyn


s/\bdouble precision\b/real (PREC)/i
s/\bd?max1\b/max/ig
s/\bd?min1\b/min/ig
s/\bd?abs\b/abs/ig

