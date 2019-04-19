
# Changes introduced in NUMFORT #

1.  Completely remove calls to the `timer` routine from `lbfgsb.f`. 
    These calls lead to extremely poor performance in multi-threaded code
    using OpenMP.
    Consequently, the elements 7-9 of the `DSAVE` array in the main routine
    `SETULB` are not updated to contain timing information.
    
2.  The `SETULB` routine and the routines in `linpack.f` were included into
    a module file.
    
    
# License #

All changes introduced in NUMFORT are licensed under the same terms as the
original code.


