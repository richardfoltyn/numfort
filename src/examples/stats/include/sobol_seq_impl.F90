
subroutine __APPEND2(test,__INTSIZE,__PREC) ()

    integer, parameter :: LENGTH = 10

    call __APPEND2(test_dim,__INTSIZE,__PREC) (1, LENGTH)
    call __APPEND2(test_dim,__INTSIZE,__PREC) (2, LENGTH)
    call __APPEND2(test_dim,__INTSIZE,__PREC) (3, LENGTH)

end subroutine

subroutine __APPEND2(test_dim,__INTSIZE,__PREC) (ndim, length)

    integer, parameter :: PREC = __PREC
    integer, intent(in) :: ndim, length

    integer :: i
    type (sobol_state) :: state
    real (PREC), dimension(:,:), allocatable :: sobol

    call sobol_init (state, ndim=ndim)

    allocate (sobol(ndim, length))
    do i = 1, length
        if (ndim == 1) then
            call sobol_next (state, sobol(1, i))
        else
            call sobol_next (state, sobol(:, i))
        end if
    end do

    call __APPEND2(print,__INTSIZE,__PREC) (sobol)

end subroutine

subroutine __APPEND2(print,__INTSIZE,__PREC) (x)
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:) :: x

    integer :: i

    print '(tr1, a, i1, a, i1, a)', 'Sobol sequence [integer(', __INTSIZE, '), real(', __PREC, ')]:'
    do i = 1, size(x, 2)
        print '(tr2, i3, tr2, *(f10.8, :, ", "))', i, x(:, i)
    end do
end subroutine
