

pure function __APPEND(norm_sup_diff_3d,__PREC) (arr1, arr2) result(res)
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:,:) :: arr1
    real (PREC), intent(in), dimension(:,:,:) :: arr2
    real (PREC) :: res

    integer :: i1, i2, i3

    res = - 1.0_PREC
    if (any(shape(arr1) /= shape(arr2))) return

    do i3 = 1, size(arr1, 3)
        do i2 = 1, size(arr1, 2)
            do i1 = 1, size(arr1, 1)
                res = max(res, abs(arr1(i1,i2,i3)-arr2(i1,i2,i3)))
            end do
        end do
    end do

end function


pure function __APPEND(norm_sup_diff_4d,__PREC) (arr1, arr2) result(res)
    integer, parameter :: PREC = __PREC

    real (PREC), intent(in), dimension(:,:,:,:) :: arr1
    real (PREC), intent(in), dimension(:,:,:,:) :: arr2
    real (PREC) :: res

    integer :: i1, i2, i3, i4
    
    res = - 1.0_PREC
    if (any(shape(arr1) /= shape(arr2))) return

    do i4 = 1, size(arr1, 4)
        do i3 = 1, size(arr1, 3)
            do i2 = 1, size(arr1, 2)
                do i1 = 1, size(arr1, 1)
                    res = max(res, abs(arr1(i1,i2,i3,i4)-arr2(i1,i2,i3,i4)))
                end do
            end do
        end do
    end do

end function
