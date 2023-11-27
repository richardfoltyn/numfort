

do k = 1, size(SIZES)
    n = SIZES(k)

    allocate (rarr(n), arr(n), arr_sorted(n))
    allocate (iorder(n), source=0_INTSIZE)

    call random_number (rarr)

    arr(:) = rarr * 1000.0

    ! ARGSORT returns rank array where the value at position i contains the index
    ! the element which would be at position i in a sorted array.
    call argsort (arr, iorder, status=status)

    forall (i=1:n) arr_sorted(i) = arr(iorder(i))

    values_ok = all (arr_sorted(1:n-1) <= arr_sorted(2:n))

    call tc%assert_true (values_ok .and. status == NF_STATUS_OK, &
        'Array size: ' // str(n))

    deallocate (rarr, arr, arr_sorted, iorder)
end do