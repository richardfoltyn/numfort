
integer (PREC), intent(in) :: n
    !*  Number set elements
integer (PREC), intent(in) :: k
    !*  Number of elements taken
logical, intent(in), optional :: repetition
    !*  If present and true, the number of combinations with repetition is
    !   computed (default: false)
integer (PREC) :: res
    !*  Number of combinations

integer (PREC) :: i

logical :: lrep
lrep = .false.

if (present(repetition)) lrep = repetition

if (lrep) then
    ! combination with repetition:
    ! this is (n + k - 1)!/(k! (n-1)!)
    ! First compute (n + k - 1)!/(n-1)! = (n+k-1) * ... * n
    res = 1
    do i = n+k-1, n, -1
        res = res * i
    end do
else
    ! combination without repetition:
    ! this is n-choose-k, ie \binomial{n}{k}
    ! compute n! / (n-k)! = n * (n-1) * ... * (n-k+1)
    res = 1
    do i = n, n-k+1, -1
        res = res * i
    end do
end if

! Adjust for k! ways to order k elements
res = res / factorial (k)
