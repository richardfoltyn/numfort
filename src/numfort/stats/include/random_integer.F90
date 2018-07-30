
integer, parameter :: PREC = real64
integer (INTSIZE) :: n, i
real (PREC), dimension(:), allocatable :: rnd
real (PREC) :: rlow, rhigh, rx, range

range = high - low + 1

rlow = real(low, PREC)
rhigh = real(high, PREC)
range = rhigh - rlow + 1.0_PREC

n = size(x)

allocate (rnd(n))
call random_number (rnd)

do i = 1, n
    rx = rnd(i) * range
    x(i) = int(rx, INTSIZE)
end do

x(:) = x + low
