real (PREC), intent(in), dimension(:) :: haystack
real (PREC), intent(in) :: needle
    
integer (INTSIZE) :: res, n
    
n = size(haystack)
    
if (n < 2) then
    res = -1
else if (needle <= haystack(1)) then
    res = 1
else if (needle >= haystack(n-1)) then
    res = n-1
else
    res = bsearch (needle, haystack)
end if