
    
integer :: lb, ub, midx
    
lb = 1
ub = size(haystack)
    
do while (ub > (lb + 1))
    midx = (ub + lb) / 2
    if (haystack(midx) > needle) then
        ub = midx
    else
        lb = midx
    end if
end do
