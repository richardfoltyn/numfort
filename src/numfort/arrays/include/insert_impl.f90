
integer :: ni, nv, no, ifrom, ito

ni = size(arr)
no = size(out)
nv = size(val)

if (idx < 1 .or. idx > (ni+1)) return

! If there is nothing to insert, copy max. feasible amount of elements from
! ARR to OUT and exit.
if (nv == 0) then
    ifrom = 1
    ito = min(ni, no)
    out(ifrom:ito) = arr(ifrom:ito)
    return
end if

! Copy over elements in ARR that lie before data to be inserted
ito = min(no, idx-1)
out(1:ito) = arr(1:ito)

! Insert elements from VAL
ifrom = ito+1
ito = min(ifrom+nv-1, no)
out(ifrom:ito) = val(1:ito-ifrom+1)

! Insert remaining elements from ARR
ifrom = ito + 1
ito = min(ifrom + (ni-idx+1), no)
out(ifrom:ito) = arr(idx+1:idx+ito-ifrom+1)
