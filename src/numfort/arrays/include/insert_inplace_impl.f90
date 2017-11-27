
integer :: ni, nv, ifrom, ito, i, n

ni = size(arr)
nv = size(val)

if (idx < 1 .or. idx > (ni+1) .or. nv == 0) return

! Copy over elements in ARR that lie before data to be inserted:
! Number of elements if ARR that will come to lie AFTER inserted sequence
n = max((ni-idx+1) - nv, 0)
! Shift at most NV elements in ARR to the right
forall (i=0:n-1) arr(ni-i) = arr(idx+n-1-i)

! Insert elements from VAL
ifrom = idx
ito = min(idx+nv-1, ni)
arr(ifrom:ito) = val(ifrom-idx+1:ito-idx+1)
