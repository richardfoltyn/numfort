! type-independent implementation for storing a scalar dataset

integer (hid_t), intent(in) :: loc_id
character (*), intent(in) :: name
logical, intent(in), optional :: deflate
    !!  If present and true, compress dataset using GZIP
    !!  Note: currently ignored
integer, intent(out), optional :: status

! gfortran required val to have target attribute in order to pass it to
! C_LOC
target :: val

integer :: lstatus
type (c_ptr) :: ptr
integer, parameter :: DIRECTION = 1
! Degenerate shape array for scalar data
integer, parameter, dimension(0) :: SHAPE_SCALAR = [ integer :: ]

integer (hid_t) :: dtype_file_id, dtype_native_id

call hdf5_file_dtype (val, dtype_file_id)
call H5Tget_native_type_f (dtype_file_id, DIRECTION, dtype_native_id, lstatus)
if (lstatus /= STATUS_OK) goto 100

ptr = c_loc(val)

call hdf5_store (loc_id, name, ptr, SHAPE_SCALAR, dtype_file_id, &
    dtype_native_id, status=lstatus)

100 continue
    if (present(status)) status = lstatus
