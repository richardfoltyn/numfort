integer (hid_t), intent(in) :: loc_id
character (*), intent(in) :: name
logical, intent(in), optional :: deflate
    !!  If present and true, compress dataset using GZIP
integer, intent(out), optional :: status

integer, dimension(NDIM) :: shp
integer :: lstatus
type (c_ptr) :: ptr
integer, parameter :: DIRECTION = 1
integer :: ssize

integer (hid_t) :: dtype_file_id, dtype_native_id

call hdf5_file_dtype (val_kind, dtype_file_id)
call H5Tget_native_type_f (dtype_file_id, DIRECTION, dtype_native_id, lstatus)
if (lstatus /= STATUS_OK) goto 100

ptr = c_loc(val)

! storage size in bits
ssize = storage_size (val_kind)
shp = shape(val)

call hdf5_store (loc_id, name, ptr, shp, dtype_file_id, dtype_native_id, &
    ssize, deflate, lstatus)

100 continue
    if (present(status)) status = lstatus
