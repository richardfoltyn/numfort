integer (hid_t), intent(in) :: loc_id
character (*), intent(in) :: name
logical, intent(in), optional :: ignore_missing
integer, intent(out), optional :: status

integer, dimension(NDIM) :: shp
integer, parameter :: DIRECTION = 1
type (c_ptr) :: ptr
integer :: lstatus
integer (hid_t) :: dtype_file_id, dtype_native_id

call hdf5_file_dtype (val_kind, dtype_file_id)
call H5Tget_native_type_f (dtype_file_id, DIRECTION, dtype_native_id, lstatus)
if (lstatus /= STATUS_OK) goto 100

call hdf5_alloc (loc_id, name, val, ignore_missing, lstatus)
if (lstatus /= STATUS_OK) goto 100

if (allocated(val)) then
    ptr = c_loc(val)
    shp = shape(val)

    call hdf5_load (loc_id, name, ptr, shp, dtype_native_id, ignore_missing, lstatus)
end if

100 continue
    if (present(status)) status = lstatus
