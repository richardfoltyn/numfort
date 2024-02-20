integer (hid_t), intent(in) :: loc_id
character (*), intent(in) :: name
logical, intent(in), optional :: ignore_missing
integer, intent(out), optional :: status

integer, dimension(NDIM) :: shp, shp_val
integer :: ndim_present
integer :: lstatus, status_ignore
logical :: exists, lignore_missing

lstatus = STATUS_OK

lignore_missing = .false.
if (present(ignore_missing)) lignore_missing = ignore_missing

call h5lexists_f (loc_id, name, exists, status_ignore)
if (exists) then
    ! HDF5_DIMS will return with error status code if SHP is not large
    ! enough to store array shape
    call hdf5_dims (loc_id, name, shp, ndim_present, lstatus)
    if (lstatus /= STATUS_OK) goto 100

    if (ndim_present /= NDIM) then
        lstatus = STATUS_FAILURE
        goto 100
    end if

    if (allocated(val)) then
        shp_val = shape(val)
        if (any(shp_val /= shp)) deallocate (val)
    end if

    if (.not. allocated(val)) then
        allocate (val(shp(1)))
    end if

else if (.not. lignore_missing) then
    lstatus = STATUS_FAILURE
end if


100 continue
    if (present(status)) status = lstatus
