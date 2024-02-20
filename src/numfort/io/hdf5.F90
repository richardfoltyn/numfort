

module numfort_hdf5

    use hdf5
    use, intrinsic :: iso_fortran_env
    use, intrinsic :: iso_c_binding

    implicit none
    private

    public :: hdf5_store
    public :: hdf5_load
    public :: hdf5_load_alloc
    public :: hdf5_dims

    integer, public, parameter :: STATUS_OK = 0
    ! HDF routines seem to return this as failure code
    integer, public, parameter :: STATUS_FAILURE = -1

    ! Default GZIP compression level
    integer, parameter :: GZIP_DEFAULT_LEVEL = 5
    ! max. chunk size for compression, in bytes
    integer, parameter :: COMPRESS_MAX_CHUNK_SIZE = 2 ** 20

    ! --------------------------------------------------------------------------
    ! Interfaces for storing HDF5 datasets

    ! Routines for storing scalars
    interface hdf5_store
        procedure hdf5_store_scalar_real64, hdf5_store_scalar_real32, &
            hdf5_store_scalar_int32, hdf5_store_scalar_bool, hdf5_store_scalar_char
    end interface

    ! Routine for storing data pointed to by c_ptr
    interface hdf5_store
        procedure hdf5_store_cptr
    end interface

    ! Routines for storing 1d arrays
    interface hdf5_store
        procedure hdf5_store_1d_real32, hdf5_store_1d_real64, &
            hdf5_store_1d_int32, hdf5_store_1d_logical, hdf5_store_1d_char
    end interface

    ! Routines for storing 2d arrays
    interface hdf5_store
        procedure hdf5_store_2d_real32, hdf5_store_2d_real64, &
            hdf5_store_2d_int32
    end interface

    ! Routines for storing 3d arrays
    interface hdf5_store
        procedure hdf5_store_3d_real32, hdf5_store_3d_real64, &
            hdf5_store_3d_int32, hdf5_store_3d_logical
    end interface

    ! Routines for storing 4d arrays
    interface hdf5_store
        procedure hdf5_store_4d_real32, hdf5_store_4d_real64, &
            hdf5_store_4d_int32
    end interface

    ! Routines for storing 5d arrays
    interface hdf5_store
        procedure hdf5_store_5d_real32, hdf5_store_5d_real64, &
            hdf5_store_5d_int8, hdf5_store_5d_int32
    end interface

    ! Routines for storing 6d arrays
    interface hdf5_store
        procedure hdf5_store_6d_real32, hdf5_store_6d_real64, &
            hdf5_store_6d_int8, hdf5_store_6d_int32
    end interface

    ! Routines for storing 7d arrays
    interface hdf5_store
        procedure hdf5_store_7d_real32, hdf5_store_7d_real64, &
            hdf5_store_7d_int8, hdf5_store_7d_int32
    end interface

    interface hdf5_store
        procedure hdf5_store_scalar_int8, hdf5_store_1d_int8, &
            hdf5_store_2d_int8, hdf5_store_3d_int8, hdf5_store_4d_int8
    end interface

    interface hdf5_store
        procedure hdf5_store_scalar_int64, hdf5_store_1d_int64, &
            hdf5_store_2d_int64, hdf5_store_3d_int64, hdf5_store_4d_int64, &
            hdf5_store_5d_int64, hdf5_store_6d_int64, hdf5_store_7d_int64
    end interface

    ! --------------------------------------------------------------------------
    ! Interfaces for loading HDF5 datasets

    ! Routines for loading scalars
    interface hdf5_load
        procedure hdf5_load_scalar_real64, hdf5_load_scalar_real32, &
            hdf5_load_scalar_int32, hdf5_load_scalar_bool
    end interface

    ! Routine for loading data pointed to by c_ptr
    interface hdf5_load
        procedure hdf5_load_cptr
    end interface

    ! Routines for loading 1d arrays
    interface hdf5_load
        procedure hdf5_load_1d_real32, hdf5_load_1d_real64, &
            hdf5_load_1d_int32, hdf5_load_1d_logical
    end interface

    ! Routines for loading 2d arrays
    interface hdf5_load
        procedure hdf5_load_2d_real32, hdf5_load_2d_real64, &
            hdf5_load_2d_int32
    end interface

    ! Routines for loading 3d arrays
    interface hdf5_load
        procedure hdf5_load_3d_real32, hdf5_load_3d_real64, &
            hdf5_load_3d_int32, hdf5_load_3d_logical
    end interface

    ! Routines for loading 4d arrays
    interface hdf5_load
        procedure hdf5_load_4d_real32, hdf5_load_4d_real64, &
            hdf5_load_4d_int32
    end interface

    ! Routines for loading 5d arrays
    interface hdf5_load
        procedure hdf5_load_5d_real32, hdf5_load_5d_real64, &
            hdf5_load_5d_int8, hdf5_load_5d_int32
    end interface

    ! Routines for loading 6d arrays
    interface hdf5_load
        procedure hdf5_load_6d_real32, hdf5_load_6d_real64, &
            hdf5_load_6d_int8, hdf5_load_6d_int32
    end interface

    ! Routines for loading 7d arrays
    interface hdf5_load
        procedure hdf5_load_7d_real32, hdf5_load_7d_real64, &
            hdf5_load_7d_int8, hdf5_load_7d_int32
    end interface

    interface hdf5_load
        procedure hdf5_load_scalar_int8, hdf5_load_1d_int8, &
            hdf5_load_2d_int8, hdf5_load_3d_int8, hdf5_load_4d_int8
    end interface

    interface hdf5_load
        procedure hdf5_load_scalar_int64, hdf5_load_1d_int64, &
            hdf5_load_2d_int64, hdf5_load_3d_int64, hdf5_load_4d_int64, &
            hdf5_load_5d_int64, hdf5_load_6d_int64, hdf5_load_7d_int64
    end interface

    ! --------------------------------------------------------------------------
    ! Interfaces for HDF5_LOAD_ALLOC

    interface hdf5_load_alloc
        procedure hdf5_load_alloc_1d_int8, hdf5_load_alloc_1d_int32, &
            hdf5_load_alloc_1d_real32, hdf5_load_alloc_1d_real64, &
            hdf5_load_alloc_1d_logical
    end interface

    interface hdf5_load_alloc
        procedure hdf5_load_alloc_2d_int8, hdf5_load_alloc_2d_int32, &
            hdf5_load_alloc_2d_real32, hdf5_load_alloc_2d_real64, &
            hdf5_load_alloc_2d_logical
    end interface

    interface hdf5_load_alloc
        procedure hdf5_load_alloc_3d_int8, hdf5_load_alloc_3d_int32, &
            hdf5_load_alloc_3d_real32, hdf5_load_alloc_3d_real64, &
            hdf5_load_alloc_3d_logical
    end interface

    interface hdf5_load_alloc
        procedure hdf5_load_alloc_4d_int8, hdf5_load_alloc_4d_int32, &
            hdf5_load_alloc_4d_real32, hdf5_load_alloc_4d_real64, &
            hdf5_load_alloc_4d_logical
    end interface

    interface hdf5_load_alloc
        procedure hdf5_load_alloc_5d_int8, hdf5_load_alloc_5d_int32, &
            hdf5_load_alloc_5d_real32, hdf5_load_alloc_5d_real64, &
            hdf5_load_alloc_5d_logical
    end interface

    interface hdf5_load_alloc
        procedure hdf5_load_alloc_6d_int8, hdf5_load_alloc_6d_int32, &
                hdf5_load_alloc_6d_real32, hdf5_load_alloc_6d_real64, &
                hdf5_load_alloc_6d_logical
    end interface

    ! --------------------------------------------------------------------------
    ! Interfaces for HDF5_ALLOC

    interface hdf5_alloc
        procedure hdf5_alloc_1d_int8, hdf5_alloc_1d_int32, hdf5_alloc_1d_real32, &
            hdf5_alloc_1d_real64
    end interface

    interface hdf5_alloc
        procedure hdf5_alloc_2d_int8, hdf5_alloc_2d_int32, hdf5_alloc_2d_real32, &
            hdf5_alloc_2d_real64
    end interface

    interface hdf5_alloc
        procedure hdf5_alloc_3d_int8, hdf5_alloc_3d_int32, hdf5_alloc_3d_real32, &
            hdf5_alloc_3d_real64
    end interface

    interface hdf5_alloc
        procedure hdf5_alloc_4d_int8, hdf5_alloc_4d_int32, hdf5_alloc_4d_real32, &
            hdf5_alloc_4d_real64
    end interface

    interface hdf5_alloc
        procedure hdf5_alloc_5d_int8, hdf5_alloc_5d_int32, hdf5_alloc_5d_real32, &
            hdf5_alloc_5d_real64
    end interface

    interface hdf5_alloc
        procedure hdf5_alloc_6d_int8, hdf5_alloc_6d_int32, hdf5_alloc_6d_real32, &
                hdf5_alloc_6d_real64
    end interface

    ! --------------------------------------------------------------------------
    ! Additional interfaces

    interface hdf5_file_dtype
        procedure hdf5_file_dtype_real32, hdf5_file_dtype_real64, &
            hdf5_file_dtype_int8, hdf5_file_dtype_int32, hdf5_file_dtype_int64
    end interface

    interface hdf5_dims
        procedure hdf5_dims_int32, hdf5_dims_int64
    end interface

contains


! ------------------------------------------------------------------------------
! Data dimension query

subroutine hdf5_dims_int64 (loc_id, name, dims, ndims, status)
    integer, parameter :: INTSIZE = int64
    integer (hid_t), intent(in) :: loc_id
        !!  Location ID of parent object which contains dataset
    character (*), intent(in) :: name
        !!  Dataset name
    integer (INTSIZE), intent(out), dimension(:), optional :: dims
        !*  On exit, contains
    integer, intent(out), optional :: ndims
        !*  If present, contains number of dimensions on exit.
    integer, intent(out), optional :: status
        !!  Status flag

    integer (hsize_t), dimension(:), allocatable  :: ldims, maxdims
    integer (hid_t) :: dspace_id, dset_id
    integer :: lstatus, rnk, lstatus_ignore

    if (present(ndims)) ndims = -1

    call h5dopen_f (loc_id, name, dset_id, lstatus)
    if (lstatus /= STATUS_OK) goto 20

    call h5dget_space_f (dset_id, dspace_id, lstatus)
    if (lstatus /= STATUS_OK) goto 10

    ! retrieve array rank
    call h5sget_simple_extent_ndims_f (dspace_id, rnk, lstatus)
    if (lstatus /= STATUS_OK) goto 10

    if (present(ndims)) ndims = rnk
    if (.not. present(dims)) goto 10

    ! make sure that stored array has the same rank as array into which it
    ! should be read
    if (rnk > size(dims)) then
        lstatus = STATUS_FAILURE
        goto 10
    end if

    allocate (ldims(rnk), maxdims(rnk))

    call h5sget_simple_extent_dims_f (dspace_id, ldims, maxdims, lstatus)
    ! above routine returns the rank of data space on successful exit and -1 otherwise
    if (lstatus /= rnk) then
        goto 10
        ldims(:) = -1
        lstatus = STATUS_FAILURE
    end if

    lstatus = STATUS_OK

10  call h5sclose_f (dspace_id, lstatus_ignore)

20  call h5dclose_f (dset_id, lstatus_ignore)

    if (present(dims)) then
        dims = -1
        if (size(dims) >= rnk) then
            dims(1:rnk) = int(ldims(1:rnk), INTSIZE)
        end if
    end if

    if (present(status)) status = lstatus

end subroutine

subroutine hdf5_dims_int32 (loc_id, name, dims, ndims, status)
    integer, parameter :: INTSIZE = int32
    integer (hid_t), intent(in) :: loc_id
        !*  Location ID of parent object which contains dataset
    character (*), intent(in) :: name
        !*  Dataset name
    integer (INTSIZE), intent(out), dimension(:) :: dims
        !*  On exit, contains
    integer, intent(out), optional :: ndims
        !*  If present, contains number of dimensions on exit.
    integer, intent(out), optional :: status
        !*  Status flag

    integer (int64), dimension(:), allocatable :: ldims

    allocate (ldims(size(dims)))
    call hdf5_dims (loc_id, name, ldims, ndims, status)
    dims = int(ldims, INTSIZE)

end subroutine

! ------------------------------------------------------------------------------
! HDF5_DTYPE overloads

pure subroutine hdf5_file_dtype_real64 (val, dtype_id)
    real (real64), intent(in) :: val
    integer (hid_t), intent(out) :: dtype_id
    dtype_id = H5T_IEEE_F64LE
end subroutine

pure subroutine hdf5_file_dtype_real32 (val, dtype_id)
    real (real32), intent(in) :: val
    integer (hid_t), intent(out) :: dtype_id
    dtype_id = H5T_IEEE_F32LE
end subroutine

pure subroutine hdf5_file_dtype_int8 (val, dtype_id)
    integer (int8), intent(in) :: val
    integer (hid_t), intent(out) :: dtype_id
    dtype_id = H5T_STD_I8LE
end subroutine

pure subroutine hdf5_file_dtype_int32 (val, dtype_id)
    integer (int32), intent(in) :: val
    integer (hid_t), intent(out) :: dtype_id
    dtype_id = H5T_STD_I32LE
end subroutine

subroutine hdf5_file_dtype_int64 (val, dtype_id)
    integer (int64), intent(in) :: val
    integer (hid_t), intent(out) :: dtype_id
    dtype_id = H5T_STD_I64LE
end subroutine

! ------------------------------------------------------------------------------
! HDF_STORE overloads for scalar arguments

subroutine hdf5_store_scalar_real64 (loc_id, name, val, deflate, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(in) :: val

#include "hdf5/hdf5_store_scalar_impl.f90"
end subroutine

subroutine hdf5_store_scalar_real32 (loc_id, name, val, deflate, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(in) :: val

#include "hdf5/hdf5_store_scalar_impl.f90"
end subroutine


subroutine hdf5_store_scalar_int8 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in) :: val

#include "hdf5/hdf5_store_scalar_impl.f90"
end subroutine

subroutine hdf5_store_scalar_int32 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in) :: val

#include "hdf5/hdf5_store_scalar_impl.f90"
end subroutine

subroutine hdf5_store_scalar_int64 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in) :: val

#include "hdf5/hdf5_store_scalar_impl.f90"
end subroutine

subroutine hdf5_store_scalar_bool (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int8
    integer (hid_t), intent(in) :: loc_id
    character (*), intent(in) :: name
    logical, intent(in) :: val
    logical, intent(in), optional :: deflate
        !!  If present and true, compress dataset using GZIP
        !!  Note: currently ignored
    integer, intent(out), optional :: status

    integer (INTSIZE) :: val_int

    val_int = 0
    if (val) val_int = 1

    call hdf5_store (loc_id, name, val_int, deflate, status)

end subroutine

! ------------------------------------------------------------------------------
! HDF_STORE overload for C pointer to array

!>  Store array in memory pointed to by ptr to HDF5 file.
subroutine hdf5_store_cptr (loc_id, name, ptr, dims, dtype_file_id, &
        dtype_native_id, ssize, deflate, status)
    integer (hid_t), intent(in) :: loc_id
        !!  Location ID of parent object which is to contains dataset
    character (*), intent(in) :: name
        !!  Dataset name
    type (c_ptr), intent(in) :: ptr
        !!  Pointer to array in memory which contains values to be stored
    integer, dimension(:) :: dims
        !!  Shape of array pointed to by ptr.
    integer (hid_t), intent(in) :: dtype_file_id
        !!  HDF5 datatype to be used to store array in file
    integer (hid_t), intent(in) :: dtype_native_id
        !!  HDF5 datatype of array in memory
    integer, intent(in), optional :: ssize
        !!  Storate size of one array element, in bits
    logical, intent(in), optional :: deflate
        !!  If true, deflate dataset using GZIP algorithm
    integer, intent(out), optional :: status
        !!  Status flag

    integer (hsize_t), dimension(size(dims)) :: ldims, chunks
    integer (hid_t) :: dspace_id, dset_id, dprop_id
    integer :: lstatus, rnk, status_ignore
    logical :: ldeflate, has_encode, has_decode, exists

    rnk = size(dims)
    ldims = dims
    ldeflate = .false.

    ! check for deflate support
    call hdf5_deflate_info (has_encode, has_decode)
    ! ssize argument is needed for compression chunk size, so enable compression
    ! only if ssize is present.
    ! Compression is enabled by default if deflate argument is not specified,
    ! compression is supported, chuck size is present and data is non-scalar.
    ! Ignore compression for scalar data, or degenerate arrays where one
    ! dimension has size 0.
    ldeflate = has_encode .and. present(ssize) .and. rnk > 0 .and. product(dims) > 0

    ! Don't compress anything that is less than 1KiB in size, compressing
    ! smaller array seems to fail and is not worth it
    if (present(ssize)) then
        ldeflate = ldeflate .and. (ssize/8 * product(dims)) >= 1024
    end if

    ! Override default behavior if explicitly specified.
    if (present(deflate)) ldeflate = ldeflate .and. deflate

    call h5screate_simple_f (rnk, ldims, dspace_id, lstatus)
    if (lstatus /= STATUS_OK) goto 20

    if (ldeflate) then
        ! compute chunk size: store stuff in column-major order; store complete
        ! column per chunk unless that is too large
        call find_chunk_size (ldims, ssize/8, COMPRESS_MAX_CHUNK_SIZE, chunks)

        ! Create dataset property to deflate data (property of dataset creation)
        call h5pcreate_f(H5P_DATASET_CREATE_F, dprop_id, lstatus)
        ! set deflate algorithm to max. compression
        call h5pset_deflate_f (dprop_id, GZIP_DEFAULT_LEVEL, lstatus)
        if (lstatus /= STATUS_OK) goto 5
        ! specify chunk size
        call h5pset_chunk_f (dprop_id, size(chunks), chunks, lstatus)
        if (lstatus /= STATUS_OK) goto 5
    end if

    ! Check whether dataset already exists
    exists = .false.
    call h5lexists_f(loc_id, name, exists, status_ignore)
    if (exists) then
        if (ldeflate) then
            call h5dopen_f (loc_id, name, dset_id, lstatus, dprop_id)
        else
            call h5dopen_f (loc_id, name, dset_id, lstatus)
        end if
    else
        ! Create new dataset and assign it deflate properties if applicable.
        if (ldeflate) then
            ! Create dataset with property list with deflate
            call h5dcreate_f (loc_id, name, dtype_file_id, dspace_id, &
                dset_id, lstatus, dprop_id)
        else
            ! create dataset without property list, w/o deflate
            call h5dcreate_f (loc_id, name, dtype_file_id, dspace_id, &
                dset_id, lstatus)
        end if
    end if

    if (lstatus /= STATUS_OK) goto 5

    call h5dwrite_f (dset_id, dtype_native_id, ptr, lstatus)

    call h5dclose_f (dset_id, status_ignore)

    ! close and release resources
5   if (ldeflate) call h5pclose_f (dprop_id, status_ignore)

20  call h5sclose_f (dspace_id, status_ignore)

    if (present(status)) status = lstatus
end subroutine


subroutine find_chunk_size (dims, ssize, max_bytes, chunks)
    !*  FIND_CHUNK_SIZE attempts to find an "optional" chunk layout with a
    !   chunk size that is close but below MAX_BYTES.
    !   Note that the values in CHUNKS are interpreted relative to the actual
    !   array to be stored, ie one chunk contains a subsection of the actual
    !   array with dimensions given in CHUNKS.
    !
    !   Implementation details: For arrays that are in total larger than
    !   MAX_BYTES, we identify the dimension with a stride larger than
    !   MAX_BYTES (in contiguous linear memory). Then we split the preceding
    !   dimension such that the array section has size smaller than MAX_BYTES.
    integer (hsize_t), intent(in), dimension(:) :: dims
    integer, intent(in) :: ssize
        !*  Element size in bytes
    integer, intent(in) :: max_bytes
    integer (hsize_t), intent(out), dimension(:) :: chunks

    integer (hsize_t), dimension(:), allocatable :: strides
    integer (hsize_t), parameter :: one = 1
    integer (hsize_t) :: ndim, i

    ! Return instantly if the whole array fits into a single chunk smaller
    ! than MAX_BYTES
    if (product(dims) * ssize <= max_bytes) then
        chunks = dims
        return
    end if

    ! At this point, array size is larger than MAX_BYTES
    ndim = size(dims)

    allocate (strides(ndim))
    strides(1) = 1
    do i = 2, ndim
        strides(i) = strides(i-1) * dims(i-1)
    end do

    ! Find largest dimension such that an increment by 1 along that dimension
    ! is less than MAX_BYTES apart in (contiguous) memory
    do i = 1, ndim
        if (strides(i) * ssize >= max_bytes) exit
    end do

    ! Loop exits with index just one past the desired dimension (or i = NDIM+1)
    ! if no such dimension is found
    i = min(i - 1, ndim)
    chunks = dims
    chunks(i+1:) = one
    ! Halve chunk size along that dimension until chunk size is below
    ! MAX_BYTES
    do while (product(chunks) * ssize > max_bytes .and. chunks(i) > 1)
        chunks(i) = chunks(i) / 2
    end do
    chunks(i) = max(one, chunks(i))

    deallocate (strides)

end subroutine


! ------------------------------------------------------------------------------
! HDF_STORE overloads for 1d-array arguments

subroutine hdf5_store_1d_real64 (loc_id, name, val, deflate, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 1

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_1d_real32 (loc_id, name, val, deflate, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 1

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_1d_int8 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in), dimension(:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 1

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_1d_int32 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 1

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_1d_int64 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 1

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_1d_logical (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int8
    logical, intent(in), dimension(:) :: val

    integer (hid_t), intent(in) :: loc_id
    character (*), intent(in) :: name
    logical, intent(in), optional :: deflate
        !!  If present and true, compress dataset using GZIP
    integer, intent(out), optional :: status

    integer (INTSIZE), dimension(:), allocatable :: ival

    allocate (ival(size(val)))

    where (val)
        ival = 1
    else where
        ival = 0
    end where

    call hdf5_store (loc_id, name, ival, deflate, status)
    deallocate (ival)

end subroutine


subroutine hdf5_store_scalar_char (loc_id, name, val, deflate, status)
    character (*), intent(in) :: val

    integer (hid_t), intent(in) :: loc_id
    character (*), intent(in) :: name
    logical, intent(in), optional :: deflate
        !!  If present and true, compress dataset using GZIP
    integer, intent(out), optional :: status

    character (:), dimension(:), allocatable :: val1d
    integer :: n

    n = len(val)
    allocate (character (n) :: val1d(1))
    val1d(1) = val

    call hdf5_store (loc_id, name, val1d, deflate=.false., status=status)

    deallocate (val1d)

end subroutine


subroutine hdf5_store_1d_char (loc_id, name, val, deflate, status)
    integer (hid_t), intent(in) :: loc_id
    character (*), intent(in) :: name
    character (*), intent(in), dimension(:), target, contiguous :: val
    logical, intent(in), optional :: deflate
        !!  If present and true, compress dataset using GZIP
    integer, intent(out), optional :: status

    integer, parameter :: NDIM = 1
    integer, dimension(NDIM) :: shp

    type (c_ptr) :: ptr
    integer :: ssize, lstatus

    integer (hid_t) :: dtype_file_id, dtype_native_id
    integer (size_t) :: char_len

    char_len = len(val(1))

    ! Create type for this specific string length
    call H5Tcopy_f(H5T_FORTRAN_S1, dtype_file_id, lstatus)
    if (lstatus /= 0) then
        if (present(status)) status = lstatus
        goto 100
    end if
    call H5Tset_size_f(dtype_file_id, char_len, lstatus)
    if (lstatus /= 0) then
        if (present(status)) status = lstatus
        goto 100
    end if
    ! Use same data type ID for both file and memory
    dtype_native_id = dtype_file_id

    ptr = c_loc(val(1)(1:1))

    ! storage size in bits
    ssize = storage_size (val)
    shp = shape(val)

    call hdf5_store (loc_id, name, ptr, shp, dtype_file_id, &
        dtype_native_id, ssize, deflate, status)

100 continue
    call H5Tclose_f(dtype_file_id, lstatus)
end subroutine

! ------------------------------------------------------------------------------
! HDF5_STORE overloads for 2d arrays

subroutine hdf5_store_2d_real64 (loc_id, name, val, deflate, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 2

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_2d_real32 (loc_id, name, val, deflate, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 2

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_2d_int8 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in), dimension(:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 2

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_2d_int32 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 2

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_2d_int64 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 2

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

! ------------------------------------------------------------------------------
! HDF5_STORE overloads for 3d arrays

subroutine hdf5_store_3d_real64 (loc_id, name, val, deflate, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 3

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_3d_real32 (loc_id, name, val, deflate, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 3

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_3d_int8 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in), dimension(:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 3

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_3d_int32 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 3

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_3d_int64 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 3

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_3d_logical (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int8
    logical, intent(in), dimension(:,:,:), contiguous :: val

    integer (hid_t), intent(in) :: loc_id
    character (*), intent(in) :: name
    logical, intent(in), optional :: deflate
    !!  If present and true, compress dataset using GZIP
    integer, intent(out), optional :: status

    integer (INTSIZE), dimension(:,:,:), allocatable :: ival

    allocate (ival(size(val,1), size(val,2), size(val,3)))

    where (val)
        ival = 1
    else where
        ival = 0
    end where

    call hdf5_store (loc_id, name, ival, deflate, status)
    deallocate (ival)

end subroutine


! ------------------------------------------------------------------------------
! HDF5_STORE overloads for 4d arrays

subroutine hdf5_store_4d_real64 (loc_id, name, val, deflate, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:,:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 4

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_4d_real32 (loc_id, name, val, deflate, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:,:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 4

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_4d_int8 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in), dimension(:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 4

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_4d_int32 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 4

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_4d_int64 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 4

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

! ------------------------------------------------------------------------------
! HDF5_STORE overloads for 5d arrays

subroutine hdf5_store_5d_real64 (loc_id, name, val, deflate, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:,:,:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 5

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_5d_real32 (loc_id, name, val, deflate, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:,:,:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 5

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_5d_int8 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 5

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_5d_int32 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 5

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_5d_int64 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 5

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

! ------------------------------------------------------------------------------
! HDF5_STORE overloads for 6d arrays

subroutine hdf5_store_6d_real64 (loc_id, name, val, deflate, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:,:,:,:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 6

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_6d_real32 (loc_id, name, val, deflate, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:,:,:,:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 6

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_6d_int8 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 6

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_6d_int32 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 6

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_6d_int64 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 6

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

! ------------------------------------------------------------------------------
! HDF5_STORE overloads for 7d arrays

subroutine hdf5_store_7d_real64 (loc_id, name, val, deflate, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(in), dimension(:,:,:,:,:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 7

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_7d_real32 (loc_id, name, val, deflate, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(in), dimension(:,:,:,:,:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 7

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_7d_int8 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 7

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_7d_int32 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 7

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

subroutine hdf5_store_7d_int64 (loc_id, name, val, deflate, status)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(in), dimension(:,:,:,:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 7

#include "hdf5/hdf5_store_array_impl.f90"
end subroutine

! ------------------------------------------------------------------------------
! HDF_LOAD overloads for scalar arguments

subroutine hdf5_load_scalar_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out) :: val

#include "hdf5/hdf5_load_scalar_impl.f90"
end subroutine

subroutine hdf5_load_scalar_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out) :: val

#include "hdf5/hdf5_load_scalar_impl.f90"
end subroutine

subroutine hdf5_load_scalar_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out) :: val

#include "hdf5/hdf5_load_scalar_impl.f90"
end subroutine

subroutine hdf5_load_scalar_int64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(out) :: val

#include "hdf5/hdf5_load_scalar_impl.f90"
end subroutine

subroutine hdf5_load_scalar_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out) :: val

#include "hdf5/hdf5_load_scalar_impl.f90"
end subroutine

subroutine hdf5_load_scalar_bool (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (hid_t), intent(in) :: loc_id
    character (*), intent(in) :: name
    logical, intent(out) :: val
    logical, intent(in), optional :: ignore_missing
    integer, intent(out), optional :: status

    integer (INTSIZE) :: val_int
    integer :: lstatus, status_ignore
    logical :: exists, lignore_missing

    lignore_missing = .false.
    if (present(ignore_missing)) lignore_missing = ignore_missing

    lstatus = STATUS_OK

    call h5lexists_f (loc_id, name, exists, status_ignore)

    if (exists) then
        call hdf5_load (loc_id, name, val_int, ignore_missing, lstatus)
        ! Alter VAL only if no error encountered!
        if (lstatus == STATUS_OK) then
            val = (val_int == 1)
        end if
    else
        if (.not. lignore_missing) lstatus = STATUS_FAILURE
    end if

    if (present(status)) status = lstatus
end subroutine

! ------------------------------------------------------------------------------
! HDF5_LOAD for c_ptr arguments

!>  Load array of any shape from HDF5 storage into memory pointed to by
!   ptr.
subroutine hdf5_load_cptr (loc_id, name, ptr, dims, dtype_native_id, &
        ignore_missing, status)
    integer (hid_t), intent(in) :: loc_id
        !!  Location ID of parent object which contains dataset
    character (*), intent(in) :: name
        !!  Dataset name
    type (c_ptr), intent(inout) :: ptr
        !!  Pointer to Fortran array into which values should be copied
    integer, dimension(:) :: dims
        !!  Shape of array pointed to by ptr. Argument will be used to
        !!  verify that data stored in HDF5 file is of the same shape as
        !!  array in memory.
    integer (hid_t), intent(in) :: dtype_native_id
        !!  HDF5 datatype of array in memory
    logical, intent(in), optional :: ignore_missing
    integer, intent(out), optional :: status
        !!  Status flag

    integer (hsize_t), dimension(:), allocatable  :: ldims, maxdims
    integer (hid_t) :: dspace_id, dset_id
    integer :: lstatus, rnk, status_ignore
    logical :: lignore_missing, exists

    lignore_missing = .false.
    if (present(ignore_missing)) lignore_missing = ignore_missing

    lstatus = STATUS_OK

    ! Check whether dataset exists
    call h5lexists_f (loc_id, name, exists, status_ignore)
    if ((.not. exists) .and. lignore_missing) goto 30

    call h5dopen_f (loc_id, name, dset_id, lstatus)
    if (lstatus /= STATUS_OK) goto 20

    call h5dget_space_f (dset_id, dspace_id, lstatus)
    if (lstatus /= STATUS_OK) goto 10

    ! retrieve array rank
    call h5sget_simple_extent_ndims_f (dspace_id, rnk, lstatus)
    if (lstatus /= STATUS_OK) goto 10

    if (rnk > 0) then
        ! make sure that stored array has the same rank as array into which it
        ! should be read
        if (rnk /= size(dims)) then
            lstatus = STATUS_FAILURE
            goto 10
        end if

        ! retrieve shape of stored array
        allocate (ldims(rnk), maxdims(rnk))
        call h5sget_simple_extent_dims_f (dspace_id, ldims, maxdims, lstatus)
        ! above routine returns the rank of data space on successful exit and -1 otherwise
        if (lstatus /= rnk) goto 10

        ! check the the dimensions of array in memory are identical to those
        ! stored in file. Not sure if larger arrays could be permitted, since
        ! a smaller file array would then not be contiguous in memory when
        ! properly copied into a larger array.
        if (any(ldims /= dims)) then
            lstatus = STATUS_FAILURE
            goto 10
        end if
    end if

    call h5dread_f (dset_id, dtype_native_id, ptr, lstatus)

10  continue

    call h5sclose_f (dspace_id, status_ignore)

20  continue

    call h5dclose_f (dset_id, status_ignore)

30  continue

    if (present(status)) status = lstatus
end subroutine

! ------------------------------------------------------------------------------
! HDF5_LOAD overloads for 1d arrays

subroutine hdf5_load_1d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 1

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine


subroutine hdf5_load_1d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 1

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_1d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out), dimension(:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 1

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_1d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 1

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_1d_int64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(out), dimension(:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 1

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_1d_logical (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    logical, intent(out), dimension(:) :: val

    integer (hid_t), intent(in) :: loc_id
    character (*), intent(in) :: name
    logical, intent(in), optional :: ignore_missing
    integer, intent(out), optional :: status

    integer (INTSIZE), dimension(:), allocatable :: ival

    integer (int64) :: i, n
    integer :: lstatus

    n = size(val)

    allocate (ival(n))
    ! Use this as an indicator whether array was modified (left unmodified
    ! if missing and IGNORE_MISSING=.TRUE.
    ival(1) = -1

    call hdf5_load (loc_id, name, ival, ignore_missing, lstatus)

    if (lstatus == STATUS_OK .and. ival(1) /= -1) then
        forall (i=1:n) val(i) = (ival(i) /= 0)
    end if

    deallocate (ival)

    if (present(status)) status = lstatus
end subroutine

! ------------------------------------------------------------------------------
! HDF5_LOAD overloads for 2d arrays

subroutine hdf5_load_2d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 2

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_2d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 2

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_2d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out), dimension(:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 2

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_2d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 2

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_2d_int64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(out), dimension(:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 2

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

! ------------------------------------------------------------------------------
! HDF5_LOAD overloads for 3d arrays

subroutine hdf5_load_3d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 3

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_3d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 3

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_3d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out), dimension(:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 3

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_3d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 3

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_3d_int64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(out), dimension(:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 3

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_3d_logical (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    logical, intent(out), dimension(:,:,:), contiguous :: val

    integer (hid_t), intent(in) :: loc_id
    character (*), intent(in) :: name
    logical, intent(in), optional :: ignore_missing
    integer, intent(out), optional :: status

    integer (INTSIZE), dimension(:,:,:), allocatable :: ival

    integer :: lstatus

    allocate (ival(size(val,1), size(val,2), size(val,3)))

    call hdf5_load (loc_id, name, ival, ignore_missing, lstatus)

    if (lstatus == STATUS_OK) then
        val = (ival /= 0)
    end if

    deallocate (ival)

    if (present(status)) status = lstatus
end subroutine



! ------------------------------------------------------------------------------
! HDF5_LOAD overloads for 4d arrays

subroutine hdf5_load_4d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:,:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 4

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_4d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:,:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 4

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_4d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out), dimension(:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 4

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_4d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 4

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_4d_int64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(out), dimension(:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 4

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine


! ------------------------------------------------------------------------------
! HDF5_LOAD overloads for 5d arrays

subroutine hdf5_load_5d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:,:,:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 5

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_5d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:,:,:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 5

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_5d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out), dimension(:,:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 5

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_5d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:,:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 5

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_5d_int64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(out), dimension(:,:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 5

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

! ------------------------------------------------------------------------------
! HDF5_LOAD overloads for 6d arrays

subroutine hdf5_load_6d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:,:,:,:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 6

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_6d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:,:,:,:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 6

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_6d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out), dimension(:,:,:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 6

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_6d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:,:,:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 6

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_6d_int64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(out), dimension(:,:,:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 6

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

! ------------------------------------------------------------------------------
! HDF5_LOAD overloads for 7d arrays

subroutine hdf5_load_7d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:,:,:,:,:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 7

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_7d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:,:,:,:,:,:,:), target, contiguous :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 7

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_7d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out), dimension(:,:,:,:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 7

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_7d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:,:,:,:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 7

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine

subroutine hdf5_load_7d_int64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int64
    integer (INTSIZE), intent(out), dimension(:,:,:,:,:,:,:), target, contiguous :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 7

#include "hdf5/hdf5_load_array_impl.f90"
end subroutine


! ------------------------------------------------------------------------------
! HDF5_LOAD_ALLOC for 1d-arrays

subroutine hdf5_load_alloc_1d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:), allocatable, target :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 1

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_1d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:), allocatable, target :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 1

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_1d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out), dimension(:), allocatable, target :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 1

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_1d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:), allocatable, target :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 1

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_1d_logical (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    logical, intent(out), dimension(:), allocatable, target :: val

    integer (hid_t), intent(in) :: loc_id
    character (*), intent(in) :: name
    logical, intent(in), optional :: ignore_missing
    integer, intent(out), optional :: status

    integer :: lstatus
    integer (INTSIZE), dimension(:), allocatable :: ival

    call hdf5_load_alloc (loc_id, name, ival, ignore_missing, lstatus)
    if (lstatus == STATUS_OK .and. allocated(ival)) then
        if (allocated(val)) then
            if (size(val) /= size(ival)) deallocate (val)
        end if

        if (.not. allocated(val)) then
            allocate (val(size(ival)))
            val(:) = (ival == 1)
        end if
    end if

    if (present(status)) status = lstatus

end subroutine



! ------------------------------------------------------------------------------
! HDF5_LOAD_ALLOC for 2d-arrays

subroutine hdf5_load_alloc_2d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:,:), allocatable, target :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 2

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_2d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:,:), allocatable, target :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 2

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_2d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out), dimension(:,:), allocatable, target :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 2

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_2d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:,:), allocatable, target :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 2

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_2d_logical (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer, parameter :: NDIM = 2
    logical, intent(out), dimension(:,:), allocatable, target :: val

    integer (hid_t), intent(in) :: loc_id
    character (*), intent(in) :: name
    logical, intent(in), optional :: ignore_missing
    integer, intent(out), optional :: status

    integer :: lstatus
    integer, dimension(NDIM) :: shp, ishp
    integer (INTSIZE), dimension(:,:), allocatable :: ival

    call hdf5_load_alloc (loc_id, name, ival, ignore_missing, lstatus)
    if (lstatus == STATUS_OK .and. allocated(ival)) then
        if (allocated(val)) then
            shp = shape(val)
            ishp = shape(ival)
            if (all(shp == ishp)) deallocate (val)
        end if

        if (.not. allocated(val)) then
            allocate (val(size(ival,1), size(ival,2)))
            val(:,:) = (ival == 1)
        end if
    end if

    if (present(status)) status = lstatus

end subroutine

! ------------------------------------------------------------------------------
! HDF5_LOAD_ALLOC for 3d-arrays

subroutine hdf5_load_alloc_3d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:,:,:), allocatable, target :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 3

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_3d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:,:,:), allocatable, target :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 3

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_3d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out), dimension(:,:,:), allocatable, target :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 3

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_3d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:,:,:), allocatable, target :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 3

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_3d_logical (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer, parameter :: NDIM = 3
    logical, intent(out), dimension(:,:,:), allocatable, target :: val

    integer (hid_t), intent(in) :: loc_id
    character (*), intent(in) :: name
    logical, intent(in), optional :: ignore_missing
    integer, intent(out), optional :: status

    integer :: lstatus
    integer, dimension(NDIM) :: shp, ishp
    integer (INTSIZE), dimension(:,:,:), allocatable :: ival

    call hdf5_load_alloc (loc_id, name, ival, ignore_missing, lstatus)
    if (lstatus == STATUS_OK .and. allocated(ival)) then
        if (allocated(val)) then
            shp = shape(val)
            ishp = shape(ival)
            if (all(shp == ishp)) deallocate (val)
        end if

        if (.not. allocated(val)) then
            allocate (val(size(ival,1), size(ival,2), size(ival,3)))
            val(:,:,:) = (ival == 1)
        end if
    end if

    if (present(status)) status = lstatus

end subroutine

! ------------------------------------------------------------------------------
! HDF5_LOAD_ALLOC for 4d-arrays

subroutine hdf5_load_alloc_4d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:,:,:,:), allocatable, target :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 4

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_4d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:,:,:,:), allocatable, target :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 4

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_4d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out), dimension(:,:,:,:), allocatable, target :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 4

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_4d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:,:,:,:), allocatable, target :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 4

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_4d_logical (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer, parameter :: NDIM = 4
    logical, intent(out), dimension(:,:,:,:), allocatable, target :: val

    integer (hid_t), intent(in) :: loc_id
    character (*), intent(in) :: name
    logical, intent(in), optional :: ignore_missing
    integer, intent(out), optional :: status

    integer :: lstatus
    integer, dimension(NDIM) :: shp, ishp
    integer (INTSIZE), dimension(:,:,:,:), allocatable :: ival

    call hdf5_load_alloc (loc_id, name, ival, ignore_missing, lstatus)
    if (lstatus == STATUS_OK .and. allocated(ival)) then
        if (allocated(val)) then
            shp = shape(val)
            ishp = shape(ival)
            if (all(shp == ishp)) deallocate (val)
        end if

        if (.not. allocated(val)) then
            allocate (val(size(ival,1), size(ival,2), size(ival,3), size(ival,4)))
            val(:,:,:,:) = (ival == 1)
        end if
    end if

    if (present(status)) status = lstatus

end subroutine

! ------------------------------------------------------------------------------
! HDF5_LOAD_ALLOC for 5d-arrays

subroutine hdf5_load_alloc_5d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:,:,:,:,:), allocatable, target :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 5

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_5d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:,:,:,:,:), allocatable, target :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 5

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_5d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out), dimension(:,:,:,:,:), allocatable, target :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 5

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_5d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:,:,:,:,:), allocatable, target :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 5

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_5d_logical (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer, parameter :: NDIM = 5
    logical, intent(out), dimension(:,:,:,:,:), allocatable, target :: val

    integer (hid_t), intent(in) :: loc_id
    character (*), intent(in) :: name
    logical, intent(in), optional :: ignore_missing
    integer, intent(out), optional :: status

    integer :: lstatus
    integer, dimension(NDIM) :: shp, ishp
    integer (INTSIZE), dimension(:,:,:,:,:), allocatable :: ival

    call hdf5_load_alloc (loc_id, name, ival, ignore_missing, lstatus)
    if (lstatus == STATUS_OK .and. allocated(ival)) then
        if (allocated(val)) then
            shp = shape(val)
            ishp = shape(ival)
            if (all(shp == ishp)) deallocate (val)
        end if

        if (.not. allocated(val)) then
            allocate (val(size(ival,1), size(ival,2), size(ival,3), size(ival,4), size(ival,5)))
            val(:,:,:,:,:) = (ival == 1)
        end if
    end if

    if (present(status)) status = lstatus

end subroutine


! ------------------------------------------------------------------------------
! HDF5_LOAD_ALLOC for 6d-arrays

subroutine hdf5_load_alloc_6d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:,:,:,:,:,:), allocatable, target :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 6

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_6d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:,:,:,:,:,:), allocatable, target :: val
    real (PREC), parameter :: val_kind = 1
    integer, parameter :: NDIM = 6

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_6d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out), dimension(:,:,:,:,:,:), allocatable, target :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 6

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_6d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:,:,:,:,:,:), allocatable, target :: val
    integer (INTSIZE), parameter :: val_kind = 1
    integer, parameter :: NDIM = 6

#include "hdf5/hdf5_load_alloc_array_impl.f90"
end subroutine

subroutine hdf5_load_alloc_6d_logical (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer, parameter :: NDIM = 6
    logical, intent(out), dimension(:,:,:,:,:,:), allocatable, target :: val

    integer (hid_t), intent(in) :: loc_id
    character (*), intent(in) :: name
    logical, intent(in), optional :: ignore_missing
    integer, intent(out), optional :: status

    integer :: lstatus
    integer, dimension(NDIM) :: shp, ishp
    integer (INTSIZE), dimension(:,:,:,:,:,:), allocatable :: ival

    call hdf5_load_alloc (loc_id, name, ival, ignore_missing, lstatus)
    if (lstatus == STATUS_OK .and. allocated(ival)) then
        if (allocated(val)) then
            shp = shape(val)
            ishp = shape(ival)
            if (all(shp == ishp)) deallocate (val)
        end if

        if (.not. allocated(val)) then
            allocate (val(size(ival,1), size(ival,2), size(ival,3), size(ival,4), size(ival,5), size(ival,6)))
            val(:,:,:,:,:,:) = (ival == 1)
        end if
    end if

    if (present(status)) status = lstatus

end subroutine

! ------------------------------------------------------------------------------
! HDF5_ALLOC for 1d-arrays

subroutine hdf5_alloc_1d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:), allocatable, target :: val
    integer, parameter :: NDIM = 1

#include "hdf5/hdf5_alloc_impl_1d.f90"
end subroutine

subroutine hdf5_alloc_1d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:), allocatable, target :: val
    integer, parameter :: NDIM = 1

#include "hdf5/hdf5_alloc_impl_1d.f90"
end subroutine

subroutine hdf5_alloc_1d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), dimension(:), allocatable, target :: val
    integer, parameter :: NDIM = 1

#include "hdf5/hdf5_alloc_impl_1d.f90"
end subroutine

subroutine hdf5_alloc_1d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:), allocatable, target :: val
    integer, parameter :: NDIM = 1

#include "hdf5/hdf5_alloc_impl_1d.f90"
end subroutine


! ------------------------------------------------------------------------------
! HDF5_ALLOC for 2d-arrays

subroutine hdf5_alloc_2d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:,:), allocatable, target :: val
    integer, parameter :: NDIM = 2

#include "hdf5/hdf5_alloc_impl_2d.f90"
end subroutine

subroutine hdf5_alloc_2d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:,:), allocatable, target :: val
    integer, parameter :: NDIM = 2

#include "hdf5/hdf5_alloc_impl_2d.f90"
end subroutine

subroutine hdf5_alloc_2d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out), dimension(:,:), allocatable, target :: val
    integer, parameter :: NDIM = 2

#include "hdf5/hdf5_alloc_impl_2d.f90"
end subroutine

subroutine hdf5_alloc_2d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:,:), allocatable, target :: val
    integer, parameter :: NDIM = 2

#include "hdf5/hdf5_alloc_impl_2d.f90"
end subroutine

! ------------------------------------------------------------------------------
! HDF5_ALLOC for 3d-arrays

subroutine hdf5_alloc_3d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:,:,:), allocatable, target :: val
    integer, parameter :: NDIM = 3

#include "hdf5/hdf5_alloc_impl_3d.f90"
end subroutine

subroutine hdf5_alloc_3d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:,:,:), allocatable, target :: val
    integer, parameter :: NDIM = 3

#include "hdf5/hdf5_alloc_impl_3d.f90"
end subroutine

subroutine hdf5_alloc_3d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out), dimension(:,:,:), allocatable, target :: val
    integer, parameter :: NDIM = 3

#include "hdf5/hdf5_alloc_impl_3d.f90"
end subroutine

subroutine hdf5_alloc_3d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:,:,:), allocatable, target :: val
    integer, parameter :: NDIM = 3

#include "hdf5/hdf5_alloc_impl_3d.f90"
end subroutine


! ------------------------------------------------------------------------------
! HDF5_ALLOC for 4d-arrays

subroutine hdf5_alloc_4d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:,:,:,:), allocatable, target :: val
    integer, parameter :: NDIM = 4

#include "hdf5/hdf5_alloc_impl_4d.f90"
end subroutine

subroutine hdf5_alloc_4d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:,:,:,:), allocatable, target :: val
    integer, parameter :: NDIM = 4

#include "hdf5/hdf5_alloc_impl_4d.f90"
end subroutine

subroutine hdf5_alloc_4d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out), dimension(:,:,:,:), allocatable, target :: val
    integer, parameter :: NDIM = 4

#include "hdf5/hdf5_alloc_impl_4d.f90"
end subroutine

subroutine hdf5_alloc_4d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:,:,:,:), allocatable, target :: val
    integer, parameter :: NDIM = 4

#include "hdf5/hdf5_alloc_impl_4d.f90"
end subroutine

! ------------------------------------------------------------------------------
! HDF5_ALLOC for 5d-arrays

subroutine hdf5_alloc_5d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:,:,:,:,:), allocatable, target :: val
    integer, parameter :: NDIM = 5

#include "hdf5/hdf5_alloc_impl_5d.f90"
end subroutine

subroutine hdf5_alloc_5d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:,:,:,:,:), allocatable, target :: val
    integer, parameter :: NDIM = 5

#include "hdf5/hdf5_alloc_impl_5d.f90"
end subroutine

subroutine hdf5_alloc_5d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out), dimension(:,:,:,:,:), allocatable, target :: val
    integer, parameter :: NDIM = 5

#include "hdf5/hdf5_alloc_impl_5d.f90"
end subroutine

subroutine hdf5_alloc_5d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:,:,:,:,:), allocatable, target :: val
    integer, parameter :: NDIM = 5

#include "hdf5/hdf5_alloc_impl_5d.f90"
end subroutine


! ------------------------------------------------------------------------------
! HDF5_ALLOC for 6d-arrays

subroutine hdf5_alloc_6d_real32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real32
    real (PREC), intent(out), dimension(:,:,:,:,:,:), allocatable, target :: val
    integer, parameter :: NDIM = 6

#include "hdf5/hdf5_alloc_impl_6d.f90"
end subroutine

subroutine hdf5_alloc_6d_real64 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: PREC = real64
    real (PREC), intent(out), dimension(:,:,:,:,:,:), allocatable, target :: val
    integer, parameter :: NDIM = 6

#include "hdf5/hdf5_alloc_impl_6d.f90"
end subroutine

subroutine hdf5_alloc_6d_int8 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int8
    integer (INTSIZE), intent(out), dimension(:,:,:,:,:,:), allocatable, target :: val
    integer, parameter :: NDIM = 6

#include "hdf5/hdf5_alloc_impl_6d.f90"
end subroutine

subroutine hdf5_alloc_6d_int32 (loc_id, name, val, ignore_missing, status)
    integer, parameter :: INTSIZE = int32
    integer (INTSIZE), intent(out), dimension(:,:,:,:,:,:), allocatable, target :: val
    integer, parameter :: NDIM = 6

#include "hdf5/hdf5_alloc_impl_6d.f90"
end subroutine

! ------------------------------------------------------------------------------
subroutine hdf5_deflate_info (has_encode, has_decode)
    logical, intent(out) :: has_encode, has_decode

    logical :: avail
    integer :: status, filter_info

    has_encode = .false.
    has_decode = .false.

    call h5zfilter_avail_f (H5Z_FILTER_DEFLATE_F, avail, status)
    if (avail) then
        call h5zget_filter_info_f (H5Z_FILTER_DEFLATE_F, filter_info, status)

        if (status < 0) return

        has_encode = iand(H5Z_FILTER_ENCODE_ENABLED_F, filter_info) == &
            H5Z_FILTER_ENCODE_ENABLED_F
        has_decode = iand(H5Z_FILTER_DECODE_ENABLED_F, filter_info) == &
            H5Z_FILTER_DECODE_ENABLED_F
    end if
end subroutine

end module
