module numfort_common_status

    use numfort_common_enums
    implicit none
    private

    public :: status_t
    public :: status_decode

    public :: status_set, status_add, status_reset
    public :: status_equals, status_contains
    public :: char, size

    integer (NF_ENUM_KIND), public, parameter :: NF_STATUS_UNDEFINED = 0
    integer (NF_ENUM_KIND), public, parameter :: NF_STATUS_OK = ishft(1, 0)
    integer (NF_ENUM_KIND), public, parameter :: NF_STATUS_INVALID_ARG = ishft(1, 10)
    integer (NF_ENUM_KIND), public, parameter :: NF_STATUS_UNKNOWN = ishft(1, 11)
    integer (NF_ENUM_KIND), public, parameter :: NF_STATUS_UNSUPPORTED_OP = ishft(1, 12)
    integer (NF_ENUM_KIND), public, parameter :: NF_STATUS_INVALID_STATE = ishft(1, 13)
    integer (NF_ENUM_KIND), public, parameter :: NF_STATUS_BOUNDS_ERROR = ishft(1, 14)
    integer (NF_ENUM_KIND), public, parameter :: NF_STATUS_MAX_ITER = ishft(1, 15)
    integer (NF_ENUM_KIND), public, parameter :: NF_STATUS_MAX_EVAL = ishft(1, 16)
    integer (NF_ENUM_KIND), public, parameter :: NF_STATUS_NOT_CONVERGED = ishft(1, 17)
    integer (NF_ENUM_KIND), public, parameter :: NF_STATUS_STORAGE_ERROR = ishft(1, 18)
    integer (NF_ENUM_KIND), public, parameter :: NF_STATUS_OTHER = ishft(1, 19)
    integer (NF_ENUM_KIND), public, parameter :: NF_STATUS_APPROX = ishft(1, 20)

    integer, public, parameter :: NF_MAX_STATUS_CODES = bit_size (NF_STATUS_OK)

    type :: status_t
        private
        integer (NF_ENUM_KIND) :: code = NF_STATUS_UNDEFINED
    end type


    interface status_add
        module procedure status_add_int, status_add_status
    end interface

    interface status_set
        module procedure status_set_int, status_set_int_int, status_set_status
    end interface

    interface status_equals
        module procedure status_equals_int, status_equals_status
    end interface

    interface status_contains
        module procedure status_contains_int, status_contains_status
    end interface

    interface char
        module procedure status_to_char
    end interface

    interface size
        module procedure status_size
    end interface
contains

! ------------------------------------------------------------------------------
! METHODS

pure subroutine status_decode (self, x, n)
    !*  DECODE disaggregates a composize status code into its
    !   components and turns their base-2 exponents.
    type (status_t), intent(in) :: self
        !!  Status container object.
    integer, dimension(:), intent(out) :: x
        !!  Array to store individual status codes. The lowest-exponent codes
        !!  up to the maximum given by size(array) are returned.
    integer, intent(out) :: n
        !!  Number of individual status codes present in composize status.
        !!  (n <= size(array)).

    integer :: i
    integer (NF_ENUM_KIND) :: pattern

    ! Since status codes are integers >= 0, initialize to something invalid.
    x = -1
    n = 0

    if ((size(x) >= 1) .and. (.not. status_equals (self, NF_STATUS_UNDEFINED))) then
        do i = 1, min(size(x), NF_MAX_STATUS_CODES)
            pattern = ishft(1, i-1)
            if (status_contains (self, pattern)) then
                n = n + 1
                x(n) = i-1
            end if
        end do
    end if
end subroutine

pure function status_size (self) result(res)
    type (status_t), intent(in) :: self
    integer :: res

    integer :: i
    integer (NF_ENUM_KIND) :: pattern
    res = 0
    do i = NF_MAX_STATUS_CODES, 1, -1
        pattern = ishft(1, i-1)
        if (status_contains (self, pattern)) then
            res = i
            return
        end if
    end do
end function

pure function status_to_char (self) result(res)
    type (status_t), intent(in) :: self
    character (:), allocatable :: res

    integer, dimension(NF_MAX_STATUS_CODES) :: b2status
    integer :: nstatus, n
    character (:), allocatable :: buf

    call status_decode (self, b2status, nstatus)

    ! compute char length including separators and parenthesis
    n = NF_MAX_STATUS_CODES * 4 + 2
    allocate (character (n) :: buf)

    write (buf, '("(", *(i0, :, ", "))') b2status(1:nstatus)
    n = len_trim (buf)
    buf(n+1:n+1) = ')'
    n = n+1
    allocate (character (n) :: res)
    res(1:n) = buf(1:n)

    deallocate (buf)
end function

! ------------------------------------------------------------------------------
! "Operators"

pure subroutine status_add_int (self, other)
    type (status_t), intent(in out) :: self
    integer (NF_ENUM_KIND), intent(in) :: other
    self%code = ior(self%code, other)
end subroutine

pure subroutine status_add_status (self, other)
    type (status_t), intent(in out) :: self
    type (status_t), intent(in) :: other
    self%code = ior(self%code, other%code)
end subroutine

pure subroutine status_set_int (self, other)
    type (status_t), intent(out) :: self
    integer (NF_ENUM_KIND), intent(in) :: other
    self%code = other
end subroutine

pure subroutine status_set_int_int (self, val1, val2)
    type (status_t), intent(out) :: self
    integer (NF_ENUM_KIND), intent(in) :: val1, val2
    self%code = ior(val1, val2)
end subroutine

pure subroutine status_set_status (self, other)
    type (status_t), intent(out) :: self
    type (status_t), intent(in) :: other
    self%code = other%code
end subroutine

pure function status_equals_status (self, other) result(res)
    type (status_t), intent(in) :: self, other
    logical :: res
    res = (self%code == other%code)
end function

pure function status_equals_int (self, other) result(res)
    type (status_t), intent(in) :: self
    integer (NF_ENUM_KIND), intent(in) :: other
    logical :: res
    res = (self%code == other)
end function

pure function status_contains_int (self, other) result(res)
    type (status_t), intent(in) :: self
    integer (NF_ENUM_KIND), intent(in) :: other
    logical :: res
    res = (iand(self%code, other) == other)
end function

pure function status_contains_status (self, other) result(res)
    type (status_t), intent(in) :: self, other
    logical :: res
    res = (iand(self%code, other%code) == other%code)
end function

pure subroutine status_reset (self)
    type (status_t), intent(out) :: self
    self%code = NF_STATUS_UNDEFINED
end subroutine

end module
