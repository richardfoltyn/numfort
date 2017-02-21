module numfort_common_status

    use numfort_common_enums
    implicit none
    private

    public :: status_t
    public :: assignment (=)
    public :: operator(+), iand, ior
    public :: operator(==), operator(/=)
    public :: operator(.in.), operator(.notin.)
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
    contains
        procedure, public, pass :: decode => status_decode
    end type

    interface operator(+)
        module procedure add_int_status, add_status_int
    end interface

    interface ior
        module procedure add_int_status, add_status_int
    end interface

    interface iand
        module procedure iand_int_status, iand_status_int
    end interface

    interface assignment (=)
        module procedure assign_int_status, assign_status_int
    end interface

    interface operator (==)
        module procedure equal_status_status, equal_status_int, equal_int_status
    end interface

    interface operator (/=)
        module procedure nequal_status_status, nequal_status_int, nequal_int_status
    end interface

    interface operator (.in.)
        module procedure operator_in_int_status
    end interface

    interface operator (.notin.)
        module procedure operator_notin_int_status
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
    class (status_t), intent(in) :: self
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

    if ((size(x) >= 1) .and. (self /= NF_STATUS_UNDEFINED)) then
        do i = 1, min(size(x), NF_MAX_STATUS_CODES)
            pattern = ishft(1, i-1)
            if (iand(self, pattern) == pattern) then
                n = n + 1
                x(n) = i-1
            end if
        end do
    end if
end subroutine

pure function status_size (self) result(res)
    class (status_t), intent(in) :: self
    integer :: res

    integer :: i
    integer (NF_ENUM_KIND) :: pattern
    res = 0
    do i = NF_MAX_STATUS_CODES, 1, -1
        pattern = ishft(1, i-1)
        if (iand(self%code, pattern) == pattern) then
            res = i
            return
        end if
    end do
end function

pure function status_to_char (self) result(res)
    class (status_t), intent(in) :: self
    character (:), allocatable :: res

    integer, dimension(NF_MAX_STATUS_CODES) :: b2status
    integer :: nstatus, n
    character (:), allocatable :: buf

    call self%decode (b2status, nstatus)

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
! Operator overloads

elemental function add_status_int (lhs, rhs) result(res)
    class (status_t), intent(in) :: lhs
    integer (NF_ENUM_KIND), intent(in) :: rhs
    type (status_t) :: res
    res%code = ior(lhs%code, rhs)
end function

elemental function add_int_status (lhs, rhs) result(res)
    integer (NF_ENUM_KIND), intent(in) :: lhs
    class (status_t), intent(in) :: rhs
    type (status_t) :: res
    res%code = ior(lhs, rhs%code)
end function

elemental function iand_status_int (lhs, rhs) result(res)
    class (status_t), intent(in) :: lhs
    integer (NF_ENUM_KIND), intent(in) :: rhs
    type (status_t) :: res
    res%code = iand(lhs%code, rhs)
end function

elemental function iand_int_status (lhs, rhs) result(res)
    integer (NF_ENUM_KIND), intent(in) :: lhs
    class (status_t), intent(in) :: rhs
    type (status_t) :: res
    res%code = iand(lhs, rhs%code)
end function

elemental subroutine assign_int_status (lhs, rhs)
    integer (NF_ENUM_KIND), intent(out) :: lhs
    class (status_t), intent(in) :: rhs
    lhs = rhs%code
end subroutine

elemental subroutine assign_status_int (lhs, rhs)
    class (status_t), intent(out) :: lhs
    integer (NF_ENUM_KIND), intent(in) :: rhs
    lhs%code = rhs
end subroutine

elemental function equal_status_status (lhs, rhs) result(res)
    class (status_t), intent(in) :: lhs, rhs
    logical :: res
    res = (lhs%code == rhs%code)
end function

elemental function equal_status_int (lhs, rhs) result(res)
    class (status_t), intent(in) :: lhs
    integer (NF_ENUM_KIND), intent(in) :: rhs
    logical :: res
    res = (lhs%code == rhs)
end function

elemental function equal_int_status (lhs, rhs) result(res)
    integer (NF_ENUM_KIND), intent(in) :: lhs
    class (status_t), intent(in) :: rhs
    logical :: res
    res = (lhs == rhs%code)
end function

elemental function nequal_status_status (lhs, rhs) result(res)
    class (status_t), intent(in) :: lhs, rhs
    logical :: res
    res = .not. (lhs == rhs)
end function

elemental function nequal_status_int (lhs, rhs) result(res)
    class (status_t), intent(in) :: lhs
    integer (NF_ENUM_KIND), intent(in) :: rhs
    logical :: res
    res = .not. (lhs == rhs)
end function

elemental function nequal_int_status (lhs, rhs) result(res)
    integer (NF_ENUM_KIND), intent(in) :: lhs
    class (status_t), intent(in) :: rhs
    logical :: res
    res = .not. (lhs == rhs)
end function

elemental function operator_in_int_status (lhs, rhs) result(res)
    integer (NF_ENUM_KIND), intent(in) :: lhs
    class (status_t), intent(in) :: rhs
    logical :: res
    res = (iand(lhs, rhs%code) == lhs)
end function

elemental function operator_notin_int_status (lhs, rhs) result(res)
    integer (NF_ENUM_KIND), intent(in) :: lhs
    class (status_t), intent(in) :: rhs
    logical :: res
    res = .not. (lhs .in. rhs)
end function

end module
