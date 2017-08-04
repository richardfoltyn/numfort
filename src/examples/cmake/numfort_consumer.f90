program numfort_consumer

    use numfort
    implicit none


    character (*), parameter :: LIB_VERSION = NUMFORT_VERSION

    print '(a, tr1, a)', "Using NUMFORT version", LIB_VERSION

end
