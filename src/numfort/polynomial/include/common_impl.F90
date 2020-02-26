
pure subroutine __APPEND(bsearch_proxy,__PREC) (xi, knots, i, cache_arr, cache, j)
    !*  BSEARCH_PROXY is a helper routine that uses the appropriate search
    !   cache depending on which arguments are present.
    integer, parameter :: PREC = __PREC
    real (PREC), intent(in) :: xi
    real (PREC), intent(in), dimension(:), contiguous :: knots
    integer, intent(in) :: i
    type (search_cache), intent(inout), dimension(:), optional, contiguous :: cache_arr
    !*  Array of interpolation-point-specific search cache objects
    type (search_cache), intent(inout) :: cache
    !*  Fallback shared search cache object
    integer, intent(out) :: j

    if (present(cache_arr)) then
        call bsearch_cached (xi, knots, j, cache_arr(i))
    else
        call bsearch_cached (xi, knots, j, cache)
    end if

end subroutine
