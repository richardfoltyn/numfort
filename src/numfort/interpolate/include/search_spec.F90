
interface bsearch
    procedure __APPEND(bsearch,__PREC)
end interface

interface bsearch_cached
    procedure __APPEND(bsearch_cached,__PREC)
end interface

interface bsearch_cached_impl
    procedure __APPEND(bsearch_cached_impl,__PREC)
end interface

interface interp_find
   procedure __APPEND(interp_find_scalar,__PREC), &
        __APPEND(interp_find_1d,__PREC)
end interface

interface interp_find
    procedure __APPEND(interp_find_cache_scalar,__PREC)
end interface

interface interp_find_impl
    procedure __APPEND(interp_find_impl_scalar,__PREC), &
        __APPEND(interp_find_impl_1d,__PREC)
end interface

interface interp_find_impl
   procedure __APPEND(interp_find_impl_cache_scalar,__PREC)
end interface

interface interp_find_impl_ext
    procedure __APPEND(interp_find_impl_ext_1d,__PREC)
end interface

interface interp_find_impl_default
    procedure __APPEND(interp_find_impl_default_1d,__PREC)
end interface

interface interp_find_check_input
    procedure __APPEND(interp_find_check_input,__PREC)
end interface

