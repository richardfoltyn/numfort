

interface bsearch
    procedure bsearch, bsearch_1d
end interface

interface bsearch_cached
    procedure bsearch_cached, bsearch_cached_1d
end interface

interface interp_find
   procedure interp_find_0d, interp_find_1d
end interface

interface interp_find
    procedure interp_find_cached_0d, interp_find_cached_1d
end interface

interface interp_find_impl
    procedure interp_find_impl_0d, interp_find_impl_1d
end interface

interface interp_find_impl
    procedure interp_find_cached_impl_0d, interp_find_cached_impl_1d
end interface

interface interp_find_impl_ext
    procedure interp_find_impl_ext_1d, interp_find_cached_impl_ext_1d
end interface

interface interp_find_impl_default
    procedure interp_find_impl_default_1d, interp_find_cached_impl_default_1d
end interface

interface interp_find_incr_impl
   procedure interp_find_incr_ext_impl
end interface

interface interp_find_decr
    procedure interp_find_decr_ext
end interface

interface interp_find_decr_impl
    procedure interp_find_decr_ext_impl
end interface
