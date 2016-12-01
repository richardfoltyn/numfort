! Note: CPP/FPP symbols __XXX__ are definitions added by build system,
! definitions __XXX are defined in the source code.

#ifdef __SUPPORTS_PDT_KIND__
#define __PDT_PARAM_DECL(name) (name)
#define __PDT_PARAM_DECL_BOTH(name1, name2) (name1,name2)
#define __PDT_KIND_DECL(name, val) integer, kind :: name = val
#else
#define __PDT_PARAM_DECL(list)
#define __PDT_PARAM_DECL_BOTH(name1, name2)
#define __PDT_KIND_DECL(name, val)
#endif
