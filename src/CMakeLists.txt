

################################################################################
# Required targets

# External third-party projects
add_subdirectory(external)

# Stand-alone interface definitions for F77 routines and F95 wrappers
add_subdirectory(interfaces)
add_subdirectory(wrappers)

# Numfort library code
add_subdirectory(numfort)