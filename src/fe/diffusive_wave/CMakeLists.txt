#=============================================================================
# AdH - Adaptive Hydraulics
# fe diffusive wave submodule
#=============================================================================
# General source files (always included)
list(APPEND fe_diffusive_srcs
  # Insert new files below here, alphabetically
  fe_diffusive.c
  fe_diffusive_solve.c
  fe_diffusive_boundary_resid.c
  fe_diffusive_body_resid.c
  fe_diffusive_conv.c
  fe_diffusive_inc.c
  fe_diffusive_init.c
  fe_diffusive_load.c
  fe_diffusive_prototypes.h
  fe_diffusive_resid.c
  fe_diffusive_update.c
  fe_diffusive_velocities.c

  # Insert new files above here, alphabetically
  
  )

############################################################################


list(REMOVE_DUPLICATES fe_diffusive_srcs)

add_library("adh_fe_diffusive" STATIC EXCLUDE_FROM_ALL ${fe_diffusive_srcs})

target_link_libraries(adh_fe_diffusive adh_friction adh_structs adh_messg adh_turbulence)

if(USE_NSM)
  target_link_libraries(adh_fe_diffusive adh_nsm)
endif(USE_NSM)

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
