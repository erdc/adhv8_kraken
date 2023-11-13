#=============================================================================
# AdH - Adaptive Hydraulics
# fe sw2 submodule
#=============================================================================
# General source files (always included)
list(APPEND fe_sw2_srcs
  # Insert new files below here, alphabetically
  fe_sw2.c
  fe_sw2_boundary_resid.c
  fe_sw2_body_resid.c
  fe_sw2_get_EEVF.c
  fe_sw2_get_roughness.c
  fe_sw2_inc.c
  fe_sw2_init.c
  fe_sw2_load.c
  fe_sw2_solve.c
  fe_sw2_prototypes.h
  fe_sw2_resid.c
  fe_sw2_update.c
  fe_sw2_wdflag.c
  fe_sw2_wet_dry_integrations.c
  # Insert new files above here, alphabetically
  
  )

############################################################################


list(REMOVE_DUPLICATES fe_sw2_srcs)

add_library("adh_fe_sw2" STATIC EXCLUDE_FROM_ALL ${fe_sw2_srcs})

target_link_libraries(adh_fe_sw2 adh_structs adh_messg adh_turbulence adh_fe)

if(USE_NSM)
  target_link_libraries(adh_fe_sw2 adh_nsm)
endif(USE_NSM)

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
