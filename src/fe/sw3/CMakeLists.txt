#=============================================================================
# AdH - Adaptive Hydraulics
# fe sw3 submodule
#=============================================================================
# General source files (always included)
list(APPEND fe_sw3_srcs
  # Insert new files below here, alphabetically
fe_hvel_body_resid.c
fe_hvel_boundary_resid.c
fe_hvel_inc.c
fe_hvel_init.c
fe_hvel_load.c
fe_hvel_resid.c
fe_hvel_update.c
fe_sw3.c
fe_sw3_solve.c
fe_sw3_get_roughness.c
fe_sw3_prototypes.h
fe_wvel_body_resid.c
fe_wvel_boundary_resid.c
fe_wvel_inc.c
fe_wvel_init.c
fe_wvel_load.c
fe_wvel_resid.c
fe_wvel_update.c

  # Insert new files above here, alphabetically
  
  )

############################################################################


list(REMOVE_DUPLICATES fe_sw3_srcs)

add_library("adh_fe_sw3" STATIC EXCLUDE_FROM_ALL ${fe_sw3_srcs})

target_link_libraries(adh_fe_sw3 adh_friction adh_structs adh_messg adh_turbulence adh_fe)

if(USE_NSM)
  target_link_libraries(adh_fe_sw3 adh_nsm)
endif(USE_NSM)

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
