#=============================================================================
# AdH - Adaptive Hydraulics
# fe sw3 submodule
#=============================================================================
# General source files (always included)
list(APPEND fe_gw_srcs
  # Insert new files below here, alphabetically
  fe_gw.c
  fe_gw_inc.c
  fe_gw_init.c
  fe_gw_update.c
  fe_gw_load.c
  fe_gw_resid.c
  fe_gw_body_resid.c
  fe_gw_boundary_resid.c
  fe_gw_solve.c
  )

############################################################################


list(REMOVE_DUPLICATES fe_gw_srcs)
add_library("adh_fe_gw" STATIC EXCLUDE_FROM_ALL ${fe_gw_srcs})

if (USE_GW)
   target_link_libraries(adh_fe_gw adh_structs adh_messg  adh_fe)
endif(USE_GW)

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
