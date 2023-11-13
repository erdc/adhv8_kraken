#=============================================================================
# AdH - Adaptive Hydraulics
# fe sw3 submodule
#=============================================================================
# General source files (always included)
list(APPEND fe_sw_hybrid_srcs
  # Insert new files below here, alphabetically
fe_sw_hybrid_inc.c
fe_sw_hybrid_init.c
fe_sw_hybrid_load.c
fe_sw_hybrid_resid.c
fe_sw_hybrid_solve.c
fe_sw_hybrid_update.c
fe_sw_hybrid_wvel_inc.c
fe_sw_hybrid_wvel_init.c
fe_sw_hybrid_wvel_load.c
fe_sw_hybrid_wvel_resid.c
fe_sw_hybrid_wvel_update.c
  # Insert new files above here, alphabetically
  )

############################################################################


list(REMOVE_DUPLICATES fe_sw_hybrid_srcs)

add_library("adh_fe_sw_hybrid" STATIC EXCLUDE_FROM_ALL ${fe_sw_hybrid_srcs})

target_link_libraries(adh_fe_sw_hybrid adh_friction adh_structs adh_messg adh_turbulence adh_fe adh_fe_sw2 adh_fe_sw3)

if(USE_NSM)
  target_link_libraries(adh_fe_sw3 adh_nsm)
endif(USE_NSM)

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
