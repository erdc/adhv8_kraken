#=============================================================================
# AdH - Adaptive Hydraulics
# fe diffusive wave and groundwater coupling submodule
#=============================================================================
# General source files (always included)
list(APPEND fe_dwgw_hybrid_srcs
  # Insert new files below here, alphabetically
    fe_dwgw_hybrid_inc.c
    fe_dwgw_hybrid_init.c
    fe_dwgw_hybrid_load.c
    fe_dwgw_hybrid_resid.c
    fe_dwgw_hybrid_solve.c
    fe_dwgw_hybrid_update.c
  # Insert new files above here, alphabetically
  )

############################################################################

list(REMOVE_DUPLICATES fe_dwgw_hybrid_srcs)

add_library("adh_fe_dwgw_hybrid" STATIC EXCLUDE_FROM_ALL ${fe_dwgw_hybrid_srcs})

target_link_libraries(adh_fe_dwgw_hybrid adh_fe adh_fe_diffusive adh_fe_gw)

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
