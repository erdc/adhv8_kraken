#=============================================================================
# AdH - Adaptive Hydraulics
# fe ns3 submodule
#=============================================================================
# General source files (always included)
list(APPEND fe_ns3_srcs
  # Insert new files below here, alphabetically
  fe_ns3.c
  fe_ns3_body_resid.c
  fe_ns3_boundary_resid.c
  fe_ns3_inc.c
  fe_ns3_init.c
  fe_ns3_load.c
  fe_ns3_resid.c
  fe_ns3_update.c
  fe_ns3_solve.c
  # Insert new files above here, alphabetically
  
  )

############################################################################


list(REMOVE_DUPLICATES fe_ns3_srcs)

add_library("adh_fe_ns3" STATIC EXCLUDE_FROM_ALL ${fe_ns3_srcs})

target_link_libraries(adh_fe_ns3 adh_friction adh_structs adh_messg adh_turbulence adh_fe)

