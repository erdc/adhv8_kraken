#=============================================================================
# AdH - Adaptive Hydraulics
# fe transport submodule
#=============================================================================
# General source files (always included)
list(APPEND fe_transport_srcs
  # Insert new files below here, alphabetically
fe_2d_transport_body_resid.c
fe_2d_transport_boundary_resid.c
fe_2d_transport_wet_dry_integrations.c
fe_3d_transport_body_resid.c
fe_3d_transport_boundary_resid.c
fe_transport.c
fe_transport_hybrid.c
fe_transport_hybrid_inc.c
fe_transport_hybrid_init.c
fe_transport_hybrid_load.c
fe_transport_hybrid_resid.c
fe_transport_hybrid_update.c
fe_transport_inc.c
fe_transport_init.c
fe_transport_load.c
fe_transport_prototypes.h
fe_transport_resid.c
fe_transport_update.c
  # Insert new files above here, alphabetically
  
  )

############################################################################


list(REMOVE_DUPLICATES fe_transport_srcs)

add_library("adh_fe_transport" STATIC EXCLUDE_FROM_ALL ${fe_transport_srcs})

target_link_libraries(adh_fe_transport adh_friction adh_structs adh_messg adh_turbulence)

if(USE_NSM)
  target_link_libraries(adh_fe_transport adh_nsm)
endif(USE_NSM)

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
