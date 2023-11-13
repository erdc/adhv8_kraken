#=============================================================================
# AdH - Adaptive Hydraulics
# fe sediment submodule
#=============================================================================

list(APPEND fe_sediment_srcs
    fe_sediment_transport.c
    )

if(USE_SEDLIB)
  list(APPEND fe_sediment_srcs
    # Insert files below here, alphabetically
    fe_bedload_boundary_resid.c
    fe_bedload_body_resid.c
    fe_bedload_inc.c
    fe_bedload_init.c
    fe_bedload_load.c
    fe_bedload_resid.c
    fe_bedload_update.c
    fe_sediment_prototypes.h
    )
endif(USE_SEDLIB)

############################################################################

list(REMOVE_DUPLICATES fe_sediment_srcs)

add_library("adh_fe_sediment" STATIC EXCLUDE_FROM_ALL ${fe_sediment_srcs})

target_link_libraries(adh_fe_sediment adh_friction adh_structs adh_messg adh_turbulence)

if(USE_NSM)
  target_link_libraries(adh_fe_sediment adh_nsm)
endif(USE_NSM)

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
