#=============================================================================
# AdH - Adaptive Hydraulics
# fe submodule
#=============================================================================
# General source files (always included)
list(APPEND fe_srcs
  # Insert new files below here, alphabetically
  fe_assign_elem_mat_db.c
  fe_assemble_matrix.c
  fe_elem1d_integrations.c
  fe_elem2d_integrations.c
  fe_elem3d_integrations.c
  fe_elem_utilities.c
  fe_get_supg_tau.c
  fe_main.c
  fe_newton.c
  fe_petsc_matrix_tools.c
  fe_sw2_wet_dry_wrapper.c
  fe_flux_couple.c
  # Insert new files above here, alphabetically
  )

############################################################################

list(REMOVE_DUPLICATES fe_srcs)

add_library("adh_fe" STATIC EXCLUDE_FROM_ALL ${fe_srcs})

target_link_libraries(adh_fe adh_fe_diffusive adh_fe_sw2 adh_fe_sw3 adh_fe_ns3 adh_fe_sw_hybrid adh_fe_transport adh_fe_sediment adh_friction adh_structs adh_messg adh_turbulence)

if(USE_ADH_STRUCTURES)
  target_link_libraries(adh_fe adh_hydro_structures)
endif(USE_ADH_STRUCTURES)

if(USE_NSM)
  target_link_libraries(adh_fe adh_nsm)
endif(USE_NSM)

if(USE_GW)
  target_link_libraries(adh_fe adh_fe_gw)
  if (USE_DWGW)
      target_link_libraries(adh_fe adh_fe_dwgw_hybrid)
  endif(USE_DWGW)
endif(USE_GW)
#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
