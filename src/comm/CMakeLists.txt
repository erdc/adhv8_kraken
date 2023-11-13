#=============================================================================
# AdH - Adaptive Hydraulics
# messg submodule
#=============================================================================
# General source files (always included)
list(APPEND messg_srcs
  # Insert new files below here, alphabetically
    comm_adj_nodes.c
    comm_check.c
    comm_elem_levels.c
    comm_edge_setup.c
    comm_flag_edges.c
    comm_new_node_nums.c
    comm_node_data.c
    comm_set_keys.c
    comm_update_con.c
    comm_update_edges.c
    comm_update_GN.c
    comm_update_TENSOR2D.c
    comm_update_VECT.c
    comm_update_VECT2D.c
    comm_update_double.c
    comm_update_int.c
    comm_update_double_surf.c
    comm_update_int_surf.c
    comm_check_int.c
    comm_update_snode.c
    comm_update_sw2.c
    comm_update_ns3.c
    comm_update_sw3.c
    comm_update_sw3_surface.c
    comm_update_swave.c

    comm_set_keys_supermodel.c # gkc for Multi-submodel runs
    comm_set_keys_supermodel_wvel.c # gkc for Multi-submodel runs

  # Insert new files above here, alphabetically
  )

##############################################################################

if(USE_SEDLIB)
  list(APPEND messg_srcs
    comm_update_sed.c
    )
endif(USE_SEDLIB)
if(USE_GW)
  list(APPEND messg_srcs
    comm_update_gw.c
    )
endif(USE_GW)
list(REMOVE_DUPLICATES messg_srcs)

add_library("adh_comm" STATIC EXCLUDE_FROM_ALL ${messg_srcs})

target_link_libraries(adh_comm adh_mesh adh_elem adh_grid adh_messg adh_structs)

#target_link_libraries()

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
