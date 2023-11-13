#=============================================================================
# AdH - Adaptive Hydraulics
# node submodule
#=============================================================================
# General source files (always included)
list(APPEND node_srcs
  # Insert new files below here, alphabetically
    node2element_map.c
    node_avg.c
    node_cmp.c
    node_get_local.c
    
    node_hash_add_entry.c
    node_hash_index.c
    node_hash_lookup.c
    node_new.c
    node_new_surface.c
    node_new_bed.c
    node_order.c
    node_renumber.c
    node_renumber_double.c
    node_renumber_gn.c
    node_renumber_int.c
    node_renumber_snode.c
    node_renumber_surface.c
    node_renumber_tensor2d.c
    node_renumber_vect.c
    node_renumber_vect2d.c 
  # Insert new files above here, alphabetically
  )
    if(USE_MPI)
      list(APPEND node_srcs
      node_pack_dcnt.c
      node_pack_dcnt_surface.c
      node_pack_dcnt_bed.c
      node_pack_icnt.c
      node_packd.c
      node_packd_surface.c
      node_packd_bed.c
      node_packi.c
      node_unpacki.c
      node_unpackd.c
      node_unpackd_surface.c
      node_unpackd_bed.c
      )
    endif(USE_MPI)

#############################################################################

list(REMOVE_DUPLICATES node_srcs)

add_library("adh_node" STATIC EXCLUDE_FROM_ALL ${node_srcs})

target_link_libraries(adh_node)

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
