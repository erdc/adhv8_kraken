#=============================================================================
# AdH - Adaptive Hydraulics
# grid submodule
#=============================================================================
# General source files (always included)
list(APPEND grid_srcs
  # Insert new files below here, alphabetically
    adpt_fix_adj.c
    adpt_fix_global_ids.c
    adpt_get_node.c
    adpt_main.c
    adpt_merge_elem.c
    adpt_rank_edges.c
    adpt_ref.c
    adpt_set_flags.c
    adpt_unref.c
    adpt_fix_parents.c
    adpt_set_node_number.c
    error_main.c
    partition_cleanup.c
    # Insert new files above here, alphabetically
  )
if(USE_MPI)
  list(APPEND grid_srcs
  findrank.c
  partition_adpt.c
  partition_form.c
  partition_main.c
  partition_form_all.c
  partition_transfer.c
  # Insert new files above here, alphabetically
  )
endif(USE_MPI)

#############################################################################

list(REMOVE_DUPLICATES grid_srcs)

add_library("adh_grid" STATIC EXCLUDE_FROM_ALL ${grid_srcs})

if(USE_MPI)
  target_link_libraries(adh_grid adh_node adh_tools adh_comm)
else(USE_MPI)
  target_link_libraries(adh_grid adh_node adh_tools)
endif(USE_MPI)
#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
