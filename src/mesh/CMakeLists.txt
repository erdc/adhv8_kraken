#=============================================================================
# AdH - Adaptive Hydraulics
# messg submodule
#=============================================================================
# General source files (always included)
list(APPEND messg_srcs
  # Insert new files below here, alphabetically
    mesh_pack_elem1d.c
    mesh_pack_elem2d.c
    mesh_pack_elem3d.c
    mesh_pack_node.c
    mesh_unpack_elem1d.c
    mesh_unpack_elem2d.c
    mesh_unpack_elem3d.c
    mesh_unpack_node.c
  # Insert new files above here, alphabetically
  )

##############################################################################

list(REMOVE_DUPLICATES messg_srcs)

add_library("adh_mesh" STATIC EXCLUDE_FROM_ALL ${messg_srcs})

target_link_libraries(adh_mesh adh_grid adh_structs adh_elem adh_tools)

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
