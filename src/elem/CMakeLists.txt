#=============================================================================
# AdH - Adaptive Hydraulics
# elem submodule
#=============================================================================
# General source files (always included)
list(APPEND elem_srcs
  # Insert new files below here, alphabetically
edge_hash_add_entry.c
edge_hash_index.c
edge_hash_lookup.c
elem1d_find_elem2d.c
elem1d_find_elem2d_init.c
elem1d_hash_add_entry.c
elem1d_hash_lookup.c
elem1d_merge.c
elem1d_new.c
elem1d_outward_nrml.c
elem1d_renumber.c
elem1d_split.c
elem2d_hash_add_entry.c
elem2d_hash_lookup.c
elem2d_merge.c
elem2d_new.c
elem2d_normal.c
elem2d_renumber.c
elem2d_split.c
elem3d_hash_add_entry.c
elem3d_hash_lookup.c
elem3d_new.c
elem3d_split.c
elem3d_merge.c
elem3d_renumber.c
elem_hash_add_entry.c
elem_hash_index.c
elem_level.c
element_djac_nrml_gradPhi.c
# Insert new files above here, alphabetically
  )

if(USE_MPI)
  list(APPEND elem_srcs
    elem1d_pack_icnt.c
    elem1d_packi.c
    elem1d_packd.c
    elem1d_unpacki.c
    elem2d_pack_icnt.c
    elem2d_packi.c
    elem2d_packd.c
    elem2d_unpacki.c
    elem3d_pack_icnt.c
    elem3d_packi.c
    elem3d_packd.c
    elem3d_unpacki.c
  )
endif(USE_MPI)
#############################################################################

list(REMOVE_DUPLICATES elem_srcs)

add_library("adh_elem" STATIC EXCLUDE_FROM_ALL ${elem_srcs})

target_link_libraries(adh_elem adh_tools)

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
