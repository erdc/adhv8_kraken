#=============================================================================
# AdH - Adaptive Hydraulics
# columns submodule
#=============================================================================

# General source files (always included)
list(APPEND columns_srcs
  # Insert new files below here, alphabetically
  column_adpt_merge_elem.c 
  column_adpt_set_flags.c  
  column_elem2d_merge.c    
  column_elem3d_merge.c    
  column_adpt_ref.c        
  column_adpt_unref.c      
  column_elem2d_split.c    
  column_elem3d_split.c
  tl_columns.c
  # Insert new files above here, alphabetically
  )


############################################################################

list(REMOVE_DUPLICATES columns_srcs)

add_library("adh_columns" STATIC EXCLUDE_FROM_ALL ${columns_srcs})

target_link_libraries(adh_columns adh_initio adh_debug adh_elem)

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
