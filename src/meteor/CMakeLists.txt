#=============================================================================
# AdH - Adaptive Hydraulics
# meteor submodule
#=============================================================================
# General source files (always included)
list(APPEND meteor_srcs
  # Insert new files below here, alphabetically
  meteor_update_series.c
  meteor_update.c
  station_node_contrib.c
  winds_proc_teeter.c
  winds_proc_wu.c 
  # Insert new files above here, alphabetically
  )
#############################################################################

list(REMOVE_DUPLICATES meteor_srcs)

add_library("adh_meteor" STATIC EXCLUDE_FROM_ALL ${meteor_srcs})

if(USE_WINDLIB)
    target_link_libraries(adh_meteor adh_wind)
else(USE_WINDLIB)
    target_link_libraries(adh_meteor)
endif(USE_WINDLIB)

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
