#=============================================================================
# AdH - Adaptive Hydraulics
# turbulence submodule
#=============================================================================
# General source files (always included)
list(APPEND turbulence_srcs
  # Insert new files below here, alphabetically
  tur_MY_2.c
  tur_smag.c
  tur_ws.c
  # Insert new files above here, alphabetically
  )


#############################################################################

list(REMOVE_DUPLICATES turbulence_srcs)

add_library("adh_turbulence" STATIC EXCLUDE_FROM_ALL ${turbulence_srcs})

target_link_libraries(adh_turbulence)

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
