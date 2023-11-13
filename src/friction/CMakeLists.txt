#=============================================================================
# AdH - Adaptive Hydraulics
# friction submodule tests
#=============================================================================
# General source files (always included)
list(APPEND friction_srcs
  # Insert new files below here, alphabetically
    fr_bedshstr_drag_coef.c
    fr_defs.h
    fr_equiv_depavg_shstr_coef.c
    fr_manningsn_to_rheight.c
    fr_rheight_to_manningsn.c
    fr_sav_drag_coef.c
    fr_stationary_ice_coef.c
    fr_urv_drag_coef.c
  # Insert new files above here, alphabetically
  )



###############################################################################


list(REMOVE_DUPLICATES friction_srcs)

add_library("adh_friction" STATIC EXCLUDE_FROM_ALL ${friction_srcs})

target_link_libraries(adh_friction adh_tools)

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
