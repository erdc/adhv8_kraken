#===========================================
# AdH - Adaptive Hydraulics
# structures submodule
#===========================================
# General source files (always included)
if(USE_ADH_STRUCTURES)
list(APPEND hydro_structures_srcs
  # Insert new files below
  flap_flow.c
  flap_compound_flux.c
  weir_flow.c
  weir_compound_flux.c
  sluice_flow.c
  sluice_compound_flux.c
  )
 
 list(REMOVE_DUPLICATES hydro_structures_srcs)
 add_library("adh_hydro_structures" STATIC EXCLUDE_FROM_ALL ${hydro_structures_srcs})
 
 target_link_libraries(adh_hydro_structures)

  
endif(USE_ADH_STRUCTURES)
