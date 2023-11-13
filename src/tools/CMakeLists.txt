#=============================================================================
# AdH - Adaptive Hydraulics
# tools submodule

# General source files (always included)
list(APPEND tools_srcs
  # Insert new files below here, alphabetically
  root_print.c
  elem2d_to_node.c
  elem3d_to_node.c
  is_my_node.c
  is_my_elem.c
  projected_node_distance.c
  tc_init.c
  tc_end.c
  tc_timeunits.c
  tc_conversion_factor.c
  tc_trap_area.c
  tc_scale.c
  tc_eval_series.c
  tl_3d_to_2d.c
  tl_bendway_correction.c
  tl_grid_mass.c
  tl_list_alloc.c
  tl_matrix.c
  tl_model_units.c
  tl_calculate_pressure.c
  tl_find_edge_mdpt_pressure.c
  tl_calculate_depthavgvel.c
  tl_find_avg_column_depth.c
  tl_find_common_int.c
  tl_vertical_adapt.c
  tl_density_calculator_metric.c
  tl_nodal_grad_avg_2d.c
  tl_long_edge.c
  tl_find_nodes_on_slice.c
  math_simple.c
  math_special_functions.c
  math_blas_wrappers.c
  printScreen_array.c
  printScreen_resid.c
  create_dbovl_elem1d.c
  isPointInTriangle.c 
  createChildMesh.c
  coriolis.c
  rotate3D.c
  doesFileExist.c
  # Insert new files above here, alphabetically
  )
if(USE_MPI)
  list(APPEND tools_srcs
    print_tools_MPI.c
    )
endif(USE_MPI)
if(USE_GW)
  list(APPEND tools_srcs
    tc_eval_slope.c
    )
endif(USE_GW)

#############################################################################


list(REMOVE_DUPLICATES tools_srcs)

add_library("adh_tools" STATIC EXCLUDE_FROM_ALL ${tools_srcs})

target_link_libraries(adh_tools adh_structs adh_messg)

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
