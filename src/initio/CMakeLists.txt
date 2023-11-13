#=============================================================================
# AdH - Adaptive Hydraulics
# initio submodule
#=============================================================================
# General source files (always included)
list(APPEND initio_srcs
  # Insert new files below here, alphabetically
  io_main.c
  parse_card.c
  print_build_info.c
  print_grid_to_file.c
  print_parents_to_file.c
  print_header.c
  print_runtime_info.c
  print_trailer.c
  print_xms_trailers.c
  open_run_files.c
  read_bc_prep.c
  read_geo.c
  read_faces.c
  read_int_field.c
  read_dbl_field.c
  read_text_field.c
  read_physics.c
  get_string_id.c
  get_material_id.c
  get_node_id.c
  get_transport_id.c
  get_bed_layer.c
  read_bc_prep.c
  read_bc.c
  read_bc_DB.c
  read_bc_FR.c
  read_bc_IP.c
  read_bc_MP.c
  read_bc_NB.c
  read_bc_NSM.c
  read_bc_OFF.c
  read_bc_OP.c
  read_bc_PC.c
  read_bc_STRINGS.c
  read_bc_TC.c
  read_bc_CN.c
  read_bc_DEBUG.c
  read_bc_FILE_OUTPUT.c
  read_bc_SCREEN_OUTPUT.c
  read_bc_NOTERM.c
  read_bc_SERIES.c
  read_bc_AWRITE.c
  read_data_set.c
  test_cases.c
  read_bc_TEST.c
  print_sw3_to_file.c
  read_interface_geo.c
  read_flux_geo.c
  # Insert new files above here, alphabetically
  )

# SEDIMENT
#if(USE_SEDLIB)
  list(APPEND initio_srcs
    # Insert new files below here, alphabetically
    read_bc_SEDLIB.c
    # Insert new files above here, alphabetically
    )
#endif(USE_SEDLIB)
if(USE_MPI)
  list(APPEND initio_srcs
    )
endif(USE_MPI)

#xdmf
if (USE_XDMF)
  list(APPEND initio_srcs 
  ps_print_xdmf.c
  ps_print_ts_xdmf.c
  ps_print_surf_ts_xdmf.c
  xdmf_init.c
  xdmf_utils.c
  xdmf_finalize.c
  xdmf_finalize_supermodel.c
  ps_print_geo_xdmf.c
  ps_print_surf_geo_xdmf.c
)
endif(USE_XDMF)

################################################################################################
# FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. #
if(USE_WINDLIB)
   list(APPEND initio_srcs
        read_bc_WINDLIB.c
       )
endif(USE_WINDLIB)
# ABOVE LINES ADDED BY GAJANAN                                                                 #
################################################################################################

if(USE_GW)
   list(APPEND initio_srcs
        read_bc_GW.c
	read_bc_GW_MP.c
	)
endif(USE_GW)
##############################################################################

# Define include directories
include_directories("${CMAKE_SOURCE_DIR}/include"
  "${CMAKE_SOURCE_DIR}/structs"
  # include the generated adh_revision.h
  ${CMAKE_CURRENT_BINARY_DIR}
  # include the version header
  "${CMAKE_SOURCE_DIR}/main"
  )
## Output file formats
if(USE_XDMF)
  include_directories(${XML_INCLUDE_DIR} ${HDF5_INCLUDE_DIR} ${SZIP_INCLUDE_DIR} ${ZLIB_INCLUDE_DIR})
endif(USE_XDMF)

##############################################################################
# Generate adh_revision.h which contains svn revision information.
# CMake will re-link every time so that correct information is passed into the header file.
# For big projects, re-linking every time might draw the ire of developers, but we do not care,
# given the speed that ADH builds.
if(WIN32 AND NOT __COMPILER_PGI)
  add_custom_target(adh_revision
    COMMAND cmd /k call "${CMAKE_CURRENT_SOURCE_DIR}/revision_script_win.bat" "${CMAKE_CURRENT_BINARY_DIR}/adh_revision.h"
    WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}
    COMMENT "Generating ${CMAKE_CURRENT_BINARY_DIR}/adh_revision.h (this will force rebuild of adh_initio)"
    )
else(WIN32 AND NOT __COMPILER_PGI)
  add_custom_target(adh_revision
    COMMAND sh "${CMAKE_CURRENT_SOURCE_DIR}/revision_script_unix.sh" "${CMAKE_CURRENT_BINARY_DIR}/adh_revision.h"
    WORKING_DIRECTORY ${CMAKE_HOME_DIRECTORY}
    COMMENT "Generating ${CMAKE_CURRENT_BINARY_DIR}/adh_revision.h (this will force rebuild of adh_initio)"
    )
endif(WIN32 AND NOT __COMPILER_PGI)

##############################################################################

set(lib_depends adh_structs adh_elem adh_debug adh_solver adh_tools adh_elem)

# MPI
if(USE_MPI)
  list(APPEND lib_depends adh_messg adh_elem)
endif(USE_MPI)

include_directories("${CMAKE_SOURCE_DIR}/friction/include")
list(APPEND lib_depends adh_friction) # pre_adh_winds)

list(REMOVE_DUPLICATES initio_srcs)

add_library("adh_initio" STATIC EXCLUDE_FROM_ALL ${initio_srcs})

target_link_libraries(adh_initio ${lib_depends})

#if(BUILD_TESTING)
#  add_subdirectory(tests)
  #endif()
