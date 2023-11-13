#=============================================================================
# AdH - Adaptive Hydraulics
# structs submodule
#=============================================================================
# General source files (always included)
list(APPEND structs_srcs
  # Insert new files below here, alphabetically
  sarray.c
  sdebug.c
  sfile_output.c
  sscreen_output.c
  stestcase.c
  sdof.c
  scon.c
  sio.c
  selem_1d.c
  selem_2d.c
  selem_3d.c
  svect.c
  svect2d.c
  stensor.c
  sgrid.c
  snode.c
  sseries.c
  smodel.c
  smodel_defaults.c
  smat.c
  smpi.c
  sstr_value.c
  sns.c
  sns_2d.c
  sns_3d.c
  ssw.c
  ssw_2d.c
  ssw_3d.c
  smeteor.c
  ssuperModel.c
  ssuperinterface.c
  squad.c
  sstructures_value.c # later only include when structures are turned on
  ssediment.c # later only include when sediment is turned on
  sinterface.c
  slist_items.c
  # Insert new files above here, alphabetically
  )

################################################################################################
# FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. #
if(USE_WINDLIB)
  list(APPEND structs_srcs
       swindlib.c
      )
endif(USE_WINDLIB)
# ABOVE LINES ADDED BY GAJANAN                                                                 #
if(USE_GW)
  list(APPEND structs_srcs
       sgw.c
       )
endif(USE_GW)
################################################################################################

#############################################################################
# Define include directories

# include sedlib header (later only when option is in)
include_directories("${CMAKE_SOURCE_DIR}/sed/include")

###############################################################################

set(lib_depends adh_initio adh_grid adh_columns adh_tools adh_elem)

###############################################################################

# MPI
if(USE_MPI)
  list(APPEND lib_depends adh_messg)
endif(USE_MPI)

################################################################################################
# FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. #
if(USE_WINDLIB)
  list(APPEND lib_depends adh_wind)
endif(USE_WINDLIB)
# ABOVE LINES ADDED BY GAJANAN                                                                 #
################################################################################################

list(REMOVE_DUPLICATES structs_srcs)

add_library("adh_structs" STATIC EXCLUDE_FROM_ALL ${structs_srcs})

target_link_libraries(adh_structs ${lib_depends} )

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
