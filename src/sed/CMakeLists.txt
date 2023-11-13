#=============================================================================
# AdH - Adaptive Hydraulics
# sedlib wrapper submodule
#=============================================================================
# General source files (always included)
list(APPEND sed_srcs
  # Insert new files below here, alphabetically
  structs/ssedlib.c
  sed_vorticity_velocity_components.c
  sedlib_link_bed.c
  sedlib_link_transport.c
  # Insert new files above here, alphabetically
  )


#############################################################################
# Define include directories
include_directories("${CMAKE_SOURCE_DIR}/structs" "${CMAKE_SOURCE_DIR}/sed/include" "${CMAKE_SOURCE_DIR}/sed/structs" "${CMAKE_SOURCE_DIR}/debug/include")

##############################################################################


list(REMOVE_DUPLICATES sed_srcs)

add_library("adh_sed" STATIC EXCLUDE_FROM_ALL ${sed_srcs})

target_link_libraries(adh_sed)

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
