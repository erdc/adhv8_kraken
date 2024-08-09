#=============================================================================
# AdH - Adaptive Hydraulics
# debug submodule
#=============================================================================
# General source files (always included)
list(APPEND debug_srcs
  # Insert new files below here, alphabetically
  tl_debug_alloc.c
  tl_debug_error.c
  tl_debug_free.c
  tl_debug_init.c
  tl_debug_realloc.c
  tag.c
  # Insert new files above here, alphabetically
  )

# Optional source files (included based on specified build options)

# DEBUG
if(BUILD_DEBUG_LEVEL GREATER "1")
  list(APPEND debug_srcs
    # Insert new files below here, alphabetically
    tl_debug.c
    # Insert new files above here, alphabetically
    )
endif(BUILD_DEBUG_LEVEL GREATER "1")


#############################################################################
# Define include directories
include_directories("${CMAKE_SOURCE_DIR}/debug/include")

##############################################################################

list(REMOVE_DUPLICATES debug_srcs)

add_library("adh_debug" STATIC EXCLUDE_FROM_ALL ${debug_srcs})

if(USE_MPI)
  target_link_libraries(adh_debug ${MPI_LIBRARIES})
endif(USE_MPI)

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
