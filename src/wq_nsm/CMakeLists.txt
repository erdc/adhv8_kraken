#===========================================
# AdH - Adaptive Hydraulics
# NSM WQ submodule
#===========================================
# General source files (always included)
if(USE_NSM)
list(APPEND nsm_srcs
  # Insert new files below
  NSMWQ.c
  )
 
 
#############################################################################
# Define include directories
include_directories("${CMAKE_SOURCE_DIR}/structs" "${CMAKE_SOURCE_DIR}/wq_nsm/include")

##############################################################################
 list(REMOVE_DUPLICATES nsm_srcs)
 add_library("adh_nsm" STATIC EXCLUDE_FROM_ALL ${nsm_srcs})
 
 target_link_libraries(adh_nsm)

  
endif(USE_NSM)
