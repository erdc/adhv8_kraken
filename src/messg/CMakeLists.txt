#=============================================================================
# AdH - Adaptive Hydraulics
# messg submodule
#=============================================================================
# General source files (always included)
list(APPEND messg_srcs
  # Insert new files below here, alphabetically
    messg_comm_rank.c
    messg_comm_size.c
    messg_dmax.c
    messg_dmin.c
    messg_dsum.c
    messg_imax.c
    messg_imin.c
    messg_isum.c
    messg_max_norm_loc.c
  # Insert new files above here, alphabetically
  )
if(USE_MPI)
list(APPEND messg_srcs
  # Insert new files below here, alphabetically
    messg_alloc.c
    messg_arecv.c
    messg_asend.c
    messg_barrier.c
    messg_buffer_alloc.c
    messg_buffer_free.c
    messg_buffer_init.c
    messg_err.c
    messg_incoming.c
    messg_pack.c
    messg_pack_alloc.c
    messg_precv.c
    messg_unpack.c
    messg_wait.c
  # Insert new files above here, alphabetically
  )
endif(USE_MPI)
##############################################################################

list(REMOVE_DUPLICATES messg_srcs)

add_library("adh_messg" STATIC EXCLUDE_FROM_ALL ${messg_srcs})

#target_link_libraries()

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()
