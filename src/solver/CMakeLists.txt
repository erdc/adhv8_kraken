#=============================================================================
# AdH - Adaptive Hydraulics
# solv submodule
#=============================================================================
# General source files (always included)
list(APPEND solver_srcs
  # Insert new files below here, alphabetically
    solv_blk_UMFPACK.c
    solv_blk_load_profile.c
    solv_blk_load_sparse.c
    solv_blk_profile.c
    solv_blk_set.c
    solv_blk_sparse.c
    solv_blk_solve.c
    solv_initialize.c
    spv_alloc.c
    spv_dot.c
    spv_init.c
    spv_jacobi.c
    spv_load.c
    spv_vscale.c
    spv_norm.c
    spv_reset.c
    bv_alloc.c
    bv_dot.c
    bv_init.c
    bv_load.c
    solv_prec.c
    solv_amult.c
    solv_bcgstab.c
    solv_return_old.c
    #solv_coarse_solve.c
    solv_linear_sys_setup.c
    solv_linear_sys_solve.c
  # Insert new files above here, alphabetically
  )


###########################################################################


list(REMOVE_DUPLICATES solver_srcs)

add_library("adh_solver" STATIC EXCLUDE_FROM_ALL ${solver_srcs})

set(lib_depends)
# MPI
if(USE_MPI)
  list(APPEND lib_depends adh_messg)
endif(USE_MPI)

if(USE_UMFPACK)
  list(APPEND lib_depends ${UMFPACK_LIBRARIES})
endif(USE_UMFPACK)


if(USE_SUPER_LIBRARY)
  list(APPEND lib_depends ${ADH_SUPER_LIBRARY})
endif(USE_SUPER_LIBRARY)

target_link_libraries(adh_solver ${lib_depends})

#if(BUILD_TESTING)
#  add_subdirectory(tests)
#endif()

