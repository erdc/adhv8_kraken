#=============================================================================
# AdH - Adaptive Hydraulics
# General executable
#=============================================================================

##############################################################################
# Build executable
# Define the source files needed to build Pre-AdH
list(APPEND adh_sources
    adh.c
    adh_func.c
    adh_cpl_func.c
    #petsc_test.c
    )

# Add AdH as build target
add_executable(adh ${adh_sources})

# Cygwin
  if(CYGWIN AND CMAKE_COMPILER_IS_GNUCC)
    # If compiling on Windows in cygwin v3.x then we need the -mno-cygwin flag
    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
      OUTPUT_VARIABLE GCC_VERSION)
    if(GCC_VERSION VERSION_LESS 4.0) # no longer supported in GCC v4
      set(compiler_flags "${compiler_flags} -mno-cygwin")
    endif()
  endif(CYGWIN AND CMAKE_COMPILER_IS_GNUCC)

  # Specify library dependencies of Pre-AdH
  target_link_libraries(adh
    adh_fe
    adh_initio
    adh_structs
    adh_columns
    adh_initio
    adh_solver
    adh_structs
    adh_fe
    adh_meteor
    adh_grid
    adh_node
    adh_elem
    adh_friction
    adh_tools
    adh_structs
    adh_debug
    adh_messg
    adh_meteor
    )
    # MPI
    if(USE_MPI)
    target_link_libraries(adh
      adh_comm
      adh_mesh
      ${MPI_LIBRARIES}
      )
    endif(USE_MPI)
  # BLAS
  if(BLAS_FOUND)
    target_link_libraries(adh ${BLAS_LIBRARIES})
  endif(BLAS_FOUND)
  #SuperLibrary
  if(USE_SUPER_LIBRARY)
    target_link_libraries(adh ${ADH_SUPER_LIBRARY})
  endif(USE_SUPER_LIBRARY)

  if(PARMETIS_FOUND)
      target_link_libraries(adh ${PARMETIS_LIBRARIES})
    endif(PARMETIS_FOUND)

  if(UMFPACK_FOUND)
      target_link_libraries(adh ${UMFPACK_LIBRARIES})
    endif(UMFPACK_FOUND)    

  if(USE_SEDLIB)
      target_link_libraries(adh adh_sed ${SED_LIBRARY})
  endif(USE_SEDLIB)
  
  if(USE_NSM)
      target_link_libraries(adh adh_nsm ${NSM_LIBRARY})
      set_property(TARGET adh PROPERTY LINKER_LANGUAGE CXX)
      set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "")
  endif(USE_NSM)

  if(USE_PETSC)
      target_link_libraries(adh ${PETSC_LIB})
  endif(USE_PETSC)

  # XDMF (Kitware/Paraview Binary)
  if(USE_XDMF)
    target_link_libraries(adh
      ${XML_LIBRARY}
      ${HDF5_LIBRARY}
      ${SZIP_LIBRARY}
      ${ZLIB_LIBRARY}
      )
     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ldl -llzma")
      message("HDF5 INCLUDE DIR: ${HDF5_INCLUDE_DIR}")
      message("HDF5 libraries: ${HDF5_LIBRARY}")
      message("XML INCLUDE DIR: ${XML_INCLUDE_DIR}")
      message("XML libraries: ${XML_LIBRARY}")
  endif(USE_XDMF)

##############################################################################

##############################################################################
# Add executable tests
if(BUILD_TESTING)
  add_subdirectory(tests)
endif(BUILD_TESTING)
