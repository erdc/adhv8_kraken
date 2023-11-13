#=============================================================================
# AdH - Adaptive Hydraulics
# General executable
#=============================================================================

##############################################################################
# Define the source files needed to build Pre-AdH
list(APPEND adh_sources
    main.c
    extrusion.c
    create_3d_geometry.c
    create_3d_bc.c
    create_3d_hot.c
    ../XDMF_generator/geo_xdmf_write.c
    )

# Build executable
add_executable(extrude ${adh_sources})

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
  target_link_libraries(extrude
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
  # BLAS
  if(BLAS_FOUND)
      target_link_libraries(extrude ${BLAS_LIBRARIES})
  endif(BLAS_FOUND)
  #SuperLibrary
  if(USE_SUPER_LIBRARY)
      target_link_libraries(extrude ${ADH_SUPER_LIBRARY})
  endif(USE_SUPER_LIBRARY)

  if(PARMETIS_FOUND)
      target_link_libraries(extrude ${PARMETIS_LIBRARIES})
    endif(PARMETIS_FOUND)

  if(UMFPACK_FOUND)
      target_link_libraries(extrude ${UMFPACK_LIBRARIES})
    endif(UMFPACK_FOUND)    

  # XDMF (Kitware/Paraview Binary)
  if(USE_XDMF)
    target_link_libraries(extrude
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

#############################################################################
# Define include directories
include_directories("${CMAKE_SOURCE_DIR}/tools/STAND_ALONE_UTILITIES/extrusion")
include_directories("${CMAKE_SOURCE_DIR}/tools/STAND_ALONE_UTILITIES/XDMF_generator")
include_directories("${CMAKE_SOURCE_DIR}/debug/include")

##############################################################################
