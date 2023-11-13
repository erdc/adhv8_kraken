#=============================================================================
# AdH - Adaptive Hydraulics
# windlib submodule
#=============================================================================
################################################################################################
# FOLLOWING FILES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. #

# General source files (always included)

# Windlib
if(USE_WINDLIB)
  list(APPEND wind_srcs
    # Insert new files below here, alphabetically
    global.F
    kdtree2.F
    owiwind.F
    sizes.F
    vortex.F
    wind.F 
    # Insert new files above here, alphabetically
    )
if(USE_MPI)
    list(APPEND wind_srcs messenger.F)
endif(USE_MPI)
############################################################################

  set(WIND_DEBUG " ")
  if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    set(WIND_FFLAGS " -O2 -mcmodel=medium -ffixed-line-length-none -march=k8 -m64 ")
    set(WIND_DP     " -DREAL8 -DLINUX -DCSCA -DCMPI ")
    set(WIND_DPRE   " -DREAL8 -DLINUX ")
    set(WIND_IMODS  " -I . ")
  elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
    set(WIND_FFLAGS " -O3 -FI -assume byterecl -132 -i-dynamic -xSSE4.2 -assume buffered_io ")
    set(WIND_DP     " -DREAL8 -DLINUX -DCSCA -DCMPI ")
    set(WIND_DPRE   " -DREAL8 -DLINUX ")
    set(WIND_IMODS  " -I . ")
    if(BUILD_DEBUG_LEVEL GREATER 0)
       set(WIND_DEBUG "${WIND_DEBUG} -g -O0 -traceback -debug -check all -i-dynamic -FI -assume byterecl -132 -DALL_TRACE -DFULL_STACK -DFLUSH_MESSAGES ")
    endif(BUILD_DEBUG_LEVEL GREATER 0)
  elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI")
    set(WIND_FFLAGS " -fastsse -Mextend -tp=istanbul -m64 ")
    set(WIND_DP     " -DREAL8 -DLINUX -DCMPI -DHAVE_MPI_MOD -DCSCA ")
    set(WIND_DPRE   " -DREAL8 -DLINUX ")
    set(WIND_IMODS  " -I . ")
    if(BUILD_DEBUG_LEVEL GREATER 0)
        set(WIND_DEBUG "${WIND_DEBUG} -g")
    endif(BUILD_DEBUG_LEVEL GREATER 0)
 else()
    message("-----------------------------------------------------")
    message("Cannot identify the compiler to compile wind library.")
    message("-----------------------------------------------------")

  endif()

  if(NOT USE_MPI)
    set(CMAKE_Fortran_FLAGS "${WIND_FFLAGS} ${WIND_DPRE} ${WIND_IMODS} ${WIND_DEBUG}")
  else(NOT USE_MPI)
    include_directories(${MPI_INCLUDE_PATH})
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${WIND_FFLAGS} ${WIND_DP} ${WIND_IMODS} ${WIND_DEBUG}")
  endif(NOT USE_MPI)
############################################################################

#list(REMOVE_DUPLICATES wind_srcs)
add_library("adh_wind" STATIC EXCLUDE_FROM_ALL ${wind_srcs})
if(USE_MPI)
  target_link_libraries(adh_wind ${MPI_Fortran_LIBRARIES})
endif(USE_MPI)

endif(USE_WINDLIB)
# ABOVE LINES ADDED BY GAJANAN                                                                 #
################################################################################################
