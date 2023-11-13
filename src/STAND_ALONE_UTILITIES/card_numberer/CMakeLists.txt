#=============================================================================
# AdH - Adaptive Hydraulics
# General executable
#=============================================================================

##############################################################################
# Define the source files needed to build Pre-AdH
list(APPEND adh_sources
    card_numberer.c
    ../../initio/parse_card.c
    )

# Build executable
add_executable(card_numberer ${adh_sources})

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
  target_link_libraries(card_numberer
      adh_initio
      adh_structs
      adh_tools
      adh_structs
      adh_debug
      )

#############################################################################
# Define include directories
include_directories("${CMAKE_SOURCE_DIR}/tools/STAND_ALONE_UTILITIES/card_numbere")
include_directories("${CMAKE_SOURCE_DIR}/initio/")

##############################################################################
