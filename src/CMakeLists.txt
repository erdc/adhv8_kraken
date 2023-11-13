#=============================================================================
# AdH - Adaptive Hydraulics
# Top level CMakeList file
#=============================================================================
cmake_minimum_required(VERSION 2.8.10)

if($ENV{HOSTNAME} MATCHES ^garnet)
  set(BOX "gar")
  set(CMAKE_Fortran_COMPILER ftn)
endif($ENV{HOSTNAME} MATCHES ^garnet)

if($ENV{HOSTNAME} MATCHES ^lightning)
  set(BOX "lightning")
endif($ENV{HOSTNAME} MATCHES ^lightning)

if($ENV{HOSTNAME} MATCHES ^spirit)
  set(BOX "spirit")
endif($ENV{HOSTNAME} MATCHES ^spirit)

if($ENV{HOSTNAME} MATCHES ^topaz)
  set(BOX "topaz")
  set(CMAKE_C_COMPILER icc)
  set(CMAKE_Fortran_COMPILER ifort)
endif($ENV{HOSTNAME} MATCHES ^topaz)

if($ENV{HOSTNAME} MATCHES ^onyx)
  set(BOX "onyx")
  set(CMAKE_C_COMPILER cc)
  set(CMAKE_Fortran_COMPILER ftn)
endif($ENV{HOSTNAME} MATCHES ^onyx)

if($ENV{HOSTNAME} MATCHES ^phoenix)
  set(BOX "phoenix")
  set(CMAKE_C_COMPILER icc)
  set(CMAKE_Fortran_COMPILER ifort) 
endif($ENV{HOSTNAME} MATCHES ^phoenix)

# Include C and Fortran language support so LAPACK/BLAS are found by the 
# FindBLAS.cmake (etc.) modules provided by cmake.
project(adh C Fortran)

# Specify location of modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/include/Cmake_Modules)

# Include CTest support
set(BUILD_TESTING OFF CACHE BOOL "Build the testing tree.")
include(CTest)

# Include Directory, needed for main/adh.c
# Only make utilities if PETSc not defined
if(NOT USE_PETSC)
include_directories("${CMAKE_SOURCE_DIR}/include" "${CMAKE_SOURCE_DIR}/structs" "${CMAKE_SOURCE_DIR}/debug/include" "${CMAKE_SOURCE_DIR}/sed/include " "${CMAKE_SOURCE_DIR}/sed/structs" "${CMAKE_SOURCE_DIR}/fe" "${CMAKE_SOURCE_DIR}/fe/sw2" "${CMAKE_SOURCE_DIR}/sed/sw3" "${CMAKE_SOURCE_DIR}/fe/transport" "${CMAKE_SOURCE_DIR}/fe/sediment" "${CMAKE_SOURCE_DIR}/sed/diffusive_wave" "${CMAKE_SOURCE_DIR}/fe/gw" "${CMAKE_SOURCE_DIR}/STAND_ALONE_UTILITIES/extrusion" "${CMAKE_SOURCE_DIR}/STAND_ALONE_UTILITIES/flumeBuilder" "${CMAKE_SOURCE_DIR}/STAND_ALONE_UTILITIES/card_numberer" "${CMAKE_SOURCE_DIR}/STAND_ALONE_UTILITIES/XDMF_generator" "${CMAKE_SOURCE_DIR}/STAND_ALONE_UTILITIES/grid_interpolator")
else()
include_directories("${CMAKE_SOURCE_DIR}/include" "${CMAKE_SOURCE_DIR}/structs" "${CMAKE_SOURCE_DIR}/debug/include" "${CMAKE_SOURCE_DIR}/sed/include " "${CMAKE_SOURCE_DIR}/sed/structs" "${CMAKE_SOURCE_DIR}/fe" "${CMAKE_SOURCE_DIR}/fe/sw2" "${CMAKE_SOURCE_DIR}/sed/sw3" "${CMAKE_SOURCE_DIR}/fe/transport" "${CMAKE_SOURCE_DIR}/fe/sediment" "${CMAKE_SOURCE_DIR}/sed/diffusive_wave" "${CMAKE_SOURCE_DIR}/fe/gw" "${CMAKE_SOURCE_DIR}/STAND_ALONE_UTILITIES/card_numberer")
endif(NOT USE_PETSC)

# Specify binary and library output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Specify common location
set(ADH_LIB_DIR "${CMAKE_SOURCE_DIR}/lib/$ENV{BOX}" CACHE PATH "Directory where AdH libraries (e.g. ICM) are located (used to automatically find libraries)")

# The WINDOWS define is needed for functions like isnan, sleep, etc.
# The _CRT_SECURE_NO_DEPRECATE is needed to hide deprecated warnings
if(WIN32) # is TRUE on Windows; no longer includes CygWin starting with CMake v2.8.4
  add_definitions(-DWINDOWS -D_CRT_SECURE_NO_DEPRECATE)
  # Extra define for Windows PGI compiler
  if(__COMPILER_PGI)
    add_definitions(-D_PC_PG)
  endif(__COMPILER_PGI) 
endif(WIN32)

if(CYGWIN AND CMAKE_COMPILER_IS_GNUCC)
  execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
    OUTPUT_VARIABLE GCC_VERSION)
  if(GCC_VERSION VERSION_LESS 4.0)
    # Set the no cygwin flag for the compiler, no longer supported in GCC v4
    set(CMAKE_C_FLAGS "-mno-cygwin")
  endif()
endif(CYGWIN AND CMAKE_COMPILER_IS_GNUCC)

### Spirit specific link libraries
if($ENV{HOSTNAME} MATCHES ^spirit)
  set(CMAKE_EXE_LINKER_FLAGS "-lpthread -lifcore -lrt -lstdc++")
endif($ENV{HOSTNAME} MATCHES ^spirit)

### Spirit specific link libraries
if($ENV{HOSTNAME} MATCHES ^topaz)
  set(CMAKE_EXE_LINKER_FLAGS "-lpthread -lifcore -mkl -lmpi -lstdc++ ")
endif($ENV{HOSTNAME} MATCHES ^topaz)

### Spirit specific link libraries
if($ENV{HOSTNAME} MATCHES ^phoenix)
  set(CMAKE_EXE_LINKER_FLAGS "-mkl -lrt")
endif($ENV{HOSTNAME} MATCHES ^phoenix)

##############################################################################
# Check build type (prevent misspelling)
if(NOT CMAKE_BUILD_TYPE MATCHES "^(|Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(FATAL_ERROR "CMAKE_BUILD_TYPE parameter should be left empty or set to Debug (-D_DEBUG -g), Release, or RelWithDebInfo (-g).")
endif(NOT CMAKE_BUILD_TYPE MATCHES "^(|Debug|Release|RelWithDebInfo|MinSizeRel)$")
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  add_definitions(-D_DEBUG -g)
elseif(CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo")
  add_definitions(-g)
endif(CMAKE_BUILD_TYPE MATCHES "Debug")

# Build options for
# Determine what level of debugging we are using
set(BUILD_DEBUG_LEVEL "0" CACHE STRING
  "0: None  1: -g  2: -g + _DEBUG  3: Level 2 + Temporary  4: Level 2 + Kitchen Sink  5: All Debugging on")
# Determine if we are using some debugging so that debuggers can show the code
if(BUILD_DEBUG_LEVEL LESS 0 OR BUILD_DEBUG_LEVEL GREATER 5)
  message(FATAL_ERROR "BUILD_DEBUG_LEVEL parameter must be an integer value of the range 0 to 5, inclusive.")
elseif(BUILD_DEBUG_LEVEL EQUAL 1)
  add_definitions(-g)
# Determine if we are using the memory checking debugging
elseif(BUILD_DEBUG_LEVEL EQUAL 2)
  add_definitions(-D_DEBUG -g)
# Determine if we are using temporary user-defined debug functions
elseif(BUILD_DEBUG_LEVEL EQUAL 3)
  add_definitions(-D_DEBUG_TEMPORARY -D_DEBUG -g)
# Determine if we are using full-check-everything debugging
elseif(BUILD_DEBUG_LEVEL EQUAL 4)
  add_definitions(-D_DEBUG_KITCHEN_SINK -D_DEBUG -g)
# Determine if we are using all debugging potential
elseif(BUILD_DEBUG_LEVEL EQUAL 5)
  add_definitions(-D_DEBUG_KITCHEN_SINK -D_DEBUG_TEMPORARY -D_DEBUG -g)
endif(BUILD_DEBUG_LEVEL LESS 0 OR BUILD_DEBUG_LEVEL GREATER 5)
# Check whether the build type and debug level are consistent
if((BUILD_DEBUG_LEVEL EQUAL 0 AND CMAKE_BUILD_TYPE MATCHES "^(Debug|RelWithDebInfo)$") 
OR (BUILD_DEBUG_LEVEL GREATER 0 AND CMAKE_BUILD_TYPE MATCHES "^(Release|MinSize)$"))
  message(FATAL_ERROR "BUILD_DEBUG_LEVEL is inconsistent with CMAKE_BUILD_TYPE; both parameters should include or exclude debug definitions.")
endif((BUILD_DEBUG_LEVEL EQUAL 0 AND CMAKE_BUILD_TYPE MATCHES "^(Debug|RelWithDebInfo)$") OR
   (BUILD_DEBUG_LEVEL GREATER 0 AND CMAKE_BUILD_TYPE MATCHES "^(Release|MinSize)$"))

##############################################################################
# Optional modules (internal to AdH)

# Determine if we are using STRUCTURES 
option(USE_ADH_STRUCTURES "Include structure processes (SW)" OFF)
if(USE_ADH_STRUCTURES)
  add_definitions(-D_ADH_STRUCTURES)
endif(USE_ADH_STRUCTURES)

# Determine if we are using BREACH 
option(USE_BREACH "Include breach processes (SW)" OFF)
if(USE_BREACH)
  add_definitions(-D_ADH_BREACH)
endif(USE_BREACH)

# Determine if we are using ICM
option(USE_ICM "Include ICM library" OFF)
if(USE_ICM)
  add_definitions(-D_ADH_ICM)
  # find ICM library
  find_library(ICM_LIBRARY libicm.a "${LIB_ICM}")
  if(ICM_LIBRARY AND NOT EXISTS "${ICM_LIBRARY}")
    message(SEND_ERROR "Missing or bad ICM_LIBRARY; Specify library or turn off PACKAGE_ICM")
  endif(ICM_LIBRARY AND NOT EXISTS "${ICM_LIBRARY}")
endif(USE_ICM)

# Determine if we are using NSM
option(USE_NSM "Include NSM library" OFF)
if(USE_NSM)
  add_definitions(-D_ADH_NSM)
  # find the NSM library
  find_library(NSM_LIBRARY libnsm.a "${ADH_LIBDIR}")
  if(NSM_LIBRARY AND NOT EXISTS "${NSM_LIBRARY}")
    message(SEND_ERROR "Missing or bad NSM_LIBRARY; Specify library or turn off PACKAGE_NSM")
  endif(NSM_LIBRARY AND NOT EXISTS "${NSM_LIBRARY}")
endif(USE_NSM)


# Determine if we are using SEDIMENT LIBRARY
option(USE_SEDLIB "Include sediment transport via SEDLIB" OFF)
if(USE_SEDLIB)
    add_definitions(-D_SEDIMENT)
    add_definitions(-D_SEDLIB)
  # find the Sediment library
  find_library(SED_LIBRARY libsed.a "${ADH_LIB_DIR}")
  if(SED_LIBRARY AND NOT EXISTS "${SED_LIBRARY}")
    message(SEND_ERROR "Missing or bad SED_LIBRARY; Specify library or turn off PACKAGE_SEDLIB")
  endif(SED_LIBRARY AND NOT EXISTS "${SED_LIBRARY}")
endif(USE_SEDLIB)

# Determine if we are using SEDIMENT LIBRARY
option(USE_PETSC "USE PETSC FOR MATRIX SOLVE" OFF)
if(USE_PETSC)
	FIND_PACKAGE(PETSc)
  add_definitions(-DUSE_PETSC)
  add_definitions(-D_PETSC)
	include_directories(${PETSC_INC})
endif(USE_PETSC)

# Determine if we are using the Sediment library
#option(USE_SEDLIB "Include sediment processes" OFF)
#if(USE_SEDLIB)
#  add_definitions(-D_SEDIMENT)
#  # find the Sediment library
#  find_library(SED_LIBRARY libsed.a "${ADH_LIB_DIR}")
#  if(SED_LIBRARY AND NOT EXISTS "${SED_LIBRARY}")
#    message(SEND_ERROR "Missing or bad SED_LIBRARY; Specify library or turn off PACKAGE_SEDLIB")
#  endif(SED_LIBRARY AND NOT EXISTS "${SED_LIBRARY}")
#endif(USE_SEDLIB)

# Determine if we are using NSM
option(USE_GW "Include Groundwater Model" OFF)
if(USE_GW)
  add_definitions(-D_ADH_GROUNDWATER)
endif(USE_GW)

##############################################################################
# Determine if we are building with XDMF support
set(CMAKE_INCLUDE_MODULES_DIR ${CMAKE_SOURCE_DIR}/include/Cmake_Modules)
option(USE_XDMF "Include XDMF I/O capability" OFF) 
if(USE_XDMF)
  add_definitions(-D_ADH_HDF5)
  include(${CMAKE_INCLUDE_MODULES_DIR}/FindHDF5_wrapper.cmake)
  include(${CMAKE_INCLUDE_MODULES_DIR}/FindLIBXML_wrapper.cmake)
  #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ldl -lz -llzma ")
  include_directories(${XML_INCLUDE_DIR} ${HDF5_INCLUDE_DIR} ${SZIP_INCLUDE_DIR} ${ZLIB_INCLUDE_DIR})
  #find_package(HDF5 REQUIRED)
  #find_package(LibXml2 REQUIRED)
endif(USE_XDMF)

##############################################################################
# Determine if we want to perform extra input file validation when reading data
set(BUILD_IO_VALIDATION ON)
#option(BUILD_IO_VALIDATION "Include extra input file validation for reading data fields" ON)
if(BUILD_IO_VALIDATION)
  add_definitions(-DIO_VALIDATION)
endif(BUILD_IO_VALIDATION)

##############################################################################
# Math libraries

# For BLAS and LAPACK, use the FindBLAS.cmake module (etc.). It should
# find these in the ${ADH_DIR}/include/Cmake_Modules directory, from where Cmake was
# compiled. A local version can be created and if found (using the
# absolute path) it will override the one provided by CMake. These
# require that Fortran be specified in the project line above.

# Determine if BLAS can be found and used
option(USE_BLAS "Include BLAS package" OFF)
mark_as_advanced(USE_BLAS)
if(USE_BLAS)
  option(BLAS_INTEL_COMPILER "If using intel compiler, extra flags are needed for BLAS linking" OFF) 
  mark_as_advanced(BLAS_INTEL_COMPILER)
  # If using the intel compiler, extra flags are required to link against BLAS
  if(BLAS_INTEL_COMPILER)
    # Intel MKL Libraries purchased for the Mac, this tells FindBLAS.cmake
    # to look for them. Inherited by LAPACK as well.
    set(BLA_VENDOR "Intel10_64lp")
  endif(BLAS_INTEL_COMPILER)
  find_package(BLAS REQUIRED)
  if(BLAS_FOUND)
    add_definitions(-D_BLAS)
  else(BLAS_FOUND) 
    message(FATAL_ERROR "BLAS NOT FOUND")
  endif(BLAS_FOUND)
endif(USE_BLAS)

# Determine if LAPACK can be found and used
option(USE_LAPACK "Include LAPACK package" OFF)
mark_as_advanced(USE_LAPACK)
if(USE_LAPACK)
  find_package(LAPACK REQUIRED)
  if(LAPACK_FOUND)
    add_definitions(-D_LAPACK)
  else(LAPACK_FOUND) 
    message(FATAL_ERROR "LAPACK NOT FOUND")
  endif(LAPACK_FOUND)
endif(USE_LAPACK)

# Determine if we have a super library, if not use individual packages
if(WIN32)
  option(USE_SUPER_LIBRARY "Use the single super library" OFF)
else(WIN32)
  option(USE_SUPER_LIBRARY "Use the single super library" ON)
endif(WIN32)

if(USE_SUPER_LIBRARY)
  if("${BOX}" MATCHES "gar")
   string(TOLOWER "$ENV{PE_ENV}" PRG_ENV)
   find_library(ADH_SUPER_LIBRARY "adh_${BOX}_${PRG_ENV}")
  else("${BOX}" MATCHES "gar")
      find_library(ADH_SUPER_LIBRARY "adh_${BOX}")
      find_library(ADH_SUPER_LIBRARY
          NAMES "adh_${BOX}" libadh_mac.a
          HINTS "/usr/local/lib /usr/lib"
          )
  endif("${BOX}" MATCHES "gar")
  set(UMFPACK_VERSION 5)
  set(UMFPACK_INT_SIZE 32)
  FIND_PATH(UMFPACK_INCLUDE_DIRS umfpack.h
     HINTS ${CMAKE_INCLUDE_PATH}/UMFPACK
     ${CMAKE_SOURCE_DIR}/include/UMFPACK)
  set(umfpack_compiler_flags "-D_UMFPACK -D_UMFPACK_INT_SIZE=${UMFPACK_INT_SIZE} -D_UMFPACK_VERSION=${UMFPACK_VERSION}")
  if(UMFPACK_INT_SIZE EQUAL 64)
    set(umfpack_compiler_flags "${umfpack_compiler_flags} -D_UMFPACK_LONG")
  endif(UMFPACK_INT_SIZE EQUAL 64)
  add_definitions(${umfpack_compiler_flags})
  include_directories(${UMFPACK_INCLUDE_DIRS})
#  set(adh_compiler_flags "${adh_compiler_flags} -D_METIS")
#  FIND_PATH(PARMETIS_INCLUDE_DIR parmetis.h
#    HINTS ${CMAKE_INCLUDE_PATH}
#    ${CMAKE_SOURCE_DIR}/include)
#  include_directories(${PARMETIS_INCLUDE_DIR})
  if("${BOX}" MATCHES "(topaz|gar|spirit|lightning|phoenix|onyx)")
  else("${BOX}" MATCHES "(topaz|gar|spirit|lightning|phoenix|onyx)")
    find_package(BLAS REQUIRED)
  endif("${BOX}" MATCHES "(topaz|gar|spirit|lightning|phoenix|onyx)")

else(USE_SUPER_LIBRARY)
# Determine if we are using UMFPACK functions
option(USE_UMFPACK "Include UMFPACK package" OFF)
if(USE_UMFPACK)
  # Require the UMFPACK library to be found
  find_package(UMFPACK REQUIRED)
  if(UMFPACK_FOUND)
    # Define the compiler flags
    set(umfpack_compiler_flags "-D_UMFPACK -D_UMFPACK_INT_SIZE=${UMFPACK_INT_SIZE} -D_UMFPACK_VERSION=${UMFPACK_VERSION}")
    if(UMFPACK_INT_SIZE EQUAL 64)
      set(umfpack_compiler_flags "${umfpack_compiler_flags} -D_UMFPACK_LONG")
    endif(UMFPACK_INT_SIZE EQUAL 64)
    include_directories(${UMFPACK_INCLUDE_DIRS})
    add_definitions(${umfpack_compiler_flags})
  else(UMFPACK_FOUND)
    message(FATAL_ERROR "UMFPACK NOT FOUND")
  endif(UMFPACK_FOUND)
endif(USE_UMFPACK)

# Determine if we are using PARMETIS functions
option(USE_PARMETIS "Include ParMETIS package" OFF)
if(USE_PARMETIS)
  if(NOT USE_MPI)
    message(FATAL_ERROR "Cannot use ParMETIS without MPI. Turn on PACKAGE_MPI")
  endif(NOT USE_MPI)
  find_package(PARMETIS REQUIRED)
  include_directories(${PARMETIS_INCLUDE_DIRS})
  # Require the PARMETIS library to be found
  #find_package(PARMETIS REQUIRED)
  #include(${CMAKE_INCLUDE_MODULES_DIR}/FindPARMETIS.cmake)
  # Define the compiler flags
  add_definitions(-D_METIS)
endif(USE_PARMETIS)
endif(USE_SUPER_LIBRARY)

##############################################################################

# Determine if we are using MPI functions
if("${BOX}" MATCHES "(topaz|gar|spirit|lightning|phoenix)")
  option(USE_MPI "Include MPI package" ON)
else("${BOX}" MATCHES "(topaz|gar|spirit|lightning|phoenix)")
  option(USE_MPI "Include MPI package" OFF)
endif("${BOX}" MATCHES "(topaz|gar|spirit|lightning|phoenix)")

# If we are using MPI then try to find and set the MPI_COMPILER
if(USE_MPI)
  # Require the MPI compiler to be found
  find_package(MPI REQUIRED) 
  # Define the compiler flags
  add_definitions(-D_MPI -D_MESSG)
  include_directories(${MPI_C_INCLUDE_PATH})
endif(USE_MPI)

################################################################################################
# FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. #
option(USE_WINDLIB "Use ADCIRC wind library via WINDLIB" OFF)
if(USE_WINDLIB)
  if(USE_MPI)
        include_directories(${MPI_INCLUDE_PATH})
        set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})       ## WARNING!!!! THIS NEEDS TO BE CHECKED!!!
        add_definitions(-DCMPI)
  endif(USE_MPI)

  if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    #add_definitions(-DHAVE_MPI_MOD)
    add_definitions(-DWIND_GNU_NM)
  elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
    add_definitions(-DWIND_INTEL_NM)
  elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "PGI")
    add_definitions(-DWIND_PGI_NM)
  else()
    message("-----------------------------------------------------")
    message("Cannot identify the compiler to compile wind library.")
    message("${CMAKE_Fortran_COMPILER_ID}")
    message("-----------------------------------------------------")
  endif()

  message("Turning WINDLIB macro on.")
  add_definitions(-DWINDLIB)
endif(USE_WINDLIB)

if(NOT WIN32)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lm")
endif(NOT WIN32)

# ABOVE LINES ADDED BY GAJANAN                                                                 #
################################################################################################
################################################################################################
# FOLLOWING LINES ADDED BY GAJANAN [ gkc May 2020 ]. These are for DW-GW coupling.             #
option(USE_DWGW "Include Diffusive wave - Groundwater coupling" OFF)
if(USE_DWGW)
    include_directories("${CMAKE_SOURCE_DIR}/fe/dwgw_hybrid")
    add_definitions(-D_DWGW_COUPLING)
endif(USE_DWGW)
# ABOVE LINES ADDED BY GAJANAN                                                                 #
################################################################################################

##############################################################################
# Testing
# Write custom information file for CTest
set(ctestcustom_file
# text of file; do not ident this block as formatting is literal
"#This file sets the maximum number of warnings and errors to display in the dashboard
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 500)
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS 500)"
  )
file(WRITE "${CMAKE_BINARY_DIR}/CTestCustom.cmake" "${ctestcustom_file}")

if(BUILD_TESTING)
  set(ADH_TEST_DIR $ENV{ADH_TEST_DIR} CACHE PATH "Directory where AdH tests are located")
  if(NOT EXISTS "${ADH_TEST_DIR}")
    message(FATAL_ERROR "Testing directory (ADH_TEST_DIR) must be specified and exist to run tests (or turn off BUILD_TESTING to avoid tests)")
  endif(NOT EXISTS "${ADH_TEST_DIR}")
  # Specify the level of tests to include.
  #   BASIC: include very low level tests that do not check output (Is executable usage
  #     note displayed when no agruments are given?).
  #   STANDARD: basic tests with simple tests that verify "toy problem" output (compare
  #     results with analytic solution).
  #   ADVANCED: basic and standard tests with large tests that verify real simulation
  #     output (these tests will require extra time to run).
  #   CUSTOM: tests that are under development or are otherwise special and should be
  #     separated from other tests.
  #   ALL: all basic, standard, advanced, and custom tests.
  set(BUILD_TEST_LEVEL "STANDARD" CACHE STRING "Specify test suite: \"BASIC\", \"STANDARD\", \"ADVANCED\", \"CUSTOM\", or \"ALL\"")
  string(TOUPPER "${BUILD_TEST_LEVEL}" BUILD_TEST_LEVEL)
  if(NOT "${BUILD_TEST_LEVEL}" MATCHES "^(BASIC|STANDARD|ADVANCED|CUSTOM|ALL)$")
    message(FATAL_ERROR "BUILD_TEST_LEVEL parameter must be set to BASIC, STANDARD, ADVANCED, CUSTOM, or ALL.")
  endif(NOT "${BUILD_TEST_LEVEL}" MATCHES "^(BASIC|STANDARD|ADVANCED|CUSTOM|ALL)$")
endif(BUILD_TESTING)

##############################################################################
# Dashboard
# Set project dashboard label
set(DASHBOARD_LABEL "Restructure" CACHE STRING "Project label for the dashboard")
mark_as_advanced(DASHBOARD_LABEL)

##############################################################################
# Add Subdirectories
# Each Subdirectory should have a rule set of its own (CMakeLists.txt file)
add_subdirectory(fe)
add_subdirectory(fe/sw2)
add_subdirectory(fe/sw3)
add_subdirectory(fe/ns3)
add_subdirectory(fe/sw_hybrid)
add_subdirectory(fe/diffusive_wave)
add_subdirectory(fe/transport)
add_subdirectory(fe/sediment)
add_subdirectory(debug)
add_subdirectory(elem)
add_subdirectory(turbulence)
add_subdirectory(initio)
add_subdirectory(solver)
add_subdirectory(tools)
add_subdirectory(messg)
add_subdirectory(friction)
add_subdirectory(grid)
add_subdirectory(columns)
add_subdirectory(structs)
add_subdirectory(meteor)
if(USE_ADH_STRUCTURES)
add_subdirectory(hydro_structures)
endif(USE_ADH_STRUCTURES)
if(USE_SEDLIB)
add_subdirectory(sed)
endif(USE_SEDLIB)
add_subdirectory(node)
if(USE_MPI)
  add_subdirectory(comm)
  add_subdirectory(mesh)
endif(USE_MPI)
if(USE_GW)
add_subdirectory(fe/gw)
endif(USE_GW)


if(NOT USE_PETSC)
add_subdirectory(STAND_ALONE_UTILITIES/extrusion)
add_subdirectory(STAND_ALONE_UTILITIES/flumeBuilder)
add_subdirectory(STAND_ALONE_UTILITIES/card_numberer)
add_subdirectory(STAND_ALONE_UTILITIES/grid_interpolator)
if (USE_XDMF)
add_subdirectory(STAND_ALONE_UTILITIES/XDMF_generator)
endif(USE_XDMF)
endif(NOT USE_PETSC)
add_subdirectory(main)

################################################################################################
# FOLLOWING LINES ADDED BY GAJANAN [ gkc July 2015 ]. These are for ADCIRC wind library usage. #
if(USE_WINDLIB)
    add_subdirectory(windlib)
endif(USE_WINDLIB)
# ABOVE LINES ADDED BY GAJANAN                                                                 #
################################################################################################
################################################################################################
# FOLLOWING LINES ADDED BY GAJANAN [ gkc May 2020 ]. These are for DW-GW coupling.             #
if(USE_DWGW)
    add_subdirectory(fe/dwgw_hybrid)
endif(USE_DWGW)
# ABOVE LINES ADDED BY GAJANAN                                                                 #
################################################################################################

