# - Find UMFPACK library
# This module finds an installed library that implements the UMFPACK 
# unsymmetric multifrontal sparse LU factorization interface 
# (see http://www.cise.ufl.edu/research/sparse/umfpack/).
#
# This module sets the following variables:
#  UMFPACK_FOUND - set to true if a library implementing the UMFPACK interface
#    is found
#  UMFPACK_LIBRARY - UMFPACK library to link against (using full path name)
#  UMFPACK_STATIC  if set on this determines what kind of linkage we do (static)
##########

# Determines whether messages are printed
set(UMFPACK_FIND_QUIETLY FALSE)

# We want to link against the static library
set(UMFPACK_STATIC TRUE)

# Query user for UMFPACK Version
set(UMFPACK_VERSION 5 CACHE STRING "Major version number of UMFPACK, i.e. 3.x would be entered as 3")
if(NOT "${UMFPACK_VERSION}" MATCHES "^(3|4|5)$")
  message(FATAL_ERROR "Only UMFPACK 3.x, 4.x, and 5.x are currently supported; UMFPACK_VERSION must be \"3\", \"4\", or \"5\"")
endif(NOT "${UMFPACK_VERSION}" MATCHES "^(3|4|5)$")
#set(UMFPACK_VERSION_LAST 0 CACHE STRING "For Updating")
set(UMFPACK_INT_SIZE 32 CACHE STRING "Select 32 for integers, 64 for longs")
if(NOT "${UMFPACK_INT_SIZE}" MATCHES "^(32|64)$")
  message(FATAL_ERROR "UMFPACK_INT_SIZE must be \"32\" or \"64\"")
endif(NOT "${UMFPACK_INT_SIZE}" MATCHES "^(32|64)$")

if(UMFPACK_VERSION EQUAL 3)
  if (NOT UMFPACK_VERSION_LAST EQUAL UMFPACK_VERSION)
    set(UMFPACK_VERSION_LAST ${UMFPACK_VERSION} CACHE INTERNAL "For Updating" FORCE)
    unset(UMFPACK_AMD_INCLUDE_DIR CACHE)
    unset(UMFPACK_AMD_LIBRARY CACHE)
    unset(UMFPACK_UFCONFIG_INCLUDE_DIR CACHE)
    unset(UMFPACK_LIBRARY CACHE)
    unset(UMFPACK_INCLUDE_DIR CACHE)
  endif(NOT UMFPACK_VERSION_LAST EQUAL UMFPACK_VERSION)
  
  FIND_PATH(UMFPACK_INCLUDE_DIR umfpack.h
    HINTS ${CMAKE_INCLUDE_PATH}/UMFPACK 
      ${CMAKE_SOURCE_DIR}/include/UMFPACK
      /usr/local/include
      /usr/include
     ~/UMFPACK
      /opt/umfpack/UMFPACK3.2
      ~/UMFPACK/UMFPACK3.2
	)
  find_library(UMFPACK_LIBRARY NAMES libumpack umfpack
      HINTS  ${CMAKE_SOURCE_DIR}/lib/
      ${CMAKE_SOURCE_DIR}/lib/${BOX}
      /usr/local/lib
      /usr/lib
      ~/UMFPACK
      /opt/umfpack/UMFPACK3.2
      ~/UMFPACK/UMFPACK3.2
      )
  elseif(UMFPACK_VERSION EQUAL 4)
    if (NOT UMFPACK_VERSION_LAST EQUAL UMFPACK_VERSION)
      set(UMFPACK_VERSION_LAST ${UMFPACK_VERSION} CACHE INTERNAL "For updating" FORCE)
      unset(UMFPACK_AMD_INCLUDE_DIR CACHE)
      unset(UMFPACK_AMD_LIBRARY CACHE)
      unset(UMFPACK_UFCONFIG_INCLUDE_DIR CACHE)
      unset(UMFPACK_LIBRARY CACHE)
      unset(UMFPACK_INCLUDE_DIR CACHE)
    endif(NOT UMFPACK_VERSION_LAST EQUAL UMFPACK_VERSION)
  FIND_PATH(UMFPACK_INCLUDE_DIR umfpack.h
    ${CMAKE_SOURCE_DIR}/include/UMFPACK
      /usr/local/include
      /usr/include
      ~/UMFPACKv4.6/UMFPACK/Include
      ~/UMFPACKv4.4/UMFPACK/Include
      ~/UMFPACK/UMFPACKv4.4/UMFPACK/Include
      )
  FIND_LIBRARY(UMFPACK_LIBRARY umfpack
      ${CMAKE_SOURCE_DIR}/lib/${BOX}
      /usr/local/lib
      /usr/lib
      ~/UMFPACKv4.6/UMFPACK/Lib
      ~/UMFPACKv4.4/UMFPACK/Lib
      ~/UMFPACK/UMFPACKv4.4/UMFPACK/Lib
      )
  FIND_PATH(UMFPACK_AMD_INCLUDE_DIR amd.h
      ${CMAKE_SOURCE_DIR}/include/UMFPACK
      /usr/local/include/
      /usr/include
      ~/UMFPACKv4.6/AMD/Include
      ~/UMFPACKv4.4/AMD/Include
      ~/UMFPACK/UMFPACKv4.4/AMD/Include
      )
  FIND_LIBRARY(UMFPACK_AMD_LIBRARY amd
      ${CMAKE_SOURCE_DIR}/lib/${BOX}
      /usr/local/lib
      /usr/lib
      ~/UMFPACKv4.6/AMD/Lib
      ~/UMFPACKv4.4/AMD/Lib
      ~/UMFPACK/UMFPACKv4.4/AMD/Lib
      )
elseif(UMFPACK_VERSION EQUAL 5)
    if (NOT UMFPACK_VERSION_LAST EQUAL UMFPACK_VERSION)
      set(UMFPACK_VERSION_LAST ${UMFPACK_VERSION} CACHE INTERNAL "For updating" FORCE)
      unset(UMFPACK_AMD_INCLUDE_DIR CACHE)
      unset(UMFPACK_AMD_LIBRARY CACHE)
      unset(UMFPACK_UFCONFIG_INCLUDE_DIR CACHE)
      unset(UMFPACK_LIBRARY CACHE)
      unset(UMFPACK_INCLUDE_DIR CACHE)
    endif(NOT UMFPACK_VERSION_LAST EQUAL UMFPACK_VERSION)
FIND_PATH(UMFPACK_INCLUDE_DIR umfpack.h
    ${CMAKE_SOURCE_DIR}/include/UMFPACK
      /usr/local/include
      /usr/include
      ~/UMFPACKv4.6/UMFPACK/Include
      ~/UMFPACKv4.4/UMFPACK/Include
      ~/UMFPACK/UMFPACKv4.4/UMFPACK/Include
      )
  FIND_LIBRARY(UMFPACK_LIBRARY umfpack
      ${CMAKE_SOURCE_DIR}/lib/${BOX}
      /usr/local/lib
      /usr/lib
      ~/UMFPACKv4.6/UMFPACK/Lib
      ~/UMFPACKv4.4/UMFPACK/Lib
      ~/UMFPACK/UMFPACKv4.4/UMFPACK/Lib
      )
  FIND_PATH(UMFPACK_AMD_INCLUDE_DIR amd.h
      ${CMAKE_SOURCE_DIR}/include/UMFPACK
      /usr/local/include/
      /usr/include
      ~/UMFPACKv4.6/AMD/Include
      ~/UMFPACKv4.4/AMD/Include
      ~/UMFPACK/UMFPACKv4.4/AMD/Include
      )
  FIND_LIBRARY(UMFPACK_AMD_LIBRARY amd
      ${CMAKE_SOURCE_DIR}/lib/${BOX}
      /usr/local/lib
      /usr/lib
      ~/UMFPACKv4.6/AMD/Lib
      ~/UMFPACKv4.4/AMD/Lib
      ~/UMFPACK/UMFPACKv4.4/AMD/Lib
      )
    FIND_PATH(UMFPACK_SUITESPARSECONFIG_INCLUDE_DIR SuiteSparse_config.h
      ${CMAKE_SOURCE_DIR}/include/UMFPACK
      /usr/local/include/
      /usr/include
      ~/UMFPACKv4.6/AMD/Include
      ~/UMFPACKv4.4/AMD/Include
      ~/UMFPACK/UMFPACKv4.4/AMD/Include
      )
  FIND_LIBRARY(UMFPACK_SUITESPARSECONFIG_LIBRARY suitesparseconfig
      ${CMAKE_SOURCE_DIR}/lib/${BOX}
      /usr/local/lib
      /usr/lib
      ~/UMFPACKv4.6/AMD/Lib
      ~/UMFPACKv4.4/AMD/Lib
      ~/UMFPACK/UMFPACKv4.4/AMD/Lib
      )
  IF(WIN32)
    OPTION(UMFPACK_HAS_CHOLMOD "Is UMFPACK Compiled with CHOLMOD" OFF)
  ELSE(WIN32)
    OPTION(UMFPACK_HAS_CHOLMOD "Is UMFPACK Compiled with CHOLMOD" ON)
  ENDIF(WIN32)
  IF(UMFPACK_HAS_CHOLMOD)
    Find_library(UMFPACK_CHOLMOD_LIB cholmod
      ${CMAKE_SOURCE_DIR}/lib/${BOX}
      /usr/local/lib
      /usr/lib
      ~/UMFPACKv4.6/AMD/Lib
      ~/UMFPACKv4.4/AMD/Lib
      ~/UMFPACK/UMFPACKv4.4/AMD/Lib
      )
  ENDIF(UMFPACK_HAS_CHOLMOD)
  IF(WIN32)
    OPTION(UMFPACK_HAS_COLAMD "Is UMFPACK Compiled with COLAMD" OFF)
  ELSE(WIN32)
    OPTION(UMFPACK_HAS_COLAMD "Is UMFPACK Compiled with COLAMD" ON)
  ENDIF(WIN32)
  IF(UMFPACK_HAS_COLAMD)
    Find_library(UMFPACK_COLAMD_LIB colamd
      ${CMAKE_SOURCE_DIR}/lib/${BOX}
      /usr/local/lib
      /usr/lib
      ~/UMFPACKv4.6/AMD/Lib
      ~/UMFPACKv4.4/AMD/Lib
      ~/UMFPACK/UMFPACKv4.4/AMD/Lib
      )
    Find_library(UMFPACK_CCOLAMD_LIB ccolamd
      ${CMAKE_SOURCE_DIR}/lib/${BOX}
      /usr/local/lib
      /usr/lib
      ~/UMFPACKv4.6/AMD/Lib
      ~/UMFPACKv4.4/AMD/Lib
      ~/UMFPACK/UMFPACKv4.4/AMD/Lib
      )
    Find_library(UMFPACK_CAMD_LIB camd
      ${CMAKE_SOURCE_DIR}/lib/${BOX}
      /usr/local/lib
      /usr/lib
      ~/UMFPACKv4.6/AMD/Lib
      ~/UMFPACKv4.4/AMD/Lib
      ~/UMFPACK/UMFPACKv4.4/AMD/Lib
      )
  ENDIF(UMFPACK_HAS_COLAMD)
  IF(WIN32)
    OPTION(UMFPACK_HAS_METIS "Is UMFPACK Compiled with METIS" OFF)
  ELSE(WIN32)
    OPTION(UMFPACK_HAS_METIS "Is UMFPACK Compiled with METIS" ON)
  ENDIF(WIN32)
  IF(UMFPACK_HAS_METIS)
    IF(NOT METIS_LIBRARY)
      find_library(METIS_LIBRARY metis
        ${CMAKE_SOURCE_DIR}/lib/${BOX}
        /usr/local/lib
        /usr/lib
        ~/UMFPACKv4.6/AMD/Lib
        ~/UMFPACKv4.4/AMD/Lib
        ~/UMFPACK/UMFPACKv4.4/AMD/Lib
      )  
    ENDIF(NOT METIS_LIBRARY)
  ENDIF(UMFPACK_HAS_METIS)
  IF(WIN32)
    OPTION(UMFPACK_HAS_BLAS "Is UMFPACK Compiled with BLAS" OFF)
  ELSE(WIN32)
    OPTION(UMFPACK_HAS_BLAS "Is UMFPACK Compiled with BLAS" ON)
  ENDIF(WIN32)
  IF(UMFPACK_HAS_BLAS)
	FIND_PACKAGE(BLAS REQUIRED)
  ENDIF(UMFPACK_HAS_BLAS)
endif(UMFPACK_VERSION EQUAL 3) 

set(UMFPACK_FOUND FALSE)
IF(UMFPACK_LIBRARY)
  IF(UMFPACK_INCLUDE_DIR)
    SET(UMFPACK_LIBRARIES ${UMFPACK_LIBRARY} ${UMFPACK_AMD_LIBRARY} ${UMFPACK_SUITESPARSECONFIG_LIBRARY} ${UMFPACK_CHOLMOD_LIB} ${UMFPACK_COLAMD_LIB} ${BLAS_LIBRARIES} ${UMFPACK_CAMD_LIB} ${UMFPACK_CCOLAMD_LIB} ${METIS_LIBRARY})
    SET(UMFPACK_INCLUDE_DIRS ${UMFPACK_INCLUDE_DIR} ${UMFPACK_SUITESPARSECONFIG_INCLUDE_DIR} ${UMFPACK_AMD_INCLUDE_DIR})
    set(UMFPACK_FOUND TRUE)
  ENDIF(UMFPACK_INCLUDE_DIR)
ENDIF(UMFPACK_LIBRARY)
 
