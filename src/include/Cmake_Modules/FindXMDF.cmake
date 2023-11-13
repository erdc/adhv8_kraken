# - Find XMDF library
# This module finds an installed library that implements the HDF5
# unsymmetric multifrontal sparse LU factorization interface 
# (see http://www.cise.ufl.edu/research/sparse/umfpack/).
#
# This module sets the following variables:
#  XMDF_FOUND - set to true if a library implementing the HDF5 interface
#    is found
#  XMDF_LIBRARY - HDF5 library to link against (using full path name)
#  XMDF_STATIC  if set on this determines what kind of linkage we do (static)
##########

# Determines whether messages are printed
set(XMDF_FIND_QUIETLY FALSE)

# We want to link against the static library
set(XMDF_STATIC TRUE)

  if(WIN32)
if(XMDF_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
else(XMDF_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".dll")
  endif(XMDF_STATIC)
  elseif(APPLE)
if(XMDF_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so")
else(XMDF_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".dylib")
  endif(XMDF_STATIC)
else(WIN32)
if(XMDF_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
else(XMDF_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".dylib;.dll")
endif(XMDF_STATIC)
endif(WIN32)

# Check to see if XMDF_LIBRARY has been set
if(NOT XMDF_LIBRARY)
    # If it hasn't try to find the XMDF library
    # Run through all of the possible suffixes and
    # break as soon as a library is found
    foreach(_suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
      find_library(XMDF_LIBRARY
				 NAMES xmdf${_suffix} libxmdf${_suffix}
				 PATHS 
				 # Generic
				 /usr/lib /usr/local/lib /opt/xmdf/intel/1.8.3/lib 
				 # Good Place To Build on your machine (hint, hint)
                                 ~/Applications/5-1.8.2-mac-intel/lib
				 /opt/XMDF/XMDF
				 # Good Place To Build if using Cygwin
				 /usr/local/XMDF
				 # Jade Path.  Needs to be generic #CWW
                 /usr/local/usp/hdf5/1.8.4-gnu/lib/
				 # Sapphire
				)
      if(XMDF_LIBRARY)
         set(XMDF_FOUND TRUE)
         # break out of the foreach loop since we have found a library
         break(_suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
      endif(XMDF_LIBRARY)
    endforeach(_suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
else(NOT XMDF_LIBRARY)
    # XMDF_LIBRARY has been set
    set(XMDF_FOUND TRUE)
endif(NOT XMDF_LIBRARY)

if(NOT XMDF_FIND_QUIETLY)
    if(XMDF_FOUND)
      # Get the path to the library
      get_filename_component(XMDF_LIBRARY_PATH ${XMDF_LIBRARY} PATH)
      # Add the directory for the include files.
      # check to see if the XMDF include files are in the same directory as the library
      #      set(XMDF_INCLUDE_PATH ${XMDF_LIBRARY_PATH}/../Include)
      FIND_PATH(XMDF_INCLUDE_DIR xmdf.h
                /usr/include/XMDF
                /usr/local/include/XMDF
                ~/Applications/5-1.8.2-mac-intel/include
                /opt/XMDF/intel/1.8.3/include/
# Jade Path.  Needs to be generic #CWW
                /usr/local/usp/XMDF/1.8.4-gnu/include
                )
      message(STATUS "A library with XMDF API found.")
    else(XMDF_FOUND)
      if(XMDF_FIND_REQUIRED)
        message(FATAL_ERROR
        "A required library with XMDF API not found. Please specify library location.")
      else(XMDF_FIND_REQUIRED)
        message(STATUS
        "A library with XMDF API not found. Please specify library location.")
      endif(XMDF_FIND_REQUIRED)
    endif(XMDF_FOUND)

endif(NOT XMDF_FIND_QUIETLY)
