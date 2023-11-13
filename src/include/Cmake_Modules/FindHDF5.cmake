# - Find HDF5 library
# This module finds an installed library that implements the HDF5
# unsymmetric multifrontal sparse LU factorization interface 
# (see http://www.cise.ufl.edu/research/sparse/umfpack/).
#
# This module sets the following variables:
#  HDF5_FOUND - set to true if a library implementing the HDF5 interface
#    is found
#  HDF5_LIBRARY - HDF5 library to link against (using full path name)
#  HDF5_STATIC  if set on this determines what kind of linkage we do (static)
##########

# This module sets the following variables:
#  SZIP_FOUND - set to true if a library implementing the HDF5 interface
#    is found
#  SZIP_LIBRARY - HDF5 library to link against (using full path name)
#  SZIP_STATIC  if set on this determines what kind of linkage we do (static)
##########

# This module sets the following variables:
#  ZLIB_FOUND - set to true if a library implementing the HDF5 interface
#    is found
#  ZLIB_LIBRARY - HDF5 library to link against (using full path name)
#  ZLIB_STATIC  if set on this determines what kind of linkage we do (static)
##########


# Determines whether messages are printed
set(HDF5_FIND_QUIETLY FALSE)

# We want to link against the static library
set(HDF5_STATIC TRUE)

  if(WIN32)
if(HDF5_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
else(HDF5_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".dll")
  endif(HDF5_STATIC)
  elseif(APPLE)
if(HDF5_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so")
else(HDF5_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".dylib")
  endif(HDF5_STATIC)
else(WIN32)
if(HDF5_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
else(HDF5_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".dylib;.dll")
endif(HDF5_STATIC)
endif(WIN32)

# Check to see if HDF5_LIBRARY has been set
if(NOT HDF5_LIBRARY)
    # If it hasn't try to find the HDF5 library
    # Run through all of the possible suffixes and
    # break as soon as a library is found
    foreach(_suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
      find_library(HDF5_LIBRARY
				 NAMES hdf5${_suffix} libhdf5${_suffix}
				 PATHS 
				 # Generic
				 /usr/lib /usr/local/lib /opt/hdf5/gcc/1.8.6/lib /opt/hdf5/intel/1.8.3/lib 
				 # Good Place To Build on your machine (hint, hint)
                                 ~/Applications/5-1.8.2-mac-intel/lib
				 /opt/HDF5/HDF5
				 # Good Place To Build if using Cygwin
				 /usr/local/HDF5
				 # Jade Path.  Needs to be generic #CWW
                 /usr/local/usp/hdf5/1.8.4-gnu/lib/
				 # Sapphire
				)
      if(HDF5_LIBRARY)
         set(HDF5_FOUND TRUE)
         # break out of the foreach loop since we have found a library
         break(_suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
      endif(HDF5_LIBRARY)
    endforeach(_suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
else(NOT HDF5_LIBRARY)
    # HDF5_LIBRARY has been set
    set(HDF5_FOUND TRUE)
endif(NOT HDF5_LIBRARY)

if(NOT HDF5_FIND_QUIETLY)
    if(HDF5_FOUND)
      # Get the path to the library
      get_filename_component(HDF5_LIBRARY_PATH ${HDF5_LIBRARY} PATH)
      # Add the directory for the include files.
      # check to see if the HDF5 include files are in the same directory as the library
      #      set(HDF5_INCLUDE_PATH ${HDF5_LIBRARY_PATH}/../Include)
      FIND_PATH(HDF5_INCLUDE_DIR hdf5.h
                /usr/include/HDF5
                /usr/local/include/HDF5
                ~/Applications/5-1.8.2-mac-intel/include
                /opt/hdf5/gcc/1.8.6/include
                /opt/hdf5/intel/1.8.3/include
# Jade Path.  Needs to be generic #CWW
                /usr/local/usp/hdf5/1.8.4-gnu/include
                )
      message(STATUS "A library with HDF5 API found.")
    else(HDF5_FOUND)
      if(HDF5_FIND_REQUIRED)
        message(FATAL_ERROR
        "A required library with HDF5 API not found. Please specify library location.")
      else(HDF5_FIND_REQUIRED)
        message(STATUS
        "A library with HDF5 API not found. Please specify library location.")
      endif(HDF5_FIND_REQUIRED)
    endif(HDF5_FOUND)


# Determines whether messages are printed
set(SZIP_FIND_QUIETLY FALSE)

# We want to link against the static library
set(SZIP_STATIC TRUE)

if(WIN32)
  if(SZIP_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  else(SZIP_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".dll")
  endif(SZIP_STATIC)
elseif(APPLE)
  if(SZIP_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so")
  else(SZIP_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".dylib")
  endif(SZIP_STATIC)
  else(WIN32)
if(SZIP_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  else(SZIP_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".dylib")
endif(SZIP_STATIC)  
endif(WIN32)

# Check to see if SZIP_LIBRARY has been set
if(NOT SZIP_LIBRARY)
    # If it hasn't try to find the SZIP library
    # Run through all of the possible suffixes and
    # break as soon as a library is found
    foreach(_suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
      find_library(SZIP_LIBRARY
				 NAMES sz${_suffix} libsz${_suffix}
				 PATHS 
				 # Generic
				 /usr/lib /usr/local/lib /opt/szip/intel/lib
				 # Good Place To Build on your machine (hint, hint)
                                 ~/Applications/szip-2.1-mac-intel/lib
				 /opt/SZIP/SZIP
				 # Good Place To Build if using Cygwin
				 /usr/local/SZIP
				 # Jade Path Needs to be Generic # CWW
                 /usr/local/usp/szip/2.1-gnu/lib
				 # Sapphire
				)
      if(SZIP_LIBRARY)
         set(SZIP_FOUND TRUE)
         # break out of the foreach loop since we have found a library
         break(_suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
      endif(SZIP_LIBRARY)
    endforeach(_suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
else(NOT SZIP_LIBRARY)
    # SZIP_LIBRARY has been set
    set(SZIP_FOUND TRUE)
endif(NOT SZIP_LIBRARY)

if(NOT SZIP_FIND_QUIETLY)
    if(SZIP_FOUND)
      # Get the path to the library
      get_filename_component(SZIP_LIBRARY_PATH ${SZIP_LIBRARY} PATH)
      # Add the directory for the include files.
      # check to see if the SZIP include files are in the same directory as the library
      #      set(SZIP_INCLUDE_PATH ${SZIP_LIBRARY_PATH}/../Include)
      FIND_PATH(SZIP_INCLUDE_DIR szlib.h
                /usr/include/SZIP
                /usr/local/include/SZIP
                ~/Applications/szip-2.1-mac-intel/include
                /opt/szip/intel/include
# Jade Path.  Needs to be generic. # CWW
                /usr/local/usp/szip/2.1-gnu/include
                )
      message(STATUS "A library with SZIP API found.")
    else(SZIP_FOUND)
      if(SZIP_FIND_REQUIRED)
        message(FATAL_ERROR
        "A required library with SZIP API not found. Please specify library location.")
      else(SZIP_FIND_REQUIRED)
        message(STATUS
        "A library with SZIP API not found. Please specify library location.")
      endif(SZIP_FIND_REQUIRED)
    endif(SZIP_FOUND)

endif(NOT SZIP_FIND_QUIETLY)

# Determines whether messages are printed
set(ZLIB_FIND_QUIETLY FALSE)

# We want to link against the static library
set(ZLIB_STATIC TRUE)

if(WIN32)
  if(ZLIB_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  else(ZLIB_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".dll")
  endif(ZLIB_STATIC)
elseif(APPLE)
  if(ZLIB_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so")
  else(ZLIB_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".dylib")
  endif(ZLIB_STATIC)
else(WIN32)
  if(ZLIB_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  else(ZLIB_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".dylib")
  endif(ZLIB_STATIC)
endif(WIN32)

# Check to see if ZLIB_LIBRARY has been set
if(NOT ZLIB_LIBRARY)
    # If it hasn't try to find the ZLIB library
    # Run through all of the possible suffixes and
    # break as soon as a library is found
    foreach(_suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
      find_library(ZLIB_LIBRARY
				 NAMES libz${_suffix} z${_suffix} 
				 PATHS 
				 # Generic
				 /usr/lib /usr/local/lib 
				 # Good Place To Build on your machine (hint, hint)
				 /opt/ZLIB/ZLIB
				 # Good Place To Build if using Cygwin
				 /usr/local/ZLIB
				 # Jade path, needs to be Generic #CWW
                 /usr/lib64
				 # Sapphire
				)
      if(ZLIB_LIBRARY)
         set(ZLIB_FOUND TRUE)
         # break out of the foreach loop since we have found a library
         break(_suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
      endif(ZLIB_LIBRARY)
    endforeach(_suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
else(NOT ZLIB_LIBRARY)
    # ZLIB_LIBRARY has been set
    set(ZLIB_FOUND TRUE)
endif(NOT ZLIB_LIBRARY)

if(NOT ZLIB_FIND_QUIETLY)
    if(ZLIB_FOUND)
      # Get the path to the library
      get_filename_component(ZLIB_LIBRARY_PATH ${ZLIB_LIBRARY} PATH)
      # Add the directory for the include files.
      # check to see if the ZLIB include files are in the same directory as the library
      #      set(ZLIB_INCLUDE_PATH ${ZLIB_LIBRARY_PATH}/../Include)
      FIND_PATH(ZLIB_INCLUDE_DIR zlib.h
                /usr/include
                /usr/local/include/ZLIB
                )
      message(STATUS "A library with ZLIB API found.")
    else(ZLIB_FOUND)
      if(ZLIB_FIND_REQUIRED)
        message(FATAL_ERROR
        "A required library with ZLIB API not found. Please specify library location.")
      else(ZLIB_FIND_REQUIRED)
        message(STATUS
        "A library with ZLIB API not found. Please specify library location.")
      endif(ZLIB_FIND_REQUIRED)
    endif(ZLIB_FOUND)

endif(NOT ZLIB_FIND_QUIETLY)
endif(NOT HDF5_FIND_QUIETLY)
