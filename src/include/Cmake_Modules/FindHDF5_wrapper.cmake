#
# Decide which module to use to try to find HDF5 and associated stuff
#

####### JADE ############
if($ENV{HOSTNAME} MATCHES jade)
  if(HDF5_LIBRARY)
    # Do nothing; use cache
  else(HDF5_LIBRARY)
    set(HDF5_LIBRARY "/usr/local/usp/hdf5/1.8.4-gnu/lib/libhdf5.a" CACHE STRING "HDF5 Library")
    set(HDF5_INCLUDE_DIR "/usr/local/usp/hdf5/1.8.4-gnu/include" CACHE STRING "HDF5 Include directory")
  endif(HDF5_LIBRARY)

  if(SZIP_LIBRARY)
    # Do nothing; use cache
  else(SZIP_LIBRARY)
    set(SZIP_LIBRARY "/usr/local/usp/szip/2.1-gnu/lib/libsz.a" CACHE STRING "SZIP Library")
    set(SZIP_INCLUDE_DIR "/usr/local/usp/szip/2.1-gnu/include" CACHE STRING "SZIP Include directory")
  endif(SZIP_LIBRARY)

  if(ZLIB_LIBRARY)
    # Do nothing; use cache
  else(ZLIB_LIBRARY)
#    set(ZLIB_LIBRARY "/usr/lib64/libz.a" CACHE STRING "ZLIB Library")
    set(ZLIB_LIBRARY "/usr/local/usp/PETtools/CE/pkgs/zlib-1.2.3/lib/libz.a" CACHE STRING "ZLIB Library")
    set(ZLIB_INCLUDE_DIR "/usr/include" CACHE STRING "ZLIB Include directory")
  endif(ZLIB_LIBRARY)


####### OTHER MACHINES ##########
elseif($ENV{HOSTNAME} MATCHES garnet)
  if(HDF5_LIBRARY)
    # Do nothing; use cache
  else(HDF5_LIBRARY)
    set(HDF5_LIBRARY "/usr/local/usp/hdf5/1.8.5p1-cle/lib/libhdf5.a" CACHE STRING "HDF5 Library")
    set(HDF5_INCLUDE_DIR "/usr/local/usp/hdf5/1.8.5p1-cle/include" CACHE STRING "HDF5 Include directory")
  endif(HDF5_LIBRARY)

  if(SZIP_LIBRARY)
    # Do nothing; use cache
  else(SZIP_LIBRARY)
    set(SZIP_LIBRARY "/usr/local/usp/szip/2.1-cle/lib/libsz.a" CACHE STRING "SZIP Library")
    set(SZIP_INCLUDE_DIR "/usr/local/usp/szip/2.1-cle/include" CACHE STRING "SZIP Include directory")
  endif(SZIP_LIBRARY)

  if(ZLIB_LIBRARY)
    # Do nothing; use cache
  else(ZLIB_LIBRARY)
    set(ZLIB_LIBRARY "/usr/lib64/libz.so" CACHE STRING "ZLIB Library")
    set(ZLIB_INCLUDE_DIR "/usr/include" CACHE STRING "ZLIB Include directory")
  endif(ZLIB_LIBRARY)

####### AJ ############
elseif($ENV{HOSTNAME} MATCHES irene.ices.utexas.edu)
  if(HDF5_LIBRARY)
    # Do nothing; use cache
  else(HDF5_LIBRARY)
    set(HDF5_LIBRARY "/workspace/libraries/CMake-hdf5-1.10.1/build/bin/libhdf5.a" CACHE STRING "HDF5 Library")
    set(HDF5_INCLUDE_DIR "/workspace/libraries/CMake-hdf5-1.10.1/hdf5-1.10.1/src" CACHE STRING "HDF5 Include directory")
  endif(HDF5_LIBRARY)

  if(SZIP_LIBRARY)
    # Do nothing; use cache
  else(SZIP_LIBRARY)
    set(SZIP_LIBRARY "/workspace/libraries/szip-2.1.1/myszip/lib/libsz.a" CACHE STRING "SZIP Library")
    set(SZIP_INCLUDE_DIR "/workspace/libraries/szip-2.1.1/myszip/include" CACHE STRING "SZIP Include directory")
  endif(SZIP_LIBRARY)

  if(ZLIB_LIBRARY)
    # Do nothing; use cache
  else(ZLIB_LIBRARY)
    set(ZLIB_LIBRARY "/workspace/libraries/zlib-1.2.11/myzlib/lib/libz.a" CACHE STRING "ZLIB Library")
    set(ZLIB_INCLUDE_DIR "/workspace/libraries/zlib-1.2.11/myzlib/include" CACHE STRING "ZLIB Include directory")
  endif(ZLIB_LIBRARY)
#endif($ENV{HOSTNAME} MATCHES irene.ices.utexas.edu)

####### Gajanan ############
elseif($ENV{HOSTNAME} MATCHES harvey.ices.utexas.edu)
  if(HDF5_LIBRARY)
    # Do nothing; use cache
  else(HDF5_LIBRARY)
    set(HDF5_LIBRARY "/workspace/gajanan/adh/libraries/CMake-hdf5-1.10.1/build/bin/libhdf5.a" CACHE STRING "HDF5 Library")
    set(HDF5_INCLUDE_DIR "/workspace/gajanan/adh/libraries/CMake-hdf5-1.10.1/hdf5-1.10.1/src" CACHE STRING "HDF5 Include directory")
  endif(HDF5_LIBRARY)

  if(SZIP_LIBRARY)
    # Do nothing; use cache
  else(SZIP_LIBRARY)
    set(SZIP_LIBRARY "/workspace/gajanan/adh/libraries/szip-2.1.1/myszip/lib/libsz.a" CACHE STRING "SZIP Library")
    set(SZIP_INCLUDE_DIR "/workspace/gajanan/adh/libraries/szip-2.1.1/myszip/include" CACHE STRING "SZIP Include directory")
  endif(SZIP_LIBRARY)

  if(ZLIB_LIBRARY)
    # Do nothing; use cache
  else(ZLIB_LIBRARY)
    set(ZLIB_LIBRARY "/workspace/gajanan/adh/libraries/zlib-1.2.11/myzlib/lib/libz.a" CACHE STRING "ZLIB Library")
    set(ZLIB_INCLUDE_DIR "/workspace/gajanan/adh/libraries/zlib-1.2.11/myzlib/include" CACHE STRING "ZLIB Include directory")
  endif(ZLIB_LIBRARY)

####### Sam ############
elseif($ENV{HOSTNAME} MATCHES katrina)
  if(HDF5_LIBRARY)
    # Do nothing; use cache
  else(HDF5_LIBRARY)
    set(HDF5_LIBRARY "/workspace/sestes/ADH/adh_libs/CMake-hdf5-1.10.1/build/bin/libhdf5.a" CACHE STRING "HDF5 Library")
    set(HDF5_INCLUDE_DIR "/workspace/sestes/ADH/adh_libs/CMake-hdf5-1.10.1/hdf5-1.10.1/src" CACHE STRING "HDF5 Include directory")
  endif(HDF5_LIBRARY)

  if(SZIP_LIBRARY)
    # Do nothing; use cache
  else(SZIP_LIBRARY)
    set(SZIP_LIBRARY "/workspace/sestes/ADH/adh_libs/szip-2.1.1/myszip/lib/libsz.a" CACHE STRING "SZIP Library")
    set(SZIP_INCLUDE_DIR "/workspace/sestes/ADH/adh_libs/szip-2.1.1/myszip/include" CACHE STRING "SZIP Include directory")
  endif(SZIP_LIBRARY)

  if(ZLIB_LIBRARY)
    # Do nothing; use cache
  else(ZLIB_LIBRARY)
    set(ZLIB_LIBRARY "/workspace/sestes/ADH/adh_libs/zlib-1.2.11/myzlib/lib/libz.a" CACHE STRING "ZLIB Library")
    set(ZLIB_INCLUDE_DIR "/workspace/sestes/ADH/adh_libs/zlib-1.2.11/myzlib/include" CACHE STRING "ZLIB Include directory")
  endif(ZLIB_LIBRARY)

else($ENV{HOSTNAME} MATCHES jade)
  include(${CMAKE_INCLUDE_MODULES_DIR}/FindHDF5.cmake)

endif($ENV{HOSTNAME} MATCHES jade)
