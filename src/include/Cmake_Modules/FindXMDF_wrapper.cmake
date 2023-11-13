#
# Decide which module to use to try to find HDF5 and associated stuff
#

####### JADE ############
#if($ENV{HOSTNAME} MATCHES jade)
#  if(HDF5_LIBRARY)
#    # Do nothing; use cache
#  else(HDF5_LIBRARY)
#    set(HDF5_LIBRARY "/usr/local/usp/hdf5/1.8.4-gnu/lib/libhdf5.a" CACHE STRING "HDF5 Library")
#    set(HDF5_INCLUDE_DIR "/usr/local/usp/hdf5/1.8.4-gnu/include" CACHE STRING "HDF5 Include directory")
#  endif(HDF5_LIBRARY)

####### OTHER MACHINES ##########
#else($ENV{HOSTNAME} MATCHES jade)
  include(${CMAKE_INCLUDE_MODULES_DIR}/FindXMDF.cmake)

#endif($ENV{HOSTNAME} MATCHES jade)
