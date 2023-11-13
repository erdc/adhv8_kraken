# - Find LIBXML library
# This module finds an installed library that implements the LIBXML 
# unsymmetric multifrontal sparse LU factorization interface 
# (see http://www.cise.ufl.edu/research/sparse/umfpack/).
#
# This module sets the following variables:
#  XML_FOUND - set to true if a library implementing the LIBXML interface
#    is found
#  XML_LIBRARY - XML library to link against (using full path name)
#  XML_STATIC  if set on this determines what kind of linkage we do (static)
##########

# Determines whether messages are printed
set(XML_FIND_QUIETLY FALSE)
set(XML_LIB "xml2.2")
set(XML_LONG_LIB "libxml2.2")
# We want to link against the static library
set(XML_STATIC TRUE)

if(WIN32)
  if(XML_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  else(XML_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".dll")
  endif(XML_STATIC)
elseif(APPLE)
  if(XML_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".2.a;.so;.2.dylib")
  else(XML_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".2.dylib")
  endif(XML_STATIC)
else(WIN32)
if(XML_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
else(XML_STATIC)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".2.dylib;.2.dll")
endif(XML_STATIC)
endif(WIN32)

# Check to see if XML_LIBRARY has been set
if(NOT XML_LIBRARY)
    # If it hasn't try to find the XML library
    # Run through all of the possible suffixes and
    # break as soon as a library is found
     foreach(_suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
      find_library(XML_LIBRARY
        NAMES xml2${_suffix} libxml2${_suffix} 

        HINTS
				 # Generic
				 /usr/lib /usr/local/lib 
				 # Good Place To Build on your machine (hint, hint)
				 /opt/XML/XML
         /opt/local/lib
				 # Good Place To Build if using Cygwin
				 /usr/local/XML
				 # Jade #CWW
         /usr/lib64
				 # Sapphire
				)
      if(XML_LIBRARY)
         set(XML_FOUND TRUE)
         # break out of the foreach loop since we have found a library
          break(_suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
      endif(XML_LIBRARY)
       endforeach(_suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
else(NOT XML_LIBRARY)
    # XML_LIBRARY has been set
    set(XML_FOUND TRUE)
endif(NOT XML_LIBRARY)

if(NOT XML_FIND_QUIETLY)
    if(XML_FOUND)
      # Get the path to the library
      get_filename_component(XML_LIBRARY_PATH ${XML_LIBRARY} PATH)
      # Add the directory for the include files.
      # check to see if the XML include files are in the same directory as the library
      #      set(XML_INCLUDE_PATH ${XML_LIBRARY_PATH}/../Include)
      FIND_PATH(XML_INCLUDE_DIR libxml/tree.h
                /usr/include/libxml2
                /usr/local/include/libxml2
                /usr/include
                )
      message(STATUS "A library with XML API found.")
    else(XML_FOUND)
      if(XML_FIND_REQUIRED)
        message(FATAL_ERROR
        "A required library with XML API not found. Please specify library location.")
      else(XML_FIND_REQUIRED)
        message(STATUS
        "A library with XML API not found. Please specify library location.")
      endif(XML_FIND_REQUIRED)
    endif(XML_FOUND)
endif(NOT XML_FIND_QUIETLY)
