#
# Decide which module to use to try to find LIBXML
#

###### JADE #######
if($ENV{HOSTNAME} MATCHES jade)
  if(XML_LIBRARY)
    # Do nothing; use cache
  else(XML_LIBRARY)
    set(XML_LIBRARY "/usr/local/usp/CTB/codes/libxml/libxml2.a" CACHE STRING "XML Library")
     # NOTE: the version of libxml2 on jade has a bug; the slightly hacked header
     #       files are located here
    set(XML_INCLUDE_DIR "/usr/local/usp/CTB/codes/libxml" CACHE STRING "XML Include directory")
  endif(XML_LIBRARY)



##### OTHER ########
# Call the traditional module
elseif($ENV{HOSTNAME} MATCHES harvey.ices.utexas.edu)
  if(XML_LIBRARY)
    # Do nothing; use cache
  else(XML_LIBRARY)
    set(XML_LIBRARY "/workspace/gajanan/adh/libraries/libxml2-2.9.4/mylibxml/lib/libxml2.a" CACHE STRING "XML Library")
     # NOTE: the version of libxml2 on jade has a bug; the slightly hacked header
     #       files are located here
    set(XML_INCLUDE_DIR "/workspace/gajanan/adh/libraries/libxml2-2.9.4/mylibxml/include/libxml2" CACHE STRING "XML Include directory")
  endif(XML_LIBRARY)

elseif($ENV{HOSTNAME} MATCHES irene.ices.utexas.edu)
  if(XML_LIBRARY)
    # Do nothing; use cache
  else(XML_LIBRARY)
    set(XML_LIBRARY "/workspace/libraries/libxml2-2.9.4/mylibxml/lib/libxml2.a" CACHE STRING "XML Library")
     # NOTE: the version of libxml2 on jade has a bug; the slightly hacked header
     #       files are located here
    set(XML_INCLUDE_DIR "/workspace/libraries/libxml2-2.9.4/mylibxml/include/libxml2" CACHE STRING "XML Include directory")
  endif(XML_LIBRARY)

elseif($ENV{HOSTNAME} MATCHES katrina)
  if(XML_LIBRARY)
    # Do nothing; use cache
  else(XML_LIBRARY)
      set(XML_LIBRARY "/workspace/sestes/ADH/adh_libs/libxml2-2.9.4/mylibxml/lib/libxml2.a" CACHE STRING "XML Library")
     # NOTE: the version of libxml2 on jade has a bug; the slightly hacked header
     #       files are located here
     set(XML_INCLUDE_DIR "/workspace/sestes/ADH/adh_libs/libxml2-2.9.4/mylibxml/include/libxml2" CACHE STRING "XML Include directory")
  endif(XML_LIBRARY)


else($ENV{HOSTNAME} MATCHES jade)
  include(${CMAKE_INCLUDE_MODULES_DIR}/FindLIBXML.cmake)

endif($ENV{HOSTNAME} MATCHES jade)
