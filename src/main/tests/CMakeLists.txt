#=============================================================================
# AdH - Adaptive Hydraulics
# General executable test definitions 
#=============================================================================
include(${CMAKE_CURRENT_LIST_DIR}/Test_Macros.cmake)

# regex for newline, includes pre and post whitespace
set(rgx_nl "[ \r]*\r?\n[ \r]*")

##############################################################################
# Check for external utilities used by the tests
set(comp_exec_name AdHDiff)
find_program(TESTING_COMPARISON_EXEC
  NAMES ${comp_exec_name}
  PATHS ${ADH_TEST_DIR}
  DOC "${comp_exec_name} executable for comparing AdH output results."
  NO_DEFAULT_PATH
  )
if(NOT TESTING_COMPARISON_EXEC)
  message(WARNING "${comp_exec_name} executable was not found (TESTING_COMPARISON_EXEC); comparison tests will not be generated.")
endif()

##############################################################################
# Include executable specific test definition files
# Basic test suite (included with all test levels except custom)  
if("${BUILD_TEST_LEVEL}" MATCHES "^(BASIC|STANDARD|ADVANCED|ALL)$")
  include(${CMAKE_CURRENT_LIST_DIR}/AdH_Basic_Tests.cmake)
endif("${BUILD_TEST_LEVEL}" MATCHES "^(BASIC|STANDARD|ADVANCED|ALL)$")

# Standard test suite (included with the "advanced" and "all" test suite)
if("${BUILD_TEST_LEVEL}" MATCHES "^(STANDARD|ADVANCED|ALL)$")
  include(${CMAKE_CURRENT_LIST_DIR}/AdH_Std_Tests.cmake)
endif("${BUILD_TEST_LEVEL}" MATCHES "^(STANDARD|ADVANCED|ALL)$")

# Advanced test suite (included with the "all" test suite)
if("${BUILD_TEST_LEVEL}" MATCHES "^(ADVANCED|ALL)$")
  include(${CMAKE_CURRENT_LIST_DIR}/AdH_Adv_Tests.cmake)
endif("${BUILD_TEST_LEVEL}" MATCHES "^(ADVANCED|ALL)$")

# Custom test suite (included with the "all" test suite)
# The Custom test suite definition files do not necessarily exist (they should not be checked into
# the SVN repository) so only include the files if found.
if("${BUILD_TEST_LEVEL}" MATCHES "^(CUSTOM|ALL)$")
  if(EXISTS ${CMAKE_CURRENT_LIST_DIR}/AdH_Custom_Tests.cmake)
    include(${CMAKE_CURRENT_LIST_DIR}/AdH_Custom_Tests.cmake)
  endif(EXISTS ${CMAKE_CURRENT_LIST_DIR}/AdH_Custom_Tests.cmake)
endif("${BUILD_TEST_LEVEL}" MATCHES "^(CUSTOM|ALL)$")
