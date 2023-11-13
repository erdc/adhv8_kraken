#=============================================================================
# AdH - Adaptive Hydraulics
# friction submodule tests
#=============================================================================
set(submod_name "Friction")
set(submod_tar "adh_friction")
set(submod_label "${DASHBOARD_LABEL}_friction")

include(${CMAKE_SOURCE_DIR}/main/tests/Test_Macros.cmake)

if("${BUILD_TEST_LEVEL}" MATCHES "^(BASIC|STANDARD|ADVANCED|ALL)$")
  # Pre-AdH
  if("${BUILD_EXEC}" MATCHES "^(PRE|ALL)$")
    # Check that submodule was built for Pre-AdH
    TestOutputCompliance(${submod_name}_Pre LIB TRUE pre_${submod_tar} pre_${submod_label})
  endif("${BUILD_EXEC}" MATCHES "^(PRE|ALL)$")
  # AdH
  if("${BUILD_EXEC}" MATCHES "^(MAIN|ALL)$")
    # Check that submodule was built for AdH
    TestOutputCompliance(${submod_name} LIB TRUE ${submod_tar} ${submod_label})
  endif("${BUILD_EXEC}" MATCHES "^(MAIN|ALL)$")
endif("${BUILD_TEST_LEVEL}" MATCHES "^(BASIC|STANDARD|ADVANCED|ALL)$")
