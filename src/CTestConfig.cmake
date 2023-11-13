set(CTEST_PROJECT_NAME "AdH")
set(CTEST_NIGHTLY_START_TIME "20:00:00 EST")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "134.164.48.118")
set(CTEST_DROP_LOCATION "/cdash/submit.php?project=AdH")
set(CTEST_DROP_SITE_CDASH TRUE)

# AdH_Common.cmake script uses CTEST_PROJECT_SUBPROJECTS
# List subprojects in order of least to most dependent
set(CTEST_PROJECT_SUBPROJECTS
  # Pre-AdH subprojects
  # independent
  "pre_${DASHBOARD_LABEL}_externals"
  # recursive dependency
  "pre_${DASHBOARD_LABEL}_elem"
  "pre_${DASHBOARD_LABEL}_node"
  "pre_${DASHBOARD_LABEL}_initio"
  "pre_${DASHBOARD_LABEL}_solv"
  "pre_${DASHBOARD_LABEL}_tools"
  # optional library; recursive dependency
  "pre_${DASHBOARD_LABEL}_columns"
  "pre_${DASHBOARD_LABEL}_friction"
  "pre_${DASHBOARD_LABEL}_winds"
  # executable (dependent on the libraries above)
  "pre_${DASHBOARD_LABEL}"

  # AdH (main executable) subprojects
  # independent
  "${DASHBOARD_LABEL}_boat"
  "${DASHBOARD_LABEL}_externals"
  "${DASHBOARD_LABEL}_messg"
  "${DASHBOARD_LABEL}_turbulence"
  # recursive dependency
  "${DASHBOARD_LABEL}_grid"
  "${DASHBOARD_LABEL}_elem"
  "${DASHBOARD_LABEL}_node"
  "${DASHBOARD_LABEL}_comm"
  "${DASHBOARD_LABEL}_fe"
  "${DASHBOARD_LABEL}_initio"
  "${DASHBOARD_LABEL}_solv"
  "${DASHBOARD_LABEL}_tools"
  # optional library; recursive dependency
  "${DASHBOARD_LABEL}_columns"
  "${DASHBOARD_LABEL}_file_io"
  "${DASHBOARD_LABEL}_friction"
  "${DASHBOARD_LABEL}_sed"
  "${DASHBOARD_LABEL}_hydro_structures"
  "${DASHBOARD_LABEL}_winds"
  # executable (dependent on the libraries above)
  "${DASHBOARD_LABEL}"
  )
