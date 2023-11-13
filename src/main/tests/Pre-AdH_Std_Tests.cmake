#=============================================================================
# AdH - Adaptive Hydraulics
# Standard test suite test definitions for the Pre-AdH Executable
#=============================================================================
set(exec_name "pre-AdH")
set(exec_tar "pre_adh")
set(exec_label "pre_${DASHBOARD_LABEL}")
set(test_level "STANDARD")

# Local macro to wrap and simpifly function call
macro(myAddSimRunTest _name _timeout)
  AddSimRunTest(PRE ${test_level} "${_name}" ${_timeout}) 
endmacro()

##############################################################################
##############################################################################
# The following tests run real Pre-AdH simulation initializations to verify
# proper adh.out output.

# Run the example simulations (these tests only runs the simulation and looks
# for successful completion; comparison of the solution is not done currently) 
myAddSimRunTest(angle_db 5)
myAddSimRunTest(angle_nb 5)
myAddSimRunTest(angle_bt 5)
myAddSimRunTest(angle_sup 5)
myAddSimRunTest(angle_wind 5)
myAddSimRunTest(riprap2d-steady 5)
myAddSimRunTest(riprap2d-vor 5)
myAddSimRunTest(vessel 5)
myAddSimRunTest(dike-adapt 5)
myAddSimRunTest(supercritical 5)
myAddSimRunTest(pool8 5)
myAddSimRunTest(basin 5)
myAddSimRunTest(theis 5)
#myAddSimRunTest(first_flume 5)
myAddSimRunTest(het_column_8 5)
#myAddSimRunTest(riprap3d 5)

if(PACKAGE_SEDIMENT)
  myAddSimRunTest(angle_sed 5)
#  myAddSimRunTest(ns_cube_snd 5)
endif()
