#=============================================================================
# AdH - Adaptive Hydraulics
# Standard test suite test definitions for the AdH Executable
#=============================================================================
set(exec_name "AdH")
set(exec_tar "adh")
set(exec_label "${DASHBOARD_LABEL}")
set(test_level "STANDARD")

# Local macro to wrap and simpifly function call
macro(myAddSimRunTest _name _timeout)
  AddSimRunTest(MAIN ${test_level} "${_name}" ${_timeout}) 
endmacro()

##############################################################################
##############################################################################
# The following tests run simple AdH "toy problem" simulations to compare
# output with analytic solutions.

# Run the example simulations (these tests only runs the simulation and looks
# for successful completion; comparison of the solution is not done currently) 
myAddSimRunTest(angle_wind 60)
myAddSimRunTest(angle_db 60)
myAddSimRunTest(angle_nb 60)
myAddSimRunTest(angle_wave 60)
myAddSimRunTest(angle_wind_3station 60)
myAddSimRunTest(angle_salt2d 60)

