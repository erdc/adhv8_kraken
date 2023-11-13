#=============================================================================
# AdH - Adaptive Hydraulics
# Advanced test suite test definitions for the Pre-AdH Executable
#=============================================================================
set(exec_name "pre-AdH")
set(exec_tar "pre_adh")
set(exec_label "pre_${DASHBOARD_LABEL}")
set(test_level "ADVANCED")

# Local macro to wrap and simpifly function call
macro(myAddSimRunTest _name _timeout)
  AddSimRunTest(PRE ${test_level} "${_name}" ${_timeout}) 
endmacro()

##############################################################################
##############################################################################
# The following tests run simple Pre-AdH "toy problem" simulation
# initializations to verify proper adh.out output.
myAddSimRunTest(SD-tide 5)
myAddSimRunTest(SD-adpt 5)
myAddSimRunTest(SD-con 5)
myAddSimRunTest(SD-salt 5)
myAddSimRunTest(backstep 5)
myAddSimRunTest(bump 5)
myAddSimRunTest(2D-convergent_wall 5)
myAddSimRunTest(2D-coriolis_X 5)
myAddSimRunTest(2D-coriolis_Y 5)
myAddSimRunTest(3D-coriolis_X 5)
myAddSimRunTest(3D-coriolis_Y 5)
myAddSimRunTest(2D-dambreak 5)
myAddSimRunTest(3D-UD_wind 5)
myAddSimRunTest(3D-mass_conservation 5)
myAddSimRunTest(3D-tide_prop 5)
myAddSimRunTest(3D-turbulence 5)
myAddSimRunTest(3D-lock_exchange 5)