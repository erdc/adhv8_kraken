#=============================================================================
# AdH - Adaptive Hydraulics
# Advanced test suite test definitions for the AdH Executable
#=============================================================================
set(exec_name "AdH")
set(exec_tar "adh")
set(exec_label "${DASHBOARD_LABEL}")
set(test_level "ADVANCED")

# Local macro to wrap and simpifly function call
macro(myAddSimRunTest _name _timeout)
  AddSimRunTest(MAIN ${test_level} "${_name}" ${_timeout}) 
endmacro()

##############################################################################
##############################################################################
# The following tests run real AdH simulations to ensure complete simulation.
myAddSimRunTest(wet-dry-mc 1500)
myAddSimRunTest(2D-coriolis_X 600)
myAddSimRunTest(2D-coriolis_Y 600)
myAddSimRunTest(3D-coriolis_X 6000)
myAddSimRunTest(3D-coriolis_Y 6000)
myAddSimRunTest(angle_nb_3d 6000)
myAddSimRunTest(angle_nb_3d_new 6000)
#myAddSimRunTest(angle_wave_3d 600)
myAddSimRunTest(angle_wind_3d 6000)
myAddSimRunTest(angle_salt3d 6000)
myAddSimRunTest(slosh 6000)
myAddSimRunTest(pressure_test 6000)

##############################################################################
##############################################################################
# The following tests run real AdH simulations to compare output with accepted
# solutions.
# Run Result Comparison Tests Based on Previously Saved "Acceptable" Results
# The run type specifies hydro (SW2, SW3, GW, NS) as well as transport or sediment (T, S)
AddSimCompareTest(SW2 angle_db 100)
AddSimCompareTest(SW2 angle_nb 100)
AddSimCompareTest(SW2 angle_wind 100)
AddSimCompareTest(SW2 angle_wave 100)
AddSimCompareTest(SW2 angle_wind_3station 100)
AddSimCompareTest(SW2T angle_salt2d 100)
AddSimCompareTest(SW2 wet-dry-mc 100)
AddSimCompareTest(SW2 2D-coriolis_X 100)
AddSimCompareTest(SW2 2D-coriolis_Y 100)
AddSimCompareTest(SW3 angle_nb_3d 600)
AddSimCompareTest(SW3 angle_nb_3d_new 600)
#AddSimCompareTest(SW3 angle_wave_3d 600)
#AddSimCompareTest(SW3 angle_wind_3d 600)
AddSimCompareTest(SW3T angle_salt3d 600)
AddSimCompareTest(SW3 slosh 600)
AddSimCompareTest(SW3 pressure_test 600)
AddSimCompareTest(SW3 3D-coriolis_X 600)
AddSimCompareTest(SW3 3D-coriolis_Y 600)
AddSimCompareAnalytic(angle_salt2d 600)
AddSimCompareAnalytic(angle_nb_3d_new 600)
AddSimCompareAnalytic(angle_wind_3d 600)
AddSimCompareAnalytic(angle_salt3d 600)
AddSimCompareAnalytic(slosh 600)
AddSimCompareAnalytic(pressure_test 600)
AddSimCompareAnalytic(3D-coriolis_X 600)
AddSimCompareAnalytic(3D-coriolis_Y 600)

