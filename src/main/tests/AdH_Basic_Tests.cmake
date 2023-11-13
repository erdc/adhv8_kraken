#=============================================================================
# AdH - Adaptive Hydraulics
# Basic test suite test definitions for the AdH Executable
#=============================================================================
set(exec_name "AdH")

# Local macro to wrap and simpifly function call
macro(myAddExecTest _name _exec_arg _timeout)
  AddExecTest(MAIN "${_name}" "" BASIC "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}" 0
    "${_exec_arg}" ${_timeout}) 
endmacro()
macro(myTestOptionEnabled _name _rgx _opt)
  TestOptionEnabled(MAIN "${exec_name}_${_name}" ${_rgx} ${_opt})
endmacro()


##############################################################################
##############################################################################
# The following tests ensure AdH runs and provides screen output (stdout).
# They do not necessarily test the results.

##############################################################################
# Was AdH built?

TestOutputCompliance(${exec_name} EXEC TRUE adh ${DASHBOARD_LABEL})

##############################################################################
# Does AdH warn about missing command line agrument?

set(test_name ${exec_name}_Missing_Command_Line_Arg)
set(test_timeout 1)

#-----------------------------------------------------------------------------
# 1. Check for uasge information
set(test_num 1)
# Create regex of expected output (done in parts for readability)
set(rgx       "Missing agrument.${rgx_nl}")
set(rgx "${rgx}Usage:${rgx_nl}")
set(rgx "${rgx}To run a simulation:${rgx_nl}")
set(rgx "${rgx}pre_adh file_base${rgx_nl}")
set(rgx "${rgx}Or,${rgx_nl}")
set(rgx "${rgx}For version information:${rgx_nl}")
set(rgx "${rgx}pre_adh -v")

myAddExecTest(${test_name}${test_num} "" ${test_timeout})
set_tests_properties(${test_name}${test_num} PROPERTIES
  PASS_REGULAR_EXPRESSION ${rgx}
  )

#-----------------------------------------------------------------------------
# 2. Check for clean exit
set(test_num 2)
myAddExecTest(${test_name}${test_num} "" ${test_timeout})

##############################################################################
# Does AdH provide version information?

set(test_name ${exec_name}_Version_Info)
set(test_timeout 1)

#-----------------------------------------------------------------------------
# 1. Check for version information
set(test_num 1)
# Create regex of expected output (done in parts for readability)
set(rgx_bi          "---*${rgx_nl}")
set(rgx_bi "${rgx_bi}${exec_name} Version${rgx_nl}")
set(rgx_bi "${rgx_bi}---*${rgx_nl}")
set(rgx_bi "${rgx_bi}Version [0-9]*\\.[0-9]*${rgx_nl}")
if(PACKAGE_SEDLIB)
  set(rgx_bi "${rgx_bi}Compatible with Sedlib v[0-9]*\\.[0-9]*${rgx_nl}")
endif()
set(rgx_bi "${rgx_bi}${rgx_nl}")
set(rgx_bi "${rgx_bi}---*${rgx_nl}")
set(rgx_bi "${rgx_bi}${exec_name} Build Information${rgx_nl}")
set(rgx_bi "${rgx_bi}---*${rgx_nl}")
myAddExecTest(${test_name}${test_num} "-v" ${test_timeout})
set_tests_properties(${test_name}${test_num} PROPERTIES
  PASS_REGULAR_EXPRESSION ${rgx_bi}
  )

#-----------------------------------------------------------------------------
# 2. Check for clean exit
set(test_num 2)
myAddExecTest(${test_name}${test_num} "-v" ${test_timeout})

##############################################################################
# Does AdH list specific build options?
set(opt_on ON)
set(opt_off OFF)

set(rgx "DEBUG enabled")
if(BUILD_DEBUG_LEVEL GREATER 1)
  myTestOptionEnabled("Debug" rgx opt_on)
else()
  myTestOptionEnabled("Debug" rgx opt_off)
endif()

set(rgx "DEBUG_KITCHEN_SINK enabled")
if(BUILD_DEBUG_LEVEL GREATER 3)
  myTestOptionEnabled("Debug_Kitchen_Sink" rgx opt_on)
else()
  myTestOptionEnabled("Debug_Kitchen_Sink" rgx opt_off)
endif()

set(rgx "DEBUG_TEMPORARY enabled")
if(BUILD_DEBUG_LEVEL EQUAL 3 OR BUILD_DEBUG_LEVEL EQUAL 5)
  myTestOptionEnabled("Debug_Temporary" rgx opt_on)
else()
  myTestOptionEnabled("Debug_Temporary" rgx opt_off)
endif()

set(rgx "GW physics enabled")
myTestOptionEnabled("GW" rgx PHYSICS_GROUNDWATER)

set(rgx "HEAT physics enabled")
myTestOptionEnabled("HT" rgx PHYSICS_HEAT)

set(rgx "NS physics enabled")
myTestOptionEnabled("NS" rgx PHYSICS_NAVIER_STOKES)

set(rgx "OVERLAND physics enabled")
myTestOptionEnabled("OL" rgx PHYSICS_OVERLAND)

set(rgx "REACT physics enabled")
myTestOptionEnabled("RCT" rgx opt_off)

set(rgx "SW2 physics enabled")
myTestOptionEnabled("SW2" rgx PHYSICS_2D_SHALLOW_WATER)

set(rgx "SW3 physics enabled")
myTestOptionEnabled("SW3" rgx PHYSICS_3D_SHALLOW_WATER)

set(rgx "XMS output file format")
myTestOptionEnabled("XMS" rgx opt_on)

set(rgx "XDMF \\(\\*\\.xml with \\*\\.h5\\) output file format")
if(BUILD_OUTPUT_FORMAT MATCHES "XDMF")
  myTestOptionEnabled("XDMF" rgx opt_on)
else()
  myTestOptionEnabled("XDMF" rgx opt_off)
endif()

set(rgx "HEC DSS file support enabled")
myTestOptionEnabled("HEC" rgx MODULE_HEC)

set(rgx "SOCKETS enabled")
myTestOptionEnabled("Sockets" rgx MODULE_SOCKETS)

set(rgx "STRUCTURES enabled")
myTestOptionEnabled("Structures" rgx MODULE_STRUCTURES)

set(rgx "UNREFINEMENT CONSERVATION enabled")
myTestOptionEnabled("Unref_Consv" rgx MODULE_UNREFINEMENT_CONSERVATION)

set(rgx "ICM library support enabled")
myTestOptionEnabled("ICM" rgx PACKAGE_ICM)

set(rgx "NSM library support enabled")
myTestOptionEnabled("NSM" rgx PACKAGE_NSM)

set(rgx "Sedlib support enabled")
myTestOptionEnabled("Sedlib" rgx PACKAGE_SEDLIB)

set(rgx "MPI enabled")
myTestOptionEnabled("MPI" rgx PACKAGE_MPI)

set(rgx "PARMETIS enabled")
myTestOptionEnabled("PARMETIS" rgx PACKAGE_PARMETIS)

set(rgx "UMFPACK (enabled|Version [345] \\(int size: (32|64) bytes)")
myTestOptionEnabled("UMFPACK" rgx PACKAGE_UMFPACK)

##############################################################################
# Does AdH provide information as screen output?

set(test_name ${exec_name}_Screen_Output)
# Initially call AdH with bad input argument
set(exec_arg bad_input)
set(test_timeout 2)

#-----------------------------------------------------------------------------
# 1. Check for build information section
set(test_num 1)
myAddExecTest(${test_name}${test_num} ${exec_arg} ${test_timeout})
# Reuse regex from above
set_tests_properties(${test_name}${test_num} PROPERTIES
  PASS_REGULAR_EXPRESSION ${rgx_bi}
  )

#-----------------------------------------------------------------------------
# 2. Check for runtime information section (default options - no super file)
set(test_num 2)
# Create regex of expected output (done in parts for readability)
set(rgx       "-------------------${rgx_nl}")
set(rgx "${rgx}Runtime Information${rgx_nl}")
set(rgx "${rgx}-------------------${rgx_nl}")
set(rgx "${rgx}${exec_name} execution Date/Time: 20[0-9][0-9]\\.[01][0-9]\\.[0-3][0-9] / [0-2][0-9]:[0-5][0-9]:[0-5][0-9]${rgx_nl}")
set(rgx "${rgx}Launching ${exec_name} with project name: ${exec_arg}${rgx_nl}")
set(rgx "${rgx}Launching ${exec_name} with run name: \\(unspecified\\)${rgx_nl}")

myAddExecTest(${test_name}${test_num} ${exec_arg} ${test_timeout})
set_tests_properties(${test_name}${test_num} PROPERTIES
  PASS_REGULAR_EXPRESSION ${rgx}
  )

#-----------------------------------------------------------------------------
# 3. Check for number of processors
#set(test_num 3)
#set(rgx "${exec_name} was launched with [0-9]* processors?")

#myAddExecTest(${test_name}${test_num} ${exec_arg} ${test_timeout})
#set_tests_properties(${test_name}${test_num} PROPERTIES
#  PASS_REGULAR_EXPRESSION ${rgx}
#  )

#-----------------------------------------------------------------------------
# 4. Check for failure from bad input
set(test_num 3)
set(rgx "ERROR: Could not find the file: ${exec_arg}")

myAddExecTest(${test_name}${test_num} ${exec_arg} ${test_timeout})
set_tests_properties(${test_name}${test_num} PROPERTIES
  PASS_REGULAR_EXPRESSION ${rgx}
  )

#-----------------------------------------------------------------------------
# 5. Check for exit error from bad input
set(test_num 4)

myAddExecTest(${test_name}${test_num} ${exec_arg} ${test_timeout})
