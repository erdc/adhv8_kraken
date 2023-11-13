#=============================================================================
# AdH - Adaptive Hydraulics
# Basic test definitions for the AdH CSTORM library
#=============================================================================
set(exec_name "AdH-CSTORM")
set(exec_tar "adh")
set(exec_label "${DASHBOARD_LABEL}")

include(CheckLibraryExists)

##############################################################################
# Was AdH CSTORM library built?

TestOutputCompliance(${exec_name} LIB TRUE ${exec_tar} ${exec_label})
