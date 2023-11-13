#=============================================================================
# AdH - Adaptive Hydraulics
# file_io submodule
#=============================================================================
include(${CMAKE_INCLUDE_MODULES_DIR}/Build_Macros.cmake)

# Define include directories
include_directories("${CMAKE_SOURCE_DIR}/include")

# This submodule is not included in Pre-AdH (and should never be).

##############################################################################
# Build AdH library

# Define additional include directories, compiler flags, and dependencies
set(compiler_flags "")
set(lib_depends)

# Create library target
Add_AdH_Library("file_io" MAIN compiler_flags lib_depends)

##############################################################################
# Add submodule tests
Add_AdH_Submodule_Tests()
