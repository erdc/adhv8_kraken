#=============================================================================
# AdH - Adaptive Hydraulics
# file_io submodule - AdH file list
#=============================================================================

# This file is included by the submodule's CMakeList.txt to get the list of
# source code files required by AdH (main executable)

##############################################################################
# Define the source files needed to build submodule

# Add new files to either the General list or to *each* applicable Specific
# build option list (add file to both SW2 and SW3 sections, if applicable to
# both). The final file list is processed to remove duplcate entrees, so there
# is no need for fancy if statements (no ANDs and ORs) in the Specific section.

# General source files (always included)
list(APPEND main_exec_srcs
  # Insert new files below here, alphabetically
  file_io.c
  # Insert new files above here, alphabetically
  )

# Optional source files (included based on specified build options)

# We currently do not need any specific cases for this submodule. The case
# template can be copied form another submodule when necessary. 
