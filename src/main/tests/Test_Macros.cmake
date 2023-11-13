#=============================================================================
# AdH - Adaptive Hydraulics
# Macros for common tests
#=============================================================================

set(mpi_serial_suffix "MPI_Serial")
set(mpi_parallel_suffix "MPI_Parallel")

##############################################################################
# Ensure AdH output target build compliance (target was built or the intent is
# that it does not exist)
# Inputs:
#   1. Output name to generate a unique test name
#   2. EXEC or LIB whether output target is an executable or library, respectively
#   3. TRUE or FALSE whether the target should exist (have been built), respectively
#   4. Name of CMake target as appearing in add_executable or add_library
#   5. Dashboard label of test (e.g. ${DASHBOARD_LABEL}_pilotpt)
macro(TestOutputCompliance _output_name _output_type _output_check _output_target _test_label)
  if(WIN32)
    set(list_cmd "dir")
  elseif(UNIX)
    set(list_cmd "ls")
  endif()

  string(TOUPPER ${_output_type} type)
  if(type STREQUAL "EXEC")
    set(work_dir ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    set(type_name "Exec")
  elseif(type STREQUAL "LIB")
    set(work_dir ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    set(type_name "Lib")
  else()
    message(FATAL_ERROR "The second parameter of TestOutputCompliance macro must be EXEC or LIB")
  endif()

  string(TOUPPER ${_output_check} check)
  if(check STREQUAL "TRUE")
    set(name ${_output_name}_${type_name}_Built)
    if(NOT TARGET ${_output_target})
      message(FATAL_ERROR "The forth parameter of TestOutputCompliance macro, the CMake target, does not exist (passing TRUE as second parameter implies that the target should exist)")
    endif()
    if(type STREQUAL "EXEC")
      add_test(NAME ${name}
        WORKING_DIRECTORY ${work_dir}
        COMMAND ${list_cmd} $<TARGET_FILE_NAME:${_output_target}>
        )
    elseif(type STREQUAL "LIB")
      add_test(NAME ${name}
        WORKING_DIRECTORY ${work_dir}
        COMMAND ${list_cmd} $<TARGET_LINKER_FILE_NAME:${_output_target}>
        )
    endif()
    set_tests_properties(${name} PROPERTIES
      FAIL_REGULAR_EXPRESSION "No such file or directory;File Not Found"
      TIMEOUT 1
      )
  elseif(check STREQUAL "FALSE")
    set(name ${_output_name}_${type_name}_Not_Included)
    add_test(NAME ${name}
      WORKING_DIRECTORY ${work_dir}
      COMMAND ${list_cmd}
      )
    # Add the NOT underscore stipulation to the regex to differentiate between pre_adh_xxx and adh_xxx targets
    set_tests_properties(${name} PROPERTIES
      WILL_FAIL ${check}
      FAIL_REGULAR_EXPRESSION "[^_]${_output_target}"
      TIMEOUT 1
      )
  else()
    message(FATAL_ERROR "The third parameter of TestOutputCompliance macro must be TRUE or FALSE")
  endif()

  #set_tests_properties(${name} PROPERTIES
  set_property(TEST ${name} PROPERTY
    LABELS "${_test_label}" basic
    )
endmacro()

##############################################################################
# Ensure that AdH reports that the physics or package is enabled
# Inputs:
#   1. Exectuable flag (PRE or MAIN to run pre-Adh or the main AdH executable,
#      respectively)
#   2. Option name to generate a unique test name
#   3. Name of variable that contains the Build Information section text for
#      this option ("Built with " is prefixed to this text within this
#      function for convenience).
#      Note: This input uses the "pass by variable name" method instead of
#      passing a literal string to permit regex syntax since functions/macros
#      consume double escapes in string arguments (see CMake bug 5389)
#   4. Name of variable that contains an expression to determine whether option
#      is enabled (expression is used in an "if" statement)
function(TestOptionEnabled _test_exec _opt_name _regex_var _flag)
  string(TOUPPER ${_test_exec} exec)
  if(exec STREQUAL "PRE")
    set(exec_tar "pre_adh")
    set(exec_label "pre_${DASHBOARD_LABEL}")
  elseif(exec STREQUAL "MAIN")
    set(exec_tar "adh")
    set(exec_label "${DASHBOARD_LABEL}")
  else()
    message(FATLA_ERROR "The first paramter of TestOptionEnabled macro must be PRE, or MAIN")
  endif()

  set(not_name "")
  set(extra_prop "")
  if(NOT ${_flag})
    set(not_name "_Not")
    set(extra_prop WILL_FAIL TRUE)
  endif()

  set(ver_name ${_opt_name}${not_name}_Enabled_Version_Build_Info)
  set(run_name ${_opt_name}${not_name}_Enabled_Executed_Build_Info)

  add_test(NAME ${ver_name}
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE_NAME:${exec_tar}> -v
    )
  add_test(NAME ${run_name}
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE_NAME:${exec_tar}> bad_input_name
    )
  set_tests_properties(${ver_name} ${run_name} PROPERTIES
    PASS_REGULAR_EXPRESSION "Built with ${${_regex_var}}"
    TIMEOUT 1
    ${extra_prop}
    )
  set_property(TEST ${ver_name} ${run_name} PROPERTY
    LABELS ${exec_label} basic
    )
endfunction()

##############################################################################
# Add test that executes AdH with given argument. This function will set
# standard test properties (labels and timeout). Additional test properties
# should be specified by a separate set_tests_properties call.
# Inputs:
#   1. Exectuable flag (PRE or MAIN to run pre-Adh or the main AdH executable,
#      respectively) 
#   2. Name of test
#   3. Extra test label (specify empty string, "", if not needed) 
#   4. Test level (BASIC, STANDARD, ADVANCED, or CUSTOM)
#   5. Working directory of test
#   6. Number of processors (0 = run serial and >0 = run parallel; used only
#      if executing AdH and built with MPI)
#   7. AdH simulation name (command line executable argument)
#   8. Timeout (how many seconds to allow for this test)
function(AddExecTest _test_exec _test_name _test_label _test_level _test_dir _num_procs _exec_arg _timeout)
  string(TOUPPER ${_test_exec} exec)
  if(exec STREQUAL "PRE")
    set(exec_tar "pre_adh")
    set(exec_label "pre_${DASHBOARD_LABEL}")
  elseif(exec STREQUAL "MAIN")
    set(exec_tar "adh")
    set(exec_label "${DASHBOARD_LABEL}")
  else()
    message(FATAL_ERROR "The first paramter of AddExecTest macro must be PRE, or MAIN")
  endif()

  string(TOLOWER ${_test_level} level) # use lowercase because this will be a added as a label
  if(NOT level MATCHES "^(basic|standard|advanced|custom)$")
    message(FATAL_ERROR "The fourth parameter of AddExecTest macro must be BASIC, STANDARD, ADVANCED, or CUSTOM")
  endif()

  if(NOT EXISTS ${_test_dir} AND NOT ${_test_dir} STREQUAL ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    message(FATAL_ERROR "The fifth parameter of AddExecText macro must be valid (directory path exists)") 
  endif()

  set(np 0)
  if(_num_procs LESS 0)
    message(FATAL_ERROR "The sixth parameter of AddExecText macro, number of processors, must be zero or greater")
  else(_num_procs GREATER 0 AND PACKAGE_MPI)
    set(np ${_num_procs})
    if(np GREATER MPIEXEC_MAX_NUMPROCS)
      set(np ${MPIEXEC_MAX_NUMPROCS})
    endif()
  endif()

  if(np EQUAL 0) # serial
    add_test(NAME ${_test_name}
      WORKING_DIRECTORY ${_test_dir}
      COMMAND $<TARGET_FILE:${exec_tar}> ${_exec_arg}
      )
  else() # parallel (launch with mpi)
    add_test(NAME ${_test_name}
      WORKING_DIRECTORY ${_test_dir}
      COMMAND ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_POSTFLAGS} $<TARGET_FILE:${exec_tar}> ${_exec_arg}
      )
  endif()

  set_tests_properties(${_test_name} PROPERTIES
    TIMEOUT ${_timeout}
    )
  set_property(TEST ${_test_name} PROPERTY
    LABELS ${exec_label} ${level} ${_test_label}
    )
endfunction()

###############################################################################
# Add test that executes AdH and checks that model successfully completes
# Inputs:
#   1. Exectuable flag (PRE or MAIN to run pre-Adh or the main AdH executable,
#      respectively) 
#   2. Test level (BASIC, STANDARD, ADVANCED, or CUSTOM)
#   3. Name of test simulation (used to find test files, as executable
#      agrument, and to generate test unique name)
#   4. Timeout (how many seconds to allow for this test)
function(AddSimRunTest _test_exec _test_level _sim_name _timeout)
  set(pre_test_name "pre-AdH_SimRun_${_sim_name}")
  string(TOUPPER ${_test_exec} exec)
  set(depends "")
  if(exec STREQUAL "PRE")
    set(test_name ${pre_test_name})
    set(rgx_end "Writing the adh file")
  elseif(exec STREQUAL "MAIN")
    set(test_name "AdH_SimRun_${_sim_name}")
    set(rgx_end "Run Time is [0-9.e+]* seconds")
    set(depends DEPENDS ${pre_test_name})
  else()
    message(FATAL_ERROR "The first paramter of AddSimRunTest macro must be PRE, or MAIN")
  endif()
  
  string(TOUPPER ${_test_level} level)
  if(NOT level MATCHES "^(BASIC|STANDARD|ADVANCED|CUSTOM)$")
    message(FATAL_ERROR "The second parameter of AddSimRunTest macro must be BASIC, STANDARD, ADVANCED, or CUSTOM")
  endif()

  set(src_dir "${ADH_TEST_DIR}/${_sim_name}")
  if(NOT EXISTS ${src_dir})
    message(FATAL_ERROR "The ${_sim_name} simulation directory does not exist in the specified ADH_TEST_DIR path (${ADH_TEST_DIR})")
  endif()

  file(GLOB input_files
    "${src_dir}/*.3dm"
    "${src_dir}/*.bc"
    "${src_dir}/*.hot"
    "${src_dir}/*.bt"
    "${src_dir}/*.sup"
    "${src_dir}/*.faces"
    "${src_dir}/*.sedlib"
    )
  list(LENGTH input_files num)
  if(num LESS 3)
    message(FATAL_ERROR "Input files are missing from simulation directory (${src_dir})")
  endif()

  set(work_dir "${src_dir}/ctest")
  file(COPY ${input_files}
    DESTINATION ${work_dir}
    )

  AddExecTest(${_test_exec} ${test_name} "" ${level} "${work_dir}" 0 "${_sim_name}" ${_timeout})
  set(extra_names)
  if(exec STREQUAL "MAIN" AND PACKAGE_MPI)
    AddExecTest(${_test_exec} ${test_name}_${mpi_serial_suffix} "" ${level} "${work_dir}" 1
      "${_sim_name};${mpi_serial_suffix}" ${_timeout})
    AddExecTest(${_test_exec} ${test_name}_${mpi_parallel_suffix} "" ${level} "${work_dir}"
      ${MPIEXEC_MAX_NUMPROCS} "${_sim_name};${mpi_parallel_suffix}" ${_timeout})
    set(extra_names ${test_name}_${mpi_serial_suffix} ${test_name}_${mpi_parallel_suffix})
  endif()
  set_tests_properties(${test_name} ${extra_names} PROPERTIES
    PASS_REGULAR_EXPRESSION "${rgx_end}"
    )
  if(exec STREQUAL "MAIN")
    set_tests_properties(${test_name} ${extra_names} PROPERTIES
      ${depends}
      )
  endif()
endfunction()

##############################################################################
# Add comparison test for result files, test on angle_nb
#   1. Run Type (SW2, SW3, NS, GW and optional T or S)
#   2. Name of test simulation (used to find test files, as executable
#      agrument, and to generate test unique name)
#   3. Timeout (how many seconds to allow for this test)

function(AddSimCompareTest _run_type _sim_name _timeout)
  # do nothing if comparison executable was not found
  if(NOT TESTING_COMPARISON_EXEC)
    return()
  endif()

  string(TOUPPER ${_run_type} type)
    if(type STREQUAL "SW2")
      set(extension "_dep" "_vel" "_error")
    elseif(type STREQUAL "SW2T")
      set(extension "_dep" "_vel" "_error" "_con1")
    elseif(type STREQUAL "SW2S")
      set(extension "_dep" "_vel" "_error" "_con1" "_dpl" "_bedload" "_alt" "_ald")
    elseif(type STREQUAL "SW3")
      set(extension "_dpl" "_dep" "_vel" "_grid_speed" "_pressure")
    elseif(type STREQUAL "SW3T")
      set(extension "_dpl" "_dep" "_vel" "_grid_speed" "_pressure" "_con1")
    elseif(type STREQUAL "GW")
      set(extension "_flx" "_hd" "_sat" "_vel")
    elseif(type STREQUAL "NS")
      set(extension "_err" "_prs" "_vel")
    elseif(type STREQUAL "NST")
      set(extension "_err" "_prs" "_vel" "_con1")
    else()
      message(FATAL_ERROR "The first paramter of AddSimCompareTest macro must be a valid run type")
    endif()


  # Loop on files (extension) for each simulation test
  foreach(ext ${extension})

    if(ext STREQUAL "_dep")
      set(_col 1)
    elseif(ext STREQUAL "_vel")
      set(_col 3)
    elseif(ext STREQUAL "_error")
      set(_col 1)
    elseif(ext STREQUAL "_con1")
      set(_col 1)
    elseif(ext STREQUAL "_dpl")
      set(_col 1)
    elseif(ext STREQUAL "_bedload")
      set(_col 3)
    elseif(ext STREQUAL "_alt")
      set(_col 1)
    elseif(ext STREQUAL "_ald")
      set(_col 1)
    elseif(ext STREQUAL "_gsp")
      set(_col 1)
    elseif(ext STREQUAL "_prs")
      set(_col 1)
    elseif(ext STREQUAL "_flx")
      set(_col 1)
    elseif(ext STREQUAL "_sat")
      set(_col 1)
    elseif(ext STREQUAL "_hd")
      set(_col 1)
    endif()

    set(test_name "AdH_SimCompare_${_sim_name}${ext}")
    set(src_dir "${ADH_TEST_DIR}/${_sim_name}")
    if(NOT EXISTS ${src_dir})
      message(FATAL_ERROR "The ${_sim_name} simulation directory does not exist in the specified ADH_TEST_DIR path (${ADH_TEST_DIR})")
    endif()
    
    set(work_dir "${src_dir}/ctest")
    set(basis_result "${src_dir}/${_sim_name}${ext}.dat.tst")
    add_test(NAME ${test_name}
      WORKING_DIRECTORY ${work_dir}
      COMMAND ${TESTING_COMPARISON_EXEC} ${basis_result} ${_sim_name}${ext}.dat ${_col}
      )
    # set specific properties
    set_tests_properties(${test_name} PROPERTIES
      DEPENDS "AdH_SimRun_${_sim_name}"
      )

    # include MPI versions
    set(extra_names)
    if(PACKAGE_MPI)
      add_test(NAME ${test_name}_${mpi_serial_suffix}
        WORKING_DIRECTORY ${work_dir}
        COMMAND ${TESTING_COMPARISON_EXEC} ${basis_result} ${_sim_name}_${mpi_serial_suffix}${ext}.dat ${_col}
        )
      # set specific properties
      set_tests_properties(${test_name}_${mpi_serial_suffix} PROPERTIES
        DEPENDS "AdH_SimRun_${_sim_name}_${mpi_serial_suffix}"
        )

      add_test(NAME ${test_name}_${mpi_parallel_suffix}
        WORKING_DIRECTORY ${work_dir}
        COMMAND ${TESTING_COMPARISON_EXEC} ${basis_result} ${_sim_name}_${mpi_parallel_suffix}${ext}.dat ${_col}
        )
      # set specific properties
      set_tests_properties(${test_name}_${mpi_parallel_suffix} PROPERTIES
        DEPENDS "AdH_SimRun_${_sim_name}_${mpi_parallel_suffix}"
        )

      add_test(NAME ${test_name}_Serial_vs_Parallel
        WORKING_DIRECTORY ${work_dir}
        COMMAND ${TESTING_COMPARISON_EXEC} ${_sim_name}${ext}.dat ${_sim_name}_${mpi_parallel_suffix}${ext}.dat ${_col}
        )
      # set specific properties
      set_tests_properties(${test_name}_Serial_vs_Parallel PROPERTIES
        DEPENDS "AdH_SimRun_${_sim_name};AdH_SimRun_${_sim_name}_${mpi_parallel_suffix}"
        )

      set(extra_names ${test_name}_${mpi_serial_suffix} ${test_name}_${mpi_parallel_suffix})
    endif()
    # set common properties
    set_tests_properties(${test_name} ${extra_names} PROPERTIES
      FAIL_REGULAR_EXPRESSION "FAIL:"
      TIMEOUT ${_timeout}
      )
    set_property(TEST ${test_name} ${extra_names} PROPERTY
      LABELS ${DASHBOARD_LABEL} advanced
      )
  endforeach()
endfunction()
##############################################################################
# Add comparison test for analytic test cases 
#   1. Name of test simulation (used to find test files, as executable
#      agrument, and to generate test unique name)
#   2. Timeout (how many seconds to allow for this test)

function(AddSimCompareAnalytic _sim_name _timeout)
    set(test_name "AdH_SimCompareAnalytic_${_sim_name}${ext}")
    set(src_dir "${ADH_TEST_DIR}/${_sim_name}")
    if(NOT EXISTS ${src_dir})
      message(FATAL_ERROR "The ${_sim_name} simulation directory does not exist in the specified ADH_TEST_DIR path (${ADH_TEST_DIR})")
    endif()
    
    set(work_dir "${src_dir}/ctest")
    set(basis_result "${src_dir}/error.out.last")
    add_test(NAME ${test_name}
      WORKING_DIRECTORY ${work_dir}
      COMMAND diff -bwt ${basis_result} "error.out" 
      )
    # set specific properties
    set_tests_properties(${test_name} PROPERTIES
      DEPENDS "AdH_SimRun_${_sim_name}"
      )

    # include MPI versions
    set(extra_names)
    if(PACKAGE_MPI)
      add_test(NAME ${test_name}_${mpi_serial_suffix}
        WORKING_DIRECTORY ${work_dir}
        COMMAND diff -bwt ${basis_result} "${_sim_name}_${mpi_serial_suffix}_error.out"
        )
      # set specific properties
      set_tests_properties(${test_name}_${mpi_serial_suffix} PROPERTIES
        DEPENDS "AdH_SimRun_${_sim_name}_${mpi_serial_suffix}"
        )

      add_test(NAME ${test_name}_${mpi_parallel_suffix}
        WORKING_DIRECTORY ${work_dir}
        COMMAND diff -bwt ${basis_result} "${_sim_name}_${mpi_parallel_suffix}_error.out" 
        )
      # set specific properties
      set_tests_properties(${test_name}_${mpi_parallel_suffix} PROPERTIES
        DEPENDS "AdH_SimRun_${_sim_name}_${mpi_parallel_suffix}"
        )

      add_test(NAME ${test_name}_Serial_vs_Parallel
        WORKING_DIRECTORY ${work_dir}
        COMMAND diff -bwt "${_sim_name}_error.out" "${_sim_name}_${mpi_parallel_suffix}_error.out"
        )
      # set specific properties
      set_tests_properties(${test_name}_Serial_vs_Parallel PROPERTIES
        DEPENDS "AdH_SimRun_${_sim_name};AdH_SimRun_${_sim_name}_${mpi_parallel_suffix}"
        )

      set(extra_names ${test_name}_${mpi_serial_suffix} ${test_name}_${mpi_parallel_suffix})
    endif()
    # set common properties
    set_tests_properties(${test_name} ${extra_names} PROPERTIES
      FAIL_REGULAR_EXPRESSION "FAIL"
      PASS_REGULAR_EXPRESSION "PASS"
      TIMEOUT ${_timeout}
      )
    set_property(TEST ${test_name} ${extra_names} PROPERTY
      LABELS ${DASHBOARD_LABEL} advanced
      )
endfunction()
