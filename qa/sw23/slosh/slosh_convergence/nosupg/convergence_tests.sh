#!/bin/bash
#########################################################################
sw="sloshfull3d."                                # Test case folder name
mesh_set="4 5"
#sw=""                                # Test case folder name
#mesh_set="3"
#mesh_set="sloshfull2d.4 slosh2d3d.4 sloshfull3d.4 sloshfull2d.5 slosh2d3d.5 sloshfull3d.5"
#mesh_set="sloshfull2d.6 slosh2d3d.6 sloshfull3d.6"

testcase="superfile"                                # Test case file name

adh_root_path="/workspace/gajanan/adh_2017/adh_01_31_2018"
TestcaseRepoPath="/workspace/gajanan/adh_2017/testcases/slosh_convergence/finalslosh3"

mailto="gajanan@utexas.edu"
umfpackonoff="ON"
windlibonoff="OFF"
wvel2donoff="ON"
xdmfonoff="ON"
debuglevel=0

# Additional inputs for parallel runs
parallel_run="mpirun -np"               # <MPI executable> <#cores flag>
compiler_module="gcc/8.2"    # This module will be loaded via module load
mpi_module="openmpi"         # This module will be loaded via module load
cores_set="6"         # Edit to change number of cores in parallel runs
computer="ICES"    # COMMENT THIS LINE IF THIS IS NOT A ICES, UT COMPUTER

#########################################################################

####################################################################################################################

#    DO NOT EDIT THIS SCRIPT BEYOND THIS LINE!!!!!!!!!!

####################################################################################################################







####################################################################################################################
clear

if [ ! -d $adh_root_path ]; then
   echo "AdH root directory path '$adh_root_path' not found. Please open the script and specify the FULL path to the variable 'adh_root_path' manually";
   echo "Terminating bash run"
   exit -1
fi
if [ ! -d $TestcaseRepoPath ]; then
   echo "Testcase repository '$TestcaseRepoPath' not found. Please open the script and specify the FULL path to the test case repository manually";
   echo "Terminating bash run"
   exit -1
fi

echo "";
echo "***********************************************";
for mesh in $mesh_set;do
  echo "Test case $sw$mesh chosen for simulation"
  if [ ! -d "$TestcaseRepoPath/$sw$mesh" ]; then
    echo "Could not find the folder $sw$mesh at path $TestcaseRepoPath/";
    echo "Terminating bash run"
    exit -1
  else
    echo "Found $sw$mesh directory: $TestcaseRepoPath/$sw$mesh";
  fi
done
sleep 0.5

runstatus=-1
savestatus=-1
cleanstatus=-1
tailoutstatus=-1
couplestatus=-1
maketarstatus=-1
mailstatus=-1
#resultdir="results"
cmake_args=""

####################################################################################################################
checkpoint(){
if [ $1 -eq -2 ]; then
  echo;echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "Error: Conflicting command line arguments supplied by user. Terminating script."
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  exit -1
fi
if [ $1 -eq -3 ]; then
  echo;echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "Error: Directory $2 not found! Terminating script."
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  exit -1
fi
if [ $1 -ne 0 ]; then
  echo;echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "Some command was prevented from execution or proceeded incorrectly. Terminating script."
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  exit -1
fi
}
shortecho() { echo ""; echo "***********************************************"; }
longecho() { echo ""; echo "*****************************************************************************************************"; }

####################################################################################################################
longecho
for arg in "$@" ;do
  shortecho
  echo "CHECKING ARGUMENT USER ARGUMENT: $arg";echo;
########################################################## For SERIAL RUNS
  if [ "$arg" = "s" ] || [ "$arg" = "serial" ]; then
    if [ $runstatus -eq 1 ]; then
      checkpoint -2
    fi
    runstatus=0
    numtimes="0"
    mpionoff="OFF"
    echo "Run type option detected: SERIAL runs chosen"
    sleep 0.5
  fi
########################################################## For PARALLEL RUNS
  if [ "$arg" = "p" ] || [ "$arg" = "parallel" ]; then
    if [ $runstatus -eq 0 ]; then
      checkpoint -2
    fi
    runstatus=1
    numtimes="$cores_set"
    mpionoff="ON"
    echo "Run type option detected: PARALLEL runs chosen"
    sleep 0.5
  fi
########################################################## For TAR ARCHIVE
  if [ "$arg" = "maketar" ] || [ "$arg" = "mt" ]; then
    maketarstatus=1
    echo "'Make tar option detected - tar archive of simulation will be created.'"
    #sleep 0.5
  fi
########################################################## For CLEANING BUILD
  if [ "$arg" = "justclean" ] || [ "$arg" = "jc" ]; then
    cleanstatus=2
    echo "'Just clean - don't build' option detected. All CMake files will be deleted."
    #sleep 0.5
  fi
########################################################## For CLEANING BUILD
  if [ "$arg" = "noclean" ] || [ "$arg" = "nc" ]; then
    cleanstatus=0
    echo "'Don't clean build' option detected. adh executable in /CMAKE/bin/ will be used as is"
    #sleep 0.5
  fi
########################################################## For TAIL OUTPUT
  if [ "$arg" = "t" ] || [ "$arg" = "tail" ] || [ "$arg" = "tailout" ]; then
    tailoutstatus=1
    echo "Tail output of all simulation runs will be printed"
    sleep 0.5
  fi
########################################################## For Debug
  if [ "$arg" = "debug" ] || [ "$arg" = "d" ] || [ "$arg" = "DEBUG" ] || [ "$arg" = "D" ]; then
    debuglevel=2
    echo "Build type DEBUG detected:"
    echo "    -- Debugging turned on"
    #sleep 0.5
  fi
########################################################## For MAILING JOB COMPLETION
  if [ "$arg" = "mail" ] || [ "$arg" = "m" ]; then
    mailstatus=1
    echo "Mail option detected - An email will be sent to $mailto upon job completion."
    #sleep 0.5
  fi
##########################################################

done
####################################################################################################################

sleep 0.5

########################################################## For CLEANING BUILD - DEFAULT
  if [ $cleanstatus -eq -1 ]; then
    cleanstatus=1
    shortecho
    echo "Build cleaning option was not specified by user"
    echo "     - Build directory /CMAKE will be cleaned by default"
    sleep 0.5
  fi
########################################################## For SERIAL RUN - DEFAULT
  if [ $runstatus -eq -1 ]; then
    shortecho
    echo "No run type (serial/parallel) was chosen by user"
    echo "     - Serial run chosen by default"
    runstatus=0
    numtimes="0"
    mpionoff="OFF"
    sleep 0.5
  fi
########################################################## FOR NOT CLEANING AFTER SIMULATION - DEFAULT
  if [ $savestatus -eq -1 ]; then
    savestatus=2
    shortecho
    echo "No options for saving simulation files specified by user"
    echo "     - All simulation files will be saved."
    sleep 0.5
  fi
########################################################## FOR NOT PRINTING TAIL OUTPUT - DEFAULT
  if [ $tailoutstatus -eq -1 ]; then
    tailoutstatus=0
    shortecho
    echo "Tail output option was not specified by user"
    echo "     - Tail output of simulation files will not be printed by default"
    sleep 0.5
  fi
##########################################################

cmake_args="-DUSE_WINDLIB=$windlibonoff -DUSE_MPI=$mpionoff -DUSE_PARMETIS=$mpionoff -DUSE_SUPER_LIBRARY=OFF -DUSE_UMFPACK=$umfpackonoff -DUSE_XDMF=$xdmfonoff -DBUILD_DEBUG_LEVEL=$debuglevel  -DINCLUDE_WVEL_STAGE_FOR_2D_MODELS=$wvel2donoff "

longecho
####################################################################################################################

##########################################################
if [ ! -d "$adh_root_path/CMAKE" ]; then  mkdir "$adh_root_path/CMAKE"; fi
if [ ! -d "$adh_root_path/CMAKE" ]; then  checkpoint -3 "$adh_root_path/CMAKE"; fi
cd "$adh_root_path/CMAKE/"


if [ $cleanstatus -eq 1 ]; then
  shortecho
  echo "Cleaning previous build, if any..."
  echo "     - Deleting folders"
  rm -r bin/adh CMakeFiles columns debug elem fe friction grid initio main messg meteor >& /dev/null
  rm -r lib node solver structs tools turbulence windlib wq_nsm comm mesh parmetis metis GKlib >& /dev/null
  echo "     - Deleting remaining files"
  rm Makefile CTestCustom.cmake cmake_install.cmake CMakeCache.txt >& /dev/null
  echo "     - Cleaning complete"
fi
if [ $cleanstatus -eq 2 ]; then
  # Just cleaning!
  exit 0;
fi
##########################################################
if [ $runstatus -eq 1 ]; then
  shortecho
  if [ "$HOSTNAME" = "harvey.ices.utexas.edu" ]; then
      echo "Loading modules gcc and openmpi"
      module load ubt18
      checkpoint $?
  fi
  echo "Loading modules $compiler_module and $mpi_module"
  module load $compiler_module
  checkpoint $?
  module load $mpi_module
  checkpoint $?
  module list
fi
##########################################################

####################################################################################################################
longecho
cmake $cmake_args ..
checkpoint $?
longecho
make
checkpoint $?
longecho
sleep 1
if [ ! -d "$TestcaseRepoPath" ]; then  checkpoint -3 "$TestcaseRepoPath"; fi
cd "$TestcaseRepoPath/"

####################################################################################################################

if [ $runstatus -eq 0 ]; then RESULTFLAG=""; else RESULTFLAG="_np"; fi

for mesh in $mesh_set;do
  TestcaseDir="$sw$mesh$RESULTFLAG"
  rm -r -f "$adh_root_path/CMAKE/bin/$TestcaseDir"* >& /dev/null


##########################################################
  for np in $numtimes;do
  
##########################################################
    if [ $runstatus -eq 0 ]; then NP=""; else NP="$np";fi
#    if [ $savestatus -ne 0 ]; then mkdir "$adh_root_path/CMAKE/bin/$RESULTDIR$NP/"; fi
  
##########################################################
    shortecho
    echo "Copying all input files to $adh_root_path/CMAKE/bin/"
    cp -r "$sw$mesh/" "$adh_root_path/CMAKE/bin/$TestcaseDir$NP"
##########################################################
   
    
    cd "$adh_root_path/CMAKE/bin/"
##########################################################
    echo "Copying adh executable to $adh_root_path/CMAKE/bin/$TestcaseDir$NP/"
    cp adh "$TestcaseDir$NP/"
####################################################################################################################
    
    longecho
    if [ $runstatus -eq 0 ]; then
      echo "Serial $sw$mesh $testcase runs starting now"
    else
      echo "Parallel $sw$mesh $testcase runs starting on $np processor(s) now"
    fi
##########################################################
    cd $TestcaseDir$NP
    longecho
    rm *.dat *.out output.* *.h5 *.xmf *.xmf_back *.h5_back >& /dev/null
    echo "$sw$mesh RUN STARTED..."
    starttime="`date +"%D %T"`";
    if [ $runstatus -eq 0 ];then ./adh $testcase > "output.$sw$mesh$RESULTFLAG$NP"; else $parallel_run $np ./adh $testcase > "output.$sw$mesh$RESULTFLAG$NP";fi
    #if [ $runstatus -eq 0 ];then ./adh $testcase 2>&1 | tee "output.$sw$mesh$RESULTFLAG$NP"; else $parallel_run $np ./adh $testcase 2>&1 | tee "output.$sw$mesh$RESULTFLAG$NP";fi
    if [ $? -eq 0 ]; then
      endtime="`date +"%D %T"`";
      startsec=`date -d "${starttime}" +%s`
      endsec=`date -d "${endtime}" +%s`
      elapsed=`expr ${endsec} - ${startsec}`
      rm "$TestcaseDir$NP/adh" >& /dev/null
      rm "$TestcaseDir$NP/*xmf_back" >& /dev/null
      echo "$sw$mesh RUN COMPLETED."
      if [ $tailoutstatus -eq 1 ]; then
        shortecho
        echo "TAIL SIMULATION OUTPUT $sw$mesh; FILE output.$sw$mesh$RESULTFLAG$NP:";
        echo ". . ."; echo ". . .";
        tail "output.$sw$mesh$RESULTFLAG$NP"
        longecho
        echo "error.out RESULTS FOR $sw$mesh $testcase:"
        echo ". . ."; echo ". . .";
        tail "$TestcaseDir$NP/error.out"
    
      fi
      if [ $mailstatus -eq 1 ]; then
        mail -s "Job completed without errors." $mailto <<< "Start time: $starttime,    End time: $endtime,    Elapsed: $elapsed s"
      fi
    else
      endtime="`date +"%D %T"`";
      startsec=`date -d "${starttime}" +%s`
      endsec=`date -d "${endtime}" +%s`
      elapsed=`expr ${endsec} - ${startsec}`
      echo "$sw$mesh RUN EXITED WITH ERRORS. FOLLOWING IS THE TAIL OUTPUT:";
      echo ". . ."; echo ". . .";
      tail "output.$sw$mesh$RESULTFLAG$NP"
      if [ $mailstatus -eq 1 ]; then
        mail -s "Job terminated with errors." $mailto <<< "Start time: $starttime,    End time: $endtime,    Elapsed: $elapsed s"
      fi
    fi
    cd ..
    if [ $maketarstatus -eq 1 ]; then
      rm "$TestcaseDir$NP/adh"  >& /dev/null
      rm "$TestcaseDir$NP/*.xmf_back"  >& /dev/null
      rm "$TestcaseDir$NP/*.h5_back"  >& /dev/null
      
      echo "Creating a Tar archive of the results"; echo ". . .";
      tar -cvf "$TestcaseDir$NP.tar"  "$TestcaseDir$NP/"
    fi

done
 
####################################################################################################################
  
##########################################################
#  longecho
#  if [ $savestatus -ne 0 ]; then
#    cp "$TestcaseDir$NP/output.$testcase""_$sw$mesh$RESULTFLAG$NP" "$RESULTDIR$NP/"
#    cp "$TestcaseDir$NP/error.out" "$RESULTDIR$NP/error.out.$testcase""_$sw$mesh$RESULTFLAG$NP"
#  fi
#  if [ $savestatus -ne 2 ]; then rm -r -f "$TestcaseDir$NP";echo "All simulation folders deleted"; fi
##########################################################
 
  cd "$TestcaseRepoPath/" 
done

####################################################################################################################

longecho
echo ""
for mesh in $mesh_set;do
  echo "$sw$mesh $testcase TEST CASE COMPLETED. PLEASE DOUBLE CHECK FOR ERRORS MANUALLY."
done
if [ $savestatus -ne 0 ]; then
  echo "Please check $adh_root_path/CMAKE/bin/ for saved simulation files"
fi
longecho
