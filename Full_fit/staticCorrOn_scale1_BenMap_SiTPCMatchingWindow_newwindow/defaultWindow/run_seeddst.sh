#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n new  # setup sPHENIX environment in the singularity container shell. Note the shell is bash by default

# Additional commands for my local environment
export SPHENIX=/sphenix/u/xyu3
export MYINSTALL=$SPHENIX/install

# Setup MYINSTALL to local directory and run sPHENIX setup local script
# to adjust PATH, LD LIBRARY PATH, ROOT INCLUDE PATH, etc
#source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL

echo "sPHENIX environment setup finished"

this_script=$BASH_SOURCE
this_script=`readlink -f $this_script`
this_dir=`dirname $this_script`
echo running: $this_script $*

nEvents=$1
InClusterDst=$2
InClusterPath=$3
InSeedDst=$4
InSeedPath=$5
OutDir=$6
OutPrefix=$7
Index=$8
StepSize=$9

root.exe -q -b Fun4All_TrackAnalysis.C\($nEvents,\"${InClusterDst}\",\"${InClusterPath}\",\"${InSeedDst}\",\"${InSeedPath}\",\"${OutDir}\",\"${OutPrefix}\",$Index,$StepSize\)
echo Script done
