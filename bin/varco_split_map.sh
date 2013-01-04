#! /bin/bash

#
# VARCO_SPLIT_MAP
#

# Copyright 2012 Joseph Tran <Joseph.Tran@versailles.inra.fr>

# This software is a computer program whose purpose is to:
# - perform reads mapping in batch mode by splitting given samples in several batches
#   and run sequentially each batch to avoid cpu overload and running out of disk space.

# This software is governed by the CeCILL license, Version 2.0 (the "License"), under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license, Version 2.0 (the "License"), as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt". 

# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 

# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 

# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license, Version 2.0 (the "License"), and that you accept its terms.

# Date: 2012-12-18
VERSION=v0.1.1

########################
# SECTION CONFIGURATION
#######################

# Include lib functions

PROD_PREFIX="/usr/local"
DEV_PREFIX="$(pwd)/.."
PREFIX=$PROD_PREFIX # TO BE CHANGED WHEN SWITCHING TO PROD
. $PREFIX/share/varco_split_map/lib/varco_split_map_lib.inc

# Set variables

ARGS=3
DATA_ROOT_DIR=$1
JOB_TAG=$2
RECIPIENT=$3

EXECUTED_COMMAND="$0 $*"

DATE=$(date '+%F_%Hh%Mm%Ss')
WORKING_DIR=$(pwd)
LOGFILE=${JOB_TAG}_${USER}_$DATE.log

NAMESPACE="VARCO"
VARCO_SPLIT_MAP_SHARED=$PREFIX/share/$(basename ${0%.*})
PROD_VARCO_SPLIT_MAP_USER_CONFIG=$WORKING_DIR/$(basename ${0%.*})_user.config
DEV_VARCO_SPLIT_MAP_USER_CONFIG=$VARCO_SPLIT_MAP_SHARED/etc/$(basename ${0%.*})_user.config
VARCO_SPLIT_MAP_USER_CONFIG=$PROD_VARCO_SPLIT_MAP_USER_CONFIG # TO BE CHANGED WHEN SWITCHING TO PROD
VARCO_SPLIT_MAP_USER_CONFIG_JOB=$WORKING_DIR/$JOB_TAG/${JOB_TAG}_$(basename ${0%.*})_user.config
SAMPLE_CONFIG=sample.config

MAX_NUMB_CORES=$(cat /proc/cpuinfo | grep processor | wc -l)
MAX_NUMB_CORES_ALLOWED=$[$MAX_NUMB_CORES/2]
MAX_BATCH_SIZE=0

CPU_CHECK_TIMEOUT=86400
CPU_CHECK_INTERVAL=5
CPU_CHECK_DELAY=5
CPU_CHECK_TIMEOUT_HR=$(echo $CPU_CHECK_TIMEOUT | gawk '{printf("%dd:%02dh:%02dm:%02ds",($1/60/60/24),($1/60/60%24),($1/60%60),($1%60))}')

CORES_REDUC_FACTOR=2
CORES_SYS_AMOUNT=2

DATA_EXPANSION_FACTOR=2

PIDS_ARR=()
WAITALL_TIMEOUT=86400
WAITALL_INTERVAL=60
WAITALL_DELAY=60
WAITALL_TIMEOUT_HR=$(echo $WAITALL_TIMEOUT | gawk '{printf("%dd:%02dh:%02dm:%02ds",($1/60/60/24),($1/60/60%24),($1/60%60),($1%60))}')

LOG_DIR=$JOB_TAG/"log"
QC_TRIM_DIR="QC_TRIM"
MAPPING_DIR="MAPPING"

ERROR_TMP="/tmp/$(basename ${0%.*})_error_${USER}_$DATE.log"

PREREQUISITES_MSG="You need to copy the user configuration file to your working directory
                and configure the options before running this script.
                To copy the user configuration template file, issue this command in your terminal:
                \$ cp ${DEV_VARCO_SPLIT_MAP_USER_CONFIG} .
                Then use any editor (emacs, gedit, nano, vi, etc.) to open it and configure the options.
                Save your changes, and you are ready to run the script."

MAPPER_VERSION=$(gsnap --version 2>&1 | awk -F"\n" 'BEGIN{info=""}; {info=(info "\n\t" $1)}; END{printf info}')

AUTHOR_INFOS="Joseph Tran\nIJPB Bioinformatics Development Team\nContact: Joseph.Tran@versailles.inra.fr"

# local functions
exit_on_error()
{
	err=$1
	msg=$2
	status=$3
	log=$4
	if [[ "$status" -ne 0 ]]; then
		echo -e "$(date '+%Y_%m_%d %T') $msg" | tee -a $err 2>&1 | tee -a $log 2>&1
		echo -e "$(date '+%Y_%m_%d %T') [Pipeline error] Exits the pipeline, with error code $status." | tee -a $err 2>&1 | tee -a $log 2>&1
		echo -e "$(date '+%Y_%m_%d %T') [Pipeline error] Tail error output:" | tee -a $log 2>&1
		echo -e "... $(tail -n 10 $err)" 2>&1 | tee -a $log 2>&1
    	echo -e "$(date '+%Y_%m_%d %T') [Pipeline error] More details can be found in $err." | tee -a $log 2>&1
		# send an email
		job_error_msg="$JOB_TAG job error: $msg\nExit status: $status\nMore details can be found in $(readlink -f $err)\nLast error output:\n... $(tail -n 10 $err)"
		[[ -n $RECIPIENT ]] && sendEmail $RECIPIENT "[$(basename ${0%.*})] $JOB_TAG job error" "$job_error_msg" 
		exit $status
	fi
}

#==============================================
# TEST if enough args else print usage message
#==============================================
[[ $# -le 1 || $# -gt "$ARGS" ]] && { printf %s "\
Program: $(basename $0)
Version: $VERSION

Copyright 2012 Joseph Tran <Joseph.Tran@versailles.inra.fr>

Licensed under the CeCILL License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

Usage: $(basename $0) samples_root_dir job_tag [recipient_email_addr]

Arguments: samples_root_dir        Path to the parent folder containing the reads samples subdirectories 
           job_tag                 <String> Prefix to attach to any output files (without space)
           [recipient_email_addr]  Valid email address to send log and error messages to (optional)

Description: This script performs reads mapping in batch mode by splitting
             given samples in several batches and run sequentially each batch 
             to avoid cpu overload and running out of disk space

Pre-requisites: ${PREREQUISITES_MSG}

User Configuration File: 
  Remarks:
  - If one line begins with a "#", it will be considered as a comment. Use comment to disable the following options.
  - Note that, when option value is boolean (TRUE|FALSE), the option can be disabled by setting its value to FALSE.

  Here is the main configuration sections and their corresponding parameters:
  [split_map] section
    - check_cpu_overload (default=TRUE)
                         This option controls if cpu overload checking has to be performed.
                         If TRUE, this script will run only if cpu average load in the last 1, 5 and 15 minutes   
                         left sufficient cores to run the job. This test is performed before running the script,
                         and before each batch. If available cores (#max_cores_limit-#cpu_load) is lower or equal
                         than maximum allowed cores (needed to perform the job), the script will wait for 5 seconds
                         before testing again cpu average load.
    - batch_size (default=4)
                         This option controls the number of samples for each batch. 
                         Its value is determined dynamically if the user value exceeds the max batch size. 
                         The max batch size equals the ratio max number of cores allowed 
                         divided by the number of threads (mapper option).
                         The max number of cores allowed is fairly set by default 
                         to 50% of max number of cores available.
                         If batch_size lower or equal to max batch size, then use batch_size value.
                         Else, batch_size equals half max batch size.
    - clean (default=TRUE)
                         This option allows to clean each sample output directory by removing sam file,
                         leaving only unsorted, sorted and indexed bam files.

  [data] section
    - include_sample_subdirs (default=^.* equivalent to all sample subdirs)
                         This option controls which sample directories to consider for the mapping.
                         Another control is performed on the list of sample subdirectories to consider
                         only sample subdirs having 2 fastq files.
    - exclude_sample_subdirs (default=NULL equivalent to none sample subdirs)
                         This option controls which sample directories to exclude from the mapping.
      
    Both options accepts extended regular expression separated by comma.
                         
  [gsnap] section
    Please refer to gsnap documentation: gsnap --help
    - k (default=15)
    - d (default=Brassica_napus_v3.0.scaffold.fa)
    - m (default=1)
    - s (default=/data/temp_projects/AZM/INDEX/Brassica_napus_v3.0.scaffold.fa/Brassica_napus_CDS_v3.0.splicesites)
    - nofails (default=TRUE)
    - A (default=sam)
    - N (default=1)
    - B (default=4)
    - t (default=2)
    - D (default=/data/temp_projects/AZM/INDEX)
	
	The v option is commented by default. Its value is currently set to ALY2H_1_Bnm1.uniq.filt.SNPv2.
	Enable the option by uncommenting the line.
	Feel free to add gsnap options, using the given syntax, 'option_name'='value'.
	When option_name do not expect any value, set the option_name value with boolean to TRUE (enable), FALSE (disable).

  [samtools_view]
    Please refer to samtools vew documentation: samtools view --help
    - b (default=TRUE)
    - S (default=TRUE)

	Feel free to add samtools view options, using the given syntax, 'option_name'='value'.
	When option_name do not expect any value, set the option_name value with boolean to TRUE (enable), FALSE (disable).

  [samtools_sort]
    Please refer to samtools vew documentation: samtools sort --help
    No current options are available.

	Feel free to add samtools sort options, using the given syntax, 'option_name'='value'.
	When option_name do not expect any value, set the option_name value with boolean to TRUE (enable), FALSE (disable).

  [samtools_index]
    Please refer to samtools vew documentation: samtools index --help
    No current options are available.

	Feel free to add samtools index options, using the given syntax, 'option_name'='value'.
	When option_name do not expect any value, set the option_name value with boolean to TRUE (enable), FALSE (disable).

Notes: 1. Current mapper is gsnap 
          $(gsnap --version 2>&1 | awk -F"\n" 'BEGIN{info=""}; {info=(info "\n\t" $1)}; END{printf info}')
       2. Current sam toolkit is samtools
          $(samtools 2>&1 | egrep 'Program|Version'| awk -F"\n" 'BEGIN{info=""}; {info=(info "\n\t" $1)}; END{printf info}')

$(echo -ne $AUTHOR_INFOS)

";
exit 1; }

#=======
# BEGIN
#=======

echo "$(date '+%Y_%m_%d %T') [$(basename $0)] Start running the pipeline (version: $VERSION)." | tee $ERROR_TMP 2>&1
echo "$(date '+%Y_%m_%d %T') [$(basename $0)] Executed command: $0 $*" | tee -a $ERROR_TMP 2>&1

#
# Create a directory named with JOB_TAG value, to save all outputs 
#
echo "$(date '+%Y_%m_%d %T') [Job directory] Creating $JOB_TAG directory ..." | tee -a $ERROR_TMP 2>&1
if [[ -d $JOB_TAG ]]; then
    echo "$(date '+%Y_%m_%d %T') [Job directory] OK $JOB_TAG directory already exists. Will output all job files in this directory." | tee -a $ERROR_TMP 2>&1
else
    mkdir $JOB_TAG 2>>$ERROR_TMP
	rtrn=$?
	job_dir_failed_msg="[Job directory] Failed Job directory, $JOB_TAG, was not created."
	exit_on_error "$ERROR_TMP" "$job_dir_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
	echo "$(date '+%Y_%m_%d %T') [Job directory] OK $JOB_TAG directory was created successfully. Will output all job files in this directory." | tee -a $ERROR_TMP 2>&1
fi

# Create log directory
echo "$(date '+%Y_%m_%d %T') [Log directory] Creating $LOG_DIR directory ..." | tee -a $ERROR_TMP 2>&1
if [[ -d $LOG_DIR ]]; then
    [[ -s $ERROR_TMP ]] && cat $ERROR_TMP > $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %T') [Log directory] OK $LOG_DIR directory already exists. Will write log files in this directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
else
    mkdir $LOG_DIR 2>>$ERROR_TMP
	rtrn=$?
	log_dir_failed_msg="[Log directory] Failed Log directory, $LOG_DIR, was not created."
	exit_on_error "$ERROR_TMP" "$log_dir_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"   
	[[ -s $ERROR_TMP ]] && cat $ERROR_TMP > $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y_%m_%d %T') [Log directory] OK $LOG_DIR directory was created sucessfully. Will write log files in this directory." | tee -a $LOG_DIR/$LOGFILE 2>&1	
fi

#
# Test for cpu average load
# 
echo -ne "$(date '+%Y_%m_%d %T') [CPU load] Checking for average cpu load every $CPU_CHECK_INTERVAL seconds until $CPU_CHECK_TIMEOUT_HR timeout ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
cpu_load_failed_msg="[CPU load] Failed checking for average cpu load."	
timeout=$CPU_CHECK_TIMEOUT
until [[ $(isCpuAvailable $CORES_REDUC_FACTOR $SYS_CORES_AMOUNT 2>${ERROR_TMP})  == "TRUE" ]]; do
	rtrn=$?
	exit_on_error "$ERROR_TMP" "$cpu_load_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"	
	echo -e "." | tee -a $LOG_DIR/$LOGFILE 2>&1
	timeout=$(exit_on_timeout $$ "$LOG_DIR/$LOGFILE" $timeout $CPU_CHECK_INTERVAL $CPU_CHECK_DELAY 2>>${ERROR_TMP})
done
rtrn=$?
exit_on_error "$ERROR_TMP" "$cpu_load_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"	
echo -e "done" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "$(date '+%Y_%m_%d %T') [CPU load] $(uptime)" | tee -a $LOG_DIR/$LOGFILE 2>&1

#
# Check for DATA_ROOT_DIR existence
#
echo "$(date '+%Y_%m_%d %T') [Data root directory] Checking $DATA_ROOT_DIR directory ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
if [[ -d $DATA_ROOT_DIR ]]; then
    echo "$(date '+%Y_%m_%d %T') [Data root directory] $DATA_ROOT_DIR exists and is a directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
else
	data_root_dir_failed_msg="[Data root directory] Failed $DATA_ROOT_DIR does not exist or is not a directory."
	echo -e $data_root_dir_failed_msg 2>&1 >$ERROR_TMP	
	exit_on_error "$ERROR_TMP" "$data_root_dir_failed_msg" 1 "$LOG_DIR/$LOGFILE"
fi

#
# Test for disk space: 
# warning only about available disk space and if enough space or not to process all samples in data root dir
#
echo -e "$(date '+%Y_%m_%d %T') [Disk space] Checking for available disk space ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
disk_space_failed_msg="[Disk space] Failed while checking for available disk space for current $JOB_TAG job."
complete_dataset=($(ls -d -1 $DATA_ROOT_DIR/*))
res=$(isDiskSpaceAvailable $PWD $DATA_EXPANSION_FACTOR "${complete_dataset[@]}" 2>${ERROR_TMP})
rtrn=$?
if [[ $res == "TRUE" ]]; then
	exit_on_error "$ERROR_TMP" "$disk_space_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"		
	echo -e "$(date '+%Y_%m_%d %T') [Disk space] OK Enough disk space to process the complete dataset for data expansion factor $DATA_EXPANSION_FACTOR." | tee -a $LOG_DIR/$LOGFILE 2>&1
else
	exit_on_error "$ERROR_TMP" "$disk_space_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
	echo -e "$(date '+%Y_%m_%d %T') [Disk space] Warning Not enough disk space to process the complete dataset for data expansion factor $DATA_EXPANSION_FACTOR." | tee -a $LOG_DIR/$LOGFILE 2>&1
fi
echo -e "$(date '+%Y_%m_%d %T') [Disk space] available disk space on working directory partition:\n$(df -h $PWD)" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "$(date '+%Y_%m_%d %T') [Disk space] dataset expected disk space expansion for data expansion factor $DATA_EXPANSION_FACTOR:\n$(du -shc "${complete_dataset[@]}" | tail -n 1 | awk '{print $1}') x $DATA_EXPANSION_FACTOR" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "$(date '+%Y_%m_%d %T') [Disk space] Checking for available disk space done" | tee -a $LOG_DIR/$LOGFILE 2>&1

#
# Test for absence of user config file
# if present ok continue else display warning message and exit
#
echo "$(date '+%Y_%m_%d %T') [Check config: user config file] Checking for $VARCO_SPLIT_MAP_USER_CONFIG user config file ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
if [[ -s $VARCO_SPLIT_MAP_USER_CONFIG ]]; then
#if [[ -s $WORKING_DIR/$(basename ${0%.*})_user.config ]]; then # for testing purpose
    echo "$(date '+%Y_%m_%d %T') [Check config: user config file] OK User config file, $VARCO_SPLIT_MAP_USER_CONFIG, exists and is not empty." | tee -a $LOG_DIR/$LOGFILE 2>&1
else
	user_config_failed_msg="[Check config: user config file] Failed User config file, $VARCO_SPLIT_MAP_USER_CONFIG, does not exist or is empty."
	echo -e "$user_config_failed_msg" 2>&1 >$ERROR_TMP
    echo "$(date '+%Y_%m_%d %T') [Check config: user config file] Warning: " | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo -e "\t\t$PREREQUISITES_MSG" | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit_on_error "$ERROR_TMP" "$user_config_failed_msg" 3 "$LOG_DIR/$LOGFILE"
fi

#
# Copy user config parameters file into job directory
# 1. Copy user config parameters file, prefixing it with the job tag
# 2. Then load config parameters from that new file, leaving the user config file in the working directory for some new job

# 1. Copy user config file 
echo "$(date '+%Y_%m_%d %T') [Check config: job user config file] Copying user config file into job directory ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
cp $VARCO_SPLIT_MAP_USER_CONFIG $VARCO_SPLIT_MAP_USER_CONFIG_JOB 2>$ERROR_TMP
rtrn=$?
cp_user_config_failed_msg="[Check config: job user config file] Failed copying user config file into job directory."
exit_on_error "$ERROR_TMP" "$cp_user_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
echo "$(date '+%Y_%m_%d %T') [Check config: job user config file] Will use copied job user config file: $VARCO_SPLIT_MAP_USER_CONFIG_JOB" | tee -a $LOG_DIR/$LOGFILE 2>&1

# 2. Load config parameters from job user config file
echo "$(date '+%Y_%m_%d %T') [Check config: job user config file] Loading job user config parameters from $VARCO_SPLIT_MAP_USER_CONFIG_JOB file ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
load_user_config_failed_msg="[Check config: job user config file] Failed loading job user config parameters from $VARCO_SPLIT_MAP_USER_CONFIG_JOB file."
for cfg in $(get_config_sections $VARCO_SPLIT_MAP_USER_CONFIG_JOB 2>$ERROR_TMP;); do
	rtrn=$?	
	exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
    echo -e "--- Config section [${cfg}] ---"
    unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${NAMESPACE}" 'BEGIN { 
          cfg = toupper(cfg);
          prefix = toupper(prefix);
       }
       /^prefix_cfg_/  { print $1 }' 2>$ERROR_TMP) $(toupper ${NAMESPACE}_${cfg}_) 2>>$ERROR_TMP
	rtrn=$?
	exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
    set_config_params $VARCO_SPLIT_MAP_USER_CONFIG_JOB ${cfg} ${NAMESPACE} 2>$ERROR_TMP
    rtrn=$?
	exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
    for params in $(set | grep ^$(toupper ${NAMESPACE}_${cfg}_) 2>$ERROR_TMP); do
		echo -e "$params"
    done
	rtrn=$?
	exit_on_error "$ERROR_TMP" "$load_user_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
done
echo "$(date '+%Y_%m_%d %T') [Check config: job user config file] OK User config file, $VARCO_SPLIT_MAP_USER_CONFIG_JOB, was loaded successfully." | tee -a $LOG_DIR/$LOGFILE 2>&1

#
# Check for parameters validity
# cf lib implement a function to check for config parameters validity: TODO
# check only internal pipeline config params not vendor config params

#
# Override  defined batch_size:
# 1. get the total number of cores
# 2. get the number of threads to use for the mapper
# 3. compute max_batch_size=#max_cores_allowed/#threads
# 4. if batch_size <= max_batch_size then ok use batch_size else batch_size=max_batch_size/2
#
#VARCO_SPLIT_MAP_batch_size=8 # for testing purpose
echo "$(date '+%Y_%m_%d %T') [Override user config: batch_size] Testing if override batch_size user defined config parameter value is needed ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
MAX_BATCH_SIZE=$[ $MAX_NUMB_CORES_ALLOWED/$VARCO_GSNAP_t ] 2>$ERROR_TMP
rtrn=$?
max_batch_size_failed_msg="[Override user config: batch_size] Failed computing max batch size."
exit_on_error "$ERROR_TMP" "$max_batch_size_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
echo -e "max_cores=$MAX_NUMB_CORES" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "max_cores_allowed=$MAX_NUMB_CORES_ALLOWED" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "threads_by_sample=$VARCO_GSNAP_t" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "max_batch_size=$MAX_BATCH_SIZE" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "batch_size=$VARCO_SPLIT_MAP_batch_size" | tee -a $LOG_DIR/$LOGFILE 2>&1
if [[ $VARCO_SPLIT_MAP_batch_size -le $MAX_BATCH_SIZE ]]; then
	echo "$(date '+%Y_%m_%d %T') [Override user config: batch_size] No need to override batch_size user defined config parameter value." | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y_%m_%d %T') [Override user config: batch_size] Keep batch_size user defined config parameter value: $VARCO_SPLIT_MAP_batch_size" | tee -a $LOG_DIR/$LOGFILE 2>&1
else
	echo "$(date '+%Y_%m_%d %T') [Override user config: batch_size] Need to override batch_size user defined config parameter value." | tee -a $LOG_DIR/$LOGFILE 2>&1
	VARCO_SPLIT_MAP_batch_size=$[$MAX_BATCH_SIZE/2] 2>$ERROR_TMP
	rtrn=$?
	batch_size_failed_msg="[Override user config: batch_size] Failed computing batch size."
	exit_on_error "$ERROR_TMP" "$batch_size_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
	echo "$(date '+%Y_%m_%d %T') [Override user config: batch_size] Override batch_size user defined config parameter value: $VARCO_SPLIT_MAP_batch_size." | tee -a $LOG_DIR/$LOGFILE 2>&1
fi

#
# Search for subdirectories with fastq files:
# 1. List all available subdirs in data root dir
# 2. Filter subdirs with include pattern(s)
# 3. Filter subdirs with exclude pattern(s)
# 4. Filter subdirs for fastq files

# 1. list all available subdirs in data root dir
echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] Searching for fastq sample subdirs ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
all_subdirs=($(find $DATA_ROOT_DIR -maxdepth 1 -type d | egrep -v ^$DATA_ROOT_DIR$ 2>$ERROR_TMP))
rtrn=$?
all_subdirs_failed_msg="[Fastq subdirs] Failed listing all sub-directories in $DATA_ROOT_DIR directory."
exit_on_error "$ERROR_TMP" "$all_subdirs_failed_msg" $rtrn "$LOG_DIR/$LOGFILE" 

#read -a all_subdirs_arr echo <<< $(echo -e $all_subdirs | tr " " "\n")
echo -e "all subdirectories count: ${#all_subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "all subdirectories list: ${all_subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1

# 2. Filter subdirs with include pattern(s)
echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] Filtering subdirs with include pattern(s) ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
#VARCO_DATA_include_sample_subdirs="^D.*,XA,XC" # for testing purpose
include_patterns=$(echo ${VARCO_DATA_include_sample_subdirs//,/|})
echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] include_patterns=($include_patterns) ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
inc_subdirs_failed_msg="[Fastq subdirs] Failed filtering subdirs with include pattern(s): $include_patterns"
inc_subdirs=($(for subdir in "${all_subdirs[@]}"; do
    res=$(echo $(basename $subdir) | egrep "($include_patterns)" 2>$ERROR_TMP)
    [[ -n $res ]] && echo -e "$subdir" 
done 2>>$ERROR_TMP))
[[ -s $ERROR_TMP ]] && rtrn=1 || rtrn=0
exit_on_error "$ERROR_TMP" "$inc_subdirs_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
echo -e "include pattern(s) subdirectories count: ${#inc_subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "include pattern(s) subdirectories list: ${inc_subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1

# 3. Filter subdirs with exclude pattern(s)
echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] Filtering subdirs with exclude pattern(s) ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
#VARCO_DATA_exclude_sample_subdirs="^X.*,DA,DI" # for testing purpose
#VARCO_DATA_exclude_sample_subdirs="" # for testing purpose
if [[ -n $VARCO_DATA_exclude_sample_subdirs ]]; then
    exclude_patterns=$(echo ${VARCO_DATA_exclude_sample_subdirs//,/|})
    echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] exclude_patterns=($exclude_patterns) ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
    subdirs=($(for subdir in "${inc_subdirs[@]}"; do
	    res=$(echo $(basename $subdir) | egrep -v "($exclude_patterns)" 2>$ERROR_TMP)
	    rtrn=$?
	    [[ -n $res ]] && echo -e "$subdir" 
	    done 2>>$ERROR_TMP))
	rtrn=$?
	excl_subdirs_failed_msg="[Fastq subdirs] Failed filtering subdirs with exclude pattern(s): $exclude_patterns"
	exit_on_error "$ERROR_TMP" "$excl_subdirs_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
else
    echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] No exclude patterns." | tee -a $LOG_DIR/$LOGFILE 2>&1
    subdirs=("${inc_subdirs[@]}")
fi
echo -e "include pattern(s) subdirectories count: ${#subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE
echo -e "include pattern(s) subdirectories list: ${subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE

# 4. Filter subdirs for fastq files
echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] Filtering subdirs for fastq files ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
fastq_forward_pattern="*_1_*.fastq"
fastq_reverse_pattern="*_2_*.fastq"
fastq_failed_message="[Fastq subdirs] Failed An error occured while filtering for subdirs with fastq files."
fastq_subdirs=($(for subdir in "${subdirs[@]}"; do
	fastq_files=($(ls "$subdir" | egrep -v "_single_" | egrep ".*.fastq$" 2>$ERROR_TMP))
	exit_on_error "$ERROR_TMP" "$fastq_failed_message" $? "$LOG_DIR/$LOGFILE"
	if [[ "${#fastq_files[@]}" == 2 ]]; then
	    echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] INFO $subdir has 2, non '_single_', fastq files." 1>&2 | tee -a $LOG_DIR/$LOGFILE 1>&2
        # get the pair if possible
	    # else warn and discard the subdir
	    sample_name_1=$(getFastqSampleName "${fastq_files[0]}" 2>>$ERROR_TMP)
		exit_on_error "$ERROR_TMP" "$fastq_failed_message" $? "$LOG_DIR/$LOGFILE"
	    sample_name_2=$(getFastqSampleName "${fastq_files[1]}" 2>>$ERROR_TMP)
	    exit_on_error "$ERROR_TMP" "$fastq_failed_message" $? "$LOG_DIR/$LOGFILE"

	    if [[ "$sample_name_1" == "$sample_name_2" ]]; then
		echo -e "$(date '+%Y_%m_%d %T') [Fastq subdirs] INFO $subdir, same fastq sample name: forward=$sample_name_1 == reverse=$sample_name_2" 1>&2 | tee -a $LOG_DIR/$LOGFILE 1>&2 
		echo -e "$subdir" 
	    else
		echo -e "$(date '+%Y_%m_%d %T') [Fastq subdirs] Warning $subdir, different fastq sample name: forward=$sample_name_1 != reverse=$sample_name_2" 1>&2 | tee -a $LOG_DIR/$LOGFILE 1>&2
		echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] Warning: $subdir will not be considered because fastq sample name are different." 1>&2| tee -a $LOG_DIR/$LOGFILE 1>&2
	    fi
	elif [[ "${#fastq_files[@]}" < 2 ]]; then
	    echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] Warning: $subdir has less than 2, non '_single_', fastq files." 1>&2 | tee -a $LOG_DIR/$LOGFILE 1>&2
	    echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] Warning: $subdir will not be considered because non '_single_' fastq files are missing." 1>&2| tee -a $LOG_DIR/$LOGFILE 1>&2
	elif [[ "${#fastq_files[@]}" > 2 ]]; then
	    # not satisfying
	    echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] Warning: $subdir has more than 2, non '_single_', fastq files." 1>&2 | tee -a $LOG_DIR/$LOGFILE 1>&2
	    echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] Warning: You need to leave in $subdir only the 2 forward (_1_) and reverse (_2_) fastq files. '_single_' fastq files are ignored." 1>&2 | tee -a $LOG_DIR/$LOGFILE 1>&2
	    echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] Warning: $subdir will not be considered, ambiguous fastq files list." 1>&2 | tee -a $LOG_DIR/$LOGFILE 1>&2
	fi
	done 2>>$LOG_DIR/$LOGFILE))

    echo -e "fastq subdirectories count: ${#fastq_subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo -e "fastq subdirectories list: ${fastq_subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1
    if [[ "${#fastq_subdirs[@]}" -ne "${#subdirs[@]}" ]]; then
	echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] Warning Some subdirectories was discarded because fastq files checking has detected divergent forward/reverse sample name or an unexpected fastq files count." | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] Warning More details about discarded subdirectories can be found in $LOG_DIR/$LOGFILE" | tee -a $LOG_DIR/$LOGFILE 2>&1
    fi

#
# Test for available disk space on fastq subdirs
#
echo -e "$(date '+%Y_%m_%d %T') [Disk space] Checking for available disk space only for fastq subdirs ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
disk_space_failed_msg="[Disk space] Failed while checking for available disk space only for fastq subdirs."
res=$(isDiskSpaceAvailable $PWD $DATA_EXPANSION_FACTOR "${fastq_subdirs[@]}" 2>${ERROR_TMP})
rtrn=$?
if [[ $res == "TRUE" ]]; then
	exit_on_error "$ERROR_TMP" "$disk_space_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"		
	echo -e "$(date '+%Y_%m_%d %T') [Disk space] OK Enough disk space to process only fastq subdirs for data expansion factor $DATA_EXPANSION_FACTOR." | tee -a $LOG_DIR/$LOGFILE 2>&1
else
	exit_on_error "$ERROR_TMP" "$disk_space_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
	echo -e "$(date '+%Y_%m_%d %T') [Disk space] Warning Not enough disk space to process only fastq subdirs for data expansion factor $DATA_EXPANSION_FACTOR." | tee -a $LOG_DIR/$LOGFILE 2>&1
fi
echo -e "$(date '+%Y_%m_%d %T') [Disk space] available disk space on working directory partition:\n$(df -h $PWD)" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "$(date '+%Y_%m_%d %T') [Disk space] dataset expected disk space expansion for data expansion factor $DATA_EXPANSION_FACTOR:\n$(du -shc "${fastq_subdirs[@]}" | tail -n 1 | awk '{print $1}') x $DATA_EXPANSION_FACTOR" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "$(date '+%Y_%m_%d %T') [Disk space] Checking for available disk space done" | tee -a $LOG_DIR/$LOGFILE 2>&1

#
# Batch mode:
# 1. Iterate over batches
# 1.1 Test for average cpu load
# 1.2 Test available disk space
# 1.3 Iterate over batch samples for getting sample infos
# 1.3.1 Create an output directory for each sample
# 1.4 Iterate over batch samples for Quality Control and Trimming
# 1.4.1 Create Quality control and trimming sub-subdir for each sample: fastqc and trimmomatic, optionnal step
# 1.4.2 Perform QC before -> waitall
# 1.4.3 Perform Trimming -> waitall
# 1.4.4 Perform QC after -> waitall
# 1.4.5 Wait for all QC+Trim processes to finish
# 1.5 Iterate over batch samples for Mapping
# 1.5.1 Create a mapping sub-subdir for each sample
# 1.5.2 Create mapping subdirs
# 1.5.3 Perform mapping -> waitall
# 1.5.4 Wait for all Mapping processes to finish
# 1.6 Iterate over commands for Conversion
# 1.6.1 Iterate over batch samples
# 1.6.2 Perform bam conversion -> waitall
# 1.6.3 Perform bam sorting -> waitall
# 1.6.4 Perform bam indexing -> waitall
# 1.6.5 Wait for all Conversion processes to finish
# 1.7 Unshifting current batch samples
# 1.8 Next batch
# 1.9 End of batches
# 1.10 Clean

echo "$(date '+%Y_%m_%d %T') [Batch mode] Running batch mode on fastq subdirs ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
unsorted_subdirs=("${fastq_subdirs[@]}")
readarray -t subdirs < <(printf '%s\0' "${unsorted_subdirs[@]}" | sort -z | xargs -0n1)

# Computing number of batches expected to be run
## if #subdirs <= batch_size then #batches=1
## if #subdirs > batch_size then #batches=ceil(#subdirs/batch_size)
echo "$(date '+%Y_%m_%d %T') [Batch mode] Computing number of batches to run ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
if [[ "${#subdirs[@]}" -le "$VARCO_SPLIT_MAP_batch_size" ]]; then
	batches=1
else
	batches=$(echo "${#subdirs[@]} $VARCO_SPLIT_MAP_batch_size" | awk '{print int( ($1/$2) + 1 )}' 2>ERROR_TMP)
	rtrn=$?
	batches_failed_msg="[Batch mode] Failed computing the number of batches expected to be run."
	exit_on_error "$ERROR_TMP" "$batches_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
fi
echo "$(date '+%Y_%m_%d %T') [Batch mode] $batches computed batche(s) expected to be run." | tee -a $LOG_DIR/$LOGFILE 2>&1

# 1. Iterate over batches
for b in $(seq 1 $[ $batches ]); do

    # 1.1 Test for average cpu load:
	if [[ $VARCO_SPLIT_MAP_check_cpu_overload == "TRUE" ]]; then
		echo -ne "$(date '+%Y_%m_%d %T') [CPU load] Checking for average cpu load before running on samples batch #$b every $CPU_CHECK_INTERVAL seconds until $CPU_CHECK_TIMEOUT_HR timeout ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
		cpu_load_failed_msg="[CPU load] Failed checking for average cpu load before running on samples batch #$b."	
		timeout=$CPU_CHECK_TIMEOUT		
		until [[ $(isCpuAvailable $CORES_REDUC_FACTOR $SYS_CORES_AMOUNT 2>${ERROR_TMP})  == "TRUE" ]]; do
			rtrn=$?
			exit_on_error "$ERROR_TMP" "$cpu_load_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"	
			echo -e "." | tee -a $LOG_DIR/$LOGFILE 2>&1
			timeout=$(exit_on_timeout $$ "$LOG_DIR/$LOGFILE" $timeout $CPU_CHECK_INTERVAL $CPU_CHECK_DELAY 2>>${ERROR_TMP})
		done
		rtrn=$?
		exit_on_error "$ERROR_TMP" "$cpu_load_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"	
		echo -e "done" | tee -a $LOG_DIR/$LOGFILE 2>&1
		echo -e "$(date '+%Y_%m_%d %T') [CPU load] $(uptime)" | tee -a $LOG_DIR/$LOGFILE 2>&1
	else
		echo -e "$(date '+%Y_%m_%d %T') [CPU load] Skipping checking for average cpu load before running on samples batch #$b." | tee -a $LOG_DIR/$LOGFILE 2>&1
		echo -e "$(date '+%Y_%m_%d %T') [CPU load] $(uptime)" | tee -a $LOG_DIR/$LOGFILE 2>&1
	fi
	
    # 1.2 Test for available disk space for current samples batch
	echo -e "$(date '+%Y_%m_%d %T') [Disk space] Checking for available disk space before running on samples batch #$b ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
	disk_space_failed_msg="[Disk space] Failed while checking for available disk space only for samples batch #$b."
	last=FALSE
	samples_batch_dirs=($(for s in $(seq 1 $VARCO_SPLIT_MAP_batch_size); do
		si=$[$s-1]
		[[ "$s" -eq "${#subdirs[@]}" ]] && last=TRUE

		echo "${subdirs[$si]}"

	# Last batch sample: have a break!
		[[ $last == "TRUE" ]] && break
	done))
	res=$(isDiskSpaceAvailable $PWD $DATA_EXPANSION_FACTOR "${samples_batch_dirs[@]}" 2>${ERROR_TMP})
	rtrn=$?
	if [[ $res == "TRUE" ]]; then
		exit_on_error "$ERROR_TMP" "$disk_space_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"		
		echo -e "$(date '+%Y_%m_%d %T') [Disk space] OK Enough disk space to process only samples batch #$b for data expansion factor $DATA_EXPANSION_FACTOR." | tee -a $LOG_DIR/$LOGFILE 2>&1
	else
		exit_on_error "$ERROR_TMP" "$disk_space_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
		disk_space_runout_msg="[Disk space] Warning Not enough disk space to process only samples batch #$b for data expansion factor $DATA_EXPANSION_FACTOR. Abort running this pipeline."
		echo -e "$(date '+%Y_%m_%d %T') [Disk space] available disk space on working directory partition:\n$(df -h $PWD)" >$ERROR_TMP
		echo -e "$(date '+%Y_%m_%d %T') [Disk space] dataset expected disk space expansion for data expansion factor $DATA_EXPANSION_FACTOR:\n$(du -shc "${samples_batch_dirs[@]}" | tail -n 1 | awk '{print $1}') x $DATA_EXPANSION_FACTOR" >>$ERROR_TMP
		exit_on_error "$ERROR_TMP" "$disk_space_runout_msg" 1 "$LOG_DIR/$LOGFILE"		
	fi
	echo -e "$(date '+%Y_%m_%d %T') [Disk space] available disk space on working directory partition:\n$(df -h $PWD)" | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo -e "$(date '+%Y_%m_%d %T') [Disk space] dataset expected disk space expansion for data expansion factor $DATA_EXPANSION_FACTOR:\n$(du -shc "${samples_batch_dirs[@]}" | tail -n 1 | awk '{print $1}') x $DATA_EXPANSION_FACTOR" | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo -e "$(date '+%Y_%m_%d %T') [Disk space] Checking for available disk space done" | tee -a $LOG_DIR/$LOGFILE 2>&1

	# 1.3 Iterate over batch samples for getting sample infos
	last=FALSE 	# for last batch sample
    echo "$(date '+%Y_%m_%d %T') [Batch mode] Running batch mode on samples batch #$b ..." | tee -a $LOG_DIR/$LOGFILE 2>&1   
    for s in $(seq 1 $VARCO_SPLIT_MAP_batch_size); do
		si=$[$s-1]
		sdi=$[$si+($b-1)*$VARCO_SPLIT_MAP_batch_size]
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode] Processing batch sample $sdi, $si: ${subdirs[$si]}"
		CURRENT_BATCH_SUBDIR=$JOB_TAG/$(basename "${subdirs[$si]}")
		sample_dir_failed_msg="[Batch mode: sample output directory] Failed Sample output directory, $CURRENT_BATCH_SUBDIR, was not created."
	# 1.3.1 Create an output directory for each sample
		if [[ "$s" -lt "${#subdirs[@]}" ]]; then
	    	echo "$(date '+%Y_%m_%d %T') [Batch mode] Creating sample output directory: $CURRENT_BATCH_SUBDIR" | tee -a $LOG_DIR/$LOGFILE 2>&1
	    	if [[ -d $CURRENT_BATCH_SUBDIR ]]; then
				echo "$(date '+%Y_%m_%d %T') [Batch mode: sample output directory] OK $CURRENT_BATCH_SUBDIR directory already exists. Will write output files in this directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
	    	else
				mkdir $CURRENT_BATCH_SUBDIR 2>$ERROR_TMP
				rtrn=$?
				exit_on_error "$ERROR_TMP" "$sample_dir_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"		
				echo "$(date '+%Y_%m_%d %T') [Batch mode: sample output directory] OK $CURRENT_BATCH_SUBDIR directory was created successfully. Will write output files in this directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
	    	fi
		elif [[ "$s" -eq "${#subdirs[@]}" ]]; then
	    	last=TRUE
	    	echo "$(date '+%Y_%m_%d %T') [Batch mode] Creating last batch sample output directory: $CURRENT_BATCH_SUBDIR" | tee -a $LOG_DIR/$LOGFILE 2>&1
	    	if [[ -d $CURRENT_BATCH_SUBDIR ]]; then
				echo "$(date '+%Y_%m_%d %T') [Batch mode: sample output directory] OK $CURRENT_BATCH_SUBDIR last batch directory already exists. Will write output files in this directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
	    	else
				mkdir $CURRENT_BATCH_SUBDIR 2>$ERROR_TMP
				rtrn=$?
				exit_on_error "$ERROR_TMP" "$sample_dir_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
	    		echo "$(date '+%Y_%m_%d %T') [Batch mode: sample output directory] OK $CURRENT_BATCH_SUBDIR last batch directory was created successfully. Will write output files in this directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
	    	fi  
		fi

	# 1.3.2 Get forward and reverse fastq files for current batch sample subdir
		echo "$(date '+%Y_%m_%d %T') [Batch mode] Listing forward and reverse fastq files in ${subdirs[$si]} current batch sample directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
		fastq_files=($(ls "${subdirs[$si]}" | egrep -v "_single_" | egrep ".*.fastq$" 2>$ERROR_TMP))
		rtrn=$?
		fastq_list_failed_msg="[Batch mode] Failed An error occured while listing forward and reverse fastq files in ${subdirs[$si]} current batch sample directory."
		exit_on_error "$ERROR_TMP" "$fastq_list_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"	
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode] fastq files count: ${#fastq_files[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode] fastq files list: ${fastq_files[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1
		getFastqFile()
		{
			direction=$1
			for f in "${fastq_files[@]}"; do 
				m=$(echo $f | egrep "_[0-9]+_${direction}_" 2>$ERROR_TMP); if [[ -n $m ]]; then echo $m; break; fi
			done
		}
		fastq_file_failed_msg="[Batch mode] Failed An error occured while getting forward/reverse fastq files in ${subdirs[$si]} current batch sample directory."		
		forward_fastq=$(getFastqFile 1)
		rtrn=$?
		exit_on_error "$ERROR_TMP" "$fastq_file_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
		reverse_fastq=$(getFastqFile 2)
		rtrn=$?
		exit_on_error "$ERROR_TMP" "$fastq_file_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
	
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode] forward fastq file: $forward_fastq" | tee -a $LOG_DIR/$LOGFILE 2>&1
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode] reverse fastq file: $reverse_fastq" | tee -a $LOG_DIR/$LOGFILE 2>&1
		
	# 1.3.3 Get current sample name
		echo "$(date '+%Y_%m_%d %T') [Batch mode] Getting sample name for ${subdirs[$si]} current batch sample directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
		sample_name=$(getFastqSampleName $forward_fastq 2>$ERROR_TMP)
		rtrn=$?
		fastq_sample_name_failed_msg="[Batch mode] Failed An error occured while getting sample name for ${subdirs[$si]} current batch sample directory."
		exit_on_error "$ERROR_TMP" "$fastq_sample_name_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode] OK Current sample name: $sample_name" | tee -a $LOG_DIR/$LOGFILE 2>&1

	# 1.3.4 Save fastq infos to config file
		echo -ne "$(date '+%Y_%m_%d %T') [Batch mode] Saving sample infos for ${subdirs[$si]} current batch sample directory ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
		sample_config_file=$CURRENT_BATCH_SUBDIR/$SAMPLE_CONFIG
		echo -e "[sample]" 2>$ERROR_TMP >$sample_config_file
		echo -e "name=${sample_name}" 2>>$ERROR_TMP >>$sample_config_file
		echo -e "fastq_forward=${forward_fastq}" 2>>$ERROR_TMP >>$sample_config_file
		echo -e "fastq_reverse=${reverse_fastq}" 2>>$ERROR_TMP >>$sample_config_file
		echo -e "done" | tee -a $LOG_DIR/$LOGFILE 2>&1
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode] sample config file:" | tee -a $LOG_DIR/$LOGFILE 2>&1
		echo -e "$(cat $sample_config_file)" | tee -a $LOG_DIR/$LOGFILE 2>&1
	# Last batch sample: have a break!
		[[ $last == "TRUE" ]] && break
	done

	# 1.4 Iterate over quality control and trimming commands: optional step TODO
	qc_cmd="fastqc"
	trimming_cmd="trimmomatic"
	qc_trim_cmds=("$qc_cmd" "$trimming_cmd" "$qc_cmd")
	trimmed=FALSE

	#VARCO_QC_TRIM_process=TRUE # for testing purpose
	if [[ $VARCO_QC_TRIM_process == "TRUE" ]]; then

	# Iterate over qc_trim_cmds
	for cmd in "${qc_trim_cmds[@]}"; do
		last=FALSE 	# for last batch sample
    	echo "$(date '+%Y_%m_%d %T') [Batch mode: quality control and trimming] Running quality control and trimming on samples batch #$b ..." | tee -a $LOG_DIR/$LOGFILE 2>&1 

	# reinitiate pids array
		PIDS_ARR=()
  
	# Iterate over batch samples for current quality control and trimming command
    	for s in $(seq 1 $VARCO_SPLIT_MAP_batch_size); do
			si=$[$s-1]
			sdi=$[$si+($b-1)*$VARCO_SPLIT_MAP_batch_size]
			echo -e "$(date '+%Y_%m_%d %T') [Batch mode: quality control and trimming] Processing batch sample $sdi, $si: ${subdirs[$si]}, for quality control and trimming ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
			CURRENT_BATCH_SUBDIR=$JOB_TAG/$(basename "${subdirs[$si]}")
			CURRENT_SAMPLE_CONFIG=$CURRENT_BATCH_SUBDIR/$SAMPLE_CONFIG
			
	# last batch sample
	[[ "$s" -eq "${#subdirs[@]}" ]] && last=TRUE

	# load sample config
			echo "$(date '+%Y_%m_%d %T') [Batch mode: quality control and trimming] Loading sample config for $CURRENT_SAMPLE_CONFIG file ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
			sample_config_failed_msg="[Batch mode: quality control and trimming] Failed loading sample config file, $CURRENT_SAMPLE_CONFIG, from $CURRENT_BATCH_SUBDIR"
			for cfg in $(get_config_sections $CURRENT_SAMPLE_CONFIG 2>$ERROR_TMP; rtrn=$?); do
				if [[ "$rtrn" -eq 0 ]]; then
    				echo -e "--- Config section [${cfg}] ---"
    				unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${NAMESPACE}" 'BEGIN { 
          				cfg = toupper(cfg);
          				prefix = toupper(prefix);
       					}
       					/^prefix_cfg_/  { print $1 }' 2>$ERROR_TMP) $(toupper ${NAMESPACE}_${cfg}_) 2>>$ERROR_TMP
					rtrn=$?
					exit_on_error "$ERROR_TMP" "$sample_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
    				set_config_params $CURRENT_SAMPLE_CONFIG ${cfg} ${NAMESPACE} 2>$ERROR_TMP
    				rtrn=$?
					exit_on_error "$ERROR_TMP" "$sample_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE" 
    				for params in $(set | grep ^$(toupper ${NAMESPACE}_${cfg}_) 2>$ERROR_TMP); do
						echo -e "$params"
    				done
				else
					exit_on_error "$ERROR_TMP" "$sample_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
				fi
			done

	# 1.4.1 Create Quality control and trimming sub-subdir: fastqc and trimmomatic, 
			echo "$(date '+%Y_%m_%d %T') [Batch mode: quality control and trimming] Creating quality control and trimming sub-directory: $CURRENT_BATCH_SUBDIR/$QC_TRIM_DIR ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
			CURRENT_QCTRIM_LOG=$CURRENT_BATCH_SUBDIR/$QC_TRIM_DIR/qc_trim.log
			CURRENT_QCTRIM_ERR=$CURRENT_BATCH_SUBDIR/$QC_TRIM_DIR/qc_trim_err.log
			if [[ -d $CURRENT_BATCH_SUBDIR/$QC_TRIM_DIR ]]; then
				echo "$(date '+%Y_%m_%d %T') [Batch mode: quality control and trimming: Quality Control and Trimming] OK $CURRENT_BATCH_SUBDIR/$QC_TRIM_DIR directory already exists. Will write output files in this directory." | tee -a $CURRENT_QCTRIM_LOG 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
			else
				mkdir $CURRENT_BATCH_SUBDIR/$QC_TRIM_DIR 2>$ERROR_TMP
				rtrn=$?
				qctrim_dir_failed_msg="[Batch mode: quality control and trimming] Failed Quality control and Trimming output directory, $CURRENT_BATCH_SUBDIR/$QC_TRIM_DIR, was not created."
				exit_on_error "$ERROR_TMP" "$qctrim_dir_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
				echo "$(date '+%Y_%m_%d %T') [Batch mode: quality control and trimming] OK $CURRENT_BATCH_SUBDIR/$QC_TRIM_DIR directory was created successfully. Will write output files in this directory." | tee -a $CURRENT_QCTRIM_LOG 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
			fi 
	
	# Build cli for the 3 cases: in fact 2 cases 1) fastqc with a) before trimming and b) after trimming 2) trimmomatic
	# 1.4.1 Quality control, 
	# 1.4.1.a before trimming (raw data): TODO
	# 1.4.1.b after trimming (trimmed data): TODO
 
	# 1.4.2 Trimming on raw data: TODO
	# ${forward_fastq_both_surviving}
	# ${reverse_fastq_both_surviving}
	# ${forward_fastq_only_surviving}
	# ${reverse_fastq_only_surviving}

	# Run the cli
	
	# Last batch sample: have a break!
			[[ $last == "TRUE" ]] && break
		done

	# waitalluntiltimeout: TODO
	# several commands to wait: fastqc (trimmed=FALSE), trimmomatic, fastqc (trimmed=TRUE)


	# check for errors: TODO
	# for trimming, 
	# set trimmed to TRUE if ok
	# set forward and reverse fastq files for mapping step: update the sample config file
			echo -ne "$(date '+%Y_%m_%d %T') [Batch mode: quality control and trimming] Saving quality control and trimming infos for ${subdirs[$si]} current batch sample directory ... " | tee -a $CURRENT_QCTRIM_LOG 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
			echo -e "fastq_forward_trimmed_both_surviving=${forward_fastq_both_surviving}" 2>>$ERROR_TMP >>$CURRENT_SAMPLE_CONFIG
			echo -e "fastq_reverse_trimmed_both_surviving=${reverse_fastq_both_surviving}" 2>>$ERROR_TMP >>$CURRENT_SAMPLE_CONFIG
			echo -e "fastq_forward_trimmed_only_surviving=${forward_fastq_only_surviving}" 2>>$ERROR_TMP >>$CURRENT_SAMPLE_CONFIG
			echo -e "fastq_reverse_trimmed_only_surviving=${reverse_fastq_only_surviving}" 2>>$ERROR_TMP >>$CURRENT_SAMPLE_CONFIG
			echo -e "done" | tee -a $CURRENT_QCTRIM_LOG 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
			echo -e "$(date '+%Y_%m_%d %T') [Batch mode: quality control and trimming] updated sample config file:" | tee -a $CURRENT_QCTRIM_LOG 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
			echo -e "$(cat $sample_config_file)" | tee -a $CURRENT_QCTRIM_LOG 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1

	done
	fi

	# 1.5 Iterate over batch samples for mapping 
	
	# reinitiate last batch sample
	last=FALSE

	# reinitiate pids array
	PIDS_ARR=()

	# iterate over batch samples
	for s in $(seq 1 $VARCO_SPLIT_MAP_batch_size); do
		si=$[$s-1]
		sdi=$[$si+($b-1)*$VARCO_SPLIT_MAP_batch_size]
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] Processing batch sample $sdi, $si: ${subdirs[$si]}, for mapping ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
		CURRENT_BATCH_SUBDIR=$JOB_TAG/$(basename "${subdirs[$si]}")
		CURRENT_SAMPLE_CONFIG=$CURRENT_BATCH_SUBDIR/$SAMPLE_CONFIG

	# last batch sample
	[[ "$s" -eq "${#subdirs[@]}" ]] && last=TRUE

	# load sample config
		echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping] Loading sample config for $CURRENT_SAMPLE_CONFIG file ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
		sample_config_failed_msg="[Batch mode: mapping] Failed loading sample config file, $CURRENT_SAMPLE_CONFIG, from $CURRENT_BATCH_SUBDIR"
		for cfg in $(get_config_sections $CURRENT_SAMPLE_CONFIG 2>$ERROR_TMP; rtrn=$?); do
			if [[ "$rtrn" -eq 0 ]]; then
				echo -e "--- Config section [${cfg}] ---"
				unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${NAMESPACE}" 'BEGIN { 
	  				cfg = toupper(cfg);
	  				prefix = toupper(prefix);
					}
					/^prefix_cfg_/  { print $1 }' 2>$ERROR_TMP) $(toupper ${NAMESPACE}_${cfg}_) 2>>$ERROR_TMP
				rtrn=$?
				exit_on_error "$ERROR_TMP" "$sample_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
				set_config_params $CURRENT_SAMPLE_CONFIG ${cfg} ${NAMESPACE} 2>$ERROR_TMP
				rtrn=$?
				exit_on_error "$ERROR_TMP" "$sample_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE" 
				for params in $(set | grep ^$(toupper ${NAMESPACE}_${cfg}_) 2>$ERROR_TMP); do
					echo -e "$params"
				done
			else
				exit_on_error "$ERROR_TMP" "$sample_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
			fi
		done

	# 1.5 Create a mapping sub-subdir
		echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping] Creating mapping sub-directory: $CURRENT_BATCH_SUBDIR/$MAPPING_DIR" | tee -a $LOG_DIR/$LOGFILE 2>&1
		if [[ -d $CURRENT_BATCH_SUBDIR/$MAPPING_DIR ]]; then
			echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping] OK $CURRENT_BATCH_SUBDIR/$MAPPING_DIR directory already exists. Will write output files in this directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
		else
			mkdir $CURRENT_BATCH_SUBDIR/$MAPPING_DIR 2>$ERROR_TMP
			rtrn=$?
			mapping_dir_failed_msg="[Batch mode: mapping] Failed Mapping output directory, $CURRENT_BATCH_SUBDIR/$MAPPING_DIR was not created."
			exit_on_error "$ERROR_TMP" "$mapping_dir_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
			echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping] OK $CURRENT_BATCH_SUBDIR/$MAPPING_DIR directory was created successfully. Will write output files in this directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
		fi

	# 1.5.1 Create mapping subdirs
		mapping_subdirs=("tmp" "log")
		for msd in "${mapping_subdirs[@]}"; do
			echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping] Creating mapping sub-sub-directory: $CURRENT_BATCH_SUBDIR/$MAPPING_DIR/$msd" | tee -a $LOG_DIR/$LOGFILE 2>&1
			if [[ -d $CURRENT_BATCH_SUBDIR/$MAPPING_DIR/$msd ]]; then
				echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping: mapping $msd directory] OK $CURRENT_BATCH_SUBDIR/$MAPPING_DIR/$msd directory already exists. Will write $msd output files in this directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
			else
				mkdir $CURRENT_BATCH_SUBDIR/$MAPPING_DIR/$msd 2>$ERROR_TMP
				rtrn=$?
				mapping_subdir_failed_msg="[Batch mode: mapping] Failed Mapping output sub-directory, $CURRENT_BATCH_SUBDIR/$MAPPING_DIR/$msd was not created."
				exit_on_error "$ERROR_TMP" "$mapping_subdir_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
		 		echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping: mapping $msd directory] OK $CURRENT_BATCH_SUBDIR/$MAPPING_DIR/$msd directory was created successfully. Will write $msd output files in this directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
			fi		
		done

	# current mapping directories
		CURRENT_MAPPING_DIR=$CURRENT_BATCH_SUBDIR/$MAPPING_DIR	
		CURRENT_MAPPING_TMP=$CURRENT_MAPPING_DIR/tmp
		CURRENT_MAPPING_LOG=$CURRENT_MAPPING_DIR/log
	# current fastq files
		if [[ "$VARCO_QC_TRIM_process" == TRUE ]]; then
			if [[ -n $VARCO_SAMPLE_fastq_forward_trimmed_both_surviving 
					&& -n $VARCO_SAMPLE_fastq_reverse_trimmed_both_surviving 
					&& -n $VARCO_SAMPLE_name ]]; then
				CURRENT_FASTQ_FORWARD=$VARCO_SAMPLE_fastq_forward_trimmed_both_surviving
				CURRENT_FASTQ_REVERSE=$VARCO_SAMPLE_fastq_reverse_trimmed_both_surviving
				CURRENT_SAMPLE_NAME=$VARCO_SAMPLE_name
			else
				fastq_files_failed_msg="[Batch mode: mapping] Failed loading sample trimmed fastq files variables."
				echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] echo \$VARCO_SAMPLE_fastq_forward_trimmed_both_surviving failed." >$ERROR_TMP
				echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] echo \$VARCO_SAMPLE_fastq_reverse_trimmed_both_surviving failed." >>$ERROR_TMP
				echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] echo \$VARCO_SAMPLE_fastq_name failed." >>$ERROR_TMP
				exit_on_error "$ERROR_TMP" "$fastq_files_failed_msg" 1 "$LOG_DIR/$LOGFILE"
			fi
		else
			if [[ -n $VARCO_SAMPLE_fastq_forward 
					&& -n $VARCO_SAMPLE_fastq_reverse 
					&& -n $VARCO_SAMPLE_name ]]; then
				CURRENT_FASTQ_FORWARD=$VARCO_SAMPLE_fastq_forward
				CURRENT_FASTQ_REVERSE=$VARCO_SAMPLE_fastq_reverse
				CURRENT_SAMPLE_NAME=$VARCO_SAMPLE_name
			else
				fastq_files_failed_msg="[Batch mode: mapping] Failed loading sample fastq files variables."
				echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] echo \$VARCO_SAMPLE_fastq_forward failed." >$ERROR_TMP
				echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] echo \$VARCO_SAMPLE_fastq_reverse failed." >>$ERROR_TMP
				echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] echo \$VARCO_SAMPLE_fastq_name failed." >>$ERROR_TMP
				exit_on_error "$ERROR_TMP" "$fastq_files_failed_msg" 1 "$LOG_DIR/$LOGFILE"
			fi
		fi
		# current mapping files
		CURRENT_MAPPING_ERROR=$CURRENT_MAPPING_LOG/$CURRENT_SAMPLE_NAME\_mapping_b$b\_s$s\_err.log
		CURRENT_MAPPING_LOGFILE=$CURRENT_MAPPING_LOG/$CURRENT_SAMPLE_NAME\_mapping_b$b\_s$s.log

		echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] Processing the following fastq sample name: $CURRENT_SAMPLE_NAME" | tee -a ${CURRENT_MAPPING_LOGFILE} 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] current sample forward fastq file: $CURRENT_FASTQ_FORWARD" | tee -a ${CURRENT_MAPPING_LOGFILE} 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] current sample reverse fastq file: $CURRENT_FASTQ_REVERSE" | tee -a ${CURRENT_MAPPING_LOGFILE} 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1

	# 1.5.2 Build gsnap mapper command
		command_name="gsnap"	
		echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping] Building $command_name mapping command line options for current batch sample $CURRENT_BATCH_SUBDIR ..." | tee -a $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
		echo -e "--- Config section [${command_name}] ---" | tee -a ${CURRENT_MAPPING_LOGFILE} 2>&1
		for params in $(set | grep ^$(toupper ${NAMESPACE}_${command_name}_)); do
	    	echo -e "params: $params" | tee -a ${CURRENT_MAPPING_LOGFILE} 2>&1
		done
		gsnap_cli_options=($(buildCommandLineOptions $command_name $NAMESPACE 2>$CURRENT_MAPPING_ERROR))
		rtrn=$?
		cli_options_failed_msg="[Batch mode: mapping] Failed An error occured while building the $command_name mapping command line for current batch sample $CURRENT_SAMPLE_NAME."
		exit_on_error "$CURRENT_MAPPING_ERROR" "$cli_options_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
		opts="${gsnap_cli_options[@]}"
		echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping] OK Successfully build mapping command line options for current batch sample $CURRENT_BATCH_SUBDIR." | tee -a $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
		echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping] $command_name options: $opts" | tee -a $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1

		echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping] Building $command_name mapping command line for current batch sample $CURRENT_BATCH_SUBDIR ..." | tee -a $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1	
		gsnap_out=$CURRENT_SAMPLE_NAME\_gsnap_out_b$b\_s$s.sam
		CURRENT_DATA_SAMPLE_DIR=$DATA_ROOT_DIR/$(basename ${subdirs[$si]})
		CURRENT_MAPPING_PID=$CURRENT_MAPPING_TMP/$CURRENT_SAMPLE_NAME.pid	
		gsnap_cli="gsnap $opts $CURRENT_DATA_SAMPLE_DIR/$CURRENT_FASTQ_FORWARD $CURRENT_DATA_SAMPLE_DIR/$CURRENT_FASTQ_REVERSE >$CURRENT_MAPPING_DIR/$gsnap_out 2>${CURRENT_MAPPING_ERROR} &"

		# append time monitoring support
		timeCmd=$(buildTimeCmd "$command_name" "$CURRENT_MAPPING_LOG" 2>${ERROR_TMP})
		rtrn=$?
		time_cmd_failed_msg="[Batch mode: mapping] Failed building time command for use with $command_name."		
		exit_on_error "$ERROR_TMP" "$time_cmd_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
		gsnap_cli="$timeCmd $gsnap_cli"

	# 1.5.3 Run the command	
		echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping] Executing $command_name mapping command line for current batch sample $CURRENT_BATCH_SUBDIR ..." | tee -a $CURRENT_MAPPING_LOGFILE 2>&1  2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] $command_name mapper version: \n$MAPPER_VERSION" | tee -a $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
		eval "$gsnap_cli" 2>$ERROR_TMP
		pid=$!
		rtrn=$?
		eval_failed_msg="[Batch mode: mapping] Failed eval gsnap cli."
		exit_on_error "$ERROR_TMP" "$eval_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] $CURRENT_SAMPLE_NAME pid: $pid" | tee -a $CURRENT_MAPPING_LOGFILE 2>&1  2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
		echo -e $pid >$CURRENT_MAPPING_PID 2>$ERROR_TMP
		PIDS_ARR=("${PIDS_ARR[@]}" "$pid")

	# Last batch sample: have a break!
		[[ $last == "TRUE" ]] && break
	done

	# wait for all gsnap processes batch to finish
	echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] gsnap pids count: ${#PIDS_ARR[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] gsnap pids list: ${PIDS_ARR[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1
	pid_list_failed_msg="[Batch mode: mapping] Failed listing process $p."	
	for p in "${PIDS_ARR[@]}"; do
		#echo -e $(ps aux | grep $p | grep $USER | grep -v grep 2>${ERROR_TMP})
		echo -e $(ps aux | grep $USER | gawk -v pid=$p '$2 ~ pid {print $0}' 2>${ERROR_TMP}) 
		rtrn=$?
		exit_on_error "$ERROR_TMP" "$pid_list_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
	done

	echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] Waiting for $command_name processes to exit until $WAITALL_TIMEOUT_HR timeout ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
	#waitall "${PIDS_ARR[@]}" 2>${ERROR_TMP}
	waitalluntiltimeout "${PIDS_ARR[@]}" 2>${ERROR_TMP}
	egrep "exited" ${ERROR_TMP} 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] All $command_name processes exited." | tee -a $LOG_DIR/$LOGFILE 2>&1

	# check for errors
	errs=0
	last=FALSE
	for s in $(seq 1 $VARCO_SPLIT_MAP_batch_size); do
		si=$[$s-1]
		sdi=$[$si+($b-1)*$VARCO_SPLIT_MAP_batch_size]
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] Checking for mapping errors for batch sample $sdi, $si: ${subdirs[$si]} ... " | tee -a $LOG_DIR/$LOGFILE 2>&1

	# last batch sample
		[[ "$s" -eq "${#subdirs[@]}" ]] && last=TRUE

	# current sample variables
		CURRENT_BATCH_SUBDIR=$JOB_TAG/$(basename ${subdirs[$si]})
		CURRENT_MAPPING_DIR=$CURRENT_BATCH_SUBDIR/$MAPPING_DIR
		CURRENT_MAPPING_LOG=$CURRENT_MAPPING_DIR/log
		CURRENT_MAPPING_TMP=$CURRENT_MAPPING_DIR/tmp
		CURRENT_SAMPLE_CONFIG=$CURRENT_BATCH_SUBDIR/$SAMPLE_CONFIG
		CURRENT_SAMPLE_CONFIG_ERR=$CURRENT_MAPPING_LOG/sample_loading_err.log	

	# load sample config
		echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping] Loading sample config for $CURRENT_SAMPLE_CONFIG file ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
		sample_config_failed_msg="[Batch mode: mapping] Failed loading sample config file, $CURRENT_SAMPLE_CONFIG, from $CURRENT_BATCH_SUBDIR"
		for cfg in $(get_config_sections $CURRENT_SAMPLE_CONFIG 2>$CURRENT_SAMPLE_CONFIG_ERR; rtrn=$?); do
			if [[ "$rtrn" -eq 0 ]]; then
				echo -e "--- Config section [${cfg}] ---"
				unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${NAMESPACE}" 'BEGIN { 
	  				cfg = toupper(cfg);
	  				prefix = toupper(prefix);
					}
					/^prefix_cfg_/  { print $1 }' 2>$CURRENT_SAMPLE_CONFIG_ERR) $(toupper ${NAMESPACE}_${cfg}_) 2>>$CURRENT_SAMPLE_CONFIG_ERR
				rtrn=$?
				exit_on_error "$CURRENT_SAMPLE_CONFIG_ERR" "$sample_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
				set_config_params $CURRENT_SAMPLE_CONFIG ${cfg} ${NAMESPACE} 2>$CURRENT_SAMPLE_CONFIG_ERR
				rtrn=$?
				exit_on_error "$CURRENT_SAMPLE_CONFIG_ERR" "$sample_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE" 
				for params in $(set | grep ^$(toupper ${NAMESPACE}_${cfg}_) 2>$CURRENT_SAMPLE_CONFIG_ERR); do
					echo -e "$params"
				done
			else
				exit_on_error "$CURRENT_SAMPLE_CONFIG_ERR" "$sample_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
			fi
		done

	# current sample variables
		CURRENT_SAMPLE_NAME=$VARCO_SAMPLE_name
		CURRENT_MAPPING_ERROR=$CURRENT_MAPPING_LOG/$CURRENT_SAMPLE_NAME\_mapping_b$b\_s$s\_err.log
		CURRENT_MAPPING_LOGFILE=$CURRENT_MAPPING_LOG/$CURRENT_SAMPLE_NAME\_mapping_b$b\_s$s.log
		CURRENT_MAPPING_PID=$CURRENT_MAPPING_TMP/$CURRENT_SAMPLE_NAME.pid

	# copy time monitoring on mapping command to logs
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] Time monitoring on mapping command, $command_name:" | tee -a  $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
		cat $CURRENT_MAPPING_LOG/${command_name}_time.log 2>&1 | tee -a  $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1

	# check for errors
		pid_status=$(egrep "exited" ${ERROR_TMP} 2>&1 | egrep $(cat $CURRENT_MAPPING_PID) 2>&1)

		if [[ -z $(echo -e $pid_status 2>&1 | egrep "zero exit status") ]]; then
			echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] $CURRENT_SAMPLE_NAME sample mapping process exited with non zero exit status:\n$pid_status" | tee -a  $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
			let errs=errs+1
			break
		else
			echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] pid status: $pid_status" 2>&1 | tee -a $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
		fi

		if [[ -z $(tail -n 1 "${CURRENT_MAPPING_ERROR}" | egrep -v "^$" | egrep "^Processed" 2>&1) ]]; then
			echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] $CURRENT_SAMPLE_NAME sample error output:" 2>&1 | tee -a $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1; 
			cat "${CURRENT_MAPPING_ERROR}" 2>&1 | tee -a  $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
			let errs=errs+1
			break
		fi
	# Last batch sample: have a break!
		[[ $last == "TRUE" ]] && break   
	done

	if [[ $errs == 0 ]]; then
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] all gsnap processes for batch #$b finished without errors." | tee -a  $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	else
		mapping_err_failed_msg="[Batch mode: mapping] some errors occured while gsnap processing samples for batch #$b."
		echo -e "$(date '+%Y_%m_%d %T') $mapping_err_failed_msg" | tee -a  $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] refer to $CURRENT_MAPPING_ERROR file for details." | tee -a  $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
		exit_on_error "$CURRENT_MAPPING_ERROR" "$mapping_err_failed_msg" 1 "$LOG_DIR/$LOGFILE"
	fi

	# 1.6 Conversion: sam to sorted bam
	echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] Processing sam conversion to sorted indexed bam ... " |  tee -a $LOG_DIR/$LOGFILE 2>&1

	# for each conversion command
	view_command="samtools view"
	view_command_ext="bam"
	sort_command="samtools sort"
	sort_command_ext="sorted.$view_command_ext"
	index_command="samtools index"
	index_command_ext="${sort_command_ext}.bai"
	conversion_commands=("$view_command" "$sort_command" "$index_command")

	for cmd in "${conversion_commands[@]}"; do
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] Processing with $cmd ... " |  tee -a $LOG_DIR/$LOGFILE 2>&1

	# reinitiate last batch sample
		last=FALSE

	# reinitiate pids array
		PIDS_ARR=()

	## for each sample
		for s in $(seq 1 $VARCO_SPLIT_MAP_batch_size); do
			si=$[$s-1]
			sdi=$[$si+($b-1)*$VARCO_SPLIT_MAP_batch_size]
			echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] Processing batch sample $sdi, $si: ${subdirs[$si]}, for sam conversion to sorted indexed bam ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
			CURRENT_BATCH_SUBDIR=$JOB_TAG/$(basename "${subdirs[$si]}")
			CURRENT_SAMPLE_CONFIG=$CURRENT_BATCH_SUBDIR/$SAMPLE_CONFIG

	# last batch sample
			[[ "$s" -eq "${#subdirs[@]}" ]] && last=TRUE

	# load sample config
			echo "$(date '+%Y_%m_%d %T') [Batch mode: conversion] Loading sample config for $CURRENT_SAMPLE_CONFIG file ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
			sample_config_failed_msg="[Batch mode: conversion] Failed loading sample config file, $CURRENT_SAMPLE_CONFIG, from $CURRENT_BATCH_SUBDIR"
			for cfg in $(get_config_sections $CURRENT_SAMPLE_CONFIG 2>$ERROR_TMP; rtrn=$?); do
				if [[ "$rtrn" -eq 0 ]]; then
					echo -e "--- Config section [${cfg}] ---"
					unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${NAMESPACE}" 'BEGIN { 
		  				cfg = toupper(cfg);
		  				prefix = toupper(prefix);
						}
						/^prefix_cfg_/  { print $1 }' 2>$ERROR_TMP) $(toupper ${NAMESPACE}_${cfg}_) 2>>$ERROR_TMP
					rtrn=$?
					exit_on_error "$ERROR_TMP" "$sample_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
					set_config_params $CURRENT_SAMPLE_CONFIG ${cfg} ${NAMESPACE} 2>$ERROR_TMP
					rtrn=$?
					exit_on_error "$ERROR_TMP" "$sample_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE" 
					for params in $(set | grep ^$(toupper ${NAMESPACE}_${cfg}_) 2>$ERROR_TMP); do
						echo -e "$params"
					done
				else
					exit_on_error "$ERROR_TMP" "$sample_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
				fi
			done

	# current mapping variables
			CURRENT_MAPPING_DIR=$CURRENT_BATCH_SUBDIR/$MAPPING_DIR
			CURRENT_MAPPING_LOG=$CURRENT_MAPPING_DIR/log
			CURRENT_MAPPING_TMP=$CURRENT_MAPPING_DIR/tmp
			CURRENT_SAMPLE_NAME=$VARCO_SAMPLE_name
			CURRENT_MAPPING_ERROR=$CURRENT_MAPPING_LOG/$CURRENT_SAMPLE_NAME\_mapping_b$b\_s$s\_err.log
			CURRENT_MAPPING_LOGFILE=$CURRENT_MAPPING_LOG/$CURRENT_SAMPLE_NAME\_mapping_b$b\_s$s.log
			CURRENT_MAPPING_PID=$CURRENT_MAPPING_TMP/$CURRENT_SAMPLE_NAME.pid			
			CURRENT_MAPPING_SAM=$CURRENT_SAMPLE_NAME\_gsnap_out_b$b\_s$s.sam
			CURRENT_MAPPING_SAM_BASE=$(basename ${CURRENT_MAPPING_SAM%.*})

	## current conversion variables 
			CURRENT_CONVERSION_PID=$CURRENT_MAPPING_TMP/${CURRENT_MAPPING_SAM_BASE}.pid
			CURRENT_CONVERSION_ERROR=$CURRENT_MAPPING_LOG/$CURRENT_MAPPING_SAM_BASE\_conversion_b$b\_s$s\_err.log
			CURRENT_CONVERSION_LOGFILE=$CURRENT_MAPPING_LOG/$CURRENT_MAPPING_SAM_BASE\_conversion_b$b\_s$s.log

	## build conversion cli
			echo -e "--- Config section [${cmd// /_}] ---" | tee -a $CURRENT_CONVERSION_LOGFILE 2>&1| tee -a $LOG_DIR/$LOGFILE 2>&1
			for params in $(set | grep ^$(toupper ${NAMESPACE}_${cmd// /_}_)); do
	    		echo -e "params: $params" | tee -a $CURRENT_CONVERSION_LOGFILE 2>&1| tee -a $LOG_DIR/$LOGFILE 2>&1
			done
			cli_options=($(buildCommandLineOptions "$cmd" $NAMESPACE 2>${ERROR_TMP}))
			rtrn=$?
			build_opts_failed_msg="[Batch mode: conversion] Failed building $cmd command line options."
			exit_on_error "$ERROR_TMP" "$build_opts_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
			opts="${cli_options[@]}"

			echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] $cmd cli options: $opts" | tee -a $CURRENT_CONVERSION_LOGFILE 2>&1| tee -a $LOG_DIR/$LOGFILE 2>&1
			case $cmd in
				"$view_command")
					if [[ -s $CURRENT_MAPPING_DIR/$CURRENT_MAPPING_SAM ]]; then
						echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] gsnap sam output file: $CURRENT_MAPPING_DIR/$CURRENT_MAPPING_SAM" | tee -a $CURRENT_CONVERSION_LOGFILE 2>&1| tee -a $LOG_DIR/$LOGFILE 2>&1		
						cmd_out=$CURRENT_MAPPING_SAM_BASE\.$view_command_ext
						echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] input: $CURRENT_MAPPING_SAM" | tee -a $CURRENT_CONVERSION_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
						echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] output: $cmd_out" | tee -a $CURRENT_CONVERSION_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
						cmd_cli="$cmd $opts $CURRENT_MAPPING_DIR/$CURRENT_MAPPING_SAM > $CURRENT_MAPPING_DIR/$cmd_out 2>>$CURRENT_CONVERSION_ERROR &"
					else
						sam_out_file_failed_msg="$CURRENT_MAPPING_DIR/$CURRENT_MAPPING_SAM file does not exist or is empty." 
						debug $sam_out_file_failed_msg | tee -a $CURRENT_CONVERSION_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1;
						exit_on_error "$CURRENT_CONVERSION_LOGFILE" "$sam_out_file_failed_msg" 3 "$LOG_DIR/$LOGFILE"
					fi
				;;
				"$sort_command")
					bam_out=$CURRENT_MAPPING_SAM_BASE\.$view_command_ext
					if [[ -s $CURRENT_MAPPING_DIR/$bam_out ]]; then
						cmd_out=$CURRENT_MAPPING_SAM_BASE\.sorted
						echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] input: $bam_out" | tee -a $CURRENT_CONVERSION_LOGFILE 2>&1| tee -a $LOG_DIR/$LOGFILE 2>&1
						echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] output: $cmd_out" | tee -a $CURRENT_CONVERSION_LOGFILE 2>&1| tee -a $LOG_DIR/$LOGFILE 2>&1
						cmd_cli="$cmd $opts $CURRENT_MAPPING_DIR/$bam_out $CURRENT_MAPPING_DIR/$cmd_out 2>>$CURRENT_CONVERSION_ERROR &"
					else
						bam_out_file_failed_msg="$CURRENT_MAPPING_DIR/$bam_out file does not exist or is empty." 
						debug $bam_out_file_failed_msg | tee -a $CURRENT_CONVERSION_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1;
						exit_on_error "$CURRENT_CONVERSION_LOGFILE" "$bam_out_file_failed_msg" 3 "$LOG_DIR/$LOGFILE"
					fi
				;;	
				"$index_command")
					sorted_bam_out=${CURRENT_MAPPING_SAM_BASE}.$sort_command_ext
					if [[ -s $CURRENT_MAPPING_DIR/$sorted_bam_out ]]; then
						cmd_out=$CURRENT_MAPPING_SAM_BASE\.$index_command_ext
						echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] input: $sorted_bam_out" | tee -a $CURRENT_CONVERSION_LOGFILE 2>&1| tee -a $LOG_DIR/$LOGFILE 2>&1
						echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] output: $cmd_out" | tee -a $CURRENT_CONVERSION_LOGFILE 2>&1| tee -a $LOG_DIR/$LOGFILE 2>&1
						cmd_cli="$cmd $opts $CURRENT_MAPPING_DIR/$sorted_bam_out 2>>${CURRENT_CONVERSION_ERROR} &"
					else
						bam_sorted_out_file_failed_msg="$CURRENT_MAPPING_DIR/$sorted_bam_out file does not exist or is empty." 
						debug $bam_sorted_out_file_failed_msg | tee -a $CURRENT_CONVERSION_LOGFILE 2>&1| tee -a $LOG_DIR/$LOGFILE 2>&1;
						exit_on_error "$CURRENT_CONVERSION_LOGFILE" "$bam_sorted_out_file_failed_msg" 3 "$LOG_DIR/$LOGFILE"
					fi
				;;
				*)
					debug "unexpected $cmd, it will not be processed." | tee -a ${CURRENT_CONVERSION_ERROR} 2>&1 | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
			esac
	## run cli
			if [[ -n $cmd_cli ]]; then
	### append time monitoring support
				timeCmd=$(buildTimeCmd "${cmd// /_}" "$CURRENT_MAPPING_LOG" 2>${ERROR_TMP})
				rtrn=$?
				time_cmd_failed_msg="[Batch mode: conversion] Failed building time command for use with $cmd."		
				exit_on_error "$ERROR_TMP" "$time_cmd_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
				cmd_cli="$timeCmd $cmd_cli"

				echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] Executing command: ${cmd_cli}" | tee -a $CURRENT_CONVERSION_LOGFILE 2>&1| tee -a $LOG_DIR/$LOGFILE 2>&1
				eval "$cmd_cli" 2>${ERROR_TMP}
				pid=$!
				rtrn=$?
				eval_conversion_failed_msg="[Batch mode] Failed eval $cmd_cli."
				exit_on_error "$ERROR_TMP" "$eval_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
				echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] $CURRENT_MAPPING_SAM_BASE pid: $pid" | tee -a $CURRENT_CONVERSION_LOGFILE 2>&1| tee -a $LOG_DIR/$LOGFILE 2>&1
				echo -e $pid >$CURRENT_CONVERSION_PID 2>${ERROR_TMP}
				PIDS_ARR=("${PIDS_ARR[@]}" "$pid")
			else
				cmd_cli_failed_msg="command name, $cmd, was not processed. cmd cli is null." 
				debug $cmd_cli_failed_msg | tee -a ${CURRENT_CONVERSION_ERROR} 2>&1 | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1;
				exit_on_error "CURRENT_CONVERSION_ERROR" "$cmd_cli_failed_msg" 1 "$LOG_DIR/$LOGFILE"
			fi

	# Last batch sample: have a break!
			[[ $last == "TRUE" ]] && break
		done

	# wait for all conversion cmd processes to finish then run the next cmd
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] $cmd pids count: ${#PIDS_ARR[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1
		echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] $cmd pids list: ${PIDS_ARR[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1		
		for p in "${PIDS_ARR[@]}"; do
			#echo -e $(ps aux | grep $p | grep $USER | grep -v grep 2>${ERROR_TMP})
			echo -e $(ps aux | grep $USER | gawk -v pid=$p '$2 ~ pid {print $0}' 2>${ERROR_TMP})
			rtrn=$?
			pid_list_failed_msg="[Batch mode: conversion] Failed listing $cmd process using its pid."
			exit_on_error "$ERROR_TMP" "$pid_list_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
		done

		echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] Waiting for $cmd processes to exit until $WAITALL_TIMEOUT_HR timeout ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
		#waitall "${PIDS_ARR[@]}" 2>${ERROR_TMP}
		waitalluntiltimeout "${PIDS_ARR[@]}" 2>${ERROR_TMP}
		egrep "exited" ${ERROR_TMP} 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] All $cmd processes exited." | tee -a $LOG_DIR/$LOGFILE 2>&1

	# check for errors
		errs=0
		last=FALSE
		for s in $(seq 1 $VARCO_SPLIT_MAP_batch_size); do
			si=$[$s-1]
			sdi=$[$si+($b-1)*$VARCO_SPLIT_MAP_batch_size]
			echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] Checking for $cmd conversion errors for batch sample $sdi, $si: ${subdirs[$si]} ... " | tee -a $LOG_DIR/$LOGFILE 2>&1

	# last batch sample
			[[ "$s" -eq "${#subdirs[@]}" ]] && last=TRUE

	# current sample variables
			CURRENT_BATCH_SUBDIR=$JOB_TAG/$(basename ${subdirs[$si]})
			CURRENT_SAMPLE_CONFIG=$CURRENT_BATCH_SUBDIR/$SAMPLE_CONFIG

	# current mapping variables
			CURRENT_MAPPING_DIR=$CURRENT_BATCH_SUBDIR/$MAPPING_DIR
			CURRENT_MAPPING_LOG=$CURRENT_MAPPING_DIR/log
			CURRENT_MAPPING_TMP=$CURRENT_MAPPING_DIR/tmp

	# load sample config
			echo "$(date '+%Y_%m_%d %T') [Batch mode: conversion] Loading sample config: $CURRENT_SAMPLE_CONFIG file ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
			sample_config_failed_msg="[Batch mode: conversion] Failed loading sample config file, $CURRENT_SAMPLE_CONFIG, from $CURRENT_BATCH_SUBDIR"
			for cfg in $(get_config_sections $CURRENT_SAMPLE_CONFIG 2>$CURRENT_SAMPLE_CONFIG_ERR; rtrn=$?); do
				if [[ "$rtrn" -eq 0 ]]; then
					echo -e "--- Config section [${cfg}] ---"
					unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${NAMESPACE}" 'BEGIN { 
		  				cfg = toupper(cfg);
		  				prefix = toupper(prefix);
						}
						/^prefix_cfg_/  { print $1 }' 2>$CURRENT_SAMPLE_CONFIG_ERR) $(toupper ${NAMESPACE}_${cfg}_) 2>>$CURRENT_SAMPLE_CONFIG_ERR
					rtrn=$?
					exit_on_error "$CURRENT_SAMPLE_CONFIG_ERR" "$sample_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
					set_config_params $CURRENT_SAMPLE_CONFIG ${cfg} ${NAMESPACE} 2>$CURRENT_SAMPLE_CONFIG_ERR
					rtrn=$?
					exit_on_error "$CURRENT_SAMPLE_CONFIG_ERR" "$sample_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE" 
					for params in $(set | grep ^$(toupper ${NAMESPACE}_${cfg}_) 2>$CURRENT_SAMPLE_CONFIG_ERR); do
						echo -e "$params"
					done
				else
					exit_on_error "$CURRENT_SAMPLE_CONFIG_ERR" "$sample_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
				fi
			done

	# current mapping variables
			CURRENT_SAMPLE_NAME=$VARCO_SAMPLE_name
			CURRENT_MAPPING_ERROR=$CURRENT_MAPPING_LOG/$CURRENT_SAMPLE_NAME\_mapping_b$b\_s$s\_err.log
			CURRENT_MAPPING_LOGFILE=$CURRENT_MAPPING_LOG/$CURRENT_SAMPLE_NAME\_mapping_b$b\_s$s.log
			CURRENT_MAPPING_PID=$CURRENT_MAPPING_TMP/$CURRENT_SAMPLE_NAME.pid			
			CURRENT_MAPPING_SAM=$CURRENT_SAMPLE_NAME\_gsnap_out_b$b\_s$s.sam
			CURRENT_MAPPING_SAM_BASE=$(basename ${CURRENT_MAPPING_SAM%.*})
			CURRENT_SAMPLE_CONFIG_ERR=$CURRENT_MAPPING_LOG/sample_loading_err.log

	## current conversion variables 
			CURRENT_CONVERSION_PID=$CURRENT_MAPPING_TMP/${CURRENT_MAPPING_SAM_BASE}.pid
			CURRENT_CONVERSION_ERROR=$CURRENT_MAPPING_LOG/$CURRENT_MAPPING_SAM_BASE\_conversion_b$b\_s$s\_err.log
			CURRENT_CONVERSION_LOGFILE=$CURRENT_MAPPING_LOG/$CURRENT_MAPPING_SAM_BASE\_conversion_b$b\_s$s.log

	# copy time monitoring on mapping command to logs
			echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] Time monitoring on conversion command, $cmd:" | tee -a  $CURRENT_CONVERSION_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
			cat $CURRENT_MAPPING_LOG/${cmd// /_}_time.log 2>&1 | tee -a  $CURRENT_CONVERSION_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1

	# check for errors
			pid_status=$(egrep "exited" ${ERROR_TMP} 2>&1 | egrep $(cat $CURRENT_CONVERSION_PID) 2>&1)

			if [[ -z $(echo -e $pid_status 2>&1 | egrep "zero exit status") ]]; then
				echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] $CURRENT_SAMPLE_NAME sample, $cmd conversion process exited with non zero exit status:\n$pid_status" | tee -a  $CURRENT_CONVERSION_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
				let errs=errs+1
				break
			else
				echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] pid status: $pid_status" 2>&1 | tee -a $CURRENT_CONVERSION_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
			fi

			if [[ -z $(tail -n 1 "${CURRENT_CONVERSION_ERROR}" | egrep -v "^$" | egrep "^\[samopen\]" 2>&1) ]]; then
				echo -e "$(date '+%Y_%m_%d %T') [Batch mode: conversion] $CURRENT_SAMPLE_NAME sample error output:" 2>&1 | tee -a $CURRENT_CONVERSION_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1; 
				cat "${CURRENT_CONVERSION_ERROR}" 2>&1 | tee -a  $CURRENT_CONVERSION_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
				let errs=errs+1
				break
			fi
		# Last batch sample: have a break!
			[[ $last == "TRUE" ]] && break   
		done

		if [[ $errs == 0 ]]; then		
			echo -e "$(date '+%Y_%m_%d %T') [Batch mode] all $cmd conversion processes for batch #$b finished without errors." | tee -a  $CURRENT_CONVERSION_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
		else
			conversion_err_failed_msg="[Batch mode] some errors occured while $cmd processing samples for batch #$b."
			echo -e "$(date '+%Y_%m_%d %T') $conversion_err_failed_msg" | tee -a  $CURRENT_CONVERSION_LOGFILE 2>&1 | tee -a  $CURRENT_CONVERSION_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
			echo -e "$(date '+%Y_%m_%d %T') [Batch mode] refer to $CURRENT_CONVERSION_ERROR file for details." | tee -a  $CURRENT_CONVERSION_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
			exit_on_error "$CURRENT_CONVERSION_LOGFILE" "$conversion_err_failed_msg" 1 "$LOG_DIR/$LOGFILE"
		fi
	done

	# find and remove all sam/unsorted unindexed bam output files in $JOB_TAG directory
	if [[ $VARCO_SPLIT_MAP_clean == "TRUE" ]]; then
		echo -e "$(date '+%Y_%m_%d %T') [Cleaning] Removing all sam output files from $WORKING_DIR/$JOB_TAG directory ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
		#find $WORKING_DIR/$JOB_TAG -type f -name "*.sam" -exec rm -f {} \;
		sam_files=($(find $WORKING_DIR/$JOB_TAG -type f -name "*.sam" 2>${ERROR_TMP}))
		rtrn=$?	
		sam_listing_failed_msg="[Cleaning] Failed listing sam files in $WORKING_DIR/$JOB_TAG."
		exit_on_error "$ERROR_TMP" "$sam_listing_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
		# remove sam files		
		for s in "${sam_files[@]}"; do
			echo -ne "$(date '+%Y_%m_%d %T') [Cleaning] removing $s file ... "
			rm -f $s 2>${ERROR_TMP}
			rtrn=$?
			sam_rm_failed_msg="[Cleaning] Failed removing $s file."
			exit_on_error "$ERROR_TMP" "$sam_rm_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
			echo -e "done"
		done
		echo -e "$(date '+%Y_%m_%d %T') [Cleaning] Removing all unsorted and unindexed bam output files from $WORKING_DIR/$JOB_TAG directory ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
		# remove unsorted and unindexed bam files
		for s in "${sam_files[@]}"; do
			uib=${s%.*}.bam
			bam_rm_failed_msg="[Cleaning] Failed removing $uib file."
			echo -ne "$(date '+%Y_%m_%d %T') [Cleaning] removing $uib file ... "
			m=$(echo $uib | egrep -v "sorted.bam$" | egrep -v "bai$" 2>${ERROR_TMP})
			rtrn=$?
			exit_on_error "$ERROR_TMP" "$bam_rm_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
			[[ -n $m ]] && rm -f "$uib" 2>${ERROR_TMP}
			rtrn=$?
			exit_on_error "$ERROR_TMP" "$bam_rm_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
			echo -e "done"
		done
		echo -e "$(date '+%Y_%m_%d %T') [Cleaning] All sam, and unsorted and unindexed bam output files was deleted from $WORKING_DIR/$JOB_TAG directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
	fi

    # Unshifting the current batch samples from fastq subdirs array
    echo -e "$(date '+%Y_%m_%d %T') [Batch mode] remaining subdirs count before unshifting: ${#subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1
    #echo -e "batch_size - 1: $[$VARCO_SPLIT_MAP_batch_size -1]" >> $LOG_DIR/$LOGFILE 2>&1 # for testing purpose
    for i in $(seq 0 $[$VARCO_SPLIT_MAP_batch_size -1]); do
		# echo -ne $i >> $LOG_DIR/$LOGFILE 2>&1 # for testing purpose
		if [[ "${#subdirs[@]}" -ge "$i" ]]; then # weird, this condition does not get the last array element to be unshifted
	    	unset subdirs[$i]
		fi
		[[ -n "${subdirs[$i]}" ]] && unset subdirs[$i] # don't know why but when get to the last array element, need to force unshifting
		#echo -e "${subdirs[$i]}" >> $LOG_DIR/$LOGFILE 2>&1 # for testing purpose
    done
    subdirs=("${subdirs[@]}")
    echo -e "$(date '+%Y_%m_%d %T') [Batch mode] remaining subdirs count after unshifting: ${#subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo -e "$(date '+%Y_%m_%d %T') [Batch mode] remaining subdirs list after unshifting: ${subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1

    # Wait until the last current batch sample finish before launching the next batch
	echo -e "$(date '+%Y_%m_%d %T') [Batch mode] Processing batch samples #$b has finished without errors." | tee -a $LOG_DIR/$LOGFILE 2>&1
	[[ "$batches" -gt 1 ]] && echo -e "$(date '+%Y_%m_%d %T') [Batch mode] Will proceed with next batch samples." | tee -a $LOG_DIR/$LOGFILE 2>&1

	# send an email
	batch_completed_msg="Batch #$b completed successfully."
	[[ -n $RECIPIENT ]] && sendEmail $RECIPIENT "[$(basename ${0%.*})] $JOB_TAG job, batch #$b/$batches completed " "$batch_completed_msg"
done


#
# Clean:
#

# unset environment variables with used namespace
echo -ne "$(date '+%Y_%m_%d %T') [Cleaning] Unsetting all environment variables using namespace: $NAMESPACE ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${NAMESPACE}" 'BEGIN { 
          cfg = toupper(cfg);
          prefix = toupper(prefix);
       }
       /^prefix_cfg_/  { print $1 }') $(toupper ${prefix}_${cfg}_)
echo -e "done" | tee -a $LOG_DIR/$LOGFILE 2>&1

#=====
# END
#=====
echo "$(date '+%Y_%m_%d %T') [$(basename $0)] Executed command: $EXECUTED_COMMAND" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -n "$(date '+%Y_%m_%d %T') [$(basename $0)] Elapsed time: " | tee -a $LOG_DIR/$LOGFILE 2>&1
echo |awk -v time="$SECONDS" '{print strftime("%Hh:%Mm:%Ss", time, 1)}' | tee -a $LOG_DIR/$LOGFILE 2>&1
echo "$(date '+%Y_%m_%d %T') [$(basename $0)] Exits the pipeline." | tee -a $LOG_DIR/$LOGFILE 2>&1
echo "$(date '+%Y_%m_%d %T') [$(basename $0)] More information about this job can be found in $LOG_DIR/$LOGFILE" | tee -a $LOG_DIR/$LOGFILE 2>&1
# send an email
job_completed_msg="$JOB_TAG job completed successfully."
[[ -n $RECIPIENT ]] && sendEmail $RECIPIENT "[$(basename ${0%.*})] $JOB_TAG job completed" "$job_completed_msg"

#exit 0


