#! /bin/bash

# TODO: add license copyright 2012

#
# VARCO_SPLIT_MAP
#
# Author: Joseph Tran <Joseph.Tran@versailles.inra.fr>
#
# date: 2012-12-18
#
# Description: This script performs reads mapping in batch mode by splitting given samples in several batches
#              and run sequentially each batch to avoid cpu overload and running out of disk space 
#

VERSION=dev

########################
# SECTION CONFIGURATION
#######################

# Include lib functions

PROD_PREFIX="/usr/local"
DEV_PREFIX="$(pwd)/.."
PREFIX=$DEV_PREFIX # TO BE CHANGED WHEN SWITCHING TO PROD
. $PREFIX/share/varco_split_map/lib/varco_split_map_lib.inc

# Set variables

ARGS=2
DATA_ROOT_DIR=$1
JOB_TAG=$2

DATE=$(date '+%Y_%m_%d_%T')
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

PIDS_ARR=()

LOG_DIR=$JOB_TAG/"log"
QC_TRIM_DIR="QC_TRIM"
MAPPING_DIR="MAPPING"

ERROR_TMP="/tmp/$(basename ${0%.*})_error_${USER}_$DATE.log"

PREREQUISITES_MSG="You need to copy the user configuration file to your working directory
                and configure the options before running this script.
                To copy the user configuration template file, issue this command in your terminal:
                \$ cp ${VARCO_SPLIT_MAP_USER_CONFIG} .
                Then use any editor (emacs, gedit, nano, vi, etc.) to open it and configure the options.
                Save your changes, and you are ready to run the script."

MAPPER_VERSION=$(gsnap --version 2>&1 | awk -F"\n" 'BEGIN{info=""}; {info=(info "\n\t" $1)}; END{printf info}')

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
		exit $status
	fi
}

#==============================================
# TEST if enough args else print usage message
#==============================================
[[ $# -ne "$ARGS" ]] && { printf %s "\
Program: $(basename $0)
Version: $VERSION
Author: Joseph Tran, IJPB Bioinformatics Dev Team
Contact: Joseph.Tran@versailles.inra.fr

TODO: add license copyright 2012

Usage: $(basename $0) samples_root_dir job_tag

Arguments: samples_root_dir Path to the parent folder containing the reads samples subdirectories 
           job_tag          <String> Prefix to attach to any output files (without space)

Description: This script performs reads mapping in batch mode by splitting
             given samples in several batches and run sequentially each batch 
             to avoid cpu overload and running out of disk space

Pre-requisites: ${PREREQUISITES_MSG}

User Configuration File: here is the main configuration sections and their corresponding parameters
  [split_map] section
    - check_cpu_overload (default=TRUE)
                         This option controls if cpu overload checking has to be performed.
                         If TRUE, this script will run only if cpu average load did not exceed 50% 
                         in the last 1, 5 and 15 minutes. This test is performed before running the script,
                         and for each batch. If cpu average load exceeds 50%, the script will wait for 1 minute
                         before testing again cpu average load. TODO: put some timeout else server will crash
    - batch_size (default=4)
                         This option controls the number of samples for each batch. 
                         Its value is determined dynamically if the user value exceeds the max batch size. 
                         The max batch size equals the ratio max number of processors allowed 
                         divided by the number of threads (mapper option).
                         The max number of processors allowed is fairly set by default 
                         to 50% of max number of processors available.
                         If batch_size lower or equal to max batch size, then use batch_size value.
                         Else, batch_size equals half max batch size.
    - clean (default=FALSE)
                         This option allows to clean each sample output directory by removing sam and bam files,
                         leaving only sorted and indexed bam files.

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

  [samtools_view]
    Please refer to samtools vew documentation: samtools view --help
    - b (default=TRUE)
    - S (default=TRUE)

  [samtools_sort]
    Please refer to samtools vew documentation: samtools sort --help
    No current options are available.

  [samtools_index]
    Please refer to samtools vew documentation: samtools index --help
    No current options are available.

Notes: 1. Current mapper is gsnap 
          $(gsnap --version 2>&1 | awk -F"\n" 'BEGIN{info=""}; {info=(info "\n\t" $1)}; END{printf info}')
       2. Current sam toolkit is samtools
          $(samtools 2>&1 | egrep 'Program|Version'| awk -F"\n" 'BEGIN{info=""}; {info=(info "\n\t" $1)}; END{printf info}')

Joseph Tran, IJPB Bioinformatics Development Team
Contact: Joseph.Tran@versailles.inra.fr

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
    if [[ $? -ne 0 ]]; then
	echo "$(date '+%Y_%m_%d %T') [Job directory] Failed Job directory, $JOB_TAG, was not created." | tee -a $ERROR_TMP 2>&1
	echo "$(date '+%Y_%m_%d %T') [Pipeline error] Exits the pipeline, with error code 126." | tee -a $ERROR_TMP 2>&1
	echo "$(date '+%Y_%m_%d %T') [Pipeline error] More details can be found in $ERROR_TMP." 2>&1
	exit 126
    else
	echo "$(date '+%Y_%m_%d %T') [Job directory] OK $JOB_TAG directory was created successfully. Will output all job files in this directory." | tee -a $ERROR_TMP 2>&1
    fi
fi

# Create log directory
echo "$(date '+%Y_%m_%d %T') [Log directory] Creating $LOG_DIR directory ..." | tee -a $ERROR_TMP 2>&1
if [[ -d $LOG_DIR ]]; then
    [[ -s $ERROR_TMP ]] && cat $ERROR_TMP > $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %T') [Log directory] OK $LOG_DIR directory already exists. Will write log files in this directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
else
    mkdir $LOG_DIR 2>>$ERROR_TMP
    if [[ $? -ne 0 ]]; then
	echo "$(date '+%Y_%m_%d %T') [Log directory] Failed Log directory, $LOG_DIR, was not created." | tee -a $ERROR_TMP 2>&1
	echo "$(date '+%Y_%m_%d %T') [Pipeline error] Exits the pipeline, with error code 126." | tee -a $ERROR_TMP 2>&1
	echo "$(date '+%Y_%m_%d %T') [Pipeline error] More details can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit 126
    else    
	[[ -s $ERROR_TMP ]] && cat $ERROR_TMP > $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y_%m_%d %T') [Log directory] OK $LOG_DIR directory was created sucessfully. Will write log files in this directory." | tee -a $LOG_DIR/$LOGFILE 2>&1	
    fi
fi

#
# Test for cpu average load: TODO
# cf lib for a function to tell if average cpu load is ok else wait a minute
echo -ne "$(date '+%Y_%m_%d %T') [Batch mode] Checking for average cpu load before running on samples batch #$[$b+1] ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
# cat /proc/loadavg # avg1 avg5 avg10 running_threads/total_threads last_running_pid
echo -e "done" | tee -a $LOG_DIR/$LOGFILE 2>&1

#
# Test for disk space: TODO
# do it after loading config parameters to evaluate the used disk space for raw data (all subdirs with fastq files) 
# if available disk space lower than raw data used disk space, then abort
#echo "$(date '+%Y_%m_%d %T') [$(basename $0)] Test for available disk space" | tee -a $LOG_DIR/$LOGFILE 2>&1
#avail_disk_space=$(df -h $PWD | tail -1 | awk '{print $4}' 2>$ERROR_TMP)
echo -ne "$(date '+%Y_%m_%d %T') [Batch mode] Checking for available disk space before running on samples batch #$[$b+1] ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
#echo "$(date '+%Y_%m_%d %T') [$(basename $0)] Test for available disk space" | tee -a $LOG_DIR/$LOGFILE 2>&1
#avail_disk_space=$(df -h $PWD | tail -1 | awk '{print $4}' 2>$ERROR_TMP)
echo -e "done" | tee -a $LOG_DIR/$LOGFILE 2>&1

#
# Check for DATA_ROOT_DIR existence
#
echo "$(date '+%Y_%m_%d %T') [Data root directory] Checking $DATA_ROOT_DIR directory ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
if [[ -d $DATA_ROOT_DIR ]]; then
    echo "$(date '+%Y_%m_%d %T') [Data root directory] $DATA_ROOT_DIR exists and is a directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
else
    echo "$(date '+%Y_%m_%d %T') [Data root directory] Failed $DATA_ROOT_DIR does not exist or is not a directory." | tee $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %T') [Pipeline error] Exits the pipeline." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %T') [Pipeline error] More details can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit 1
fi

#
# Test for absence of user config file
# if present ok continue else display warning message and exit
#
echo "$(date '+%Y_%m_%d %T') [Check config: user config file] Checking for $VARCO_SPLIT_MAP_USER_CONFIG user config file ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
if [[ -s $VARCO_SPLIT_MAP_USER_CONFIG ]]; then
#if [[ -s $WORKING_DIR/$(basename ${0%.*})_user.config ]]; then # for testing purpose
    echo "$(date '+%Y_%m_%d %T') [Check config: user config file] OK User config file, $VARCO_SPLIT_MAP_USER_CONFIG, exists and is not empty." | tee -a $LOG_DIR/$LOGFILE 2>&1
else
    echo "$(date '+%Y_%m_%d %T') [Check config: user config file] Failed User config file, $VARCO_SPLIT_MAP_USER_CONFIG, does not exist or is empty" | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %T') [Check config: user config file] Warning: " | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo -e "\t\t$PREREQUISITES_MSG" | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %T') [Pipeline error] Exits the pipeline, with error code 3." | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit 3
fi

#
# Copy user config parameters file into job directory
# 1. Copy user config parameters file, prefixing it with the job tag
# 2. Then load config parameters from that new file, leaving the user config file in the working directory for some new job

# 1. Copy user config file 
echo "$(date '+%Y_%m_%d %T') [Check config: job user config file] Copying user config file into job directory ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
cp $VARCO_SPLIT_MAP_USER_CONFIG $VARCO_SPLIT_MAP_USER_CONFIG_JOB 
echo "$(date '+%Y_%m_%d %T') [Check config: job user config file] Will use copied job user config file: $VARCO_SPLIT_MAP_USER_CONFIG_JOB" | tee -a $LOG_DIR/$LOGFILE 2>&1

# 2. Load config parameters from job user config file
echo "$(date '+%Y_%m_%d %T') [Check config: job user config file] Loading job user config parameters from $VARCO_SPLIT_MAP_USER_CONFIG_JOB file ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
for cfg in $(get_config_sections $VARCO_SPLIT_MAP_USER_CONFIG_JOB 2>$ERROR_TMP; rtrn=$?); do
    echo -e "--- Config section [${cfg}] ---"
    unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${NAMESPACE}" 'BEGIN { 
          cfg = toupper(cfg);
          prefix = toupper(prefix);
       }
       /^prefix_cfg_/  { print $1 }' 2>>$ERROR_TMP) $(toupper ${NAMESPACE}_${cfg}_) 2>>$ERROR_TMP
    set_config_params $VARCO_SPLIT_MAP_USER_CONFIG_JOB ${cfg} ${NAMESPACE} 2>>$ERROR_TMP
    rtrn=$?
    for params in $(set | grep ^$(toupper ${NAMESPACE}_${cfg}_) 2>>$ERROR_TMP); do
	echo -e "$params"
    done
done
if [[ ! -s $ERROR_TMP ]]; then
    echo "$(date '+%Y_%m_%d %T') [Check config: job user config file] OK User config file, $VARCO_SPLIT_MAP_USER_CONFIG_JOB, was loaded successfully." | tee -a $LOG_DIR/$LOGFILE 2>&1
else
    echo "$(date '+%Y_%m_%d %T') [Check config: job user config file] Failed loading user config file, $VARCO_SPLIT_MAP_USER_CONFIG_JOB" | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %T') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %T') [Pipeline error] More details can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit $rtrn
fi

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
MAX_BATCH_SIZE=$[ $MAX_NUMB_CORES_ALLOWED/$VARCO_GSNAP_t ]
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
	VARCO_SPLIT_MAP_batch_size=$[$MAX_BATCH_SIZE/2]
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
if [[ -s $ERROR_TMP ]]; then 
    echo "$(date '+%Y_%m_%d %T') [Pipeline error] More details can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1; 
    exit 1
fi
#read -a all_subdirs_arr echo <<< $(echo -e $all_subdirs | tr " " "\n")
echo -e "all subdirectories count: ${#all_subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "all subdirectories list: ${all_subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1

# 2. Filter subdirs with include pattern(s)
echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] Filtering subdirs with include pattern(s) ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
#VARCO_DATA_include_sample_subdirs="^D.*,XA,XC" # for testing purpose
include_patterns=$(echo ${VARCO_DATA_include_sample_subdirs//,/|})
echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] include_patterns=($include_patterns) ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
inc_subdirs=($(for subdir in "${all_subdirs[@]}"; do
    res=$(echo $(basename $subdir) | egrep "($include_patterns)" 2>$ERROR_TMP)
    rtrn=$?
    [[ -n $res ]] && echo -e "$subdir" 
done 2>>$ERROR_TMP))
if [[ -s $ERROR_TMP ]]; then
        echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] Failed An error occured while filtering for subdirs to include." | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %T') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %T') [Pipeline error] More details can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit $rtrn
else
    echo -e "include pattern(s) subdirectories count: ${#inc_subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo -e "include pattern(s) subdirectories list: ${inc_subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1
fi

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
else
    echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] No exclude patterns." | tee -a $LOG_DIR/$LOGFILE 2>&1
    subdirs=("${inc_subdirs[@]}")
fi
if [[ -s $ERROR_TMP ]]; then
    echo "$(date '+%Y_%m_%d %T') [Fastq subdirs] Failed An error occured while filtering for subdirs to exclude." | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %T') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %T') [Pipeline error] More details can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit $rtrn
else
    echo -e "include pattern(s) subdirectories count: ${#subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE
    echo -e "include pattern(s) subdirectories list: ${subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE
fi

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
# 1.5 Wait for all QC+Trim processes to finish
# 1.6 Iterate over batch samples for Mapping
# 1.6.1 Create a mapping sub-subdir for each sample
# 1.6.2 Create mapping subdirs
# 1.6.3 Perform mapping -> waitall
# 1.7 Wait for all Mapping processes to finish
# 1.8 Iterate over commands for Conversion
# 1.8.1 Iterate over batch samples
# 1.8.2 Perform bam conversion -> waitall
# 1.8.3 Perform bam sorting -> waitall
# 1.8.4 Perform bam indexing -> waitall
# 1.9 Wait for all Conversion processes to finish
# 1.10 Unshifting current batch samples
# 1.11 Next batch
# 1.12 End of batches
# 1.13 Clean

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
for b in $(seq 0 $[ $batches-1 ]); do

    # 1.1 Test for average cpu load: TODO
    echo -ne "$(date '+%Y_%m_%d %T') [Batch mode] Checking for average cpu load before running on samples batch #$[$b+1] ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
	# cat /proc/loadavg # avg1 avg5 avg10 running_threads/total_threads last_running_pid
	echo -e "done" | tee -a $LOG_DIR/$LOGFILE 2>&1
	
    # 1.2 Test for available disk space: TODO
	echo -ne "$(date '+%Y_%m_%d %T') [Batch mode] Checking for available disk space before running on samples batch #$[$b+1] ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
	#echo "$(date '+%Y_%m_%d %T') [$(basename $0)] Test for available disk space" | tee -a $LOG_DIR/$LOGFILE 2>&1
	#avail_disk_space=$(df -h $PWD | tail -1 | awk '{print $4}' 2>$ERROR_TMP)
	echo -e "done" | tee -a $LOG_DIR/$LOGFILE 2>&1
    
	# 1.3 Iterate over batch samples for getting sample infos
	last=FALSE 	# for last batch sample
    echo "$(date '+%Y_%m_%d %T') [Batch mode] Running batch mode on samples batch #$[$b+1] ..." | tee -a $LOG_DIR/$LOGFILE 2>&1   
    for s in $(seq 1 $VARCO_SPLIT_MAP_batch_size); do
		si=$[$s-1]
		sdi=$[$si+$b*$VARCO_SPLIT_MAP_batch_size]
		echo -e "$sdi, $si: ${subdirs[$si]}"
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
    	echo "$(date '+%Y_%m_%d %T') [Batch mode] Running quality control and trimming on samples batch #$[$b+1] ..." | tee -a $LOG_DIR/$LOGFILE 2>&1 

	# reinitiate pids array
		PIDS_ARR=()
  
	# Iterate over batch samples for current quality control and trimming command
    	for s in $(seq 1 $VARCO_SPLIT_MAP_batch_size); do
			si=$[$s-1]
			sdi=$[$si+$b*$VARCO_SPLIT_MAP_batch_size]
			echo -e "$(date '+%Y_%m_%d %T') [Batch mode] Processing batch sample $sdi, $si: ${subdirs[$si]}, for quality control and trimming ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
			CURRENT_BATCH_SUBDIR=$JOB_TAG/$(basename "${subdirs[$si]}")
			CURRENT_SAMPLE_CONFIG=$CURRENT_BATCH_SUBDIR/$SAMPLE_CONFIG

	# load sample config
			echo "$(date '+%Y_%m_%d %T') [Batch mode] Loading sample config for $CURRENT_SAMPLE_CONFIG file ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
			sample_config_failed_msg="[Batch mode] Failed loading sample config file, $CURRENT_SAMPLE_CONFIG, from $CURRENT_BATCH_SUBDIR"
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
			echo "$(date '+%Y_%m_%d %T') [Batch mode] Creating quality control and trimming sub-directory: $CURRENT_BATCH_SUBDIR/$QC_TRIM_DIR ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
			CURRENT_QCTRIM_LOG=$CURRENT_BATCH_SUBDIR/$QC_TRIM_DIR/qc_trim.log
			CURRENT_QCTRIM_ERR=$CURRENT_BATCH_SUBDIR/$QC_TRIM_DIR/qc_trim_err.log
			if [[ -d $CURRENT_BATCH_SUBDIR/$QC_TRIM_DIR ]]; then
				echo "$(date '+%Y_%m_%d %T') [Batch mode: Quality Control and Trimming] OK $CURRENT_BATCH_SUBDIR/$QC_TRIM_DIR directory already exists. Will write output files in this directory." | tee -a $CURRENT_QCTRIM_LOG 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
			else
				mkdir $CURRENT_BATCH_SUBDIR/$QC_TRIM_DIR 2>$ERROR_TMP
				rtrn=$?
				qctrim_dir_failed_msg="[Batch mode] Failed Quality control and Trimming output directory, $CURRENT_BATCH_SUBDIR/$QC_TRIM_DIR, was not created."
				exit_on_error "$ERROR_TMP" "$qctrim_dir_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
				echo "$(date '+%Y_%m_%d %T') [Batch mode] OK $CURRENT_BATCH_SUBDIR/$QC_TRIM_DIR directory was created successfully. Will write output files in this directory." | tee -a $CURRENT_QCTRIM_LOG 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
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

	# waitall: TODO
	# several commands to wait: fastqc (trimmed=FALSE), trimmomatic, fastqc (trimmed=TRUE)


	# check for errors: TODO
	# for trimming, 
	# set trimmed to TRUE if ok
	# set forward and reverse fastq files for mapping step: update the sample config file
			echo -ne "$(date '+%Y_%m_%d %T') [Batch mode] Saving quality control and trimming infos for ${subdirs[$si]} current batch sample directory ... " | tee -a $CURRENT_QCTRIM_LOG 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
			echo -e "fastq_forward_trimmed_both_surviving=${forward_fastq_both_surviving}" 2>>$ERROR_TMP >>$CURRENT_SAMPLE_CONFIG
			echo -e "fastq_reverse_trimmed_both_surviving=${reverse_fastq_both_surviving}" 2>>$ERROR_TMP >>$CURRENT_SAMPLE_CONFIG
			echo -e "fastq_forward_trimmed_only_surviving=${forward_fastq_only_surviving}" 2>>$ERROR_TMP >>$CURRENT_SAMPLE_CONFIG
			echo -e "fastq_reverse_trimmed_only_surviving=${reverse_fastq_only_surviving}" 2>>$ERROR_TMP >>$CURRENT_SAMPLE_CONFIG
			echo -e "done" | tee -a $CURRENT_QCTRIM_LOG 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
			echo -e "$(date '+%Y_%m_%d %T') [Batch mode] updated sample config file:" | tee -a $CURRENT_QCTRIM_LOG 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
			echo -e "$(cat $sample_config_file)" | tee -a $CURRENT_QCTRIM_LOG 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1

	done
	fi

#	# 1.5 Iterate over batch samples for mapping 
#	
#	# reinitiate pids array
#	PIDS_ARR=()

#	# iterate over batch samples
#	for s in $(seq 1 $VARCO_SPLIT_MAP_batch_size); do
#		si=$[$s-1]
#		sdi=$[$si+$b*$VARCO_SPLIT_MAP_batch_size]
#		echo -e "$(date '+%Y_%m_%d %T') [Batch mode] Processing batch sample $sdi, $si: ${subdirs[$si]}, for mapping ... " | tee -a $LOG_DIR/$LOGFILE 2>&1
#		CURRENT_BATCH_SUBDIR=$JOB_TAG/$(basename "${subdirs[$si]}")
#		CURRENT_SAMPLE_CONFIG=$CURRENT_BATCH_SUBDIR/$SAMPLE_CONFIG

#	# load sample config
#		echo "$(date '+%Y_%m_%d %T') [Batch mode] Loading sample config for $CURRENT_SAMPLE_CONFIG file ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
#		sample_config_failed_msg="[Batch mode] Failed loading sample config file, $CURRENT_SAMPLE_CONFIG, from $CURRENT_BATCH_SUBDIR"
#		for cfg in $(get_config_sections $CURRENT_SAMPLE_CONFIG 2>$ERROR_TMP; rtrn=$?); do
#			if [[ "$rtrn" -eq 0 ]]; then
#				echo -e "--- Config section [${cfg}] ---"
#				unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${NAMESPACE}" 'BEGIN { 
#	  				cfg = toupper(cfg);
#	  				prefix = toupper(prefix);
#					}
#					/^prefix_cfg_/  { print $1 }' 2>$ERROR_TMP) $(toupper ${NAMESPACE}_${cfg}_) 2>>$ERROR_TMP
#				rtrn=$?
#				exit_on_error "$ERROR_TMP" "$sample_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
#				set_config_params $CURRENT_SAMPLE_CONFIG ${cfg} ${NAMESPACE} 2>$ERROR_TMP
#				rtrn=$?
#				exit_on_error "$ERROR_TMP" "$sample_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE" 
#				for params in $(set | grep ^$(toupper ${NAMESPACE}_${cfg}_) 2>$ERROR_TMP); do
#					echo -e "$params"
#				done
#			else
#				exit_on_error "$ERROR_TMP" "$sample_config_failed_msg" $rtrn "$LOG_DIR/$LOGFILE"
#			fi
#		done

#	# 1.5 Create a mapping sub-subdir
#	echo "$(date '+%Y_%m_%d %T') [Batch mode] Creating mapping sub-directory: $CURRENT_BATCH_SUBDIR/$MAPPING_DIR" | tee -a $LOG_DIR/$LOGFILE 2>&1
#	if [[ -d $CURRENT_BATCH_SUBDIR/$MAPPING_DIR ]]; then
#		echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping output directory] OK $CURRENT_BATCH_SUBDIR/$MAPPING_DIR directory already exists. Will write output files in this directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
#	else
#		mkdir $CURRENT_BATCH_SUBDIR/$MAPPING_DIR 2>$ERROR_TMP
#		
#		if [[ $? -ne 0 ]]; then
#		 	echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping output directory] Failed Mapping output directory, $CURRENT_BATCH_SUBDIR/$MAPPING_DIR was not created." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
#			echo "$(date '+%Y_%m_%d %T') [Pipeline error] Exits the pipeline, with error code 126." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
#			echo "$(date '+%Y_%m_%d %T') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
#			exit 126
#		else
#	 		echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping output directory] OK $CURRENT_BATCH_SUBDIR/$MAPPING_DIR directory was created successfully. Will write output files in this directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
#		fi
#	fi

#	# 1.5.1 Create mapping subdirs
#	mapping_subdirs=("tmp" "log")
#	for msd in "${mapping_subdirs[@]}"; do
#	echo "$(date '+%Y_%m_%d %T') [Batch mode] Creating mapping sub-sub-directory: $CURRENT_BATCH_SUBDIR/$MAPPING_DIR/$msd" | tee -a $LOG_DIR/$LOGFILE 2>&1
#	if [[ -d $CURRENT_BATCH_SUBDIR/$MAPPING_DIR ]]; then
#		echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping $msd directory] OK $CURRENT_BATCH_SUBDIR/$MAPPING_DIR/$msd directory already exists. Will write $msd output files in this directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
#	else
#		mkdir $CURRENT_BATCH_SUBDIR/$MAPPING_DIR/$msd 2>$ERROR_TMP
#		if [[ $? -ne 0 ]]; then
#		 	echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping $msd directory] Failed Mapping $msd output directory, $CURRENT_BATCH_SUBDIR/$MAPPING_DIR/$msd was not created." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
#			echo "$(date '+%Y_%m_%d %T') [Pipeline error] Exits the pipeline, with error code 126." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
#			echo "$(date '+%Y_%m_%d %T') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
#			exit 126
#		else
#	 		echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping $msd directory] OK $CURRENT_BATCH_SUBDIR/$MAPPING_DIR/$msd directory was created successfully. Will write $msd output files in this directory." | tee -a $LOG_DIR/$LOGFILE 2>&1
#		fi
#	fi		
#	done

#	# current mapping directories
#	CURRENT_MAPPING_DIR=$CURRENT_BATCH_SUBDIR/$MAPPING_DIR	
#	CURRENT_MAPPING_TMP=$CURRENT_MAPPING_DIR/tmp
#	CURRENT_MAPPING_LOG=$CURRENT_MAPPING_DIR/log
#	# current mapping files
#	CURRENT_MAPPING_ERROR=$CURRENT_MAPPING_TMP/error_mapping_b$b\_s$s\_$LOGFILE
#	CURRENT_MAPPING_LOGFILE=$CURRENT_MAPPING_LOG/mapping_b$b\_s$s\_$LOGFILE

#	# 1.5.2 Build gsnap mapper command
#	echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping] Building mapping command line options for current batch sample $CURRENT_BATCH_SUBDIR ..." | tee -a $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1	
#	gsnap_cli_options=($(buildCommandLineOptions gsnap $NAMESPACE 2>>$CURRENT_MAPPING_ERROR))
#	rtrn=$?
#	if [[ $rtrn -ne 0 ]]; then
#		echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping] Failed An error occured while building the mapping command line for current batch sample $CURRENT_BATCH_SUBDIR." | tee -a $ERROR_TMP 2>&1 | tee -a $CURRENT_MAPPING_ERROR 2>&1 | tee -a $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
#		echo "$(date '+%Y_%m_%d %T') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $CURRENT_MAPPING_ERROR 2>&1 | tee -a $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
#		echo "$(date '+%Y_%m_%d %T') [Pipeline error] More information can be found in $CURRENT_MAPPING_ERROR." | tee -a $ERROR_TMP 2>&1 | tee -a $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
#		exit $rtrn			
#	else
#		echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping] OK Successfully build mapping command line options for current batch sample $CURRENT_BATCH_SUBDIR." | tee -a $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
#	fi
#	gsnap_out=$CURRENT_SAMPLE_NAME\_gsnap_out_b$b\_s$s.sam	
#	gsnap_cli="gsnap $gsnap_cli_options $forward_fastq $reverse_fastq > $gsnap_out"
#	
#	# 1.5.2 Run the command	
#	echo "$(date '+%Y_%m_%d %T') [Batch mode: mapping] Executing mapping command line for current batch sample  $CURRENT_BATCH_SUBDIR ..." | tee -a $CURRENT_MAPPING_LOGFILE 2>&1  2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
#	echo -e "$(date '+%Y_%m_%d %T') [Batch mode: mapping] Mapper \n$MAPPER_VERSION" | tee -a $CURRENT_MAPPING_LOGFILE 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1 
#	eval "$gsnap_cli" 2>$CURRENT_MAPPING_ERROR &


#	# get the pid and push it to the pids array	
#	pid=$!

#	# Last batch sample: have a break!
#	[[ $last == "TRUE" ]] && break
#    done
#	done

#	# Run all mapping jobs then wait for these jobs to finish before running sam conversion to sorted bam files
#	# Check for errors in tmp directories, if exist if not empty    

#	# 1.6  Convert sam to sorted bam
#	last=FALSE
#    for s in $(seq 1 $VARCO_SPLIT_MAP_batch_size); do
#	si=$[$s-1]
#	sdi=$[$si+$b*$VARCO_SPLIT_MAP_batch_size]
#	echo -e "$sdi, $si: ${subdirs[$si]}"
#	CURRENT_BATCH_SUBDIR=$JOB_TAG/$(basename "${subdirs[$si]}")
#	if [[ "$s" -eq "${#subdirs[@]}" ]]; then
#	    last=TRUE
#	fi

#	# 1.6.1 Build samtools view command
# 
#	#.1.6.2 Run the command

#	# 1.6.3 Build samtools sort command

#	# 1.6.4 Run the command

#	# 1.6.5 Build samtools index command
#	
#	# 1.6.6 Run the command



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




done




















#
# Clean: TODO
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
echo "$(date '+%Y_%m_%d %T') [$(basename $0)] Executed command: $0 $*" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -n "$(date '+%Y_%m_%d %T') [$(basename $0)] Elapsed time: " | tee -a $LOG_DIR/$LOGFILE 2>&1
echo |awk -v time="$SECONDS" '{print strftime("%Hh:%Mm:%Ss", time, 1)}' | tee -a $LOG_DIR/$LOGFILE 2>&1
echo "$(date '+%Y_%m_%d %T') [$(basename $0)] Exits the pipeline." | tee -a $LOG_DIR/$LOGFILE 2>&1
echo "$(date '+%Y_%m_%d %T') [$(basename $0)] More information about this job can be found in $LOG_DIR/$LOGFILE" | tee -a $LOG_DIR/$LOGFILE 2>&1

#exit 0


