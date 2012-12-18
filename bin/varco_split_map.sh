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

# Inclusion de la librairie de fonctions

PROD_PREFIX="/usr/local"
DEV_PREFIX="$(pwd)/.."
PREFIX=$DEV_PREFIX # TO BE CHANGED WHEN SWITCHING TO PROD
. $PREFIX/share/varco_split_map/lib/varco_split_map_lib.inc

# Positionnement des variables

ARGS=2
DATA_ROOT_DIR=$1
JOB_TAG=$2

DATE=$(date '+%Y_%m_%d_%H_%M_%S')
WORKING_DIR=$(pwd)
LOGFILE=${JOB_TAG}_$DATE.log

NAMESPACE="VARCO"
VARCO_SPLIT_MAP_SHARED=$PREFIX/share/$(basename ${0%.*})
PROD_VARCO_SPLIT_MAP_USER_CONFIG=$WORKING_DIR/$(basename ${0%.*})_user.config
DEV_VARCO_SPLIT_MAP_USER_CONFIG=$VARCO_SPLIT_MAP_SHARED/etc/$(basename ${0%.*})_user.config
VARCO_SPLIT_MAP_USER_CONFIG=$DEV_VARCO_SPLIT_MAP_USER_CONFIG # TO BE CHANGED WHEN SWITCHING TO PROD

MAX_NUMB_CORES=$(cat /proc/cpuinfo | grep processor | wc -l)
MAX_NUMB_CORES_ALLOWED=$[$MAX_NUMB_CORES/2]
MAX_BATCH_SIZE=0

LOG_DIR="log"
TRIMMING_DIR="01_Trimming"
MAPPING_DIR="02_Mapping"

ERROR_TMP="/tmp/$(basename ${0%.*})_error_${USER}_$DATE.log"

PREREQUISITES_MSG="You need to copy the user configuration file to your working directory
                and configure the options before running this script.
                To copy the user configuration template file, issue this command in your terminal:
                \$ cp ${VARCO_SPLIT_MAP_USER_CONFIG} .
                Then use any editor (emacs, gedit, nano, vi, etc.) to open it and configure the options.
                Save your changes, and you are ready to run the script."

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
           job_tag          Prefix to attach to any output files       

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

# Create log directory

if [[ -d $LOG_DIR ]]; then
    echo "$(date '+%Y_%m_%d %R') [$(basename $0)] Start running the pipeline (version: $VERSION)." | tee $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %R') [$(basename $0)] Executed command: $0 $*" | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %R') [Log directory] OK $LOG_DIR directory already exists. Will write log files in this directory." >> $LOG_DIR/$LOGFILE
else
    mkdir $LOG_DIR 2>$ERROR_TMP
    if [[ $? -ne 0 ]]; then
	echo "$(date '+%Y_%m_%d %R') [Log directory] Failed Log directory, $LOG_DIR, was not created." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y_%m_%d %R') [Pipeline error] Exits the pipeline, with error code 126." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y_%m_%d %R') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit 126
    else
	echo "$(date '+%Y_%m_%d %R') [$(basename $0)] Start running the pipeline." | tee $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y_%m_%d %R') [$(basename $0)] Executed command: $0 $*" | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y_%m_%d %R') [Log directory] OK $LOG_DIR directory was created sucessfully. Will write log files in this directory." >> $LOG_DIR/$LOGFILE	
    fi
fi

#
# Test for cpu average load: TODO
# cf lib for a function to tell if average cpu load is ok else wait a minute

#
# Check for DATA_ROOT_DIR existence
#
if [[ -d $DATA_ROOT_DIR ]]; then
    echo "$(date '+%Y_%m_%d %R') [$(basename $0)] $DATA_ROOT_DIR exists and is a directory." | tee $LOG_DIR/$LOGFILE 2>&1
else
    echo "$(date '+%Y_%m_%d %R') [$(basename $0)] Failed $DATA_ROOT_DIR does not exist or is not a directory." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %R') [Pipeline error] Exits the pipeline." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %R') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit 1
fi

#
# Test for absence of user config file
# if present ok continue else display warning message and exit
#
if [[ -s $VARCO_SPLIT_MAP_USER_CONFIG ]]; then
#if [[ -s $WORKING_DIR/$(basename ${0%.*})_user.config ]]; then # for testing purpose
    echo "$(date '+%Y_%m_%d %R') [Check config: user config file] OK User config file, $VARCO_SPLIT_MAP_USER_CONFIG, exists and is not empty." | tee -a $LOG_DIR/$LOGFILE 2>&1
else
    echo "$(date '+%Y_%m_%d %R') [Check config: user config file ] Failed User config file, $VARCO_SPLIT_MAP_USER_CONFIG, does not exist or is empty" | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %R') [Check config: user config file ] Warning: " | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo -e "\t\t$PREREQUISITES_MSG" | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %R') [Pipeline error] Exits the pipeline, with error code 3." | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit 3
fi

#
# Load config parameters
# 
echo "$(date '+%Y_%m_%d %R') [Check config: user config file] Loading user config parameters from $VARCO_SPLIT_MAP_USER_CONFIG file ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
for cfg in $(get_config_sections $VARCO_SPLIT_MAP_USER_CONFIG 2>$ERROR_TMP; rtrn=$?); do
    echo -e "--- Config section [${cfg}] ---"
    unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${NAMESPACE}" 'BEGIN { 
          cfg = toupper(cfg);
          prefix = toupper(prefix);
       }
       /^prefix_cfg_/  { print $1 }') $(toupper ${NAMESPACE}_${cfg}_)
    set_config_params $VARCO_SPLIT_MAP_USER_CONFIG ${cfg} ${NAMESPACE} 2>$ERROR_TMP
    rtrn=$?
    for params in $(set | grep ^$(toupper ${NAMESPACE}_${cfg}_) 2>$ERROR_TMP); do
	echo -e "$params"
    done
done
if [[ ! -s $ERROR_TMP ]]; then
    echo "$(date '+%Y_%m_%d %R') [Check config: user config file] OK User config file, $VARCO_SPLIT_MAP_USER_CONFIG, was loaded successfully." | tee -a $LOG_DIR/$LOGFILE 2>&1
else
    echo "$(date '+%Y_%m_%d %R') [Check config: user config file ] Failed loading user config file, $VARCO_SPLIT_MAP_USER_CONFIG" | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %R') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %R') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit $rtrn
fi

#
# Check for parameters validity
# cf lib implement a function to check for config parameters validity: TODO
#

#
# Override  defined batch_size:
# 1. get the number of cores
# 2. get the number of threads to use for gsnap
# 3. compute max_batch_size=#max_cores_allowed/#threads
# 4. if batch_size <= max_batch_size then ok use batch_size else batch_size=max_batch_size/2
#
#VARCO_SPLIT_MAP_batch_size=8 # for testing purpose
echo "$(date '+%Y_%m_%d %R') [Override user config: batch_size] Testing if override batch_size user defined config parameter value is needed ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
MAX_BATCH_SIZE=$[ $MAX_NUMB_CORES_ALLOWED/$VARCO_GSNAP_t ]
echo -e "max_cores=$MAX_NUMB_CORES" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "max_cores_allowed=$MAX_NUMB_CORES_ALLOWED" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "threads_by_sample=$VARCO_GSNAP_t" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "max_batch_size=$MAX_BATCH_SIZE" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "batch_size=$VARCO_SPLIT_MAP_batch_size" | tee -a $LOG_DIR/$LOGFILE 2>&1
if [[ $VARCO_SPLIT_MAP_batch_size -le $MAX_BATCH_SIZE ]]; then
	echo "$(date '+%Y_%m_%d %R') [Override user config: batch_size] No need to override batch_size user defined config parameter value." | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y_%m_%d %R') [Override user config: batch_size] Keep batch_size user defined config parameter value: $VARCO_SPLIT_MAP_batch_size" | tee -a $LOG_DIR/$LOGFILE 2>&1
else
	echo "$(date '+%Y_%m_%d %R') [Override user config: batch_size] Need to override batch_size user defined config parameter value." | tee -a $LOG_DIR/$LOGFILE 2>&1
	VARCO_SPLIT_MAP_batch_size=$[$MAX_BATCH_SIZE/2]
	echo "$(date '+%Y_%m_%d %R') [Override user config: batch_size] Override batch_size user defined config parameter value: $VARCO_SPLIT_MAP_batch_size." | tee -a $LOG_DIR/$LOGFILE 2>&1
fi

#
# Search for subdirectories with fastq files: TODO
# 

#
# Test for disk space: TODO
# do it after loading config parameters to evaluate the used disk space for raw data (all subdirs with fastq files) 
# if available disk space lower than raw data used disk space, then abort
#echo "$(date '+%Y_%m_%d %R') [$(basename $0)] Test for available disk space" | tee -a $LOG_DIR/$LOGFILE 2>&1
#avail_disk_space=$(df -h $PWD | tail -1 | awk '{print $4}' 2>$ERROR_TMP)


#
# Mapping: TODO
#

# create a directory named with JOB_TAG value, to save all outputs 

# iterate over batches
## test for average cpu load
## create an output directory for each sample
## create Quality control and trimming subdir: fastqc and trimmomatic, optionnal step
## create a mapping subdir for each sample
## map in parallel all samples in one batch
## test for disk space
## convert, sort and index
## if clean true, clean each sample subdir



#
# Clean: TODO
#

# unset environment variables with used namespace
unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${NAMESPACE}" 'BEGIN { 
          cfg = toupper(cfg);
          prefix = toupper(prefix);
       }
       /^prefix_cfg_/  { print $1 }') $(toupper ${prefix}_${cfg}_)

#=====
# END
#=====
echo "$(date '+%Y_%m_%d %R') [$(basename $0)] Executed command: $0 $*" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -n "$(date '+%Y_%m_%d %R') [$(basename $0)] Elapsed time: " | tee -a $LOG_DIR/$LOGFILE 2>&1
echo |awk -v time="$SECONDS" '{print strftime("%Hh:%Mm:%Ss", time, 1)}' | tee -a $LOG_DIR/$LOGFILE 2>&1
echo "$(date '+%Y_%m_%d %R') [$(basename $0)] Exits the pipeline." | tee -a $LOG_DIR/$LOGFILE 2>&1
echo "$(date '+%Y_%m_%d %R') [$(basename $0)] More information about the analysis can be found in $LOG_DIR/$LOGFILE" | tee -a $LOG_DIR/$LOGFILE 2>&1

#exit 0


