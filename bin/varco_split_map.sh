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

VARCO_SPLIT_MAP_SHARED=$PREFIX/share/$(basename ${0%.*})
PROD_VARCO_SPLIT_MAP_USER_CONFIG=$WORKING_DIR/$(basename ${0%.*})_user.config
DEV_VARCO_SPLIT_MAP_USER_CONFIG=$VARCO_SPLIT_MAP_SHARED/etc/$(basename ${0%.*})_user.config
VARCO_SPLIT_MAP_USER_CONFIG=$DEV_VARCO_SPLIT_MAP_USER_CONFIG # TO BE CHANGED WHEN SWITCHING TO PROD

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
                         The max batch size equals the ratio number of processors 
                         divided by the number of threads (mapper option).
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
    echo "$(date '+%Y_%m_%d %r') [$(basename $0)] Start running the pipeline (version: $VERSION)." | tee $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %r') [$(basename $0)] Executed command: $0 $*" | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %r') [Log directory] OK $LOG_DIR directory already exists. Will write log files in this directory." >> $LOG_DIR/$LOGFILE
else
    mkdir $LOG_DIR 2>$ERROR_TMP
    if [[ $? -ne 0 ]]; then
	echo "$(date '+%Y_%m_%d %r') [Log directory] Failed Log directory, $LOG_DIR, was not created." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y_%m_%d %r') [Pipeline error] Exits the pipeline, with error code 126." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y_%m_%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit 126
    else
	echo "$(date '+%Y_%m_%d %r') [$(basename $0)] Start running the pipeline." | tee $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y_%m_%d %r') [$(basename $0)] Executed command: $0 $*" | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y_%m_%d %r') [Log directory] OK $LOG_DIR directory was created sucessfully. Will write log files in this directory." >> $LOG_DIR/$LOGFILE	
    fi
fi

#
# Test for absence of user config file
# if present ok continue else display warning message and exit
#
if [[ -s $VARCO_SPLIT_MAP_USER_CONFIG ]]; then
#if [[ -s $WORKING_DIR/$(basename ${0%.*})_user.config ]]; then # for testing purpose
    echo "$(date '+%Y_%m_%d %r') [Check config: user config file] OK User config file, $VARCO_SPLIT_MAP_USER_CONFIG, exists and is not empty." | tee -a $LOG_DIR/$LOGFILE 2>&1
else
    echo "$(date '+%Y_%m_%d %r') [Check config: user config file ] Failed User config file, $VARCO_SPLIT_MAP_USER_CONFIG, does not exist or is empty" | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %r') [Check config: user config file ] Warning: " | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo -e "\t\t$PREREQUISITES_MSG" | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y_%m_%d %r') [Pipeline error] Exits the pipeline, with error code 3." | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit 3
fi

#
# Test for cpu average load: TODO
# cf lib for a function to tell if average cpu load is ok else wait a minute


#
# Load config parameters: TODO
# 
# 

#
# Check for parameters validity: TODO
# 1. check DATA_ROOT_DIR path existence
# 2. cf lib implement a function to check for config parameters validity

#
# Override batch_size: TODO
# 1. get the number of cores
# 2. get the number of threads to use for gsnap
# 3. compute max_batch_size=#cores/#threads
# 4. if batch_size <= max_batch_size then ok use batch_size else batch_size=max_batch_size/2
#

#
# Search for subdirectories with fastq files: TODO
# 

#
# Test for disk space: TODO
# do it after loading config parameters to evaluate the used disk space for raw data (all subdirs with fastq files) 
# if available disk space lower than raw data used disk space, then abort
#echo "$(date '+%Y_%m_%d %r') [$(basename $0)] Test for available disk space" | tee -a $LOG_DIR/$LOGFILE 2>&1
#avail_disk_space=$(df -h $PWD | tail -1 | awk '{print $4}' 2>$ERROR_TMP)



#
# Quality control: TODO, fastqc and trimmomatic, optionnal step
#


#
# Mapping: TODO
#

# iterate over batches
## test for average cpu load
## map in parallel all samples in one batch
## test for disk space
## convert, sort and index
## if clean true, clean each sample subdir



#
# Clean: TODO
#

# unset environment variables with used namespace

#=====
# END
#=====
echo "$(date '+%Y_%m_%d %r') [$(basename $0)] Executed command: $0 $*" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -n "$(date '+%Y_%m_%d %r') [$(basename $0)] Elapsed time: " | tee -a $LOG_DIR/$LOGFILE 2>&1
echo |awk -v time="$SECONDS" '{print strftime("%Hh:%Mm:%Ss", time, 1)}' | tee -a $LOG_DIR/$LOGFILE 2>&1
echo "$(date '+%Y_%m_%d %r') [$(basename $0)] Exits the pipeline." | tee -a $LOG_DIR/$LOGFILE 2>&1
echo "$(date '+%Y_%m_%d %r') [$(basename $0)] More information about the analysis can be found in $LOG_DIR/$LOGFILE" | tee -a $LOG_DIR/$LOGFILE 2>&1

#exit 0


