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
LOGFILE=${JOB_TAG}_${USER}_$DATE.log

NAMESPACE="VARCO"
VARCO_SPLIT_MAP_SHARED=$PREFIX/share/$(basename ${0%.*})
PROD_VARCO_SPLIT_MAP_USER_CONFIG=$WORKING_DIR/$(basename ${0%.*})_user.config
DEV_VARCO_SPLIT_MAP_USER_CONFIG=$VARCO_SPLIT_MAP_SHARED/etc/$(basename ${0%.*})_user.config
VARCO_SPLIT_MAP_USER_CONFIG=$DEV_VARCO_SPLIT_MAP_USER_CONFIG # TO BE CHANGED WHEN SWITCHING TO PROD
JOB_VARCO_SPLIT_MAP_USER_CONFIG=$JOB_TAG/${JOB_TAG}_$(basename ${0%.*})_user.config

MAX_NUMB_CORES=$(cat /proc/cpuinfo | grep processor | wc -l)
MAX_NUMB_CORES_ALLOWED=$[$MAX_NUMB_CORES/2]
MAX_BATCH_SIZE=0

LOG_DIR=$JOB_TAG/"log"
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

#
# Create a directory named with JOB_TAG value, to save all outputs 
#
if [[ -d $JOB_TAG ]]; then
    echo "$(date '+%Y_%m_%d %R') [Job directory] OK $JOB_TAG directory already exists. Will output all job files in this directory." | tee -a $LOG_DIR/$LOGFILE
else
    mkdir $JOB_TAG 2>$ERROR_TMP
    if [[ $? -ne 0 ]]; then
	echo "$(date '+%Y_%m_%d %R') [Job directory] Failed Job directory, $JOB_TAG, was not created." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y_%m_%d %R') [Pipeline error] Exits the pipeline, with error code 126." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y_%m_%d %R') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit 126
    else
    echo "$(date '+%Y_%m_%d %R') [Job directory] OK $JOB_TAG directory already exists. Will output all job files in this directory." | tee -a $LOG_DIR/$LOGFILE
    fi
fi

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
# Copy user config parameters file into job directory: TODO
# 1. Copy user config parameters file, prefixing it with the job tag
# 2. Then load config parameters from that new file, leaving the user config file in the working directory for some new job

# 1. Copy user config file 
echo "$(date '+%Y_%m_%d %R') [Check config: job user config file] Copying user config file into job directory ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
cp $VARCO_SPLIT_MAP_USER_CONFIG $JOB_VARCO_SPLIT_MAP_USER_CONFIG 
echo "$(date '+%Y_%m_%d %R') [Check config: job user config file] Will use copied job user config file: $JOB_VARCO_SPLIT_MAP_USER_CONFIG" | tee -a $LOG_DIR/$LOGFILE 2>&1

# 2. Load config parameters from job user config file
echo "$(date '+%Y_%m_%d %R') [Check config: job user config file] Loading job user config parameters from $JOB_VARCO_SPLIT_MAP_USER_CONFIG file ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
for cfg in $(get_config_sections $JOB_VARCO_SPLIT_MAP_USER_CONFIG 2>$ERROR_TMP; rtrn=$?); do
    echo -e "--- Config section [${cfg}] ---"
    unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${NAMESPACE}" 'BEGIN { 
          cfg = toupper(cfg);
          prefix = toupper(prefix);
       }
       /^prefix_cfg_/  { print $1 }') $(toupper ${NAMESPACE}_${cfg}_)
    set_config_params $JOB_VARCO_SPLIT_MAP_USER_CONFIG ${cfg} ${NAMESPACE} 2>$ERROR_TMP
    rtrn=$?
    for params in $(set | grep ^$(toupper ${NAMESPACE}_${cfg}_) 2>$ERROR_TMP); do
	echo -e "$params"
    done
done
if [[ ! -s $ERROR_TMP ]]; then
    echo "$(date '+%Y_%m_%d %R') [Check config: job user config file] OK User config file, $JOB_VARCO_SPLIT_MAP_USER_CONFIG, was loaded successfully." | tee -a $LOG_DIR/$LOGFILE 2>&1
else
    echo "$(date '+%Y_%m_%d %R') [Check config: job user config file ] Failed loading user config file, $JOB_VARCO_SPLIT_MAP_USER_CONFIG" | tee -a $LOG_DIR/$LOGFILE 2>&1
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
# Search for subdirectories with fastq files:
# 1. list all available subdirs in data root dir
# 2. Filter subdirs with include pattern(s)
# 3. Filter subdirs with exclude pattern(s)
# 4. Filter subdirs for fastq files

# 1. list all available subdirs in data root dir
echo "$(date '+%Y_%m_%d %R') [Fastq subdirs] Searching for fastq sample subdirs ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
all_subdirs=($(find $DATA_ROOT_DIR -maxdepth 1 -type d | egrep -v ^$DATA_ROOT_DIR$ 2>$ERROR_TMP))
if [[ -s $ERROR_TMP ]]; then 
    echo "$(date '+%Y_%m_%d %R') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1; 
    exit 1
fi
#read -a all_subdirs_arr echo <<< $(echo -e $all_subdirs | tr " " "\n")
echo -e "all subdirectories count: ${#all_subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "all subdirectories list: ${all_subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1

# 2. Filter subdirs with include pattern(s)
echo "$(date '+%Y_%m_%d %R') [Fastq subdirs] Filtering subdirs with include pattern(s) ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
#VARCO_DATA_include_sample_subdirs="^D.*,XA,XC" # for testing purpose
include_patterns=$(echo ${VARCO_DATA_include_sample_subdirs//,/|})
echo "$(date '+%Y_%m_%d %R') [Fastq subdirs] include_patterns=($include_patterns) ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
inc_subdirs=($(for subdir in "${all_subdirs[@]}"; do
    res=$(echo $(basename $subdir) | egrep "($include_patterns)" 2>>$LOG_DIR/$LOGFILE)
    [[ -n $res ]] && echo -e "$subdir" 
done 2>>$LOG_DIR/$LOGFILE))
echo -e "include pattern(s) subdirectories count: ${#inc_subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "include pattern(s) subdirectories list: ${inc_subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1

# 3. Filter subdirs with exclude pattern(s)
echo "$(date '+%Y_%m_%d %R') [Fastq subdirs] Filtering subdirs with exclude pattern(s) ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
#VARCO_DATA_exclude_sample_subdirs="^X.*,DA,DI" # for testing purpose
#VARCO_DATA_exclude_sample_subdirs="" # for testing purpose
if [[ -n $VARCO_DATA_exclude_sample_subdirs ]]; then
    exclude_patterns=$(echo ${VARCO_DATA_exclude_sample_subdirs//,/|})
    echo "$(date '+%Y_%m_%d %R') [Fastq subdirs] exclude_patterns=($exclude_patterns) ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
    subdirs=($(for subdir in "${inc_subdirs[@]}"; do
	    res=$(echo $(basename $subdir) | egrep -v "($exclude_patterns)" 2>>$LOG_DIR/$LOGFILE)
	    [[ -n $res ]] && echo -e "$subdir" 
	    done 2>>$LOG_DIR/$LOGFILE))
else
    echo "$(date '+%Y_%m_%d %R') [Fastq subdirs] No exclude patterns." | tee -a $LOG_DIR/$LOGFILE 2>&1
    subdirs=("${inc_subdirs[@]}")
fi
echo -e "include pattern(s) subdirectories count: ${#subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE
echo -e "include pattern(s) subdirectories list: ${subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE

# 4. Filter subdirs for fastq files
echo "$(date '+%Y_%m_%d %R') [Fastq subdirs] Filtering subdirs for fastq files ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
fastq_forward_pattern="*_1_*.fastq"
fastq_reverse_pattern="*_2_*.fastq"
fastq_subdirs=($(for subdir in "${subdirs[@]}"; do
	fastq_files=($(ls "$subdir" | egrep -v "_single_" | egrep ".*.fastq$"))
	if [[ "${#fastq_files[@]}" == 2 ]]; then
	    echo "$(date '+%Y_%m_%d %R') [Fastq subdirs] INFO $subdir has 2, non '_single_', fastq files." 1>&2 | tee -a $LOG_DIR/$LOGFILE 1>&2
            # get the pair if possible
	    # else warn and discard the subdir
	    sample_name_1=$(echo "${fastq_files[0]}" | gawk '
	    function getSampleName(str) {
              sample="";
              while(match(str, /^(.*_[0-9])_[1-2]_(.*).fastq$/, a)) 
              {
                sample=a[1]"_"a[2];
                str = substr(str, RSTART+RLENGTH)
              }
              print sample;
            }
            {
              getSampleName($0)
            }' 2>>$LOG_DIR/$LOGFILE
	    )
	    sample_name_2=$(echo "${fastq_files[1]}" | gawk '
	    function getSampleName(str) {
              sample="";
              while(match(str, /^(.*_[0-9])_[1-2]_(.*).fastq$/, a)) 
              {
                sample=a[1]"_"a[2];
                str = substr(str, RSTART+RLENGTH)
              }
              print sample;
            }
            {
              getSampleName($0) 
            }' 2>>$LOG_DIR/$LOGFILE
	    )
	    echo -e "$(date '+%Y_%m_%d %R') [Fastq subdirs] INFO $subdir, sample_name: forward=$sample_name_1 == reverse=$sample_name_2" 1>&2 | tee -a $LOG_DIR/$LOGFILE 1>&2
	    [[ $sample_name_1 == $sample_name_2 ]] && echo -e "$subdir"
	elif [[ "${#fastq_files[@]}" < 2 ]]; then
	    echo "$(date '+%Y_%m_%d %R') [Fastq subdirs] Warning: $subdir has less than 2, non '_single_', fastq files." 1>&2 | tee -a $LOG_DIR/$LOGFILE 1>&2
	    echo "$(date '+%Y_%m_%d %R') [Fastq subdirs] Warning: $subdir will not be considered because non '_single_' fastq files are missing." 1>&2| tee -a $LOG_DIR/$LOGFILE 1>&2
	elif [[ "${#fastq_files[@]}" > 2 ]]; then
	    # not satisfying
	    echo "$(date '+%Y_%m_%d %R') [Fastq subdirs] Warning: $subdir has more than 2, non '_single_', fastq files." 1>&2 | tee -a $LOG_DIR/$LOGFILE 1>&2
	    echo "$(date '+%Y_%m_%d %R') [Fastq subdirs] Warning: You need to leave in $subdir only the 2 forward (_1_) and reverse (_2_) fastq files. '_single_' fastq files are ignored." 1>&2| tee -a $LOG_DIR/$LOGFILE 1>&2
	    echo "$(date '+%Y_%m_%d %R') [Fastq subdirs] Warning: $subdir will not be considered, ambiguous fastq files list." 1>&2 | tee -a $LOG_DIR/$LOGFILE 1>&2
	fi
	done 2>>$LOG_DIR/$LOGFILE))
echo -e "fastq subdirectories count: ${#fastq_subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -e "fastq subdirectories list: ${fastq_subdirs[@]}" | tee -a $LOG_DIR/$LOGFILE 2>&1

#
# Test for disk space: TODO
# do it after loading config parameters to evaluate the used disk space for raw data (all subdirs with fastq files) 
# if available disk space lower than raw data used disk space, then abort
#echo "$(date '+%Y_%m_%d %R') [$(basename $0)] Test for available disk space" | tee -a $LOG_DIR/$LOGFILE 2>&1
#avail_disk_space=$(df -h $PWD | tail -1 | awk '{print $4}' 2>$ERROR_TMP)


#
# Batch mode:
# 1. Iterate over batches

# iterate over batches
## test for average cpu load
## create an output directory for each sample
## create Quality control and trimming subdir: fastqc and trimmomatic, optionnal step
## create a mapping subdir for each sample
## map in parallel all samples in one batch
## test for disk space
## convert, sort and index
## if clean true, clean each sample subdir

# 1. Iterate over batches
for b in $(seq 1 $[ ${#fastq_subdirs[@]}/$VARCO_SPLIT_MAP_batch_size ]); do
    echo $b
done




















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


