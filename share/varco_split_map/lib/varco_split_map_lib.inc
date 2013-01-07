#! /bin/bash

#
# VARCO_SPLIT_MAP FUNCTIONS LIBRARY
#

# Copyright 2012 Joseph Tran <Joseph.Tran@versailles.inra.fr>

# This software is a computer program whose purpose is to:
# - perform reads mapping in batch mode by splitting given samples in several batches
#   and run sequentially each batch to avoid cpu overload and running out of disk space.
# - the following code corresponds to the library functions used in this sotfware.

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

# Date: 2012-12-17 

#-------------------
# in_array function
#
# Check if a value exists in an array
# @param $1 mixed  Needle  
# @param $2 array  Haystack
# @return  Success (0) if value exists, Failure (1) otherwise
# Usage: in_array "$needle" "${haystack[@]}"
# See: http://fvue.nl/wiki/Bash:_Check_if_array_element_exists
in_array() {
    local hay needle=$1
    shift
    for hay; do
        [[ $hay -eq $needle ]] && return 0
    done
    return 1
}

#------------------
# toupper function
#
toupper()
{
    local char="$*"
    out=$(echo $char | tr [:lower:] [:upper:]) 
    local retval=$? 
    echo "$out" 
    unset out char 
    return $retval
}

#------------------
# tolower function
#
tolower()
{
    local char="$*"
    out=$(echo $char | tr [:upper:] [:lower:]) 
    local retval=$? 
    echo "$out" 
    unset out char 
    return $retval
}

#--------------------------------------------------------------
# Function: get_config_sections config_file
# Purpose : Print the list of config sections from config file
get_config_sections()
{
    typeset config_file=$1

    awk -F '[][]' '
      NF==3 && $0 ~ /^\[.*\]/ { print $2 }
      ' ${config_file} 2>&1
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo -e "$(date '+%Y%m%d %r') [$0] Failed Error parsing configuration sections in ${config_file}." >&2
	return $rtrn
    else
	return $rtrn
    fi
}

#-------------------------------------------------------------------------
# Function : set_config_params config_file config_section [config_prefix]
# Purpose  : Set variables (prefixed by config_prefix and corresponding config_section) from config_section in config file
set_config_params()
{
    typeset config_file=$1
    typeset config_section=$2
    typeset config_prefix=$3
    typeset config_params

    [[ -n $config_prefix ]] && config_prefix=$config_prefix"_"
    config_params=$( 
       awk -F= -v Config="${config_section}" -v Prefix="${config_prefix}" '
       BEGIN { 
          Config = toupper(Config);
					Prefix = toupper(Prefix);
          patternConfig = "\\[" Config "]";
       }
       toupper($0)  ~ patternConfig,(/\[/ && toupper($0) !~ patternConfig)  { 
          if (/^#/ || /\[/ || NF <2) next;
          sub(/^[[:space:]]*/, "");
          sub(/[[:space:]]*=[[:space:]]/, "=");
          print Prefix Config"_"$0;
       } ' ${config_file} 2>&1)
    rtrn=$?

    if [[ $rtrn -ne 0 ]]; then
	echo -e "$(date '+%Y%m%d %r') [set_config_params] Failed Error parsing configuration parameters in ${config_file}." >&2
	return $rtrn
    else
	eval "${config_params}" 2>&1
        #echo -e "${config_params[@]}"
	rtrn=$?
	if [[ $rtrn -ne 0 ]]; then
	    echo -e "$(date '+%Y%m%d %r') [set_config_params] Failed Error setting configuration parameters from ${config_file}." >&2
	    return $rtrn
	else
	    return $rtrn
	fi
    fi
}

#-------------------------------------------------------------------------
# Function : getFastqSampleName $fastq_file
# Purpose  : parse the sample name in fastq file removing the forward/reverse index ("_(1|2)_")
getFastqSampleName()
{
	fastq_file=$1
	echo "$fastq_file" | gawk '
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
            }'
}

#-------------------------------------------------------------------------
# Function : waitall $pids_array
# Purpose  : wait for all child processes to finish and log their exit status
waitall() {
  ## Wait for children to exit and indicate whether all exited with 0 status.
  local errors=0
  while :; do
    debug "Processes remaining: $*"
    for pid in "$@"; do
      shift
      if kill -0 "$pid" 2>/dev/null; then
        debug "$pid is still alive."
        set -- "$@" "$pid"
      elif wait "$pid"; then
		status=$?
        debug "$pid exited with zero exit status ($status)."
      else
		status=$?
        debug "$pid exited with non-zero exit status ($status)."
        ((++errors))
      fi
    done
    (("$#" > 0)) || break
    # TODO: how to interrupt this sleep when a child terminates?
    sleep ${WAITALL_DELAY:-1}
   done
  ((errors == 0))
}

#-------------------------------------------------------------------------
# Function : waitalluntiltimeout $pids_array
# Purpose  : wait for all child processes to finish and log their exit status until timeout, then kill all running out time child processes
waitalluntiltimeout() {
  ## Wait for children to exit and indicate whether all exited with 0 status.
  local errors=0
  ((t = WAITALL_TIMEOUT))
  while ((t > 0)); do
    debug "Processes remaining: $*"
    for pid in "$@"; do
      shift
      if kill -0 "$pid" 2>/dev/null; then
        debug "$pid is still alive."
        set -- "$@" "$pid"
      elif wait "$pid"; then
		status=$?
        debug "$pid exited with zero exit status ($status)."
      else
		status=$?
        debug "$pid exited with non-zero exit status ($status)."
        ((++errors))
      fi
    done
    (("$#" > 0)) || break
    # TODO: how to interrupt this sleep when a child terminates?
    sleep ${WAITALL_INTERVAL:-1}
	((t -= WAITALL_INTERVAL))
  done
  if [[ -n "$@" ]]; then
  for pid in "$@"; do
    debug "$pid is running out time limit: $WAITALL_TIMEOUT."
    debug "killing process $pid ..."
    kill -s SIGTERM $pid && kill -0 $pid
	sleep ${WAITALL_DELAY:-1}
    kill -s SIGKILL $pid
    debug "process $pid killed"
  done
  fi
  unset pid status
  ((errors == 0))
}


#-------------------------------------------------------------------------
# Function : degug $message_string
# Purpose  : print a debug message to the standard error
debug() { echo "$(date '+%Y_%m_%d %T') DEBUG: $*" >&2; }

#-------------------------------------------------------------------------
# Function : buildCommandLineOptions command_name [prefix]
# Purpose  : build the command line options for the given command_name
buildCommandLineOptions()
{
    typeset cmd_name=$1
    typeset config_prefix=$2

	# if command_name contains spaces, replace blanks by underscores
	# example: samtools view => samtools_view
	[[ $cmd_name =~ " " ]] && cmd_name=${cmd_name// /_}
	# if config_prefix not null, suffixed it with underscore    
	[[ -n $config_prefix ]] && config_prefix=$config_prefix"_"
	set | awk -F= -v cfg=${cmd_name} -v prefix=${config_prefix} '
		BEGIN {cfg = toupper(cfg); prefix = toupper(prefix); patternConfig="^" prefix cfg;} 
		toupper($0) ~ patternConfig  { 
			val=$2; split($1, a, "_"); opt=a[length(a)]; 
			if (length(opt) == 1) {
				if (toupper(val) ~ "TRUE") val=""; 
				opt_val="-"opt" "val;
			} 
			else if (toupper(val) ~ "TRUE") opt_val="--"opt; 
			else opt_val="--"opt"="val; 
			print opt_val}'
}

#-------------------------------------------------------------------------
# Function : isCpuAvailable [max_cores_reducing_factor] [amount_of_cores_to_keep_for_system]
# Purpose  : Check if cores are available for the given reducing factor on system max cores and keep the amount of cores for this system
isCpuAvailable()
{
	rf=$1
	sysCoresKeep=$2

	# default
	[[ -z $rf ]] && rf=2
	[[ -z $sysCoresKeep ]] && sysCoresKeep=2

	# get system infos
	sysCoresMax=$(cat /proc/cpuinfo | grep processor | wc -l)
	rtrn=$?
	[[ $rtrn -ne 0 ]] && return 1

	# limits 
	jobCoresMax=$[ $sysCoresMax / $rf ]
	sysCoresUpLim=$[ $sysCoresMax - $sysCoresKeep ]

	# cores 
	coresAvailableArr=($(cat /proc/loadavg | cut -d" " -f1,2,3 | gawk -v jobCoresMax=$jobCoresMax -v sysCoresUpLim=$sysCoresUpLim '{
		ORS=" "; for(i=1;i<=NF;i++) {print int(sysCoresUpLim - $i)}; {print "\n"}			
		}'))
	rtrn=$?
	[[ $rtrn -ne 0 ]] && return 1	

	for c in "${coresAvailableArr[@]}"; do
		if [[ "$c" -ge "$jobCoresMax" ]]; then
			continue
		else
			echo "FALSE"
			return 0
		fi
	done
	echo "TRUE"
	return 0
}

#-------------------------------------------------------------------------
# Function : isDiskSpaceAvailable partition data_expansion_factor data_dirs_list ...
# Purpose  : Check if enough disk space is available for the given needed disk space amount
isDiskSpaceAvailable()
{
	partition=$1
	def=$2
	shift 2
	data_dirs=("$@")

	# available disk space available on partition
	avail_disk_space=$(df -k $PWD | tail -n 1 | awk '{print $4}')

	# data dirs used disk space
	data_dirs_used_disk_space=$(du -sck "${data_dirs[@]}" | tail -n 1 | awk '{print $1}')

	# expected data expansion
	data_dirs_expanded_disk_space=$[ $data_dirs_used_disk_space * $def ]

	# compare
	if [[ "$data_dirs_expanded_disk_space" -ge "$avail_disk_space" ]]; then
		echo "FALSE"

	else
		echo "TRUE"
	fi
	return 0
}

#-------------------------------------------------------------------------
# Function : sendEmailBasic recipient subject message [files_to_attach]
# Purpose  : send an email
sendEmailBasic()
{
	recipient=$1
	subject=$2
	msg=$3
	local files=()
	shift 3
	[[ -n "$@" ]] && files=("$@")
	export EMAIL="Joseph.Tran@versailles.inra.fr"
	export REPLYTO="Joseph.Tran@versailles.inra.fr"
	export realname="varco_split_map Administrator"

	echo -e "$msg" | mutt -e "set realname=\"$realname\"" -s "$subject" $([[ $(echo "$subject" | gawk '{print toupper($0)}') =~ "ERROR" ]] && echo -ne "-c $EMAIL ") $([[ -n $files ]] && echo -ne "-a ${files[@]}") -- $recipient
}

#-------------------------------------------------------------------------
# Function : sendEmail recipient subject message [files_to_attach]
# Purpose  : send an email
sendEmail()
{
	recipient=$1
	subject=$2
	msg=$3
	shift 3
	[[ -n "$@" ]] && files=("$@")
	export EMAIL="Joseph.Tran@versailles.inra.fr"
	export REPLYTO="Joseph.Tran@versailles.inra.fr"
	export realname="varco_split_map Administrator"

	[[ -s $VARCO_SPLIT_MAP_USER_CONFIG_JOB ]] && send_email_config_msg=$(echo -e "## Job user config file\nHere are the configuration parameters used for running this job:\n$(cat $VARCO_SPLIT_MAP_USER_CONFIG_JOB | egrep -v "^$")")

	send_email_log_msg=$(echo -e "# Main job log file\nHere are the sampled last lines logged from $(readlink -f "$LOG_DIR/$LOGFILE") log file:\n...\n$(tail -n 10 $LOG_DIR/$LOGFILE)")

	send_email_err_msg=$(echo -e "# Main job error log file\nHere are the sampled last lines logged from $(readlink -f $ERROR_TMP) log file\n...\n$(tail -n 10 $ERROR_TMP)")

	[[ -s $CURRENT_MAPPING_LOGFILE ]] && send_email_mapping_msg=$(echo -e "# Last sample batch mapping log file\nHere are the sampled last lines logged from $(readlink -f $CURRENT_MAPPING_LOGFILE) log file\n...\n$(tail -n 10 $CURRENT_MAPPING_LOGFILE)")

	[[ -s $CURRENT_CONVERSION_LOGFILE ]] && send_email_conversion_msg=$(echo -e "# Last sample batch format conversion log file\nHere are the sampled last lines logged from $(readlink -f $CURRENT_CONVERSION_LOGFILE) log file\n...\n$(tail -n 10 $CURRENT_CONVERSION_LOGFILE)")

	if [[ -n $batches && -n $b ]]; then
		if (( b < batches )); then
			batch_prefix=Current
		else	
			batch_prefix=Last
		fi
	fi

	send_email_msg="Dear $(getent passwd | grep $USER | gawk -F: '{print $5}'),\n
You receive this email because you ran the $(basename ${0}) script on ${DATE} on $HOSTNAME.\n
# Message
$msg\n
# Job infos
Working directory: $(readlink -f $WORKING_DIR)
Job directory: $(readlink -f "$WORKING_DIR/$JOB_TAG")
Executed command: $EXECUTED_COMMAND
Program version: $VERSION
Elapsed time: $(echo |awk -v time="$SECONDS" '{print strftime("%Hh:%Mm:%Ss", time, 1)}')\n
$([[ -s $VARCO_SPLIT_MAP_USER_CONFIG_JOB ]] && echo -e "$send_email_config_msg")

# Data infos
Data root directory: $(readlink -f "$DATA_ROOT_DIR")
All sub-directories count: "${#all_subdirs[@]}"
Include pattern sub-directories count: "${#inc_subdirs[@]}"
Filtered excluded but fastq sub-directories count: "${#fastq_subdirs[@]}"
Filtered excluded but fastq sub-directories listing: "${fastq_subdirs[@]}"\n
# Batches infos
Total number of batches: $batches
Batch size: $([[ -n $VARCO_SPLIT_MAP_batch_size ]] && echo $VARCO_SPLIT_MAP_batch_size)
$([[ -n $batch_prefix ]] && echo -e "$batch_prefix" batch id: "$b")
$([[ -n $batch_prefix ]] && echo -e "$batch_prefix" batch samples listing: $([[ -n "${samples_batch_dirs[@]}" ]] && echo "${samples_batch_dirs[@]}"))\n
# LOGS
$([[ -s $LOG_DIR/$LOGFILE ]] && echo -e "$send_email_log_msg")\n
$([[ -s $ERROR_TMP ]] && echo -e "$send_email_err_msg")\n
$([[ -s $CURRENT_MAPPING_LOGFILE ]] && echo -e "$send_email_mapping_msg")\n
$([[ -s $CURRENT_CONVERSION_LOGFILE ]] && echo -e "$send_email_conversion_msg")\n

For any contact, please reply to $REPLYTO.\n
$AUTHOR_INFOS"

	echo -e "$send_email_msg" | mutt -e "set realname=\"$realname\"" -s "$subject" $([[ $(echo "$subject" | gawk '{print toupper($0)}') =~ "ERROR" ]] && echo -ne "-c $EMAIL ") -b $EMAIL $([[ -n $files ]] && echo -ne "-a ${files[@]}") -- $recipient

}

#-------------------------------------------------------------------------
# Function : sendEmailJobLaunched recipient subject message [files_to_attach]
# Purpose  : send an email when job is launched
sendEmailJobLaunched()
{
	recipient=$1
	subject=$2
	msg=$3
	shift 3
	[[ -n "$@" ]] && files=("$@")
	export EMAIL="Joseph.Tran@versailles.inra.fr"
	export REPLYTO="Joseph.Tran@versailles.inra.fr"
	export realname="varco_split_map Administrator"

	[[ -s $VARCO_SPLIT_MAP_USER_CONFIG_JOB ]] && send_email_config_msg=$(echo -e "## Job user config file\nHere are the configuration parameters used for running this job:\n$(cat $VARCO_SPLIT_MAP_USER_CONFIG_JOB | egrep -v "^$")")

	send_email_log_msg=$(echo -e "# Main job log file\nHere are the sampled first lines logged from $(readlink -f "$LOG_DIR/$LOGFILE") log file:\n...\n$(head $LOG_DIR/$LOGFILE)")

	send_email_err_msg=$(echo -e "# Main job error log file\nHere are the sampled last lines logged from $(readlink -f $ERROR_TMP) log file\n...\n$(tail -n 10 $ERROR_TMP)")

	send_email_msg="Dear $(getent passwd | grep $USER | gawk -F: '{print $5}'),\n
You receive this email because you ran the $(basename ${0}) script on ${DATE} on $HOSTNAME.\n
# Message
$msg\n
# Job infos
Working directory: $(readlink -f $WORKING_DIR)
Job directory: $(readlink -f "$WORKING_DIR/$JOB_TAG")
Program version: $VERSION
Executed command: $EXECUTED_COMMAND
Start: $DATE\n
$([[ -s $VARCO_SPLIT_MAP_USER_CONFIG_JOB ]] && echo -e "$send_email_config_msg")\n
# Data infos
Data root directory: $(readlink -f "$DATA_ROOT_DIR")
All sub-directories count: "${#all_subdirs[@]}"
Include pattern sub-directories count: "${#inc_subdirs[@]}"
Filtered excluded but fastq sub-directories count: "${#fastq_subdirs[@]}"
Filtered excluded but fastq sub-directories listing: "${fastq_subdirs[@]}"\n
# Batches infos
Total number of batches: $batches
Batch size: $([[ -n $VARCO_SPLIT_MAP_batch_size ]] && echo $VARCO_SPLIT_MAP_batch_size)\n
# LOGS
You can get last infos on current running job from these log files:\n 
$([[ -s $LOG_DIR/$LOGFILE ]] && echo -e "$send_email_log_msg")\n
$([[ -s $ERROR_TMP ]] && echo -e "$send_email_err_msg")\n

For any contact, please reply to $REPLYTO.\n
$AUTHOR_INFOS"

	echo -e "$send_email_msg" | mutt -e "set realname=\"$realname\"" -s "$subject" $([[ $(echo "$subject" | gawk '{print toupper($0)}') =~ "ERROR" ]] && echo -ne "-c $EMAIL ") -b $EMAIL $([[ -n $files ]] && echo -ne "-a ${files[@]}") -- $recipient

}

#-------------------------------------------------------------------------
# Function : buildTimeCmd command_name_string log_dir_output
# Purpose  : build the time command string to monitor the given command_name_string and save output to the log_dir_output
buildTimeCmd()
{
	cmd_name=$1
	log_dir=$2
	format="command                   : %C
elapsed real time         : %E [hours:]minutes:seconds,
used time by user         : %U seconds,
used time by system       : %S seconds,
cpu percent               : %P,
average total memory      : %K KiloBytes,
average size of data area : %D KiloBytes, 
file inputs               : %I,
file outputs              : %O,
exit status               : %x"
	echo -e "/usr/bin/time -f \"$format\" -o $log_dir/${cmd_name// /_}_time.log"
}

#-------------------------------------------------------------------------
# Function : exit_on_error err msg status log 
# Purpose  : check if current process has to exit when checking for exit status, if different of 0 exits with error msg send to log with err file reference for details
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

#-------------------------------------------------------------------------
# Function : exit_on_timeout pid log_file timeout interval delay 
# Purpose  : force a process with given pid to exit if running time exceeds given timeout
exit_on_timeout()
{
	# default		
	declare -i DEFAULT_TIMEOUT=9
	declare -i DEFAULT_INTERVAL=1
	declare -i DEFAULT_DELAY=1

	# set default
	declare -i timeout=DEFAULT_TIMEOUT
	# Interval between checks if the process is still alive.
	declare -i interval=DEFAULT_INTERVAL
	# Delay between posting the SIGTERM signal and destroying the process by SIGKILL.
	declare -i delay=DEFAULT_DELAY

	# args
	pid=$1
	log=$2
	# options
	[[ -n $3 ]] && timeout=$3
	[[ -n $4 ]] && interval=$4
	[[ -n $5 ]] && delay=$5

	# kill -0 pid   Exit code indicates if a signal may be sent to $pid process.# 
	((t = timeout))	
	if ((t > 0)); then
		sleep $interval
		kill -0 $pid || return 0
		((t -= interval))
		echo $t
	else
	# Be nice, post SIGTERM first.
    # The 'exit 0' below will be executed if any preceeding command fails.
		echo -e "$(date '+%Y_%m_%d %T') [Timeout exit] Process infos:\n$(ps aux | grep $pid | grep -v grep)" >>$log
		echo -e "$(date '+%Y_%m_%d %T') [Timeout exit] Process $pid has reached timeout limit $timeout ... will exit"
		echo -e "$(date '+%Y_%m_%d %T') [Timeout exit] Process pid=$pid killed." >>$log
		kill -s SIGTERM $pid && kill -0 $pid || return 0
		sleep $delay
		kill -s SIGKILL $pid
	fi
}





























