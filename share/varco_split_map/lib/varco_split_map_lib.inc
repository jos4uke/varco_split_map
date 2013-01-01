#! /bin/bash

# TODO: add license copyright 2012

#
# VARCO_SPLIT_MAP FUNCTIONS LIBRARY
#
# Author: Joseph Tran <Joseph.Tran@versailles.inra.fr>
#
# Date: 2012-12-17 
#


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
          if (/\[/ || NF <2) next;
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
	
	# limits 
	jobCoresMax=$[ $sysCoresMax / $rf ]
	sysCoresUpLim=$[ $sysCoresMax - $sysCoresKeep ]

	# cores 
	coresAvailableArr=($(cat /proc/loadavg | cut -d" " -f1,2,3 | gawk -v jobCoresMax=$jobCoresMax -v sysCoresUpLim=$sysCoresUpLim '{
		ORS=" "; for(i=1;i<=NF;i++) {print int(sysCoresUpLim - $i)}; {print "\n"}			
		}'))
	
	for c in "${coresAvailableArr[@]}"; do
		if [[ "$c" -ge "$jobCoresMax" ]]; then
			continue
		else
			echo FALSE	
			return 1
		fi
	done
	echo TRUE
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
		echo FALSE
	else
		echo TRUE
	fi
}


































