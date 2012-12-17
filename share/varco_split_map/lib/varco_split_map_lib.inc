#! /bin/bash

#
# VARCO_SPLIT_MAP FUNCTIONS LIBRARY
#

# author: Joseph Tran <Joseph.Tran@versailles.inra.fr>
# TODO: add license copyright 2012

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
   ' ${config_file}
}

#--------------------------------------------------------------
# Function : get_config_params config_file config
# Purpose  : Get variables (prefixed by corresponding config) from config in config file
get_config_params()
{
	typeset config_file=$1
  typeset config=$2
  typeset config_params

  config_params=$( 
       awk -F= -v Config="${config}" -v Prefix="${config}" '
       BEGIN { 
          Config = toupper(Config);
          patternConfig = "\\[" Config "]";
       }
       toupper($0)  ~ patternConfig,(/\[/ && toupper($0) !~ patternConfig)  { 
          if (/\[/ || NF <2) next;
          sub(/^[[:space:]]*/, "");
          sub(/[[:space:]]*=[[:space:]]/, "=");
          print Prefix"_"$0;
       } ' ${config_file} )

	# eval "${config_params}"
  echo -e "${config_params[@]}"
}















