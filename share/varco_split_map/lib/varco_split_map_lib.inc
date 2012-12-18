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
# Function : set_config_params config_file config_section [config_prefix]
# Purpose  : Set variables (prefixed by config_prefix and corresponding config_section) from config_section in config file
set_config_params()
{
	typeset config_file=$1
  typeset config_section=$2
	typeset config_prefix=$3
  typeset config_params

  config_params=$( 
       awk -F= -v Config="${config_section}" -v Prefix="${config_prefix}_" '
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
       } ' ${config_file} )

	eval "${config_params}"
  #echo -e "${config_params[@]}"
}















