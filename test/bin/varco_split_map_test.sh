#! /bin/bash

# TODO: add license copyright 2012

#
# VARCO_SPLIT_MAP TESTS SUITE
#
# Author: Joseph Tran <Joseph.Tran@versailles.inra.fr>
#
# Date: 2012-12-17
#


#--------------------------------------
# testFailedPrintingConfigSectionsList
#
testFailedPrintingConfigSectionsList()
{
    res=$(get_config_sections $USER_CONFIG_FILE 2>${stderrF})
    read -a array echo <<< $(echo -e $res | tr " " "\n")
    echo -e "config sections list: "${array[@]}
    assertTrue 'Unexpected array size, should be equal to 7' "[ ${#array[@]} == 7 ]"
    assertFalse 'Unexpected output to stderr' "[ -s ${stderrF} ]"
}

#-----------------------------------
# testFailedReadingUserConfigParams
#
testFailedReadingUserConfigParams()
{
    for cfg in $(get_config_sections $USER_CONFIG_FILE); do
	echo -e "--- Config section [${cfg}] ---"
	unset $(set | awk -F= -v cfg="${cfg}" 'BEGIN { 
          cfg = toupper(cfg);
       }
			/^cfg_/  { print $1 }') $(toupper ${cfg}_)
	set_config_params $USER_CONFIG_FILE ${cfg} 2>${stderrF}
	set | grep ^$(toupper ${cfg}_)
    done	
    
    [[ -s ${stderrF} ]] && echo -e "stderr output:"; cat ${stderrF}
    assertFalse 'Unexpected output to stderr' "[ -s ${stderrF} ]"
}

#-----------------------------------
# testFailedLoadingUserConfigParams
#
testFailedLoadingUserConfigParams()
{
    prefix="varco"
    for cfg in $(get_config_sections $USER_CONFIG_FILE); do
	echo -e "--- Config section [${cfg}] ---"
	unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${prefix}" 'BEGIN { 
          cfg = toupper(cfg);
          prefix = toupper(prefix);
       }
       /^prefix_cfg_/  { print $1 }') $(toupper ${prefix}_${cfg}_)
	set_config_params $USER_CONFIG_FILE ${cfg} ${prefix} 2>${stderrF}
	for params in $(set | grep ^$(toupper ${prefix}_${cfg}_)); do
	    echo -e "params: $params"
	    delimiter="="
	    declare -a arr
	    arr=($(echo ${params//$delimiter/ }))
	    echo -e "variable: ${arr[0]}"
	    echo -e "value: ${arr[1]}"
	done
	unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${prefix}" 'BEGIN { 
          cfg = toupper(cfg);
          prefix = toupper(prefix);
       }
       /^prefix_cfg_/  { print $1 }') $(toupper ${prefix}_${cfg}_)
    done	
    
    [[ -s ${stderrF} ]] && echo -e "stderr output:"; cat ${stderrF}
    assertFalse 'Unexpected output to stderr' "[ -s ${stderrF} ]"
}

#-------------------------------------------
# testFailedBuildingGsnapCommandLineOptions
#
testFailedBuildingCommandLineOptions()
{
	# get and set options for gsnap
    prefix="mytest"
    for cfg in $(get_config_sections $USER_CONFIG_FILE); do
		unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${prefix}" 'BEGIN { 
          cfg = toupper(cfg);
          prefix = toupper(prefix);
       }
       /^prefix_cfg_/  { print $1 }') $(toupper ${prefix}_${cfg}_)
		set_config_params $USER_CONFIG_FILE ${cfg} ${prefix} 2>${stderrF}
	done

	# build options command line
	command_name="gsnap"
	for params in $(set | grep ^$(toupper ${prefix}_${command_name}_)); do
	    echo -e "params: $params"
	done

	cli_options=($(buildCommandLineOptions $command_name $prefix 2>>${stderrF}))
	res="${cli_options[@]}"
	echo -e $res
	expected="-A sam -B 4 -D /data/temp_projects/AZM/INDEX -N 1 -d Brassica_napus_v3.0.scaffold.fa -k 15 -m 1 --nofails -s /data/temp_projects/AZM/INDEX/Brassica_napus_v3.0.scaffold.fa/Brassica_napus_CDS_v3.0.splicesites -t 2"

	unset $(set | awk -F= -v cfg="${command_name}" -v prefix="${prefix}" 'BEGIN { 
          cfg = toupper(cfg);
          prefix = toupper(prefix);
       }
       /^prefix_cfg_/  { print $1 }') $(toupper ${prefix}_${cfg}_)

	assertTrue "Unexpected options command line length, should be equivalent to ${#expected}" "[ ${#expected} == ${#res} ]"
	assertTrue "Unexpected options command line, should be equivalent to $expected" '[[ "${expected}" == "${res}" ]]'

    if [[ -s ${stderrF} ]]; then
		echo -e "stderr output:"; cat ${stderrF}
	fi    
	assertFalse 'Unexpected output to stderr' "[ -s ${stderrF} ]"
}



















#================

#
# Configuration
#

oneTimeSetUp()
{
    test_start=`date +%H:%M:%S`
    . ../../share/varco_split_map/lib/varco_split_map_lib.inc
    
    TEST_OUTPUT_DIR="../output"

    USER_CONFIG_FILE="../etc/varco_split_map_user.config"
}

setUp()
{
    OUTPUT_DIR="${SHUNIT_TMPDIR}/OUTPUT"
    stdoutF="${OUTPUT_DIR}/stdoutF"
    stderrF="${OUTPUT_DIR}/stderrF"
    mkdir $OUTPUT_DIR
}

tearDown()
{  
    echo "Test starts ${test_start}"
    test_end=`date +%H:%M:%S`
    echo "Test ends ${test_end}"
    #exec_start_time=`date +%s -d ${test_start}`
    #exec_end_time=`date +%s -d ${test_end}`
    #exec_time=$[${exec_end_time}-${exec_start_time}]
    #echo |awk -v time="$exec_time" '{print "execution time: " strftime("%Hh:%Mm:%Ss", time, 1)}'
    rm -rf $OUTPUT_DIR
    echo "------------"
}

. /usr/local/share/shunit2-2.1.6/src/shunit2
