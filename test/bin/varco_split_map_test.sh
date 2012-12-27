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

#-------------------------------
# testFailedWaitingCliProcesses
#
testFailedWaitingCliProcesses()
{
	# test variables
	samples=("DG" "DI")
	pids_arr=()
	# test functions
	exit_on_error()
	{
		[[ -s ${stderrF} ]] && exit_on_error
	}

	# get and set cli options
	echo -ne "Getting and setting cli options ... " | tee -a ${stdoutF} 2>&1
    prefix="mytest"
    for cfg in $(get_config_sections $USER_CONFIG_FILE); do
		unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${prefix}" 'BEGIN { 
          cfg = toupper(cfg);
          prefix = toupper(prefix);
       }
       /^prefix_cfg_/  { print $1 }') $(toupper ${prefix}_${cfg}_)
		set_config_params $USER_CONFIG_FILE ${cfg} ${prefix} 2>${stderrF}
	done
	[[ -s ${stderrF} ]] && exit_on_error
	echo -e "done" | tee -a ${stdoutF} 2>&1

	# for each sample
	for s in "${samples[@]}"; do
		echo -e "Processing sample $s ... " | tee -a ${stdoutF} 2>&1

	## create sample output dir
		CURRENT_SAMPLE_DIR=$TEST_OUTPUT_DIR/$s	
		if [[ ! -d $CURRENT_SAMPLE_DIR ]]; then
			echo -ne "Creating sample output directory $CURRENT_SAMPLE_DIR ... " | tee -a ${stdoutF} 2>&1
			mkdir $CURRENT_SAMPLE_DIR 2>${stderrF}
			echo -e "done" | tee -a ${stdoutF} 2>&1
		else
			echo -e "Sample output directory $CURRENT_SAMPLE_DIR already exists." | tee -a ${stdoutF} 2>&1
		fi

	## get fastq files
		fastq_files=($(ls "$TEST_DATA_ROOT_DIR/$s" | egrep ".*.fastq$" 2>${stderrF}))
		if [[ ! -s ${stderrF} ]]; then 
			echo -e "fastq files count: ${#fastq_files[@]}" | tee -a ${stdoutF} 2>&1
			echo -e "fastq files list: ${fastq_files[@]}" | tee -a ${stdoutF} 2>&1
		else
			exit_on_error	
		fi
		forward_fastq=$(for f in "${fastq_files[@]}"; do 
		m=$(echo $f | egrep "_[0-9]+_1_" 2>${stderrF}); if [[ -n $m ]]; then echo $m; break; fi
		done 2>${stderrF})
		exit_on_error
		echo -e "forward fastq: $forward_fastq"
		reverse_fastq=$(for f in "${fastq_files[@]}"; do 
		m=$(echo $f | egrep "_[0-9]+_2_" 2>${stderrF}); if [[ -n $m ]]; then echo $m; break; fi
		done 2>${stderrF})
		exit_on_error
		echo -e "reverse fastq: $reverse_fastq"
		CURRENT_SAMPLE_NAME=$(echo "${fastq_files[0]}" | gawk '
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
           	 }' 2>${stderrF})
		exit_on_error
		echo -e "current sample name: $CURRENT_SAMPLE_NAME"

	## build mapping cli
		command_name="gsnap"
		echo -e "--- Config section [${command_name}] ---"
		for params in $(set | grep ^$(toupper ${prefix}_${command_name}_)); do
	    	echo -e "params: $params"
		done

		cli_options=($(buildCommandLineOptions $command_name $prefix 2>${stderrF}))
		if [[ ! -s ${stderrF} ]]; then
			res="${cli_options[@]}"	
			gsnap_out=$CURRENT_SAMPLE_NAME\_gsnap_out_s$s.sam
			CURRENT_SAMPLE_ROOT_DIR=$TEST_DATA_ROOT_DIR/$s
			CURRENT_SAMPLE_LOG=$CURRENT_SAMPLE_DIR/$CURRENT_SAMPLE_NAME.log
			CURRENT_SAMPLE_ERR=$CURRENT_SAMPLE_DIR/$CURRENT_SAMPLE_NAME\_err.log	
			gsnap_cli="gsnap $res $CURRENT_SAMPLE_ROOT_DIR/$forward_fastq $CURRENT_SAMPLE_ROOT_DIR/$reverse_fastq > $CURRENT_SAMPLE_DIR/$gsnap_out 2>${CURRENT_SAMPLE_ERR} &"
	## run cli
			echo -e "Executing command: ${gsnap_cli}"
			eval "$gsnap_cli" 2>&1  
			pid=$!
			echo -e "$CURRENT_SAMPLE_NAME pid: $pid"
			pids_arr=("${pids_arr[@]}" "$pid")
			echo -e "${pids_arr[@]}"
		else
			exit_on_error
		fi
	done

	# wait for all gsnap processes to finish
	echo -e "pids array: ${pids_arr[@]}"
	echo -e "gsnap processes: ${#pids_arr[@]}"
	for p in "${pids_arr[@]}"; do
		echo -e $(ps aux | grep $p | grep -v grep )
	done
	WAITALL_DELAY=60
	waitall "${pids_arr[@]}"
	echo -e "${pids_arr[@]}"

	# check for errors
	errs=0
	for s in "${samples[@]}"; do
		CURRENT_SAMPLE_DIR=$TEST_OUTPUT_DIR/$s
		CURRENT_SAMPLE_ERR=$CURRENT_SAMPLE_DIR/$(ls "$CURRENT_SAMPLE_DIR" | egrep ".*_err.log$" | grep -v egrep 2>${stderrF})
		if [[ -z $(tail -1 ${CURRENT_SAMPLE_ERR} | egrep "^Processed") ]]; then
			echo -e "sample $s stderr output:"; cat ${CURRENT_SAMPLE_ERR}
			let errs=errs+1
		fi    
		assertTrue "Expected output to stderr for sample $s" "[ -s ${CURRENT_SAMPLE_ERR} ]"
	done

	if [[ $errs == 0 ]]; then
		echo -e "Processing all samples without errors." | tee -a ${stdoutF} 2>&1
	else
		echo -e "Some errors occured while processing samples." | tee -a ${stdoutF} 2>&1
	fi

	# for each sample
	## build conversion cli
	## run cli
	# wait for all conversion processes to finish

	# clean
	# test for assertions
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

	TEST_DATA_ROOT_DIR="../data/SEQUENCES"
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
