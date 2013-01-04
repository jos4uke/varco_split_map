#! /bin/bash

#
# VARCO_SPLIT_MAP TESTS SUITE
#

# Copyright 2012 Joseph Tran <Joseph.Tran@versailles.inra.fr>

# This software is a computer program whose purpose is to:
# - perform reads mapping in batch mode by splitting given samples in several batches
#   and run sequentially each batch to avoid cpu overload and running out of disk space.
# - the following code is a tests suite intended to test the library functions of this software.

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
	WAITALL_TIMEOUT=86400
	WAITALL_INTERVAL=60
	WAITALL_DELAY=60
	
	# test functions
	exit_on_error()
	{
		[[ -s ${stderrF} ]] && $(cat ${stderrF} 2>&1 | tee -a $TEST_PREFIX_ERR 2>&1; exit 1)
	}

	# main prefix files
	TEST_PREFIX_DIR=$TEST_OUTPUT_DIR/$prefix
	TEST_PREFIX_LOG=$TEST_PREFIX_DIR/$prefix.log
	TEST_PREFIX_ERR=$TEST_PREFIX_DIR/$prefix\_err.log

	# processing samples
	echo -e "$(date '+%Y_%m_%d %T') [Init] Processing the following samples: ${samples[@]}" | tee -a ${stdoutF} 2>&1

	# get and set cli options
	echo -ne "$(date '+%Y_%m_%d %T') [Config] Getting and setting cli options ... " | tee -a ${stdoutF} 2>&1
    prefix="mytest"
    for cfg in $(get_config_sections $USER_CONFIG_FILE); do
		unset $(set | awk -F= -v cfg="${cfg}" -v prefix="${prefix}" 'BEGIN { 
          cfg = toupper(cfg);
          prefix = toupper(prefix);
       }
       /^prefix_cfg_/  { print $1 }') $(toupper ${prefix}_${cfg}_)
		set_config_params $USER_CONFIG_FILE ${cfg} ${prefix} 2>${stderrF}
	done
	exit_on_error
	echo -e "done" | tee -a ${stdoutF} 2>&1

	# create the prefix test output directory
	if [[ ! -d $TEST_PREFIX_DIR ]]; then
		echo -ne "$(date '+%Y_%m_%d %T') [Directory] Creating test prefix $prefix output directory $TEST_PREFIX_DIR ... " | tee -a ${stdoutF} 2>&1
		mkdir $TEST_PREFIX_DIR 2>${stderrF}
		exit_on_error
		echo -e "done" | tee -a ${stdoutF} 2>&1
	else
		echo -e "$(date '+%Y_%m_%d %T') [Directory] Test prefix output directory $TEST_PREFIX_DIR already exists." | tee -a ${stdoutF} 2>&1
	fi

	# copy std out to main prefix log
	cat ${stdoutF} 2>&1 > $TEST_PREFIX_LOG

	# for each sample
	for s in "${samples[@]}"; do
		echo -e "$(date '+%Y_%m_%d %T') [Mapping] Processing sample $s ... " | tee ${stdoutF} 2>&1

	## create sample output dir
		CURRENT_SAMPLE_DIR=$TEST_PREFIX_DIR/$s	
		if [[ ! -d $CURRENT_SAMPLE_DIR ]]; then
			echo -ne "$(date '+%Y_%m_%d %T') [Mapping] Creating sample output directory $CURRENT_SAMPLE_DIR ... " | tee -a ${stdoutF} 2>&1
			mkdir $CURRENT_SAMPLE_DIR 2>${stderrF}
			exit_on_error
			echo -e "done" | tee -a ${stdoutF} 2>&1
		else
			echo -e "$(date '+%Y_%m_%d %T') [Mapping] Sample output directory $CURRENT_SAMPLE_DIR already exists." | tee -a ${stdoutF} 2>&1
		fi

	## get fastq files
		fastq_files=($(ls "$TEST_DATA_ROOT_DIR/$s" | egrep ".*.fastq$" 2>${stderrF}))
		if [[ ! -s ${stderrF} ]]; then 
			echo -e "$(date '+%Y_%m_%d %T') [Mapping] fastq files count: ${#fastq_files[@]}" | tee -a ${stdoutF} 2>&1
			echo -e "$(date '+%Y_%m_%d %T') [Mapping] fastq files list: ${fastq_files[@]}" | tee -a ${stdoutF} 2>&1
		else
			exit_on_error	
		fi
		forward_fastq=$(for f in "${fastq_files[@]}"; do 
		m=$(echo $f | egrep "_[0-9]+_1_" 2>${stderrF}); if [[ -n $m ]]; then echo $m; break; fi
		done 2>${stderrF})
		exit_on_error
		echo -e "$(date '+%Y_%m_%d %T') [Mapping] forward fastq: $forward_fastq" | tee -a ${stdoutF} 2>&1
		reverse_fastq=$(for f in "${fastq_files[@]}"; do 
		m=$(echo $f | egrep "_[0-9]+_2_" 2>${stderrF}); if [[ -n $m ]]; then echo $m; break; fi
		done 2>${stderrF})
		exit_on_error
		echo -e "$(date '+%Y_%m_%d %T') [Mapping] reverse fastq: $reverse_fastq" | tee -a ${stdoutF} 2>&1
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
		echo -e "$(date '+%Y_%m_%d %T') [Mapping] current sample name: $CURRENT_SAMPLE_NAME" | tee -a ${stdoutF} 2>&1

	## build mapping cli
		command_name="gsnap"
		echo -e "--- Config section [${command_name}] ---" | tee -a ${stdoutF} 2>&1
		for params in $(set | grep ^$(toupper ${prefix}_${command_name}_)); do
	    	echo -e "params: $params" | tee -a ${stdoutF} 2>&1
		done

		cli_options=($(buildCommandLineOptions $command_name $prefix 2>${stderrF}))
		if [[ ! -s ${stderrF} ]]; then
			res="${cli_options[@]}"	
			gsnap_out=$CURRENT_SAMPLE_NAME\_gsnap_out_s$s.sam
			CURRENT_SAMPLE_ROOT_DIR=$TEST_DATA_ROOT_DIR/$s
			CURRENT_SAMPLE_LOG=$CURRENT_SAMPLE_DIR/$CURRENT_SAMPLE_NAME.log
			CURRENT_SAMPLE_ERR=$CURRENT_SAMPLE_DIR/$CURRENT_SAMPLE_NAME\_err.log
			CURRENT_SAMPLE_PID=$CURRENT_SAMPLE_DIR/$CURRENT_SAMPLE_NAME.pid	
			gsnap_cli="gsnap $res $CURRENT_SAMPLE_ROOT_DIR/$forward_fastq $CURRENT_SAMPLE_ROOT_DIR/$reverse_fastq > $CURRENT_SAMPLE_DIR/$gsnap_out 2>${CURRENT_SAMPLE_ERR} &"
	## run cli
			echo -e "$(date '+%Y_%m_%d %T') [Mapping] Executing mapping command: ${gsnap_cli}" | tee -a ${stdoutF} 2>&1
			eval "$gsnap_cli" 2>${stderrF}
			pid=$!
			echo -e "$(date '+%Y_%m_%d %T') [Mapping] $CURRENT_SAMPLE_NAME pid: $pid" | tee -a ${stdoutF} 2>&1
			echo -e $pid >$CURRENT_SAMPLE_PID 2>>${stderrF}
			pids_arr=("${pids_arr[@]}" "$pid")
			#echo -e "${pids_arr[@]}" # for testing purpose
	
	## copy stdout to log
			cat ${stdoutF} 2>&1 | tee -a $CURRENT_SAMPLE_LOG 2>&1 | tee -a $TEST_PREFIX_LOG >/dev/null
			exit_on_error
		else
			exit_on_error
		fi
	done

	# wait for all gsnap processes to finish
	echo -e "$(date '+%Y_%m_%d %T') [Mapping] gsnap pids count: ${#pids_arr[@]}" | tee ${stdoutF} 2>&1
	echo -e "$(date '+%Y_%m_%d %T') [Mapping] gsnap pids list: ${pids_arr[@]}" | tee -a ${stdoutF} 2>&1
	for p in "${pids_arr[@]}"; do
		echo -e $(ps aux | grep $p | grep $USER | grep -v grep 2>${stderrF})
	done
	exit_on_error
	#waitall "${pids_arr[@]}" 2>${stderrF}
	waitalluntiltimeout "${pids_arr[@]}" 2>${stderrF}
	egrep "exited" ${stderrF} 2>&1 | tee -a ${stdoutF} 2>&1
	cat ${stdoutF} 2>&1 | tee -a $TEST_PREFIX_LOG >/dev/null

	# check for errors
	errs=0
	for s in "${samples[@]}"; do
		CURRENT_SAMPLE_DIR=$TEST_PREFIX_DIR/$s
		CURRENT_SAMPLE_ERR=$CURRENT_SAMPLE_DIR/$(ls "$CURRENT_SAMPLE_DIR" | egrep ".*_err.log$" | grep -v egrep 2>>${stderrF})
		CURRENT_SAMPLE_LOG=$CURRENT_SAMPLE_DIR/$(ls "$CURRENT_SAMPLE_DIR" | egrep -v ".*_err.log$" | egrep ".*.log$" | grep -v egrep 2>>${stderrF})
		CURRENT_SAMPLE_PID=$CURRENT_SAMPLE_DIR/$(ls "$CURRENT_SAMPLE_DIR" | egrep ".*.pid$" | grep -v egrep 2>>${stderrF})

		pid_status=$(egrep "exited" ${stderrF} 2>&1 | egrep $(cat $CURRENT_SAMPLE_PID) 2>&1)
		echo -e "$(date '+%Y_%m_%d %T') [Mapping] pid status: $pid_status" 2>&1 | tee -a $CURRENT_SAMPLE_LOG 2>&1 | tee -a $TEST_PREFIX_ERR 2>&1

		if [[ -z $(echo -e $pid_status 2>&1 | egrep "zero exit status") ]]; then
			echo -e "$(date '+%Y_%m_%d %T') [Mapping] sample $s, non zero exit status:\n$pid_status" | tee -a  $CURRENT_SAMPLE_LOG 2>&1 | tee -a $TEST_PREFIX_ERR 2>&1
			let errs=errs+1
		fi

		if [[ -z $(tail -n 1 "${CURRENT_SAMPLE_ERR}" | egrep -v "^$" | egrep "^Processed" 2>&1) ]]; then
			echo -e "$(date '+%Y_%m_%d %T') [Mapping] sample $s stderr output:" 2>&1 | tee -a $CURRENT_SAMPLE_LOG 2>&1 | tee -a $TEST_PREFIX_ERR 2>&1; 
			cat "${CURRENT_SAMPLE_ERR}" 2>&1 | tee -a $CURRENT_SAMPLE_LOG 2>&1 | tee -a $TEST_PREFIX_ERR 2>&1
			let errs=errs+1
		fi    
		assertTrue "Expected output to stderr for sample $s" "[ -s ${CURRENT_SAMPLE_ERR} ]"
	done

	if [[ $errs == 0 ]]; then
		echo -e "$(date '+%Y_%m_%d %T') [Mapping] gsnap processing all samples without errors." | tee -a ${TEST_PREFIX_LOG} 2>&1
	else
		echo -e "$(date '+%Y_%m_%d %T') [Mapping] some errors occured while gsnap processing samples." | tee -a ${TEST_PREFIX_LOG} 2>&1
		echo -e "$(date '+%Y_%m_%d %T') [Mapping] refer to $TEST_PREFIX_ERR file for details." | tee -a ${TEST_PREFIX_LOG} 2>&1
		exit 1
	fi

	# for each conversion command
	view_command="samtools view"
	view_command_ext="bam"
	sort_command="samtools sort"
	sort_command_ext="sorted.$view_command_ext"
	index_command="samtools index"
	index_command_ext="${sort_command_ext}.bai"
	conversion_commands=("$view_command" "$sort_command" "$index_command")
	for cmd in "${conversion_commands[@]}"; do
		echo -e "$(date '+%Y_%m_%d %T') [Conversion] Processing with $cmd ... " | tee ${stdoutF} 2>&1 | tee -a $CURRENT_SAMPLE_LOG 2>&1 | tee -a $TEST_PREFIX_LOG 2>&1

	# reinitiate pids array 
		pids_arr=()

	## for each sample
		for s in "${samples[@]}"; do
			echo -e "$(date '+%Y_%m_%d %T') [Conversion] Processing sample $s ... " | tee ${stdoutF} 2>&1
			
	## get the sam file
			CURRENT_SAMPLE_DIR=$TEST_PREFIX_DIR/$s
			sam_out=$(ls "$CURRENT_SAMPLE_DIR" | egrep ".*_gsnap_out_s$s.sam$" | grep -v egrep 2>${stderrF})
			if [[ ! -s ${stderrF} ]]; then
				#echo -e "$(date '+%Y_%m_%d %T') [Conversion] gsnap sam output file: $sam_out" | tee -a ${stdoutF} 2>&1 # for testing purpose
				CURRENT_SAMPLE_NAME=$(basename ${sam_out%.*})
				CURRENT_ORI_SAMPLE_NAME=$(echo ${CURRENT_SAMPLE_NAME%_gsnap_out_*})
				CURRENT_SAMPLE_LOG=$CURRENT_SAMPLE_DIR/${CURRENT_ORI_SAMPLE_NAME}.log
				CURRENT_SAMPLE_ERR=$CURRENT_SAMPLE_DIR/${CURRENT_ORI_SAMPLE_NAME}_err.log
				CURRENT_SAMPLE_PID=$CURRENT_SAMPLE_DIR/${CURRENT_ORI_SAMPLE_NAME}.pid
				#echo -e ${CURRENT_ORI_SAMPLE_NAME};
				#echo -e ${CURRENT_SAMPLE_ERR};
			else
				exit_on_error
			fi
	## build conversion cli
			echo -e "--- Config section [${cmd// /_}] ---" | tee -a ${stdoutF} 2>&1
			for params in $(set | grep ^$(toupper ${prefix}_${cmd// /_}_)); do
	    		echo -e "params: $params" | tee -a ${stdoutF} 2>&1
			done
			cli_options=($(buildCommandLineOptions "$cmd" $prefix 2>${stderrF}))
			if [[ ! -s ${stderrF} ]]; then
				res="${cli_options[@]}"
				echo -e "$(date '+%Y_%m_%d %T') [Conversion] $cmd cli options: $res" | tee -a ${stdoutF} 2>&1
				case $cmd in
					"$view_command")
						if [[ -s $CURRENT_SAMPLE_DIR/$sam_out ]]; then
							echo -e "$(date '+%Y_%m_%d %T') [Conversion] gsnap sam output file: $CURRENT_SAMPLE_DIR/$sam_out" | tee -a ${stdoutF} 2>&1		
							cmd_out=$CURRENT_SAMPLE_NAME\.$view_command_ext
							echo -e "$(date '+%Y_%m_%d %T') [Conversion] input: $sam_out" | tee -a ${stdoutF} 2>&1
							echo -e "$(date '+%Y_%m_%d %T') [Conversion] output: $cmd_out" | tee -a ${stdoutF} 2>&1
							cmd_cli="$cmd $res $CURRENT_SAMPLE_DIR/$sam_out > $CURRENT_SAMPLE_DIR/$cmd_out 2>>${CURRENT_SAMPLE_ERR} &"
						else
							debug "$CURRENT_SAMPLE_DIR/$sam_out file does not exist or is empty." | tee -a${CURRENT_SAMPLE_ERR} 2>&1 | tee -a $TEST_PREFIX_ERR 2>&1; exit 1
						fi
					;;
					"$sort_command")
						bam_out=$CURRENT_SAMPLE_NAME\.$view_command_ext
						if [[ -s $CURRENT_SAMPLE_DIR/$bam_out ]]; then
							cmd_out=$CURRENT_SAMPLE_NAME\.sorted
							echo -e "$(date '+%Y_%m_%d %T') [Conversion] input: $bam_out" | tee -a ${stdoutF} 2>&1
							echo -e "$(date '+%Y_%m_%d %T') [Conversion] output: $cmd_out" | tee -a ${stdoutF} 2>&1
							cmd_cli="$cmd $res $CURRENT_SAMPLE_DIR/$bam_out $CURRENT_SAMPLE_DIR/$cmd_out 2>>${CURRENT_SAMPLE_ERR} &"
						else
							debug "$CURRENT_SAMPLE_DIR/$bam_out file does not exist or is empty." | tee -a${CURRENT_SAMPLE_ERR} 2>&1 | tee -a $TEST_PREFIX_ERR 2>&1; exit 1
						fi
					;;	
					"$index_command")
						sorted_bam_out=${CURRENT_SAMPLE_NAME}.$sort_command_ext
						if [[ -s $CURRENT_SAMPLE_DIR/$sorted_bam_out ]]; then
							cmd_out=$CURRENT_SAMPLE_NAME\.$index_command_ext
							echo -e "$(date '+%Y_%m_%d %T') [Conversion] input: $sorted_bam_out" | tee -a ${stdoutF} 2>&1
							echo -e "$(date '+%Y_%m_%d %T') [Conversion] output: $cmd_out" | tee -a ${stdoutF} 2>&1
							cmd_cli="$cmd $res $CURRENT_SAMPLE_DIR/$sorted_bam_out 2>>${CURRENT_SAMPLE_ERR} &"
						else
							debug "$CURRENT_SAMPLE_DIR/$sorted_bam_out file does not exist or is empty." | tee -a${CURRENT_SAMPLE_ERR} 2>&1 | tee -a $TEST_PREFIX_ERR 2>&1; exit 1
						fi
					;;
					*)
						debug "unexpected $cmd, it will not be processed." | tee -a${CURRENT_SAMPLE_ERR} 2>&1 | tee -a $TEST_PREFIX_ERR 2>&1
				esac
	## run cli
				if [[ -n $cmd_cli ]]; then
					echo -e "$(date '+%Y_%m_%d %T') [Conversion] Executing command: ${cmd_cli}" | tee -a ${stdoutF} 2>&1
					eval "$cmd_cli" 2>&1  
					pid=$!
					echo -e "$(date '+%Y_%m_%d %T') [Conversion] $CURRENT_SAMPLE_NAME pid: $pid" | tee -a ${stdoutF} 2>&1
					echo -e $pid >$CURRENT_SAMPLE_PID 2>${stderrF}
					pids_arr=("${pids_arr[@]}" "$pid")
					#echo -e "$(date '+%Y_%m_%d %T') [Conversion] $cmd pids count: ${#pids_arr[@]}" | tee -a ${stdoutF} 2>&1 # for testing purpose
					#echo -e "$(date '+%Y_%m_%d %T') [Conversion] $cmd pids list: ${pids_arr[@]}" | tee -a ${stdoutF} 2>&1 # for testing purpose

				else
					debug "cmd cli is null. cmd name, $cmd, was not processed" | tee -a ${CURRENT_SAMPLE_ERR} 2>&1 | tee -a $TEST_PREFIX_ERR 2>&1; exit 1
				fi

	## copy stdout to logs
				cat ${stdoutF} 2>&1 | tee -a $CURRENT_SAMPLE_LOG 2>&1 | tee -a $TEST_PREFIX_LOG >/dev/null			
			else
				exit_on_error
			fi
		done
	
	# wait for all conversion cmd processes to finish then run the next cmd
		echo -e "$(date '+%Y_%m_%d %T') [Conversion] $cmd pids count: ${#pids_arr[@]}" | tee ${stdoutF} 2>&1
		echo -e "$(date '+%Y_%m_%d %T') [Conversion] $cmd pids list: ${pids_arr[@]}" | tee -a ${stdoutF} 2>&1		
		for p in "${pids_arr[@]}"; do
			echo -e $(ps aux | grep $p | grep $USER | grep -v grep 2>${stderrF})
		done
		exit_on_error
		#waitall "${pids_arr[@]}" 2>${stderrF}
		waitalluntiltimeout "${pids_arr[@]}" 2>${stderrF}
		egrep "exited" ${stderrF} 2>&1 | tee -a ${stdoutF} 2>&1
		cat ${stdoutF} 2>&1 >>${TEST_PREFIX_LOG}
		#cat ${stderrF} 2>&1 >>${TEST_PREFIX_ERR}

	# check for errors
		errs=0
		for s in "${samples[@]}"; do
			CURRENT_SAMPLE_DIR=$TEST_PREFIX_DIR/$s
			CURRENT_SAMPLE_ERR=$CURRENT_SAMPLE_DIR/$(ls "$CURRENT_SAMPLE_DIR" | egrep ".*_err.log$" | grep -v egrep 2>>${stderrF})
			CURRENT_SAMPLE_LOG=$CURRENT_SAMPLE_DIR/$(ls "$CURRENT_SAMPLE_DIR" | egrep -v ".*_err.log$" | egrep ".*.log$" | grep -v egrep 2>>${stderrF})
			CURRENT_SAMPLE_PID=$CURRENT_SAMPLE_DIR/$(ls "$CURRENT_SAMPLE_DIR" | egrep ".*pid$" | grep -v egrep 2>>${stderrF})

			pid_status=$(egrep "exited" ${stderrF} 2>&1 | egrep $(cat $CURRENT_SAMPLE_PID) 2>>${stderrF})
			echo -e "$(date '+%Y_%m_%d %T') [Conversion] pid status: $pid_status" 2>&1 | tee -a $CURRENT_SAMPLE_LOG 2>&1 | tee -a $TEST_PREFIX_ERR 2>&1
	
			if [[ -z $(echo -e $pid_status 2>&1 | egrep "zero exit status") ]]; then
				echo -e "$(date '+%Y_%m_%d %T') [Conversion] sample $s, non zero exit status:\n$pid_status" | tee -a  $CURRENT_SAMPLE_LOG 2>&1 | tee -a $TEST_PREFIX_ERR 2>&1
				let errs=errs+1
			fi		

			if [[ -z $(tail -n 1 ${CURRENT_SAMPLE_ERR} | egrep -v "^$" | egrep "^\[samopen\]") ]]; then
				echo -e "$(date '+%Y_%m_%d %T') [Conversion] sample $s stderr output:" 2>&1 | tee -a ${CURRENT_SAMPLE_LOG} 2>&1 | tee -a $TEST_PREFIX_ERR 2>&1; cat ${CURRENT_SAMPLE_ERR} 2>&1 | tee -a $CURRENT_SAMPLE_LOG 2>&1 | tee -a $TEST_PREFIX_ERR 2>&1
				let errs=errs+1
			fi    
			#assertFalse "Unexpected output to stderr for sample $s" "[ -s ${CURRENT_SAMPLE_ERR} ]"
		done

		if [[ $errs == 0 ]]; then
			echo -e "$(date '+%Y_%m_%d %T') [Conversion] $cmd processing all samples without errors." | tee -a ${TEST_PREFIX_LOG} 2>&1
		else
			echo -e "$(date '+%Y_%m_%d %T') [Conversion] some errors occured while $cmd processing samples." | tee -a ${TEST_PREFIX_LOG} 2>&1
			echo -e "$(date '+%Y_%m_%d %T') [Conversion] refer to $TEST_PREFIX_ERR file for details." | tee -a ${TEST_PREFIX_LOG} 2>&1
			exit 1
		fi
	done

	# clean
	## unset environment variables with used namespace
	echo -ne "$(date '+%Y_%m_%d %T') [Cleaning] Unsetting all environment variables using namespace: $NAMESPACE ... " | tee -a ${TEST_PREFIX_LOG} 2>&1
	unset $(set | awk -F= -v prefix="${prefix}" 'BEGIN { 
	          prefix = toupper(prefix);
	       }
	       /^prefix_/  { print $1 }') $(toupper ${prefix}_)
	echo -e "done" | tee -a ${TEST_PREFIX_LOG} 2>&1

	# test for assertions: todo
}

#-------------------------------
# testFailedIsCpuAvailable
#
testFailedIsCpuAvailable()
{
	if [[ $(isCpuAvailable 2 2 2>${stderrF}) == "TRUE" ]]; then
		if [[ "$?" == 0 ]]; then
			echo -e "isCpuAvailable: TRUE was returned" > ${stdoutF}
		else
			echo -e "isCpuAvailable: non zero exit status" >> ${stderrF}
		fi
	else
		echo -e "isCpuAvailable: FALSE was returned." >> ${stderrF}
	fi

	assertTrue "Expected output to stdout." "[ -s ${stdoutF} ]"
	assertFalse "Unexpected output to stderr." "[ -s ${stderrF} ]"
	[[ -s ${stderrF} ]] && (echo -e "stderr output:"; cat ${stderrF} 2>&1) 
}

#-------------------------------
# testFailedIsDiskSpaceAvailable
#
testFailedIsDiskSpaceAvailable()
{
	if [[ $(isDiskSpaceAvailable $PWD 2 $TEST_DATA_ROOT_DIR/DG $TEST_DATA_ROOT_DIR/DI 2>${stderrF}) == "TRUE" ]]; then
		if [[ "$?" == 0 ]]; then
			echo -e "isDiskSpaceAvailable: TRUE was returned" > ${stdoutF}
		else
			echo -e "isDiskSpaceAvailable: non zero exit status" >> ${stderrF}
		fi
	else
		echo -e "isDiskSpaceAvailable: FALSE was returned." >> ${stderrF}
	fi

	assertTrue "Expected output to stdout." "[ -s ${stdoutF} ]"
	assertFalse "Unexpected output to stderr." "[ -s ${stderrF} ]"
	[[ -s ${stderrF} ]] && (echo -e "stderr output:"; cat ${stderrF} 2>&1)
}

#----------------------------------------------
# testFailedIsDiskSpaceAvailableForSubdirsList
#
testFailedIsDiskSpaceAvailableForSubdirsList()
{
	complete_dataset=($(ls -d -1 $COMPLETE_DATA_ROOT_DIR/*))
	echo -e "partition: $(df -k $PWD | tail -n 1 | awk '{print $4}')"
	echo -e "data: $(du -sck "${complete_dataset[@]}" | tail -n 1 | awk '{print $1}')"
	echo -e "expanded data: $[$(du -sck "${complete_dataset[@]}" | tail -n 1 | awk '{print $1}')*2]"
	res=$(isDiskSpaceAvailable $PWD 2 "${complete_dataset[@]}" 2>${stderrF})
	rtrn=$?
	echo "*** isDiskSpaceAvailable exit status: $rtrn ***"

	if [[ $res == "TRUE" ]]; then
		if [[ "$rtrn" == 0 ]]; then
			echo -e "isDiskSpaceAvailable: TRUE was returned" > ${stdoutF}
		else
			echo -e "isDiskSpaceAvailable: TRUE and non zero exit status" >> ${stderrF}
		fi
	else
		if [[ "$rtrn" == 0 ]]; then
			echo -e "isDiskSpaceAvailable: FALSE was returned." > ${stdoutF}
		else
			echo -e "isDiskSpaceAvailable: FALSE and non zero exit status" >> ${stderrF}		
		fi
	fi

	assertTrue "Expected output to stdout." "[ -s ${stdoutF} ]"
	assertFalse "Unexpected output to stderr." "[ -s ${stderrF} ]"
	[[ -s ${stderrF} ]] && (echo -e "stderr output:"; cat ${stderrF} 2>&1)
}

#-------------------------------
# testFailedSendEmail
#
testFailedSendEmail()
{
	recipient="Joseph.Tran@versailles.inra.fr"
	subject_noattach="Test sendEmail without attachment files"
	subject_attach="Test sendEmail with attachment files"
	subject_error="Test sendEmail with error in subject"
	message="This is a test message."
	files=($(ls -d -1 $TEST_DATA_ROOT_DIR/DG/*))

	# without attachment files
	sendEmailBasic $recipient "$subject_noattach" "$message"

	# with attachment files
	sendEmailBasic $recipient "$subject_attach" "$message" "${files[@]}"

	# with error in subject
	sendEmailBasic $recipient "$subject_error" "$message"

	# assertions: TODO
}

#-------------------------------
# testFailedBuildTimeCmd
#
testFailedBuildTimeCmd()
{
	cmd="date"
	buildTimeCmd $cmd $TEST_OUTPUT_DIR >${stdoutF} 2>${stderrF}

	assertTrue "Expected output to stdout." "[ -s ${stdoutF} ]"
	[[ -s ${stdoutF} ]] && (echo -e "stdout output:"; cat ${stdoutF} 2>&1)
	assertFalse "Unexpected output to stderr." "[ -s ${stderrF} ]"
	[[ -s ${stderrF} ]] && (echo -e "stderr output:"; cat ${stderrF} 2>&1)

	# execute
	timeCmd=$(cat ${stdoutF} 2>${stderrF})
	eval "$timeCmd date" 2>${stderrF}

	assertTrue "Expected output to log file: $TEST_OUTPUT_DIR/${cmd// /_}_time.log." "[ -s $TEST_OUTPUT_DIR/${cmd// /_}_time.log ]"
	[[ -s $TEST_OUTPUT_DIR/${cmd// /_}_time.log ]] && (echo -e "log output:"; cat $TEST_OUTPUT_DIR/${cmd// /_}_time.log 2>&1)
	assertFalse "Unexpected output to stderr." "[ -s ${stderrF} ]"
	[[ -s ${stderrF} ]] && (echo -e "stderr output:"; cat ${stderrF} 2>&1)
}

#-------------------------------
# testFailedExitOnTimeout
#
testFailedExitOnTimeout()
{
	sleep 20 &
	pid=$!
	timeout=3
	while kill -0 $pid 2>${stderrF}; do
		echo "before timeout: $timeout"
		timeout=$(exit_on_timeout $pid "$TEST_OUTPUT_DIR/kill_$pid.log" $timeout 1 1 2>>${stderrF})
		echo "after timeout: $timeout"	
	done
	
	assertTrue "Expected output to stdout." "[ -s $TEST_OUTPUT_DIR/kill_$pid.log ]"
	[[ -s $TEST_OUTPUT_DIR/kill_$pid.log ]] && (echo -e "kill log output:"; cat $TEST_OUTPUT_DIR/kill_$pid.log 2>&1)
	assertTrue "Expected output to stderr." "[ -s ${stderrF} ]"
	[[ -s ${stderrF} ]] && (echo -e "stderr output:"; cat ${stderrF} 2>&1)
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

	COMPLETE_DATA_ROOT_DIR="/data/temp_projects/AZM/SEQUENCES"
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
