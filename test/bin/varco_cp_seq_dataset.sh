#! /bin/bash

#
# copy a subset of AZM/SEQUENCES to test/data
# selecting only the 10000 first reads /fastq
#

DATA_ROOT_DIR=/data/temp_projects/AZM/SEQUENCES
ERROR_TMP=../output/error_cp_seq_dataset.log
DEST_DIR=../data/SEQUENCES

exit_on_error()
{
	[[ -s $ERROR_TMP ]] && $(cat "error log: " $ERROR_TMP; exit 1)
}

# get all subdirs from data root dir
all_subdirs=($(find $DATA_ROOT_DIR -maxdepth 1 -type d | egrep -v ^$DATA_ROOT_DIR$ 2>$ERROR_TMP))
exit_on_error

# create the destination dir
if [[ ! -d $DEST_DIR ]]; then
	mkdir $DEST_DIR 2>$ERROR_TMP
	exit_on_error
fi

# Iterate over the data root dir subdirs
# for each subdir
# create the sample subdir in destination dir
# get the fastq files in the data root subdir
# copy the first 10000 reads/fastq to the destination sample subdir

for subdir in "${all_subdirs[@]}"; do
# create sample subdir in destination dir
CURRENT_SUBDIR=$DEST_DIR/$(basename $subdir)
if [[ ! -d $CURRENT_SUBDIR ]]; then
	mkdir $CURRENT_SUBDIR 2>$ERROR_TMP
	exit_on_error
fi

# get the fastq files in the data root subdir
fastq_files=($(ls "$subdir" | egrep -v "_single_" | egrep ".*.fastq$" 2>$ERROR_TMP))
exit_on_error
if [[ "${#fastq_files[@]}" == 2 ]]; then
	echo -e "$(date '+%Y_%m_%d %T') [Fastq files] INFO $subdir, has 2 fastq files: ${fastq_files[@]}" 2>&1
	for fastq in "${fastq_files[@]}"; do
		# copy the first 10000 reads (40000 first lines)/fastq to the destination sample subdir
		head -40000 $subdir/$fastq > $CURRENT_SUBDIR/${fastq%.*}_10000.fastq 2>$ERROR_TMP
		exit_on_error
	done
	fastq_files_subset=($(ls "$CURRENT_SUBDIR" | egrep -v "_single_" | egrep ".*.fastq$" 2>$ERROR_TMP))
	echo -e "fastq files subset: ${fastq_files_subset[@]}"
elif [[ "${#fastq_files[@]}" < 2 ]]; then
	echo -e "$(date '+%Y_%m_%d %T') [Fastq files] WARNING $subdir, has less than 2 fastq files: ${fastq_files[@]}" 2>&1 | tee -a $ERROR_TMP 2>&1
	exit_on_error
elif [[ "${#fastq_files[@]}" > 2 ]]; then
	echo -e "$(date '+%Y_%m_%d %T') [Fastq files] WARNING $subdir, has more than 2 fastq files: ${fastq_files[@]}" 2>&1 | tee -a $ERROR_TMP 2>&1
	exit_on_error
fi
done





