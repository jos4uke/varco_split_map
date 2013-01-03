#! /bin/bash

#
# VARCO_CP_SEQ_DATASET.SH
#

# Copyright 2012 Joseph Tran <Joseph.Tran@versailles.inra.fr>

# This software is a computer program whose purpose is to:
# - perform reads mapping in batch mode by splitting given samples in several batches
#   and run sequentially each batch to avoid cpu overload and running out of disk space.
# - the following code is a intended to generate test dataset from complete RNA-Seq dataset:
#   it copies a subset of AZM/SEQUENCES to test/data selecting only the 10000 first reads /fastq.

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





