#!/usr/bin/env bash

#
# INSTALL.SH
#

# Copyright 2012 Joseph Tran <Joseph.Tran@versailles.inra.fr>

# This software is a computer program whose purpose is to:
# - perform reads mapping in batch mode by splitting given samples in several batches
#   and run sequentially each batch to avoid cpu overload and running out of disk space.
# - the following code is intended to install and deploy this software to production server.

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


ARGV=("$@")
ARGC=$#

PROGRAM="varco_split_map.sh"
ERROR_LOG=/tmp/${PROGRAM%.*}_install_error.log

if [[ $ARGC -ne 2 ]]; then
    echo -e "Incorrect number of required positional arguments. Should be equal to 2: \n- the first argument corresponding to the tag version of $PROGRAM pipeline in git repository, \n- and the second and last argument to the prefix used to deploy the pipeline." | tee $ERROR_LOG 2>&1
    #usage;
    exit 1
else
	echo
	echo
	echo "************************************************************"
	echo "*             STARTING  INSTALLATION                       *"
	echo "************************************************************"
	echo
	echo
    echo -e "$(date '+%F %T') Starting the $PROGRAM pipeline deployment" | tee $ERROR_LOG 2>&1
    PROGRAM_VERSION=${ARGV[0]}
    echo "$(date '+%F %T') Will deploy version tag: $PROGRAM_VERSION" | tee -a $ERROR_LOG 2>&1
    PREFIX=${ARGV[1]}
	echo "$(date '+%F %T') Will deploy $PROGRAM pipeline in: $PREFIX" | tee -a $ERROR_LOG 2>&1
fi 

# export $PROGRAM pipeline from git repo 
PROGRAM_VERSION_DIR=/usr/local/archives/${PROGRAM%.*}/$PROGRAM_VERSION
mkdir -p $PROGRAM_VERSION_DIR 2>$ERROR_LOG
rtrn=$?
if [[ rtrn -ne 0 ]]; then
    echo -e "$(date '+%F %T') $PROGRAM_VERSION_DIR directory creation failed with exit error code $rtrn." | tee -a $ERROR_LOG 2>&1
    echo -e "$(date '+%F %T') You can get more information about $PROGRAM pipeline installation in $ERROR_LOG." 2>&1
    exit $rtrn
else
    echo -e "$(date '+%F %T') $PROGRAM_VERSION_DIR directory was created successfully." | tee -a $ERROR_LOG 2>&1
    echo -e "$(date '+%F %T') Will export $PROGRAM pipeline project from git repository." | tee -a $ERROR_LOG 2>&1
fi
git archive --prefix=${PROGRAM%.*}_$PROGRAM_VERSION$PREFIX/ $PROGRAM_VERSION | (cd $PROGRAM_VERSION_DIR && tar xf -)
rtrn=$?
if [[ rtrn -ne 0 ]]; then
    echo -e "$(date '+%F %T') git archive failed to export $PROGRAM pipeline project with exit error code $rtrn." | tee -a $ERROR_LOG 2>&1
    echo -e "$(date '+%F %T') You can get more information about $PROGRAM pipeline installation in $ERROR_LOG." 2>&1
    exit $rtrn
else
    echo -e "$(date '+%F %T') $PROGRAM_VERSION_DIR directory was created successfully." | tee -a $ERROR_LOG 2>&1
    echo -e "$(date '+%F %T') Will deploy $PROGRAM pipeline files to $PREFIX" | tee -a $ERROR_LOG 2>&1
fi

cd  $PROGRAM_VERSION_DIR/${PROGRAM%.*}_$PROGRAM_VERSION/ 2>$ERROR_LOG
rm -rf usr/local/test usr/local/install.sh 2>>$ERROR_LOG
for f in $(find . -name \* -print 2>>$ERROR_LOG); do
    echo -e "$(date '+%F %T') cp $f in /$f" | tee -a $ERROR_LOG 2>&1
    cp --parents $f / 2>>$ERROR_LOG
done

# set privileges
echo -ne "$(date '+%F %T') Setting executable privileges on $PROGRAM script ..." | tee -a $ERROR_LOG 2>&1
chmod 755 $PREFIX/bin/$PROGRAM 2>>$ERROR_LOG
echo -e "OK" | tee -a $ERROR_LOG 2>&1

cat $ERROR_LOG 2>&1

echo
echo "NOTE: You can get more information about $PROGRAM pipeline installation in $ERROR_LOG."
echo
echo
echo "************************************************************"
echo "*             INSTALLATION IS NOW COMPLETE                 *"
echo "************************************************************"
echo
echo
echo "You can now type $PROGRAM to get usage help for the pipeline."
echo
