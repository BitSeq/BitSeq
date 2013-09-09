#!/bin/bash

# Script that automizes creating new BitSeq release.
# Copies relevant files listed in releaseList, uses _release_Makefile as new Makefile
# (make sure it's correct) and copies directories boost and samtools.

if [ $# -ne 1 ]
then
   echo "reselaseDo.sh [dirName]"
   exit
fi

DIR=$1

if [ -d $DIR ]
then
   echo "Direcotry $DIR already exists!";
   exit
fi

mkdir $DIR

# Cleanup:
make clean-all

#svn export asa103 $DIR/asa103
if [[ -d .svn ]]
then
   svn export asa103 $DIR/asa103
   svn export boost $DIR/boost
   svn export samtools $DIR/samtools
else
   echo "Copying boost to '$DIR'."
   cp -r boost $DIR
   cp -rv asa103 boost samtools $DIR
fi

cp -v _release_Makefile $DIR/Makefile

cp -v $( cat releaseList ) $DIR

echo "==================" >> $DIR/README
date >> $DIR/README
if [[ -d .svn ]]
then
   svn info | grep -e "^Revision:" >> $DIR/README
else
   git log -1 | grep "commit" >> $DIR/README
fi

echo "REMINDERs:"
echo "File Makefile contains current version of BitSeq, please update if haven't done already."
