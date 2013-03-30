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

#svn export asa103 $DIR/asa103
svn export boost $DIR/boost
svn export samtools $DIR/samtools

cp -v _release_Makefile $DIR/Makefile

cp -v $( cat releaseList ) $DIR

echo "==================" >> $DIR/README
date >> $DIR/README
svn info | grep -e "^Revision:" >> $DIR/README

echo "REMINDERs:"
echo "Delete licensing line about asa libs from README."
echo "File Makefile contains current version of BitSeq, please update if haven't done already."
