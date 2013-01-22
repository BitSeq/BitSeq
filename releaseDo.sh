#!/bin/bash

if [ $# -ne 1 ]
then
   echo "reselaseDo.sh [dirName]"
   exit
fi

DIR="/localhome/work/$1"

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
