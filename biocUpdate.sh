#!/bin/bash

# Script that automizes copying into Bioconductor repository.
# THE ASSUMPTION IS that there were no updates made to Bioc sources, thus can be just replaced by sources from C++ version.


if [ $# -ne 1 ]
then
   echo "Usage: provide path to bioc sources directory (.../devel/src)"
   echo "  biocUpdate.sh [dirPath]"
   exit
fi

for i in `ls $1`
do
   if [  -e $i  -a  -f $i  ]
   then
      other="${1}/${i}"
      if ! diff -q $i $other > /dev/null
      then
         cp -v $i $other;
      fi
   fi
done
