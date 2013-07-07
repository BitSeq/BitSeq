#!/usr/bin/python

from __future__ import print_function

"""
Due to an error in our code, GENE EXPRESSION and WITHIN GENE EXPRESSION results
produced with BitSeq versions up to 0.5.3 might be wrong for some genes.

This only applies to SOME genes (or transcripts of that genes) IF the transcripts
in .tr file were not grouped by gene.
This program can check exactly which genes have wrong expression estimates.

The results of other genes and their transcripts are correct.
"""

import sys

def checkFile(fileName):
   try:
      inf = open(fileName);
   except:
      sys.exit("Can't open file "+str(fileName));
   genesSeen = set();

   giName = ""
   prevName = "";

   genesWrong = list();
   g2ts ={};

   for line in inf:
      if line[0]=='#':continue;
      try:
         gn = line.split()[0];
         tr = line.split()[1];
      except:
         sys.exit("The TR does not seem to have valid format in line:\n"+line);

      if gn in g2ts:
         g2ts[gn].append(tr);
      else:
         g2ts[gn] = [tr];
         
      if gn == prevName:
         if prevName != giName:
            if giName not in genesWrong: genesWrong.append(giName);
            if prevName not in genesWrong: genesWrong.append(prevName);
      else:
         if gn not in genesSeen:
            prevName=gn;
            giName=gn;
            genesSeen.add(gn);
         else:
            giName=gn;
   if len(genesWrong) == 0:
      print("Everything seems to be fine.")
   else:
      print("These",len(genesWrong),"( out of",len(genesSeen),") have wrong GENE EXPRESSION results:")
      trCount = 0;
      for it in genesWrong:
         print(it, end=' ')
         trCount+=len(g2ts[it]);
      print("\nThese:",trCount,"transcripts have wrong WITHIN GENE EXPRESSION results:");
      for it in genesWrong:
         for trit in g2ts[it]:
            print(trit, end=' ')
      print() ;

if __name__ == "__main__":
   if len(sys.argv) <2:
      sys.exit("Please provide file name as argument.");
   print("Checking file",sys.argv[1]);
   checkFile(sys.argv[1]);


