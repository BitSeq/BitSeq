#!/usr/bin/python
# Initialization {{{
import sys
from optparse import OptionParser
parser = OptionParser(usage="%prog [options] <inputFile> <outputFile>\n\n\
      This program extracts information about transcripts from reference Fasta file.
      This is partially replaced by using SAM header, which however does not include information about transcript-gene grouping.\n\
      Current version of parseAlignment extracts this information from a reference sequence file (making this script obsolete).\
")
parser.add_option("-v", "--verbose", default=False, dest="verbose",  action="store_true", help="Verbose output")
parser.add_option("-t","--type",dest="type", type="string",help="Type of file to parse: ensembl, other");

(options, args) = parser.parse_args()
def verbose(str):
   if options.verbose:
      print str;

if len(args)<2: 
   sys.exit("Missing arguments");

try:
   inF = open(args[0],"r");
except:
   sys.exit("Unable to open input file: "+args[0]+" .");


try:
   outF = open(args[1],"w");
except:
   sys.exit("Unable to open output file: "+args[1]+" .");
#}}}

seqName="";
geneName="";
seqLen=0;
seqCount=0;

result = [];

if options.type and options.type=="ensembl":
   print "Expecting header line format:\n>[tr Name] .* gene:[gene Name] .*";
   for line in inF:
      if line[0] == '>':
         if seqName!="":
            result.append([geneName,seqName,str(seqLen)]);
         seqLen=0;
         seqName = line.split()[0][1:];
         for it in line.split():
            if "gene:" in it:
               geneName=it[5:];
         seqCount+=1;
      else: seqLen+=len(line)-1;
   if seqName!="":
      result.append([geneName,seqName,str(seqLen)]);
else:
   print "Expecting header line format:\n>[tr Name] .*\n -> using \"none\" as gene names";
   for line in inF:
      if line[0] == '>':
         if seqName!="":
            result.append(["none",seqName,str(seqLen)]);
         seqLen=0;
         seqName = line.split()[0][1:];
         seqCount+=1;
      else: seqLen+=len(line)-1;
   if seqName!="":
      result.append(["none",seqName,str(seqLen)]);
inF.close();

verbose(str(seqCount)+" sequences processed.");

outF.write("# M "+str(seqCount)+"\n");
for it in result:
   outF.write(it[0]+" "+it[1]+" "+it[2]+"\n");

outF.close();

