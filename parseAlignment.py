#!/usr/bin/python
# Initialization {{{
import sys
import numpy as np
def normpdf(x,m,s):
   return 1./(s*2.5066282746310002)*np.exp(-1./(2.0*s*s)*(x-m)**2.)
import os, time # needed for this:
time_str = time.strftime("%b %e %Y %H:%M:%S", time.gmtime(os.lstat(sys.argv[0]).st_mtime));
print "###",os.path.basename(sys.argv[0]),"build:",time_str;
# {{{ parse arguments and set filenames
from optparse import OptionParser
parser = OptionParser(usage="%prog [options]\n -a -t are necessary\n -e is adviced")
parser.add_option("-T", "--transcriptPrefix", dest="tPref", help="Prefix of transcript names within MAP file (e.g. hg19_ensGene_ for ensembl genes from UCSC)", type="string")
parser.add_option("-p", "--prefix", dest="pref", help="Experiment prefix, use same prefix for all files (.map, .tr, .prob)", type="string")
parser.add_option("-a", "--alignmentFile", dest="aFile", help="Alignments file name", type="string")
parser.add_option("-A", "--alignmentFileType", dest="aType", default="bowtie", help="Alignments file type", type="string")
parser.add_option("-t", "--transcriptFile", dest="tFile", help="File with with list of transcripts (second column) and their lengths (third column, used later).", type="string")
parser.add_option("-o", "--out", dest="oFile", help="Output name (should end with .prob).", type="string")
parser.add_option("-N", "--totalN", dest = "totalN", help="Total number of reads. If <name>.map.bowtieLog does not exist this number has to be provided", type="int")
parser.add_option("-i", "--inputType", dest = "inputType", help="Input file type determines the assignemnt of probability for each read (fastq, fastq33, fasta)", default="fastq");
parser.add_option("-v", "--verbose", default=False, dest="verbose",  action="store_true", help="Verbose output")
parser.add_option("--vv", default=False, dest="veryVerbose",  action="store_true", help="Very verbose output")
parser.add_option("--paired", default=False, dest="paired", action="store_true", help="Flag fo paired alignemnts")
parser.add_option("--IamSure", default=False, dest="amSure",  action="store_true", help="I am sure I want to use this.")


(options, args) = parser.parse_args()

if not options.amSure:
   sys.exit("Please use new implementation of parsing algorithm \"parseAlignment\". If you really want to use this program use the option --IamSure.");


if options.tPref !=None:
   prefixL = len(options.tPref);
else:
   prefixL = 0;

if options.pref :
   aFileName=options.pref+".map"
   oFileName=options.pref+".prob"
   tFileName=options.pref+".tr"
else:
   if not options.aFile:
      sys.exit("Need alignemnt file name.");
   if not options.oFile:
      sys.exit("Need output file name.");
   if not options.tFile:
      sys.exit("Need transcript file name.");
if options.aFile:
   aFileName=options.aFile
if options.oFile:
   oFileName=options.oFile
if options.tFile:
   tFileName=options.tFile;
#}}}
#{{{ get total number of reads, possibly from <file>.map.bowtieLog
Ntotal = 0
if options.totalN :
   Ntotal = options.totalN;
else:
   try:
      bLog = open(aFileName+".bowtieLog");
      for line in bLog:
         if line.find("# reads processed:")>-1:
            Ntotal = int( line[line.find("# reads processed:")+18:].split()[0] ); 
            # in other words take first wor after "reads processed:" and convert it to Ntotal
            break;
      bLog.close();
      if Ntotal <= 0: 
         sys.exit("File read, but Ntotal was "+str(Ntotal));
   except:
      sys.exit( "Was not able to read file "+aFileName+".bowtieLog . Please provide number of reads (-N atribute) or the log file.")
#}}}
def nuc2i(str):#{{{
   if str.lower() == "a": return 0;
   if str.lower() == "c": return 1;
   if str.lower() == "g": return 2;
   if str.lower() == "t": return 3;
   return 4;
#}}}
def verbose(str):#{{{
   if options.verbose:
      print str;
#}}}
verbose("Using files:\n   "+aFileName+" for reading alignments\n   "+oFileName+" for writing probabilities\n   "+tFileName+" for writing transcript info");
# {{{ reading transcript info
try:
   tFile = open(tFileName,"r")
except:
   sys.exit("Unable to open transcript file: "+tFileName+" .");

trMap=dict()
i=0;
for line in tFile:
   if line[0] == '#': continue;
   trMap[line.split()[1]]=i+1;
   #trMap[line.split()[1][prefixL:]]=i+1;
   i+=1;
trN=i;
tFile.close();
#}}}
# {{{ open output file
try:
   oFile = open(oFileName,"w");
except:
   sys.exit("Unable to open output file: "+oFileName+" .");
#}}}
#{{{ open alignment file and check number of columns
if options.aType != "bowtie":
   sys.exit("Unrecognized alignment type.");
try:
   aFile = open(aFileName,"r")
except:
   sys.exit("Unable to open alignments file: "+aFileName+" .");

alignment=aFile.readline().rstrip().split("\t");
columnN=len(alignment)+1; # expect no mismatch info
try:
   x = int(alignment[columnN-2]); # this works if last column is NOT mismatch info
except:
   columnN -= 1; # otherwise decrease number of columns
colS = columnN - 8; # if 8 columns, no shift necessary 
verbose("columns: "+str(columnN));
aFile.seek(0);
#}}}
# }}}

if options.inputType=="fasta": #{{{
   minReadLength=25;
   pseudoCount = 1.0;
   nucProb = [[[pseudoCount for i in range(5)] for k in range(5)] for j in range(minReadLength)];
   noiseProb = [pseudoCount for i in range(5)];

   verbose("Estimating mismatch probability.");  # {{{
   readId=""
   mismatch=""
   hadMismatches=True;
   readN = 0;
   verbose("Use all reads, not only unique.");
   for line in aFile:
      alignment=line.rstrip().split("\t");

      readSeq=alignment[4+colS]
      if alignment[1+colS]=="-":
         readSeq = readSeq[::-1];
      
      if alignment[0] != readId or readSeq != seq:
         readId=alignment[0];
         readN+=1;
         if not hadMismatches:
            while len(seq) > len(nucProb):
                  nucProb.append([[pseudoCount for i in range(5)] for k in range(5)]);
            for i in range(len(seq)):
               nuc1 = nuc2i(seq[i]);
               nucProb[i][nuc1][nuc1]+=1;

         hadMismatches=False;
         seq = readSeq;
         mismatch=""
         for nuc in seq:
            noiseProb[nuc2i(nuc)]+=1

      if len(alignment)==columnN:
         if alignment[columnN-1] != mismatch:
            while len(seq) > len(nucProb):
                  nucProb.append([[pseudoCount for i in range(5)] for k in range(5)]);
            for i in range(len(seq)):
               nuc1 = nuc2i(seq[i]);
               nucProb[i][nuc1][nuc1]+=1;
            hadMismatches = True;

            mismatch=alignment[columnN-1]
            mismatchArray = mismatch.split(",");
            for mis in mismatchArray:
               pos = int( mis.split(":")[0] );
               nuc1 = nuc2i( mis.split(":")[1].split(">")[0] );
               nuc2 = nuc2i( mis.split(":")[1].split(">")[1] );
   #            while pos >= len(nucProb):
   #               nucProb.append([[pseudoCount for i in range(5)] for k in range(5)]);
               nucProb[pos][ nuc2 ][ nuc2 ]-=1;
               if nucProb[pos][nuc2][nuc2]<1 : print pos,nuc2,seq,mismatch;
               nucProb[pos][ nuc1 ][ nuc2 ]+=1;
   # }}}
   """verbose("Using only unique reads");#{{{ 
for line in aFile:
   alignment=line.split();
   if alignment[4] != seq:
      seq=alignment[4]

      if mismatch != "":
         mismatchArray = mismatch.split(",");
         for mis in mismatchArray:
            pos = int( mis.split(":")[0] );
            nuc1 = mis.split(":")[1].split(">")[0];
            nuc2 = mis.split(":")[1].split(">")[1];
            while pos <= len(nucProb):
               nusProb.append([[pseudoCount for i in range(5)] for k in range(5)]);
            nucProb[pos][ nuc2i(nuc1) ][ nuc2i(nuc2) ]+=1;
      if len(alignment>7):
         mismatch=alignment[7];
      
      for nuc in seq:
         noiseProb[nuc2i(nuc)]+=1
   else:
      mismatch=""
#}}}"""
   verbose("Estimating probability of noise from aligned reads.") #{{{
   total=sum(noiseProb);
   for i in range(5):
      noiseProb[i] /= total;

   verbose("Estimating nucleotide mismatch matrix.");
   for i in range(len(nucProb)):
      for j in range(5):
         total = sum( nucProb[i][j] );
         for k in range(5):
            nucProb[i][j][k] /= total;

   if options.veryVerbose:
      print "Noise probabilities: ";
      print "   ",;
      print noiseProb;
      print "Nucleotide mismatch matrix:";
      for i in range(len(nucProb)):
         print "Position ",i,":\n   ",;
         print nucProb[i];
   #}}}
   verbose("Writing alignment probabilities"); # {{{
   aFile.seek(0);

   alignment=aFile.readline().rstrip().split("\t");
   readId=alignment[0];
   if alignment[1+colS] == "+":
      seq=alignment[4+colS];
   else:
      seq=alignment[4+colS][::-1];
   prob = 1.0;
   for nuc in seq:
      prob *= noiseProb[nuc2i(nuc)];
   alignments=[(0,alignment[1+colS],prob)];

   aFile.seek(0);
   alN = 0;
   oFile.write("# Ntotal "+str(Ntotal)+"\n");
   oFile.write("# Nmap "+str(readN)+"\n");

   for line in aFile:
      alignment=line.rstrip().split("\t");
      alN+=1;
      
      readSeq=alignment[4+colS]
      if alignment[1+colS]=="-":
         readSeq = readSeq[::-1];
      
      # write old and init new reads
      if readId!=alignment[0] or readSeq!=seq:
         readId = readId.replace(" ","_");
         oFile.write(readId+" "+str(len(alignments))+" alignments:");
         for align in alignments:
            oFile.write(" " + str(align[0]) + " " + align[1] + " " + str(align[2]));

         oFile.write("\n");

         readId=alignment[0];
         seq = readSeq;
         del alignments[:]
         prob = 1.0;
         for nuc in seq:
            prob *= noiseProb[nuc2i(nuc)];
         alignments.append((0,alignment[1+colS],prob));

      # set transcript id
      if alignment[2+colS][prefixL:] in trMap:
         trans = trMap[ alignment[2+colS][prefixL:] ];
      else:
         trans = 0;
         print "Transcript '"+alignment[2+colS]+"' or '"+alignment[2+colS][prefixL:]+"' was not found in the transcript file.";
         #print alignment;
      # calculate probabilities
      prob=1.0;
      for i in range(len(seq)):
         nuc1 = nuc2i(seq[i]);
         prob *= nucProb[i][nuc1][nuc1];
      
      if len(alignment)==columnN:
            mismatch=alignment[columnN-1]
            mismatchArray = mismatch.split(",");
            for mis in mismatchArray:
               pos = int( mis.split(":")[0] );
               nuc1 = nuc2i( mis.split(":")[1].split(">")[0] );
               nuc2 = nuc2i( mis.split(":")[1].split(">")[1] );
               prob /= nucProb[pos][ nuc2 ][ nuc2 ];
               prob *= nucProb[pos][ nuc1 ][ nuc2 ];
      # add new alignment to list
      alignments.append( (trans, alignment[1+colS], prob) );
   #   if len(alignments)>2 and alignments[len(alignments)-1][2]!=alignments[len(alignments)-2][2]:
   #      print readId;

      
   readId = readId.replace(" ","_");
   oFile.write(readId+" "+str(len(alignments))+" alignments:");
   for align in alignments:
      oFile.write(" " + str(align[0]) + " " + str(align[1]) + " " + str(align[2]));
   oFile.write("\n");
   # }}}
# end if options.inputType=="fasta" }}}
else:
   # {{{ qTOp functions
   if options.inputType=="fastq": Qshift=64;
   if options.inputType=="fastq33": Qshift=33;
   phredWarning = False;
   def qTOp(Q):
      phredS = float(ord(Q)-Qshift);
      if phredS<0:
         if not phredWarning:
            print "WARNING: Phred score too low (",int(phredS),") perhpas use --inputType fastq33.";
            phredWarning=True;
      elif phredS>65:
         if not phredWarning:
            print "NOTE: Phred score unnaturally high (",int(phredS),") check your input type and perhaps set --inputType fastq.";
            phredWarning=True;
      return 1-10**( phredS / -10);
   def qTOpInvert(Q):
      p = 1-10**(float(ord(Q)-Qshift) / -10);
      if p==0: return 1;
      return (1-p)/p;
   #}}}
   # {{{ counting reads
   readN = 0
   rId = "";
   seq = "";
   phread = "";
   aFile.seek(0);
   frags=[]
   while True:
      line = aFile.readline();
      if line == "": break; # empty line means end of file
      if options.paired: 
         line2=aFile.readline();
      
      alignment=line.rstrip().split("\t");
      readId=alignment[0];
      readSeq=alignment[4+colS]
      readPhread=alignment[5+colS];
      
      if readId != rId or readSeq != seq or readPhread != phread:
         readN+=1;
         rId=readId;
         seq=readSeq;
         phread=readPhread;
         if options.paired:
            frags.append( int(line2.rstrip().split("\t")[3+colS]) - int(alignment[3+colS]) );
   if options.paired:
      fragMu = np.mean(frags)
      fragStD = np.std(frags)
   # }}}
   verbose("Writing alignment probabilities");
   aFile.seek(0);
   #{{{ read first read identificators
   alignment=aFile.readline().rstrip().split("\t");
   readId=alignment[0];
   if alignment[1+colS] == "+":
      seq=alignment[4+colS];
      phread=alignment[5+colS]
   else:
      seq=alignment[4+colS][::-1];
      phread=alignment[5+colS][::-1]
   prob=1.0;
   for Q in phread:
      prob *= qTOp(Q);

   if options.paired: #secon pair
      align2 = aFile.readline().rstrip().split("\t");
      fragL = int( align2[3+colS]) - int(alignment[3+colS]);
      prob *= normpdf(fragL,fragMu,fragStD);
      if align2[1+colS] == "+":
         phread2=align2[5+colS]
      else:
         phread2=align2[5+colS][::-1]
      for Q in phread2:
         prob *= qTOp(Q);

   alignments=[]
   aFile.seek(0);
   #}}}
   alN = 0;
   oFile.write("# Ntotal "+str(Ntotal)+"\n");
   oFile.write("# Nmap "+str(readN)+"\n");

   while True:
      line=aFile.readline();
      if line == "": break; # empty line means end of file
      alignment=line.rstrip().split("\t");

      alN+=1;
      
      readSeq=alignment[4+colS]
      readPhread=alignment[5+colS]
      if alignment[1+colS]=="-":
         readPhread = readPhread[::-1]
         readSeq = readSeq[::-1];
      if options.paired:
         align2 = aFile.readline().rstrip().split("\t")
         r2Phread = align2[5+colS]
         if align2[1+colS]=="-":
            r2Phread = r2Phread[::-1];
      else: r2Phread = "";
      
      # write old and init new reads
      if readId!=alignment[0] or readSeq!=seq or readPhread!=phread:
         readId = readId.replace(" ","_");
         oFile.write(readId+" "+str(len(alignments)+1)+" alignments:");
         minProb = 1;
         for align in alignments:
            if minProb > align[2]: minProb=align[2];
            oFile.write(" " + str(align[0]) + " " + align[1] + " " + str(align[2]));
         oFile.write(" 0 + " + str(minProb*qTOpInvert(phread[0])*qTOpInvert(phread[1])*qTOpInvert(phread[2])));
         # add noise alignment with 3 extra mismatches on first bases
         oFile.write("\n");

         readId=alignment[0];
         seq = readSeq;
         phread=readPhread;
         del alignments[:]
         prob=1.0;
         for Q in phread:
            prob *= qTOp(Q);
         if options.paired:
            fragL = int(align2[3+colS])-int(alignment[3+colS]);
            prob *= normpdf(fragL, fragMu, fragStD);
            phread2=r2Phread;
            for Q in phread2:
               prob *= qTOp(Q);
      # set transcript id
      if alignment[2+colS][prefixL:] in trMap:
         trans = trMap[ alignment[2+colS][prefixL:] ];
      else:
         trans = 0;
         print "Transcript '"+alignment[2+colS]+"' or '"+alignment[2+colS][prefixL:]+"' was not found in the transcript file.";
         #print alignment;
      # calculate probabilities
      probLoc = prob;
      if len(alignment)==columnN:
            mismatch=alignment[columnN-1]
            mismatchArray = mismatch.split(",");
            for mis in mismatchArray:
               try:
                  pos = int( mis.split(":")[0] );
               except:
                  pos=0;
                  print 'X',mis,'X',alignment;
               probLoc = probLoc * qTOpInvert(phread[pos]);
      if options.paired and len(align2)==columnN:
            mismatch=align2[columnN-1]
            mismatchArray = mismatch.split(",");
            for mis in mismatchArray:
               try:
                  pos = int( mis.split(":")[0] );
               except:
                  pos=0;
                  print mis
               probLoc = probLoc * qTOpInvert(phread2[pos]);

      # add new alignment to list
      alignments.append( (trans, alignment[1+colS], probLoc) );
   #   if len(alignments)>2 and alignments[len(alignments)-1][2]!=alignments[len(alignments)-2][2]:
   #      print readId;

   readId = readId.replace(" ","_");
   oFile.write(readId+" "+str(len(alignments)+1)+" alignments:");
   minProb = 1;
   for align in alignments:
      if minProb > align[2]: minProb=align[2];
      oFile.write(" " + str(align[0]) + " " + str(align[1]) + " " + str(align[2]));
   oFile.write(" 0 + " + str(minProb*qTOpInvert(phread[0])*qTOpInvert(phread[1])*qTOpInvert(phread[2])));
   # add noise alignment with 1 extra mismatch on first base
   oFile.write("\n");


print "Processed:\n  ",alN,"alignments + (",readN,"noise alignments)\n  ",readN,"reads\n  ",trN,"transcripts\nTotal reads: ",Ntotal,"\n";
aFile.close();
oFile.close();
