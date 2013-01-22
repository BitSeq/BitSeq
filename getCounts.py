#!/usr/bin/python
# Initialization {{{
import sys
import numpy as np
#import os, time # needed for this:
#time_str = time.strftime("%b %e %Y %H:%M:%S", time.gmtime(os.lstat(sys.argv[0]).st_mtime));
#print "###",os.path.basename(sys.argv[0]),"build:",time_str;

from optparse import OptionParser
parser = OptionParser(usage="%prog [options] [<inputFile.thetaMeans>]+\n\n\
      This program reads supplied .thetaMeans files and using either information from .prob files or Nmap option generates read counts for each input file provided.")
parser.add_option("-o", "--outFile", dest="out", help="Output file", type="string")
parser.add_option("-v", "--verbose", default=False, dest="verbose",  action="store_true", help="Verbose output")
parser.add_option("-p", "--probDir", dest="probDir", help="Directory with .prob files. The program will look in here for files with same name except fot extension .prob in order to find out total-aligned-read counts for each experiment.", type="string")
parser.add_option("-n", "--Nmap", dest="Nmap", help = "Comma separated list of total aligned-read-counts for each experiment.",type="string");
def verbose(str):
   if options.verbose:
      print str;
(options, args) = parser.parse_args()

if len(args)==0:
   sys.exit("Please supply .thetaMeans filenames as arguments.");
if not options.out:
   sys.exit("Please supply output file");
if (not options.probDir) and (not options.Nmap):
   sys.exit("Please use either --Nmap or --probDir.");
#}}}

if options.Nmap:
   try:
      N = [ float(it) for it in options.Nmap.split(",")]
      if len(N) != len(args):
         raise;
   except:
      sys.exit("Unable to turn '"+options.Nmap+"' into "+str(len(args))+" numbers.");
else:
   N = []
   for arg in args:
      fn = arg.split("/")[-1];
      if fn[-11:] == '.thetaMeans':
         fn = options.probDir +"/"+fn[:-11]+".prob";
      else:
         fn = options.probDir +"/"+fn+".prob";
      try:
         inF = open(fn);
      except:
         sys.exit("Unable to open file: "+fn);
      print "Reading file: ",fn;
      Nmap = 0;
      for line in inF:
         if line[0]!="#": break;
         ls=line.split();
         for i in xrange(len(ls)-1): 
            if ls[i] == "Nmap": Nmap = int(ls[i+1]);
      inF.close();
      if Nmap <= 0:
         sys.exit("Unable to find valid Nmap in: "+fn);
      N.append(Nmap);


means = [np.transpose(np.loadtxt(arg))[1] for arg in args];
print "Files:";
for j in xrange(len(args)):
   print "  ",args[j],N[j];

try:
   outF = open(options.out,"w");
except:
   sys.exit("Unable to open output file: ",options.out);

for i in xrange(len(means[0])):
   for j in xrange(len(means)):
      outF.write(str(long(round(means[j][i]*N[j])))+" ");
   outF.write("\n");

outF.close();


