// DECLARATIONS: {{{
#include<cmath>
#include<set>

using namespace std;

#include "ArgumentParser.h"
#include "misc.h"
#include "MyTimer.h"
#include "ReadDistribution.h"
#include "TranscriptExpression.h"
#include "TranscriptInfo.h"
#include "TranscriptSequence.h"

#include "common.h"
//}}}

//#define DEBUG_AT(x) message(x) 
#define DEBUG_AT(x)

namespace ns_parseAlignment {
class TagAlignment{//{{{
   protected:
      int_least32_t trId;
//      bool strand; // true = forward; false = reverse
      double prob,lowProb;
   public:
      TagAlignment(long t=0,double p = 0,double lp = 0){
         trId=(int_least32_t)t;
//         strand=s;
         prob=p;
         lowProb=lp;
      }
      long getTrId()const {return trId;}
      double getProb()const {return prob;}
      double getLowProb()const {return lowProb;}
      void setProb(double p){prob=p;}
}; //}}}

// Check if next fragment is different.
bool nextFragDiffers(const ns_rD::fragmentP curF, const ns_rD::fragmentP nextF, bool mateNamesDiffer);
// String comparison allowing last cmpEPS bases different as long as length
// is the same.
long readNameCmp(const char *str1, const char *str2);
// Read Fragment from SAM file.
// Copies data from 'next' fragment into 'cur' fragment and reads new fragment information into 'next'.
// Fragment is either both paired-ends or just single read.
bool readNextFragment(samfile_t* samData, ns_rD::fragmentP &cur, ns_rD::fragmentP &next);

// Determine input format base either on --format flag or on the file extension.
// Sets format to bam/sam and returns true, or returns false if format is unknown.
bool setInputFormat(const ArgumentParser &args, string *format);

bool openSamFile(const string &name, const string &inFormat, samfile_t **samFile);

bool initializeInfoFile(const ArgumentParser &args, samfile_t *samFile, TranscriptInfo **trInfo, long *M);
} // namespace ns_parseAlignment

extern "C" int parseAlignment(int *argc,char* argv[]){
string programDescription =
"Pre-computes probabilities of (observed) reads' alignments.\n\
   [alignment file] should be in either SAM or BAM format.\n";
   TranscriptInfo *trInfo=NULL;
   TranscriptSequence *trSeq=NULL;
   TranscriptExpression *trExp=NULL;
   MyTimer timer;
   timer.start();
   timer.start(7);
   long Ntotal = 0, Nmap = 0, M=0, i;
   string inFormat;
   samfile_t *samData=NULL;
   ReadDistribution readD;
   ns_rD::fragmentP curF = new ns_rD::fragmentT, nextF = new ns_rD::fragmentT, validAF = new ns_rD::fragmentT;
   // This could be changed to either GNU's hash_set or C++11's unsorted_set, once it's safe.
   set<string> ignoredReads;
   long ignoredMaxAlignments = 0, ignoredSingletons = 0;
   // Intro: {{{
   // Set options {{{
   ArgumentParser args(programDescription,"[alignment file]",1);
   args.addOptionS("o","outFile","outFileName",1,"Name of the output file.");
   args.addOptionS("f","format","format",0,"Input format: either SAM, BAM.");
   args.addOptionS("t","trInfoFile","trInfoFileName",0,"File to save transcript information extracted from [BS]AM file and reference.");
   //args.addOptionS("t","trInfoFile","trInfoFileName",0,"If transcript(reference sequence) information is contained within SAM file, program will write this information into <trInfoFile>, otherwise it will look for this information in <trInfoFile>.");
   args.addOptionS("s","trSeqFile","trSeqFileName",1,"Transcript sequence in FASTA format --- for non-uniform read distribution estimation.");
   args.addOptionS("","trSeqHeader","trSeqHeader",0,"Transcript sequence header format enables gene name extraction (standard/gencode).","standard");
   args.addOptionS("e","expressionFile","expFileName",0,"Transcript relative expression estimates --- for better non-uniform read distribution estimation.");
   args.addOptionL("N","readsN","readsN",0,"Total number of reads. This is not necessary if [SB]AM contains also reads with no valid alignments.");
   args.addOptionS("","failed","failed",0,"File name where to save names of reads that failed to align.");
   args.addOptionB("","uniform","uniform",0,"Use uniform read distribution.");
   args.addOptionD("","lenMu","lenMu",0,"Set mean of log fragment length distribution. (l_frag ~ LogNormal(mu,sigma^2))");
   args.addOptionD("","lenSigma","lenSigma",0,"Set sigma^2 (or variance) of log fragment length distribution. (l_frag ~ LogNormal(mu,sigma^2))");
   args.addOptionS("","distributionFile","distributionFileName",0,"Name of file to which read-distribution should be saved.");
   args.addOptionL("P","procN","procN",0,"Maximum number of threads to be used. This provides speedup mostly when using non-uniform read distribution model (i.e. no --uniform flag).",4);
   args.addOptionB("V","veryVerbose","veryVerbose",0,"Very verbose output.");
   args.addOptionL("","noiseMismatches","numNoiseMismatches",0,"Number of mismatches to be considered as noise.",ns_rD::LOW_PROB_MISSES);
   args.addOptionL("l","limitA","maxAlignments",0,"Limit maximum number of alignments per read. (Reads with more alignments are skipped.)");
   args.addOptionB("","unstranded","unstranded",0,"Paired read are not strand specific.");
   args.addOptionB("","show1warning","show1warning",0,"Show first alignments that are considered wrong (TID unknown, TID mismatch, wrong strand).");
   args.addOptionB("","excludeSingletons","excludeSingletons",0,"Exclude single mate alignments for paired-end reads.");
   args.addOptionB("","mateNamesDiffer","mateNamesDiffer",0,"Mates from paired-end reads have different names.");
   if(!args.parse(*argc,argv))return 0;
   if(args.verbose)buildTime(argv[0],__DATE__,__TIME__);
   readD.setProcN(args.getL("procN"));
   if(args.flag("show1warning"))readD.showFirstWarnings();
   // }}}
   if(!ns_parseAlignment::setInputFormat(args, &inFormat))return 1;
   if(!ns_parseAlignment::openSamFile(args.args()[0], inFormat, &samData))return 1;
   if(!ns_parseAlignment::initializeInfoFile(args, samData, &trInfo, &M))return 1;
   // Read expression and initialize transcript sequence {{{
   if(args.verbose)message("Initializing fasta sequence reader.\n");
   // Initialize fasta sequence reader.
   trSeq = new TranscriptSequence();
   if(args.getLowerS("trSeqHeader") == "gencode"){
      trSeq->readSequence(args.getS("trSeqFileName"), GENCODE); 
   }else{
      trSeq->readSequence(args.getS("trSeqFileName"), STANDARD); 
   }
   // Check numbers for transcripts match.
   if(trSeq->getM() != M){
      error("Main: Number of transcripts in the alignment file and the sequence file are different: %ld vs %ld\n",M,trSeq->getM());
      return 1;
   }
   // Check that length of each transcript matches.
   for(i=0;i<M;i++){
      if(trInfo->L(i) != (long)(trSeq->getTr(i).size())){
         error("Main: Transcript info length and sequence length of transcript %ld DO NOT MATCH! (%ld %d)\n",i,trInfo->L(i),(int)(trSeq->getTr(i).size()));
         return 1;
      }
   }
   // If there were gene names in transcript sequence, assign them to transcript info.
   if(trSeq->hasGeneNames() && (trSeq->getG()>1)){
      if(trInfo->getG() == 1){
         // If just one gene present, then assign gene names.
         if(args.verbose)message("Found gene names in sequence file, updating transcript information.\n");
         trInfo->updateGeneNames(trSeq->getGeneNames());
      }else{
         // If there is more than one gene name already, don't fix.
         if(trInfo->getG() != trSeq->getG()){
            warning("Main: Different number of genes detected in transcript information and sequence file (%ld %ld).\n   You might want to check your data.\n", trInfo->getG(), trSeq->getG());
         }
      }
   }
   // If format is GENCODE and transcript names were extracted, update.
   if((args.getLowerS("trSeqHeader") == "gencode")&&(trSeq->hasTrNames())){
      if(args.flag("veryVerbose"))message("Updating transcript names.\n");
      if(!trInfo->updateTrNames(trSeq->getTrNames())){
         if(args.flag("veryVerbose"))warning("Transcript names update failed.\n");
      }
   }
   if(!args.flag("uniform")){
      // Try loading expression file from previous estimation for non-uniform read model.
      if(args.isSet("expFileName")){
         if(args.verbose)message("Loading transcript initial expression data.\n");
         trExp = new TranscriptExpression(args.getS("expFileName"), GUESS);
         if(trExp->getM() != M){
            error("Main: Number of transcripts in the alignment file and the expression file are different: %ld vs %ld\n",M,trExp->getM());
            return 1;
         }
      }
   }
   //}}}
   timer.split(0,'m');
   //}}}

   // Estimating probabilities {{{
   bool analyzeReads = false;

   if(args.isSet("lenMu") && args.isSet("lenSigma")){
      readD.setLength(args.getD("lenMu"),args.getD("lenSigma"));
   }else{
      analyzeReads = true;
   }
   if(args.flag("uniform")){
      if(args.verbose)message("Using uniform read distribution.\n");
      readD.initUniform(M,trInfo,trSeq,args.flag("veryVerbose"));
   }else{
      if(args.verbose)message("Estimating non-uniform read distribution.\n");
      readD.init(M,trInfo,trSeq,trExp,args.flag("unstranded"),args.flag("veryVerbose"));
      if(args.flag("veryVerbose"))message(" ReadDistribution initialization done.\n");
      analyzeReads = true;
   }
   if(args.isSet("numNoiseMismatches")){
      readD.setLowProbMismatches(args.getL("numNoiseMismatches"));
   }
   // fill in "next" fragment:
   // Counters for all, Good Alignments; and weird alignments
   long observeN, pairedGA, firstGA, secondGA, singleGA, weirdGA, allGA, pairedBad;
   bool storedValidA = false;
   long RE_noEndInfo, RE_weirdPairdInfo, RE_nameMismatch;
   long maxAlignments = 0;
   if(args.isSet("maxAlignments") && (args.getL("maxAlignments")>0))
      maxAlignments = args.getL("maxAlignments");
   // start counting (and possibly estimating):
   observeN = pairedGA = firstGA = secondGA = singleGA = weirdGA = pairedBad = 0;
   RE_noEndInfo = RE_weirdPairdInfo = RE_nameMismatch = 0;
   ns_parseAlignment::readNextFragment(samData, curF, nextF);
   while(ns_parseAlignment::readNextFragment(samData,curF,nextF)){
      R_INTERUPT;
      if( !(curF->first->core.flag & BAM_FUNMAP) ){
         // (at least) The first read was mapped.
         if( curF->paired ) {
            // Fragment's both reads are mapped as a pair.
            // Check mates' names.
            if((ns_parseAlignment::readNameCmp(bam1_qname(curF->first), bam1_qname(curF->second))==0) || 
               (args.flag("mateNamesDiffer"))){
               pairedGA++;
            }else{
               pairedBad++;
               if(RE_nameMismatch == 0){
                  warning("Paired read name mismatch: %s %s\n",bam1_qname(curF->first), bam1_qname(curF->second));
               }
               RE_nameMismatch++;
               if(RE_nameMismatch>10)break;
            }
         }else {
            if (curF->first->core.flag & BAM_FPAIRED) {
               // Read was part of pair (meaning that the other is unmapped).
               if (curF->first->core.flag & BAM_FREAD1) {
                  firstGA++;
               } else if (curF->first->core.flag & BAM_FREAD2) {
                  secondGA++;
               } else weirdGA ++;
            } else {
               // Read is single end, with valid alignment.
               singleGA++;
            }
         }
         // Unless pairedBad>0 the alignment is valid.
         // If excludeSingletons is set, only use paired alignment and alignments of single-end reads.
         if((!storedValidA) && 
            (((!args.flag("excludeSingletons")) && (pairedBad == 0)) ||
             (pairedBad + firstGA + secondGA + weirdGA == 0))){
            validAF->copyFragment(curF);
            storedValidA=true;
         }
      }
      // Next fragment is different.
      if(ns_parseAlignment::nextFragDiffers(curF, nextF, args.flag("mateNamesDiffer"))){
         Ntotal++;
         allGA = singleGA + pairedGA + firstGA +secondGA+ weirdGA;
         if( allGA == 0 ){ // No good alignment.
            // Just in case:
            storedValidA=false;
            pairedBad = 0;
            continue;
         }
         Nmap ++;
         if(weirdGA)RE_noEndInfo++;
         if((singleGA>0) && (pairedGA>0)) RE_weirdPairdInfo++;
         // If it's good uniquely aligned fragment/read, add it to the observation.
         if(( allGA == 1) && analyzeReads && (pairedBad == 0) && storedValidA){
            if(readD.observed(validAF))observeN++;
         }else if(maxAlignments && (allGA>maxAlignments)) {
            // This read will be ignored.
            ignoredReads.insert(bam1_qname(curF->first));
            ignoredMaxAlignments++;
            Nmap --;
         }else if(args.flag("excludeSingletons") && (pairedGA + singleGA == 0)){
            // When excluding singletons only alignments of full pair or single-end read count.
            ignoredReads.insert(bam1_qname(curF->first));
            ignoredSingletons++;
            Nmap --;
         }
         pairedGA = firstGA = secondGA = singleGA = weirdGA = pairedBad = 0;
         storedValidA = false;
      }
   }
   if(RE_nameMismatch>10){
      error("Names of paired mates didn't match at least 10 times.\n"
            "   Something is possibly wrong with your data or the reads have to be renamed.\n");
      return 1;
   }
   message("Reads: all(Ntotal): %ld  mapped(Nmap): %ld\n",Ntotal,Nmap);
   if(args.verbose)message("  %ld reads were used to estimate empirical distributions.\n",observeN);
   if(ignoredMaxAlignments>0)message("  %ld reads are skipped due to having more than %ld alignments.\n",ignoredMaxAlignments, maxAlignments);
   if(ignoredSingletons>0)message("  %ld reads skipped due to having just single mate alignments.\n",ignoredSingletons);
   if(RE_noEndInfo)warning("  %ld reads that were paired, but do not have \"end\" information.\n  (is your alignment file valid?)", RE_noEndInfo);
   if(RE_weirdPairdInfo)warning("  %ld reads that were reported as both paired and single end.\n  (is your alignment file valid?)", RE_weirdPairdInfo);
   readD.writeWarnings();
   if(args.flag("veryVerbose"))timer.split(0,'m');
   // Normalize read distribution:
   if(args.flag("veryVerbose"))message("Normalizing read distribution.\n");
   readD.normalize();
   if(args.isSet("distributionFileName")){
      readD.logProfiles(args.getS("distributionFileName"));
   }
   timer.split(0,'m');
   // }}}

   // Writing probabilities: {{{
   // Re-opening alignment file 
   if(!ns_parseAlignment::openSamFile(args.args()[0], inFormat, &samData))return 1;
   if(args.verbose)message("Writing alignment probabilities.\n");
   double prob,probNoise,minProb;
   prob = probNoise = 0;
   set<string> failedReads;
   vector<ns_parseAlignment::TagAlignment> alignments;
   // Open and initialize output file {{{
   ofstream outF(args.getS("outFileName").c_str());
   if(!outF.is_open()){
      error("Main: Unable to open output file.\n");
      return 1;
   }
   outF<<"# Ntotal "<<Ntotal<<"\n# Nmap "<<Nmap<<"\n# M "<<M<<endl;
   outF<<"# LOGFORMAT (probabilities saved on log scale.)\n# r_name num_alignments (tr_id prob )^*{num_alignments}"<<endl;
   outF.precision(9);
   outF<<scientific;
   // }}}
   
   // start reading:
   timer.start(1);
   bool invalidAlignment = false;
   long readC, pairedN, singleN, firstN, secondN, weirdN, invalidN, noN;
   readC = pairedN = singleN = firstN = secondN = weirdN = invalidN = noN = 0;
   RE_nameMismatch = 0 ;
   // fill in "next" fragment:
   ns_parseAlignment::readNextFragment(samData, curF, nextF);
   while(ns_parseAlignment::readNextFragment(samData,curF,nextF)){
      R_INTERUPT;
      // Skip all alignments of this read.
      if(ignoredReads.count(bam1_qname(curF->first))>0){
         DEBUG_AT(" ignore\n");
         // Read reads while the name is the same.
         while(ns_parseAlignment::readNextFragment(samData,curF,nextF)){
            DEBUG_AT(" ignore\n");
            if(ns_parseAlignment::nextFragDiffers(curF, nextF, args.flag("mateNamesDiffer")))
               break;
         }
         readC++;
         if(args.verbose){ if(progressLog(readC,Ntotal,10,' '))timer.split(1,'m');}
         continue;
      }
      if( !(curF->first->core.flag & BAM_FUNMAP) ){
         DEBUG_AT("M");
         // (at least) The first read was mapped.
         // Check mates' names.
         if(curF->paired && (ns_parseAlignment::readNameCmp(bam1_qname(curF->first), bam1_qname(curF->second))!=0) && (!args.flag("mateNamesDiffer"))){
            if(RE_nameMismatch == 0){
               warning("Paired read name mismatch: %s %s\n",bam1_qname(curF->first), bam1_qname(curF->second));
            }
            RE_nameMismatch++;
            if(RE_nameMismatch>10)break;
            invalidAlignment = true;
         }else if((!args.flag("excludeSingletons")) || curF->paired || (! (curF->first->core.flag & BAM_FPAIRED))){
            // We only calculate probabilties and add alignments if: 
            // (singletons are not exlucded) OR  (it is a proper paired alignments) OR (it is single-end read)
            if(readD.getP(curF, prob, probNoise)){
               // We calculated valid probabilities for this alignment.   
               // Add alignment:
               alignments.push_back(ns_parseAlignment::TagAlignment(curF->first->core.tid+1, prob, probNoise));
               // Update counters:
               if( curF->paired ) {
                  // Fragment's both reads are mapped as a pair.
                  pairedN++;
                  DEBUG_AT(" P\n");
               }else {
                  if (curF->first->core.flag & BAM_FPAIRED) {
                     // Read was part of pair (meaning that the other is unmapped).
                     if (curF->first->core.flag & BAM_FREAD1) {
                        firstN++;
                        DEBUG_AT(" 1\n");
                     } else if (curF->first->core.flag & BAM_FREAD2) {
                        secondN++;
                        DEBUG_AT(" 2\n");
                     } else {
                        weirdN ++;
                        DEBUG_AT(" W\n");
                     }
                  } else {
                     // Read is single end, with valid alignment.
                     singleN++;
                     DEBUG_AT(" S\n");
                  }
               }
            } else {
               // Calculation of alignment probabilities failed.
               invalidAlignment = true;
            }
         }
      }else DEBUG_AT("UNMAP\n");
      // next fragment has different name
      if(ns_parseAlignment::nextFragDiffers(curF, nextF, args.flag("mateNamesDiffer"))){
         DEBUG_AT("  last\n");
         readC++;
         if(args.verbose){ if(progressLog(readC,Ntotal,10,' '))timer.split(1,'m');}
         if(!alignments.empty()){
            outF<<bam1_qname(curF->first)<<" "<<alignments.size()+1;
            minProb = 1;
            for(i=0;i<(long)alignments.size();i++){
               if(minProb>alignments[i].getLowProb())minProb = alignments[i].getLowProb();
               outF<<" "<<alignments[i].getTrId()
//                   <<" "<<getStrandC(alignments[i].getStrand())
                   <<" "<<alignments[i].getProb();
            }
            outF<<" 0 "<<minProb<<endl;
            alignments.clear();
         }else{
            // read has no valid alignments:
            if(invalidAlignment){
               // If there were invalid alignments, write a mock record in order to keep Nmap consistent.
               invalidN++;
               outF<<bam1_qname(curF->first)<<" 1 0 0"<<endl;
            }else {
               noN++;
            }
            if(args.isSet("failed")){
               // Save failed reads.
               failedReads.insert(bam1_qname(curF->first));
               if(curF->paired)failedReads.insert(bam1_qname(curF->second));
            }
         }
         invalidAlignment = false;
      }
   }
   if(RE_nameMismatch>10){
      error("Names of paired mates didn't match at least 10 times.\n"
            "   Something is possibly wrong with your data or the reads have to be renamed.\n");
      return 1;
   }
   outF.close();
   timer.split(0,'m');
   if(args.verbose){
      message("Analyzed %ld reads:\n",readC);
      if(ignoredMaxAlignments>0)message(" %ld ignored due to --limitA flag\n",ignoredMaxAlignments);
      if(invalidN>0)message(" %ld had only invalid alignments (see warnings)\n",invalidN);
      if(noN>0)message(" %ld had no alignments\n",noN);
      message("The rest had %ld alignments:\n",pairedN+singleN+firstN+secondN+weirdN);
      if(pairedN>0)message(" %ld paired alignments\n",pairedN);
      if(firstN+secondN+weirdN>0)
         message(" %ld half alignments (paired-end mates aligned independently)\n",firstN+secondN+weirdN);
      if(singleN>0)message(" %ld single-read alignments\n",singleN);
      //flushStdout();
      messageFlush();
   }else {
      messageF("Alignments: %ld.\n",pairedN+singleN+firstN+secondN+weirdN);
   }
   readD.writeWarnings();
   if(args.flag("veryVerbose")){
      message("Number of weights cached: %ld\n",readD.getWeightNormCount());
   }
   // Deal with reads that failed to align {{{
   if(args.isSet("failed")){
      outF.open(args.getS("failed").c_str());
      if(outF.is_open()){
         for(set<string>::iterator setIt=failedReads.begin(); setIt!=failedReads.end();setIt++)
            outF<<*setIt<<endl;
         outF.close();
      }
   } //}}}
   // Compute effective length and save transcript info {{{
   if(args.isSet("trInfoFileName")){
      if(args.verbose)messageF("Computing effective lengths.\n");
      trInfo->setEffectiveLength(readD.getEffectiveLengths());
      if(! trInfo->writeInfo(args.getS("trInfoFileName"))){
         warning("Main: File %s probably already exists.\n"
                 "   Will save new transcript info into %s-NEW.\n",(args.getS("trInfoFileName")).c_str(),(args.getS("trInfoFileName")).c_str());
         if(! trInfo->writeInfo(args.getS("trInfoFileName")+"-NEW", true)){ // DO OVERWRITE
            warning("Main: Writing into %s failed!.",(args.getS("trInfoFileName")+"-NEW").c_str());
         }
      }else {
         if(args.verbose)message("Transcript information saved into %s.\n",(args.getS("trInfoFileName")).c_str());
      }
      if(args.verbose)timer.split(0,'m');
   } //}}}
   // Close, free and write failed reads if filename provided {{{
   delete curF;
   delete nextF;
   delete validAF;
   delete trInfo;
   delete trSeq;
   delete trExp;
   samclose(samData);
   // }}}
   // }}}
   if(args.verbose)message("DONE. ");
   timer.split(7,'m');
   return 0;
}

#ifndef BIOC_BUILD
int main(int argc,char* argv[]){
   return parseAlignment(&argc,argv);
}
#endif

namespace ns_parseAlignment {

bool nextFragDiffers(const ns_rD::fragmentP curF, const ns_rD::fragmentP nextF, bool mateNamesDiffer){//{{{
   if(readNameCmp(bam1_qname(curF->first), bam1_qname(nextF->first))==0) return false;
   if(mateNamesDiffer && (readNameCmp(bam1_qname(curF->first), bam1_qname(nextF->second))==0)) return false;
   return true;
}//}}}

long readNameCmp(const char *str1, const char *str2){//{{{
   // Check first character(so that we can look back later).
   if(*str1 != *str2)return *str1 - *str2;
   while(*str1 || *str2){
      if(*str1 != *str2){
         // They can differ in last character if its preceeeded by /:_.
         if(*str1 && *str2 && (!*(str1+1)) && (!*(str2+1)) && 
            ((*(str1-1) == '/') || (*(str1-1) == ':') || (*(str1-1) == '_'))){
            return 0;
         }
         return *str1 - *str2;
      }
      str1++;
      str2++;
   }
   return 0;
}//}}}

bool readNextFragment(samfile_t* samData, ns_rD::fragmentP &cur, ns_rD::fragmentP &next){//{{{
   static ns_rD::fragmentP tmpF = NULL;
   bool currentOK = true;
   // switch current to next:
   tmpF = cur;
   cur = next;
   next = tmpF;
   // check if current fragment is valid
   if( !cur->first->data || ( *(cur->first->data) == '\0')){
      // current fragment is invalid
      currentOK = false;
   }
   // try reading next fragment:
   if(samread(samData,next->first)<0){
      // read failed: set next reads name to empty string
      *(next->first->data) = '\0';
      return currentOK;
   }
   // if paired then try reading second pair
   if( !(next->first->core.flag & BAM_FPROPER_PAIR) || 
       (samread(samData,next->second)<0)){
      next->paired = false;
   }else{
      /* look for last read 
      ignored -- expecting always only singles and pairs
      THIS WOULD NOT WORK 
        - because SAM FLAG does not indicate whether it's first or last read of a pair, but which "file" is the read from 
        - this was observed from bowtie alignment data
      while( !(next->second->core.flag & BAM_FREAD2) &&
             (samread(samData,next->second)>=0)) ;
      */
      next->paired = true;
   }
   return currentOK;
}//}}}

bool setInputFormat(const ArgumentParser &args, string *format){//{{{
   if(args.isSet("format")){
      *format = args.getLowerS("format");
      if((*format =="sam")||(*format == "bam")){
         return true;
      }
      warning("Unknown format '%s'.\n",format->c_str());
   }
   string fileName = args.args()[0];
   string extension = fileName.substr(fileName.rfind(".")+1);
   *format = ns_misc::toLower(extension);
   if((*format =="sam")||(*format == "bam")){
      if(args.verb())message("Assuming alignment file in '%s' format.\n",format->c_str());
      return true;
   }
   message("Unknown extension '%s'.\n",extension.c_str());
   error("Couldn't determine the type of input file, please use --format and check your input.\n");
   return false;
}//}}}

bool openSamFile(const string &name, const string &inFormat, samfile_t **samFile){//{{{
   if(*samFile != NULL)samclose(*samFile);
   if(inFormat=="bam") *samFile = samopen(name.c_str(), "rb" , NULL);
   else *samFile = samopen(name.c_str(), "r" , NULL);
   if(*samFile == NULL){
      error("Failed re-reading alignments.\n");
      return false;
   }
   return true;
}//}}}

bool initializeInfoFile(const ArgumentParser &args, samfile_t *samFile, TranscriptInfo **trInfo, long *M){//{{{
   if((samFile->header == NULL)||(samFile->header->n_targets == 0)){
      if(! args.isSet("trInfoFileName")){
         error("Main: alignment file does not contain header, or the header is empty.\n"
               "  Please either include header in alignment file or provide transcript information file.\n"
               "  (option --trInfoFile, file should contain lines with <gene name> <transcript name> <transcript length>.\n");
         return false;
      }else{
         if(args.verb())message("Using %s for transcript information.\n",(args.getS("trInfoFileName")).c_str());
         if((*trInfo = new TranscriptInfo(args.getS("trInfoFileName"))) && (*trInfo)->isOK()){
            *M=(*trInfo)->getM();
         }else {
            error("Main: Can't get transcript information.\n");
            return false;
         }
      }
   }else{
      if(args.verbose)message("Using alignments' header for transcript information.\n");
      *M = samFile->header->n_targets;
      vector<string> trNames(*M);
      vector<long> trLengths(*M);
      for(long i=0;i<*M;i++){
         trNames[i] = samFile->header->target_name[i];
         trLengths[i] = samFile->header->target_len[i]; 
      }
      *trInfo = new TranscriptInfo();
      if(! (*trInfo)->setInfo(vector<string>(*M,"none"), trNames, trLengths)){
         error("TranscriptInfo not initialized.\n");
         return false;
      }
   }
   return true;
}//}}}

} // namespace ns_parseAlignment
