// DECLARATIONS: {{{
#include<cmath>
#include<omp.h>
#include<set>

using namespace std;

#include "ArgumentParser.h"
#include "common.h"
#include "MyTimer.h"
#include "ReadDistribution.h"
#include "TranscriptExpression.h"
#include "TranscriptInfo.h"
#include "TranscriptSequence.h"

//}}}
#define Sof(x) (long)x.size()

// Read is not marked as not mapped, and is either not marked as paired, or it is marked as proper pair. 
#define FRAG_IS_ALIGNED(x) \
   ( !(x->first->core.flag & BAM_FUNMAP) && \
     ( !(x->first->core.flag & BAM_FPAIRED) || \
       (x->first->core.flag & BAM_FPROPER_PAIR) ) )

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

string toLower(string str);

bool readNextFragment(samfile_t* samData, fragmentP &cur, fragmentP &next);

string getInputFormat(const ArgumentParser &args);
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
   long M=0,i;
   samfile_t *samData=NULL;
   // Intro: {{{
   // Set options {{{
   ArgumentParser args(programDescription,"[alignment file]",1);
   args.addOptionS("o","outFile","outFileName",1,"Name of the output file.");
   args.addOptionS("f","format","format",0,"Input format: either SAM, BAM.");
   args.addOptionS("t","trInfoFile","trInfoFileName",0,"If transcript(reference sequence) information is contained within SAM file, program will write this information into <trInfoFile>, otherwise it will look for this information in <trInfoFile>.");
   args.addOptionS("s","trSeqFile","trSeqFileName",1,"Transcript sequence in FASTA format --- for non-uniform read distribution estimation.");
   args.addOptionS("e","expressionFile","expFileName",0,"Transcript relative expression estimates --- for better non-uniform read distribution estimation.");
   args.addOptionL("N","readsN","readsN",0,"Total number of reads. This is not necessary if [SB]AM contains also reads with no valid alignments.");
   args.addOptionS("","failed","failed",0,"File name where to save names of reads that failed to align as pair.");
   args.addOptionB("","uniform","uniform",0,"Use uniform read distribution.");
   args.addOptionD("","lenMu","lenMu",0,"Set mean of log fragment length distribution. (l_frag ~ LogNormal(mu,sigma^2))");
   args.addOptionD("","lenSigma","lenSigma",0,"Set sigma^2 (or variance) of log fragment length distribution. (l_frag ~ LogNormal(mu,sigma^2))");
   args.addOptionS("","distributionFile","distributionFileName",0,"Name of file to which read-distribution should be saved.");
   args.addOptionL("P","procN","procN",0,"Maximum number of threads to be used. This provides parallelization only when computing non-uniform read distribution (i.e. runs without --uniform flag).",3);
   args.addOptionB("V","veryVerbose","veryVerbose",0,"Very verbose output.");
   args.addOptionL("","noiseMismatches","numNoiseMismatches",0,"Number of mismatches to be considered as noise.",LOW_PROB_MISSES);
   if(!args.parse(*argc,argv))return 0;
   if(args.verbose)buildTime(argv[0],__DATE__,__TIME__);
#ifdef SUPPORT_OPENMP
   omp_set_num_threads(args.getL("procN"));
#endif
   /*if((! args.flag("uniform"))&&(! args.isSet("trSeqFileName"))){
      error("Please provide transcript sequence file in fasta format (option --trSeqFile) for non-uniform read distribution estimation.\n");
      return 1;
   }*/
   // }}}
   // Read transcriptInfo and initialize alignment file {{{
   string inFormat=ns_parseAlignment::getInputFormat(args);
   if(inFormat=="bam")
      samData = samopen(args.args()[0].c_str(), "rb" , NULL);
   else 
      samData = samopen(args.args()[0].c_str(), "r" , NULL);
   if(samData == NULL){
      error("Failed reading alignments from %s.\n",args.args()[0].c_str());
      return 1;
   }
   if(samData->header == NULL){
      if(! args.isSet("trInfoFileName")){
         //error("Main: %s file does not contain header.\n   Need transcript information file containing lines with <gene name> <transcript name> <transcript length>.\n   Use option --trInfoFile\n",(args.getS("format")).c_str());
         //return 1;
      }else{
         if(args.verbose)message("Using %s for transcript information.\n",(args.getS("trInfoFileName")).c_str());
         if(! ( trInfo = new TranscriptInfo(args.getS("trInfoFileName")))){
            error("Main: Can't get transcript information\n");
            return 1;
         }
         M=trInfo->getM();
      }
   }else{
      if(args.verbose)message("Using %s header for transcript information.\n",(args.getS("format")).c_str());
      M = samData->header->n_targets;
      vector<string> trNames(M);
      vector<long> trLengths(M);
      for(i=0;i<M;i++){
         trNames[i] = samData->header->target_name[i];
         trLengths[i] = samData->header->target_len[i]; 
      }
      trInfo = new TranscriptInfo();
      if(! trInfo->setInfo(vector<string>(M,"none"), trNames, trLengths)){
         error("TranscriptInfo not initialized.\n");
         return 1;
      }
   }//}}}
   // Read expression and initialize transcript sequence {{{
   if(args.verbose)message("Initializing fasta sequence reader.\n");
   // Initialize fasta sequence reader.
   trSeq = new TranscriptSequence();
   trSeq->readSequence(args.getS("trSeqFileName")); 
   // Check numbers for transcripts match.
   if(trSeq->getM() != M){
      error("Main: Number of transcripts in the alignment(%s) file and the sequence file are different: %ld vs %ld\n",args.getS("format").c_str(),M,trSeq->getM());
      return 1;
   }
   // Check that length of each transcript matches.
   for(i=0;i<M;i++){
      if(trInfo->L(i) != (long)(trSeq->getTr(i))->size()){
         error("Main: Transcript info length and sequence length of transcript %ld DO NOT MATCH! (%ld %d)\n",i,trInfo->L(i),(int)((trSeq->getTr(i))->size()));
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
   if(!args.flag("uniform")){
      // Try loading expression file from previous estimation for non-uniform read model.
      if(args.isSet("expFileName")){
         if(args.verbose)message("Loading transcript initial expression data.\n");
         trExp = new TranscriptExpression(args.getS("expFileName"));
         if(trExp->getM() != M){
            error("Main: Number of transcripts in the alignment(%s) file and the expression are different: %ld vs %ld\n",args.getS("format").c_str(),M,trExp->getM());
            return 1;
         }
      }
   }
   //}}}
   //}}}
   timer.split(0,'m');

   fragmentP curF = new fragmentT, nextF = new fragmentT;
   long Ntotal=0, Nmap=0;
   ReadDistribution readD(M);
   // Estimating probabilities {{{
   bool analyzeReads = false;

   if(args.isSet("lenMu") && args.isSet("lenSigma")){
      readD.setLength(args.getD("lenMu"),args.getD("lenSigma"));
   }else{
      analyzeReads = true;
   }
   if(args.flag("uniform")){
      if(args.verbose)message("Using uniform read distribution.\n");
      readD.initUniform(trInfo,trSeq,args.flag("veryVerbose"));
   }else{
      if(args.verbose)message("Estimating non-uniform read distribution.\n");
      readD.init(trInfo,trSeq,trExp,args.flag("veryVerbose"));
      message("Init done\n");
      analyzeReads = true;
   }
   if(args.isSet("numNoiseMismatches")){
      readD.setLowProbMismatches(args.getL("numNoiseMismatches"));
   }
   // fill in "next" fragment:
   ns_parseAlignment::readNextFragment(samData, curF, nextF);
   long alN = 0, alGoodN = 0;
   // start counting (and possibly estimating):
   while(ns_parseAlignment::readNextFragment(samData,curF,nextF)){
      alN ++;
      if( FRAG_IS_ALIGNED(curF) ) alGoodN++;
      // Next read is different.
      if(strcmp(bam1_qname(curF->first), bam1_qname(nextF->first))!=0){
         Ntotal++;
         if( alGoodN > 0 ) Nmap ++;
         // If it's good uniquely aligned read, add it to the observation.
         if(( alGoodN == 1) && ( alN == 1 )  && analyzeReads)
               readD.observed(curF);
         alN = 0;
         alGoodN = 0;
      }
   }
   if(args.verbose)message("Ntotal: %ld  Nmap: %ld\n",Ntotal,Nmap);
   // Normalize read distribution:
   timer.split(0,'m');
   if(args.verbose)message("Normalizing read distribution.\n");
   readD.normalize();
   timer.split(0,'m');
   if(args.isSet("distributionFileName")){
      readD.logProfiles(args.getS("distributionFileName"));
   }

   // Re-opening alignment file {{{
   samclose(samData);
   if(inFormat=="bam")
      samData = samopen(args.args()[0].c_str(), "rb" , NULL);
   else 
      samData = samopen(args.args()[0].c_str(), "r" , NULL);
   if(samData == NULL){
      error("Failed re-reading alignments.\n");
      return 1;
   }//}}}
   // }}}

   // Writing probabilities: {{{
   if(args.verbose)message("Writing alignment probabilities\n");
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
   outF<<"# Ntotal "<<Ntotal<<"\n# Nmap "<<Nmap<<endl;
   outF<<"# NEWFORMAT \n# r_name num_alignments (tr_id prob )^*{num_alignments}"<<endl;
   outF.precision(9);
   outF<<scientific;
   // }}}
   
   // fill in "next" fragment:
   ns_parseAlignment::readNextFragment(samData, curF, nextF);
   // start reading:
   long readC = 0;
   timer.start(1);
   long pairedN = 0;
   long singleN = 0;
   while(ns_parseAlignment::readNextFragment(samData,curF,nextF)){
//      if(curF->paired)message("P %s %s\n",bam1_qname(curF->first),bam1_qname(curF->second));
      if( FRAG_IS_ALIGNED(curF) ){
         if(curF->paired)pairedN++;
         else singleN++;
         readD.getP(curF, prob, probNoise);
         alignments.push_back(ns_parseAlignment::TagAlignment(curF->first->core.tid+1, prob, probNoise));
      }
      // next fragment has different name
      if(strcmp(bam1_qname(curF->first), bam1_qname(nextF->first))!=0){
         readC++;
         if(args.verbose){ if(progressLog(readC,Ntotal,10))timer.split(1,'m');}
         if(Sof(alignments)>0){
            outF<<bam1_qname(curF->first)<<" "<<Sof(alignments)+1;
            minProb = 1;
            for(i=0;i<Sof(alignments);i++){
               if(minProb>alignments[i].getLowProb())minProb = alignments[i].getLowProb();
               outF<<" "<<alignments[i].getTrId()
//                   <<" "<<getStrandC(alignments[i].getStrand())
                   <<" "<<alignments[i].getProb();
            }
            outF<<" 0 "<<minProb<<endl;
            alignments.clear();
         }else{
            // read has no valid alignments:
            failedReads.insert(bam1_qname(curF->first));
            if(curF->paired)failedReads.insert(bam1_qname(curF->second));
         }
      }
      R_INTERUPT;
   }
   if(args.verbose)message("Analyzed %ld single reads and %ld paired-end reads.\n",singleN,pairedN);
   outF.close();
   if(args.verbose)message("Found %ld single reads without proper alignment.\n",Sof(failedReads));
   // Deal with reads that failed to align
   if(args.isSet("failed")){
      outF.open(args.getS("failed").c_str());
      if(outF.is_open()){
         for(set<string>::iterator setIt=failedReads.begin(); setIt!=failedReads.end();setIt++)
            outF<<*setIt<<endl;
         outF.close();
      }
   }
   timer.split(0,'m');
   // Compute effective length and save transcript info
   if(args.isSet("trInfoFileName")){
      if(args.verbose)message("Computing effective lengths.\n");
      trInfo->setEffectiveLength(readD.getEffectiveLengths());
      if(args.verbose)message("Writing transcript information into %s.\n",(args.getS("trInfoFileName")).c_str());
      if(! trInfo->writeInfo(args.getS("trInfoFileName"))){
         warning("Main: writing to %s failed.\nWill try %s-NEW but you should rename it afterwards if you're planning to use it.\n",(args.getS("trInfoFileName")).c_str(),(args.getS("trInfoFileName")).c_str());
         trInfo->writeInfo(args.getS("trInfoFileName")+"-NEW", true); // DO OVERWRITE
      }
      timer.split(0,'m');
   }
   // Close, free and write failed reads if filename provided {{{
   delete curF;
   delete nextF;
   delete trInfo;
   delete trSeq;
   delete trExp;
   samclose(samData);
   // }}}
   // }}}
   if(args.verbose)message("DONE.\n");
   timer.split(0,'m');
   return 0;
}

#ifndef BIOC_BUILD
int main(int argc,char* argv[]){
   return parseAlignment(&argc,argv);
}
#endif

namespace ns_parseAlignment {

string toLower(string str){//{{{
   for(long i=0;i<Sof(str);i++)
      if((str[i]>='A')&&(str[i]<='Z'))str[i]=str[i]-'A'+'a';
   return str;
}//}}}

bool readNextFragment(samfile_t* samData, fragmentP &cur, fragmentP &next){//{{{
   static fragmentP tmpF = NULL;
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
   if( !( next->first->core.flag & BAM_FPAIRED) || 
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

string getInputFormat(const ArgumentParser &args){
   return "sam";
}

} // namespace ns_parseAlignment
