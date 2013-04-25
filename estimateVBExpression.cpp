
#include "ArgumentParser.h"
#include "common.h"
#include "FileHeader.h"
#include "misc.h"
#include "MyTimer.h"
#include "SimpleSparse.h"
#include "TagAlignments.h"
#include "transposeFiles.h"
#include "VariationalBayes.h"


SimpleSparse* readData(ArgumentParser &args){//{{{
/*
 As parse(filename,maxreads=None) in python
 Python diferece:
  - missing maxreads check 
    (abort if more than maxreads reads were processed)
*/
   long i,j,num,tid;
   double prb;
   long Ntotal=0,Nmap=0, M=0;
   string readId,strand,blank;
   ifstream inFile;
   MyTimer timer;
   TagAlignments *alignments = new TagAlignments();

   // Read alignment probabilities {{{
   inFile.open(args.args()[0].c_str());
   FileHeader fh(&inFile);
   ns_fileHeader::AlignmentFileType format;
   if((!fh.probHeader(&Nmap,&Ntotal,&format)) || (Nmap ==0)){//{{{
      error("Prob file header read failed.\n");
      return NULL;
   }//}}}
   if(format == ns_fileHeader::OLD_FORMAT){
      error("Please use new/log format of probfile.");
      return NULL;
   }
   message("Mappings: %ld\n",Nmap);
   message("Ntotal: %ld\n",Ntotal);
   alignments->init(Nmap,0,M);
   long mod=10000;
   long bad = 0;
   timer.start();
   for(i = 0; i < Nmap; i++) {
      inFile>>readId>>num;
      if(!inFile.good())break;
     //    message("%s %ld\n",(readId).c_str(),num);
      for(j = 0; j < num; j++) {
         inFile>>tid>>prb;
         if(inFile.fail()){
            inFile.clear();
            // ignore rest of line
            j=num;
            // this read goes to noise
            tid=0;
            // 10 means either 10 or exp(10), but should be still be large enough
            prb=10;
            bad++;
         }
         switch(format){
            case ns_fileHeader::NEW_FORMAT:
               alignments->pushAlignment(tid, prb);
               break;
            case ns_fileHeader::LOG_FORMAT:
               alignments->pushAlignmentL(tid, prb);
               break;
            default:;
         } 
      }
      // ignore rest of line
      inFile.ignore(10000000,'\n');

      alignments->pushRead();
      
      R_INTERUPT;
      if((i % mod == 0)&&(i>0)){
         message("  %ld ",i);
         timer.split();
         mod*=10;
      }
   }
   if(bad>0)warning("Main: %ld reads' alignment information were corrupted.\n",bad);
   inFile.close();
   long Nhits,NreadsReal;
   alignments->finalizeRead(&M, &NreadsReal, &Nhits);
   //}}}
   if(i<Nmap)message("Read only %ld reads.\n",NreadsReal);
   message("Finished Reading!\nTotal hits = %ld\n",Nhits);
   message("Isoforms: %ld\n",M);
   Nmap = NreadsReal;

   SimpleSparse *beta = new SimpleSparse(Nmap, M, Nhits);

   for(i=0;i<=Nmap;i++)beta->rowStart[i]=alignments->getReadsI(i);
   for(i=0;i<Nhits;i++){
      beta->val[i]=alignments->getProb(i);
      beta->col[i]=alignments->getTrId(i);
   }

   delete alignments;
   return beta;
}//}}}

extern "C" int estimateVBExpression(int *argc, char* argv[]) {//{{{
string programDescription =
"Estimates expression given precomputed probabilities of (observed) reads' alignments.\n\
   Uses Variational Bayes algorithm to produce parameters for distribution of relative abundances.\n";
   // Set options {{{
   ArgumentParser args;
   args.init(programDescription,"[prob file]",1);
   args.addOptionS("o","outPrefix","outFilePrefix",1,"Prefix for the output files.");
   //args.addOptionS("O","outType","outputType",0,"Output type (theta, RPKM, counts, tau).","counts");
   //args.addOptionS("p","parFile","parFileName",0,"File containing parameters for the sampler, which can be otherwise specified by --MCMC* options. As the file is checked after every MCMC iteration, the parameters can be adjusted while running.");
   //args.addOptionS("t","trInfoFile","trInfoFileName",0,"File containing transcript information. (Necessary for RPKM)");
   args.addOptionS("m","method","optMethod",0,"Optimalization method (steepest, PR, FR, HS).","FR");
   args.addOptionL("s","seed","seed",0,"Random initialization seed.");
   args.addOptionL("","maxIter","maxIter",0,"Maximum number of iterations.");
   args.addOptionL("P","procN","procN",0,"Limit the maximum number of threads to be used.",1);
   args.addOptionL("","samples","samples",0,"Number of samples to be sampled from the distribution.");
   if(!args.parse(*argc,argv))return 0;
   if(args.verbose)buildTime(argv[0],__DATE__,__TIME__);
   OPT_TYPE optM;
   if(args.isSet("optMethod")){
      if(args.getS("optMethod")=="steepest")optM = OPTT_STEEPEST;
      else if(args.getS("optMethod")=="PR")optM = OPTT_PR;
      else if(args.getS("optMethod")=="FR")optM = OPTT_FR;
      else if(args.getS("optMethod")=="HS")optM = OPTT_HS;
      else optM = OPTT_FR;
   }else  optM = OPTT_FR;
   // }}}
   MyTimer timer;
   timer.start(2);
   long M; 
   SimpleSparse *beta;

   // {{{ Read transcriptInfo and .prob file 
   if(args.verbose)message("Reading data.\n");
/*   if((!args.isSet("trInfoFileName"))||(!trInfo.readInfo(args.getS("trInfoFileName")))){
      if(outTypeI==RPKM){
         error("Main: Missing transcript info file. This will cause problems if producing RPKM.");
         return 1;
      }
      M = 0;
   }else{
      M = trInfo.getM()+1;
   }*/
   beta = readData(args);
   if(! beta){
      error("Main: Reading probabilitites failed.");
      return 1;
   }
   M = beta->M;
   if(M<=0){
      error("Main: Invalid number of transcripts in .prob file.\n");
      return 1;
   }
   // }}}

   if(args.verbose)timer.split();

   if(args.verbose)message("Initializing VB.\n");

   VariationalBayes varB(beta,NULL,ns_misc::getSeed(args),args.getL("procN"));
   
   if(args.verbose)timer.split();
   if(args.verbose)message("Starting VB optimization.\n");
   
#ifdef LOG_CONV
   varB.setLog(args.getS("outFilePrefix")+".convLog",&timer);
#endif

   if(args.isSet("maxIter")) varB.optimize(false,optM,args.getL("maxIter"),1e-7,1e-7);
   else varB.optimize(false,optM);

   if(args.verbose){timer.split(0,'m');}
   double *alpha = varB.getAlphas();
   double alphaSum = 0 ;
   long i;
   for(i=0;i<M;i++)alphaSum+=alpha[i];
   ofstream outF;
   if(! ns_misc::openOutput((args.getS("outFilePrefix")+".m_alphas"), &outF)){
      return 1;
   }
   outF<<"# "<<args.args()[0]<<endl;
   outF<<"# M "<<M<<endl;
   outF<<"# mean alpha beta"<<endl;
   outF<<scientific;
   outF.precision(9);
   for(i=0;i<M;i++){
      outF<<alpha[i]/alphaSum<<" "<<alpha[i]<<" "<<alphaSum-alpha[i]<<endl;
   }
   outF.close();
   // free memory
   delete beta;
   delete[] alpha;
   if(args.isSet("samples") && (args.getL("samples")>0)){
      string samplesFName = args.getS("outFilePrefix")+".VBtheta";
      string samplesTmpName = args.getS("outFilePrefix")+".VBthetaTMP"; 
      timer.start(0);
      if(args.verbose)messageF("Generating samples into temporary file %s. ",samplesTmpName.c_str());
      if(!ns_misc::openOutput(samplesTmpName, &outF)) return 1;
      outF<<"# M "<<M<<" N "<<args.getL("samples")<<endl;
      varB.generateSamples(args.getL("samples"), &outF);
      outF.close();
      if(args.verbose)timer.split(0);
      if(transposeFiles(vector<string>(1, samplesTmpName), samplesFName, args.verbose, "")){
         if(args.verbose)message("Removing temporary file %s.\n", samplesTmpName.c_str());
         remove(samplesTmpName.c_str());
      }else {
         error("Main: Transposing samples failed.");
         return 1;
      }
   }
   if(args.verbose){message("DONE. "); timer.split(2,'m');}
   return 0;
}//}}}

#ifndef BIOC_BUILD
int main(int argc, char* argv[]) {
   return estimateVBExpression(&argc,argv);
}
#endif
