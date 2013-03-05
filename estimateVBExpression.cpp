
#include "FileHeader.h"
#include "MyTimer.h"
#include "ArgumentParser.h"
#include "SimpleSparse.h"
#include "VariationalBayes.h"
#include "TagAlignments.h"
#include "common.h"


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
   bool newformat=true;
   if((!fh.probHeader(Nmap,Ntotal,newformat)) || (Nmap ==0)){//{{{
      error("Prob file header read failed.\n");
      return NULL;
   }//}}}
   if(! newformat){
      error("Please use new format of probfile.");
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
            //inFile.ignore(10000000,'\n');
            // ignore other read's alignments
            j=num;
            // this read goes to noise
            tid=0;
            prb=100;
            bad++;
         }
         alignments->pushAlignment(tid, prb);         
      }
      // ignore rest of line
      inFile.ignore(10000000,'\n');

      alignments->pushRead();

#ifdef BIOC_BUILD
      R_CheckUserInterrupt();
#endif
      if((i % mod == 0)&&(i>0)){
         message("  %ld ",i);
         timer.split();
         mod*=10;
      }
   }
   //message("Bad: %ld\n",bad);
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
   Uses Variational Bayes algorithm to produce parameters for distribution  of relative abundances.\n";
   buildTime(argv[0],__DATE__,__TIME__);
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
   if(!args.parse(*argc,argv))return 0;
   OPT_TYPE optM;
   if(args.isSet("optMethod")){
      if(args.getS("optMethod")=="steepest")optM = OPTT_STEEPEST;
      else if(args.getS("optMethod")=="PR")optM = OPTT_PR;
      else if(args.getS("optMethod")=="FR")optM = OPTT_FR;
      else if(args.getS("optMethod")=="HS")optM = OPTT_HS;
      else optM = OPTT_FR;
   }else  optM = OPTT_FR;
   long seed=0;
   if(args.isSet("seed"))seed=args.getL("seed");
   // }}}
   MyTimer timer;
   timer.start(17);
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

   VariationalBayes varB(beta,NULL,seed,args.getL("procN"));
   
   if(args.verbose)timer.split();
   if(args.verbose)message("Starting VB optimization.\n");
   
#ifdef LOG_CONV
   varB.setLog(args.getS("outFilePrefix")+".convLog",&timer);
#endif

   if(args.isSet("maxIter")) varB.optimize(false,optM,args.getL("maxIter"),1e-7,1e-7);
   else varB.optimize(false,optM);

   if(args.verbose){message("DONE. "); timer.split(0,'m');}
   timer.split(17,'m');
   double *alpha = varB.getAlphas();
   double alphaSum = 0 ;
   long i;
   for(i=0;i<M;i++)alphaSum+=alpha[i];
   ofstream outF((args.getS("outFilePrefix")+".m_alphas").c_str());
   if(!outF.is_open()){
      error("Main: Unable to open output file.");
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
   return 0;
}//}}}

#ifndef BIOC_BUILD
int main(int argc, char* argv[]) {
   return estimateVBExpression(&argc,argv);
}
#endif
