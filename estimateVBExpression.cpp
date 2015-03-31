#include "ArgumentParser.h"
#include "FileHeader.h"
#include "misc.h"
#include "MyTimer.h"
#include "SimpleSparse.h"
#include "TagAlignments.h"
#include "transposeFiles.h"
#include "VariationalBayes.h"

#include "common.h"

SimpleSparse* readData(const ArgumentParser &args, long trM){//{{{
/*
 As parse(filename,maxreads=None) in python
 Python difference:
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
   if((!fh.probHeader(&Nmap,&Ntotal,&M,&format)) || (Nmap ==0)){//{{{
      error("Prob file header read failed.\n");
      return NULL;
   }//}}}
   if(format == ns_fileHeader::OLD_FORMAT){
      error("Please use new/log format of Prob file.");
      return NULL;
   }
   message("N mapped: %ld\n",Nmap);
   messageF("N total:  %ld\n",Ntotal);
   if(args.verb())message("Reading alignments.\n");
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
      if(args.verb() && (i % mod == 0) && (i>0)){
         message("  %ld ",i);
         timer.split();
         mod*=10;
      }
   }
   if(bad>0)warning("Main: %ld reads' alignment information were corrupted.\n",bad);
   inFile.close();
   long Nhits,NreadsReal;
   alignments->finalizeRead(&M, &NreadsReal, &Nhits);
   // Increase M based on number of transcripts in trInfo file.
   if(M<trM)M = trM;
   //}}}
   if(i<Nmap)message("Read only %ld reads.\n",NreadsReal);
   message("All alignments: %ld\n",Nhits);
   messageF("Isoforms: %ld\n",M);
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
   args.addOptionS("O","outType","outputType",0,"Output type (theta, RPKM, TPM, counts) of the samples sampled from the distribution.","theta");
   args.addOptionS("t","trInfoFile","trInfoFileName",0,"File containing transcript information. (Necessary for RPKM and TPM samples)");
   args.addOptionL("P","procN","procN",0,"Limit the maximum number of threads to be used.",4);
   args.addOptionS("m","method","optMethod",0,"Optimization method (steepest, PR, FR, HS).","FR");
   args.addOptionL("s","seed","seed",0,"Random initialization seed.");
   args.addOptionL("","maxIter","maxIter",0,"Maximum number of iterations.",(long)1e4);
   args.addOptionD("","optLimit","limit",0,"Optimisation limit in terms of minimal gradient or change of bound.",1e-5); 
   args.addOptionL("","samples","samples",0,"Number of samples to be sampled from the distribution.");
   args.addOptionB("V","veryVerbose","veryVerbose",0,"More verbose output, better if output forwarded into file.");
   if(!args.parse(*argc,argv))return 0;
   if(args.verbose)buildTime(argv[0],__DATE__,__TIME__);
   OPT_TYPE optM;
   if(args.isSet("optMethod")){
      if((args.getLowerS("optMethod")=="steepest")||
         (args.getLowerS("optMethod")=="vbem"))optM = OPTT_STEEPEST;
      else if(args.getLowerS("optMethod")=="pr")optM = OPTT_PR;
      else if(args.getLowerS("optMethod")=="fr")optM = OPTT_FR;
      else if(args.getLowerS("optMethod")=="hs")optM = OPTT_HS;
      else optM = OPTT_FR;
   }else  optM = OPTT_FR;
   args.updateS("outputType", ns_expression::getOutputType(args, "theta"));
   if(args.getS("outputType") == "tau"){
      error("Main: 'tau' is not valid output type.\n");
      return 1;
   }
   // }}}
   MyTimer timer;
   timer.start(2);
   long M = 0; 
   SimpleSparse *beta;
   TranscriptInfo trInfo;

   // {{{ Read transcriptInfo and .prob file 
   if((!args.isSet("trInfoFileName"))||(!trInfo.readInfo(args.getS("trInfoFileName")))){
     if(args.isSet("samples") && (args.getL("samples")>0) && ((args.getS("outputType") == "rpkm") || (args.getS("outputType") == "tpm"))){
         error("Main: Missing transcript info file. The file is necessary for producing RPKM or TPM samples.\n");
         return 1;
      }
   }else{
      M = trInfo.getM()+1;
   }
   beta = readData(args,M);
   if(! beta){
      error("Main: Reading probabilities failed.\n");
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

   // Optimize:
   if(!args.verbose)varB.beQuiet();
   varB.optimize(args.flag("veryVerbose"),optM,args.getL("maxIter"),args.getD("limit"),args.getD("limit"));

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
   outF<<"# M "<<M<<"\n"
         "# List includes also 'noise' transcript (first line)\n"
         "# <alpha> - parameter of Dirichlet distribution\n"
         "# <alpha> <beta> - parameters of the marginal Gamma distribution\n"
         "# columns: <mean theta> <alpha> <beta>"<<endl;
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
      string outTypeS = args.getS("outputType");
      string samplesFName = args.getS("outFilePrefix")+".VB" + outTypeS;
      string samplesTmpName = args.getS("outFilePrefix")+".VB"+outTypeS+"TMP"; 
      timer.start(0);
      if(args.verbose)messageF("Generating samples into temporary file %s. ",samplesTmpName.c_str());
      if(!ns_misc::openOutput(samplesTmpName, &outF)) return 1;
      // Samples are generated without the "noise transcript".
      outF<<"# M "<<M-1<<" N "<<args.getL("samples")<<endl;
      varB.generateSamples(args.getL("samples"), outTypeS, trInfo.getShiftedLengths(), &outF);
      outF.close();
      if(args.verbose)timer.split(0);
      if(transposeFiles(vector<string>(1, samplesTmpName), samplesFName, args.verbose, "")){
         if(args.verbose)message("Removing temporary file %s.\n", samplesTmpName.c_str());
         remove(samplesTmpName.c_str());
      }else {
         error("Main: Transposing samples failed.\n");
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
