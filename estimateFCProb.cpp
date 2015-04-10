/*
 *
 * Compute posterior variance of samples.
 *
 *
 */
#include<cmath>

using namespace std;

#include "ArgumentParser.h"
#include "misc.h"
#include "PosteriorSamples.h"

#include "common.h"

extern "C" int estimateFCProb(int *argc,char* argv[]){
   string programDescription=
"Compute a probability of a specified fold-change from MCMC samples of 2 conditions.\n\
   [sample Files] should contain MCMC samples from estimateExpression.";   
   // Set options {{{
   ArgumentParser args(programDescription,"[sampleFiles]",1);
   args.addOptionS("o","outFile","outFileName",1,"Name of the output file.");
   args.addOptionB("l","log","log",0,"Use logged values.",0);
   args.addOptionD("t","logFCThreshold","logFCThreshold",0,"(Logged) Fold-change cut-off.", 2);
   args.addOptionS("","norm","normalization",0,"Normalization constants for each input file provided as comma separated list of doubles (e.g. 1.0017,1.0,0.9999 ).");
   if(!args.parse(*argc,argv)){return 0;}
   if(args.verbose)buildTime(argv[0],__DATE__,__TIME__);
   bool doLog=args.flag("log"),logged=false;
   double fcThreshold = args.getD("logFCThreshold");
   // }}}
   
   long i,j,r,N,RN,M=0;

   Conditions cond;
   if(! (cond.init("NONE", args.args(), &M, &N))){
      error("Main: Failed loading MCMC samples.\n");
      return 1;
   }
   if(doLog){ 
      logged = true;
      if(cond.logged()){
         doLog=false;
         if(args.verbose)message("Samples are already logged.\n");
      }else{
         if(args.verbose)message("Using logged values.\n");
      }
   }else{
      if(args.verbose)message("NOT using logged values.\n");
      if(cond.logged())logged=true;
   }
   
   if(args.isSet("normalization")){
      if(! cond.setNorm(args.getTokenizedS2D("normalization"))){
         error("Main: Applying normalization constants failed.\n");
         return 1;
      }
   }

   RN=cond.getRN();  

   if(args.verbose)message("replicates: %ld samples: %ld transcripts: %ld\n",RN,N,M);
   
   ofstream outFile(args.getS("outFileName").c_str());
   if(! outFile.is_open()){
      error("Main: File write failed!\n");
      return 1;
   }
   vector<double> prob(M);
   // To summarise: [0] median, [1] mean, [2] std.
   vector< vector<double> > condStat1(M, vector<double>(3));
   vector< vector<double> > condStat2(M, vector<double>(3));
   vector<double> tr(N),tr2(N); 
   
      for(j=0;j<M;j++){
         if((j%10000==0)&&(j>0)&&args.verbose)message("%ld\n",j);
         
         // Prepare 1st set:
         if(cond.getTranscript(0,j,tr,N)){
           if(doLog){
               for(i=0;i<N;i++)
                  tr[i]= log2(tr[i]);
            }
            cond.transcriptStat(tr, j, condStat1);
         }else{
            warning("Error at %ld %ld\n",j,0L);
         }
         // Prepare 2nd set:
         if(cond.getTranscript(1,j,tr2,N)){
            if(doLog){
               for(i=0;i<N;i++)
                  tr2[i]= log2(tr2[i]);
            }
            cond.transcriptStat(tr2, j, condStat2);
         }else{
            warning("Error at %ld %ld\n",j,1L);
         }

         if((tr.size()==0) || (tr2.size()==0)){
            warning("no samples for transcript: %ld.\n",j);
            prob[j] = -47;
         }else{
            // Compute probability of fold-change:
            prob[j] = cond.probFC(tr, tr2, fcThreshold);
         }
      }//}}} 
   cond.close();
   
   outFile<<"# Probability of a fold-change ( logFCThreshold: "<<fcThreshold<<" )"<<endl;
   outFile<<"# files: ";
   for(r=0;r<RN;r++)outFile<<args.args()[r]<<" ";
   outFile<<endl;
   if(logged)outFile<<"# L -> values logged"<<endl;
   outFile<<"# M "<<M<<endl;
   (outFile<<scientific).precision(9);
   for(i=0;i<M;i++){
      if(prob[i]==-47)outFile<<"NaN "<<endl;
      else {
         outFile<<prob[i]<<" ";
         for(int s=0;s<3;s++) outFile<<condStat1[i][s]<<" ";
         for(int s=0;s<3;s++) outFile<<condStat2[i][s]<<" ";
         outFile<<endl;
      }
   }
   outFile.close();
   if(args.verbose)message("DONE\n");
   return 0;
}
