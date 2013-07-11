/*
 *
 * Compute posterior variance of samples.
 *
 *
 */
#include<cmath>

using namespace std;

#include "PosteriorSamples.h"
#include "ArgumentParser.h"
#include "common.h"
#include "misc.h"

extern "C" int getVariance(int *argc,char* argv[]){
   string programDescription=
"Estimates variance of MCMC samples from 1 or multiple replicates.\n\
   [sample Files] should contain transposed MCMC samples from replicates.";   
   // Set options {{{
   ArgumentParser args(programDescription,"[sampleFiles]",1);
   args.addOptionS("o","outFile","outFileName",1,"Name of the output file.");
   args.addOptionB("l","log","log",0,"Use logged values.");
   args.addOptionS("t","type","type",0,"Type of variance, possible values: [sample,sqDif] for sample variance or squared difference.","sample");
   args.addOptionS("","norm","normalization",0,"Normalization constants for each input file provided as comma separated list of doubles (e.g. 1.0017,1.0,0.9999 ).");
   if(!args.parse(*argc,argv)){return 0;}
   if(args.verbose)buildTime(argv[0],__DATE__,__TIME__);
   bool doLog=args.flag("log"),logged=false;
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
         if(args.verbose)message("Samples are already logged, computing mean.\n");
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
   if((args.getS("type")=="sqDif")&&(RN>2)&&(args.verbose)){//{{{
      i=0;
      while(args.args()[i]=="C")i++;
      message("using only: %s ",(args.args()[i]).c_str());
      i++;
      while(args.args()[i]=="C")i++;
      message("%s\n",(args.args()[i]).c_str());
   }//}}}
   if(args.verbose)message("replicates: %ld samples: %ld transcripts: %ld\n",RN,N,M);
   
   ofstream outFile(args.getS("outFileName").c_str());
   if(! outFile.is_open()){
      error("Main: File write failed!\n");
      return 1;
   }
   vector<double> mean(M),var(M);
   vector<double> tr,tr2; 
   double m,mSq,count,sqDif;
   bool good=true;
   if(args.getS("type")=="sample"){ //{{{
      for(j=0;j<M;j++){
         if((j%10000==0)&&(j>0)&&args.verbose)message("%ld\n",j);
         
         m = mSq = count = 0;
         for(r=0;r<RN;r++){
            if(cond.getTranscript(r,j,tr,N/RN)){
               for(i=0;i<N/RN;i++){
                  if(doLog){
                     tr[i]=tr[i]<=0?ns_misc::LOG_ZERO:log(tr[i]);
                  }
                  m+=tr[i];
                  mSq += tr[i]*tr[i];
               }
               count+=N/RN;
            }else{
               warning("Error at %ld %ld\n",j,r);
            }
   //         message("%ld  %ld\n",m,count);
         }
         if(count==0){
            warning("no samples for transcript: %ld.\n",j);
            //for(i,Sof(tr))message("%lf "=0;i,Sof(tr))message("%lf "<tr[i];i,Sof(tr))message("%lf "++);
            mean[j] = -47;
            var[j] = -47;
         }else{
            mean[j] = m / count;
            var[j] = mSq/count - m*m/(count*count);
         }
      }//}}}
   }else{ // "sqDif" {{{
      for(j=0;j<M;j++){
         if((j%10000==0)&&(j>0)&&args.verbose)message("%ld\n",j);
         m = sqDif = 0;
         if(RN==1){
            if(! cond.getTranscript(0,j,tr,N)){
               mean[j] = -47;
               var[j] = -47;
               good=false;
               continue;
            }
            tr2.resize(N/2);
            for(i=0;i<N/2;i++)
               tr2[i]=tr[i+N/2];
         }else{
            if(! (cond.getTranscript(0,j,tr,N/2)&&
                  cond.getTranscript(1,j,tr2,N/2))){
               mean[j] = -47;
               var[j] = -47;
               good=false;
               continue;
            }
         }
         if(good){
            for(i=0;i<N/2;i++){
               if(doLog){
                  tr[i]=tr[i]<=0?ns_misc::LOG_ZERO:log(tr[i]);
                  tr2[i]=tr2[i]<=0?ns_misc::LOG_ZERO:log(tr2[i]);
               }
               m+=tr[i]+tr2[i];
               sqDif+=(tr[i]-tr2[i])*(tr[i]-tr2[i]);
            }
            mean[j] = m / N;
            var[j] = sqDif / N; // == ( sqDif / (N/2) ) / 2
         }
      }
   } //}}}
   cond.close();
   
   outFile<<"# Transcripts mean expression and "<<args.getS("type")<<" variance."<<endl;
   outFile<<"# files: ";
   for(r=0;r<RN;r++)outFile<<args.args()[r]<<" ";
   outFile<<endl;
   if(logged)outFile<<"# L -> values logged"<<endl;
   outFile<<"# M "<<M<<endl;
   (outFile<<scientific).precision(9);
   for(i=0;i<M;i++){
      if((mean[i]==-47)&&(var[i]==-47))outFile<<"NaN 0 "<<endl;
      else outFile<<mean[i]<<" "<<var[i]<<endl;
   }
   outFile.close();
   if(args.verbose)message("DONE\n");
   return 0;
}


#ifndef BIOC_BUILD
int main(int argc,char* argv[]){
   return getVariance(&argc,argv);
}
#endif
