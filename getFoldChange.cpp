/*
 *
 * Compute Fold Change between expression samples.
 *
 *
 */
#include <cmath>
#include <iostream>

using namespace std;

#include "PosteriorSamples.h"
#include "ArgumentParser.h"
#include "common.h"


int main(int argc,char* argv[]){
   string programDescription=
"Computes log_2 Fold Change from MCMC expression samples.\n\
   [sampleFiles] should contain transposed MCMC samples from replicates.\n\
                  (use --log option if they are not logged)";   
   // Set options {{{
   ArgumentParser args(programDescription,"[sampleFiles]",1);
   args.addOptionS("o","outFile","outFileName",1,"Name of the output file.");
   args.addOptionB("l","log","log",0,"Use logged values.");
//   args.addOptionS("t","type","type",0,"Type of variance, possible values: [sample,sqDif] for sample variance or sqared difference.","sample");
   if(!args.parse(argc,argv))return 0;
   if(args.verbose)buildTime(argv[0],__DATE__,__TIME__);
   // }}}
   bool doLog=args.flag("log");
   if(doLog){
      if(args.verbose)cout<<"Will log expression samples to produce log_2 Fold Chnage."<<endl;
   }else{
      if(args.verbose)cout<<"Assuming samples are logged, producing log_2 Fold Change."<<endl;
   }
   
   long i,j,r,N,RN,M=0,C;

   Conditions cond;
   if(! (cond.init("NONE", args.args(), &C, &M, &N))){
      cerr<<"ERROR: Main: Failed loading MCMC samples."<<endl;
      return 0;
   }
   RN=cond.getRN();  
   if((RN>2)&&(C!=2)){//{{{
      cout<<"Please specify exactly 2 conditions when using more than two sample files.\n";
      cout<<"  such as: [sample Files from first condition] C [sample files from second condition]"<<endl;
      return 0;
   }//}}}
   if(args.verbose)cout<<"Samples: "<<N<<" transcripts: "<<M<<endl;
   
   ofstream outFile(args.getS("outFileName").c_str());
   if(! outFile.is_open()){
      cerr<<"ERROR: Main: File write failed!"<<endl;
      return 0;
   }
   outFile<<"# log_2 Fold Change in expression."<<endl;
   outFile<<"# files: ";
   for(r=0;r<2;r++)outFile<<args.args()[r]<<" ";
   outFile<<endl;
   outFile<<"# T (M rows,N cols)"<<endl; 
   outFile<<"# M "<<M<<endl;
   outFile<<"# N "<<N<<endl;
   vector<double> tr,tr2,res(N); 
   double l2=log(2.0);
   long RC;
   for(j=0;j<M;j++){
      if(args.verbose)progressLog(j,M);
      if(RN==2){
         if(cond.getTranscript(0,j,tr)&&cond.getTranscript(1,j,tr2)){
            for(i=0;i<N;i++){
               if(doLog)outFile<<log(tr2[i]/tr[i])/l2;
               outFile<<(tr2[i]-tr[i])/l2<<" ";
            }
            outFile<<endl;
         }else{
            cerr<<"Failed loading "<<j<<" transcript."<<endl;
         }
      }else{
         // Comparing arithmetic means of log samples which are geometric means of samples
         res.assign(N,0);
         RC = cond.getRC(1);
         for(r=0;r< RC;r++){
            if(cond.getTranscript(1,r,j,tr)){
               for(i=0;i<N;i++)
                  if(doLog)res[i]+=log(tr[i])/RC;
                  else res[i]+=tr[i]/RC;
            }else{
               cerr<<"Failed loading "<<j<<" transcript from condition 1 replicate "<<r<<endl;
            }
         }
         RC = cond.getRC(0);
         for(r=0;r<RC;r++){
            if(cond.getTranscript(0,r,j,tr)){
               for(i=0;i<N;i++)
                  if(doLog)res[i]-=log(tr[i])/RC;
                  else res[i]-=tr[i]/RC;
            }else{
               cerr<<"Failed loading "<<j<<" transcript from condition 0 replicate "<<r<<endl;
            }
         }
         for(i=0;i<N;i++)
            outFile<<res[i]/l2<<" ";
         outFile<<endl;
      }
   }
   cond.close();
   
   outFile.close();
   if(args.verbose)cout<<"DONE"<<endl;
   return 0;
}
