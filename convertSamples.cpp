#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

#include "ArgumentParser.h"
#include "common.h"
#include "FileHeader.h"
#include "misc.h"
#include "TranscriptInfo.h"

namespace ns_convertS {
double r2c(double sample, double norm, double len){
   return sample * norm * len;
}
double c2r(double sample, double norm, double len){
   return sample / norm / len;
}
double t2rl(double sample, double Lnorm, double len){
   return log(sample / len) + Lnorm;
}
double norm(double sample, double norm, double len = 1){
   return sample * norm;
}
double logNorm(double sample, double Lnorm, double len = 1){
   return log(sample) + Lnorm;
}
}

int main(int argc,char* argv[]){
   string programDescription=
"Converts or normalizes MCMC expression samples.\n\
   [sampleFile] should contain transposed MCMC samples.";
   // Set options {{{
   ArgumentParser args(programDescription,"[sampleFile]",1);
   args.addOptionS("o","outFile","outFileName",1,"Name of the output file.");
   string actionDesc =
"Action to perform options:\n\
      T2R - theta to rpkm\n\
      R2T - rpkm to theta\n\
      T2RL - theta to log-rpkm\n\
      C2R - counts to rpkm\n\
      R2C - rpkm 2 counts\n\
      NORM - normalize (samples are multiplied by Nmap)\n\
      LOGNORM - log+normalize (samples are multiplied by Nmap and logged).";
   args.addOptionS("a","action","action",1,actionDesc);
   args.addOptionD("","Nmap","Nmap",0,"Total number of aligned reads. Or a normalization constant, when normalizing.");
   args.addOptionS("t","trInfoFile","trInfoFileName",0,"File containing transcript information.");
   if(!args.parse(argc,argv))return 0;
   if(args.verbose)buildTime(argv[0],__DATE__,__TIME__);
   string action = args.getS("action");
   if(! ((action=="T2R")||(action=="T2RL")||(action=="R2T")||(action=="C2R")||
         (action=="R2C")||(action=="NORM")||(action=="LOGNORM"))){
      error("Main: Unknown action: %s.\n",action.c_str());
      return 1;
   }

   // }}}
   
   long M=0,i,j,m,N;
   double Nmap=0;
   // Check Nmap //{{{
   if(args.isSet("Nmap")){
      Nmap=args.getD("Nmap");
      if((action=="T2R")||(action=="T2RL")||(action=="R2T")){
         warning("Main: Using %lf as normalization constant, are you sure about this?\n",Nmap);
      }
   }else{
      if((action=="C2R")||(action=="R2C")){
         error("Main: Need Nmap (total number of mapped reads) for converting from/to counts.\n");
         return 1;
      }
      if((action=="NORM")||(action=="LOGNORM")){
         error("Main: Need Nmap (normalization constant) for normalization.\n");
         return 1;
      }
   }
   //}}}
   // T2R is just C2R with Nmap = 1. //{{{
   if(action=="T2R"){
      action="C2R";
      if(!args.isSet("Nmap"))Nmap = 1;
   }
   if(action=="R2T"){
      action="R2C";
      if(!args.isSet("Nmap"))Nmap = 1;
   }
   //}}}
   bool trans;
   ifstream inFile;
   FileHeader fh;
   string geName,trName;
   TranscriptInfo trInfo;

   // Load TR file if necessary {{{
   if(!((action=="NORM")||(action=="LOGNORM"))){
      if((! args.isSet("trInfoFileName")) || (! trInfo.readInfo(args.getS("trInfoFileName")))){
         error("Main: Transcript info file read failed. Please provide valid file with --trInfoFile option.\n");
         return 1;
      }
      M=trInfo.getM();
   } //}}}

   ofstream outFile;
   if(!ns_misc::openOutput(args,&outFile))return 1;

   inFile.open(args.args()[0].c_str());
   fh.setFile(&inFile);
   if(!fh.samplesHeader(&N,&m,&trans)){//{{{
      error("Main: Unable to open samples file.\n");
      return 1;
/*   }else if((trans)&&(! ((action=="--RPKMtoCOVERAGE")||(action=="-R2C")) )){
      error("File should not be transposed");
      return 0;*/ //}}}
   }else if((m==0)||((M!=0)&&(M!=m))){ //{{{
      error("Main: Wrong number of transcripts %ld vs %ld.\n",M,m);
      return 1;
   }//}}}
   M=m;
   outFile<<"# "<<args.args()[0];
   if((action=="LOGNORM")||(action=="T2RL"))outFile<<"\n# L ";
   if(trans) outFile<<"\n# T (M rows,N cols)";
   else outFile<<"\n# (N rows,M cols)";
   outFile<<"\n# M "<<M<<"\n# N "<<N<<endl;
   outFile.precision(9);
   outFile<<scientific;

   double sample;
   double (*comp)(double a, double b, double c)=NULL;
   double normC=1;
   if(action=="R2C"){
      normC = 1e-9*Nmap;
      comp = &ns_convertS::r2c;
   } else if(action=="C2R"){
      normC = 1e-9*Nmap;
      comp = &ns_convertS::c2r;
   } else if(action=="T2RL"){
      if(args.isSet("Nmap")) normC = log(Nmap * 1e9);
      else normC = log(1e9);
      comp = &ns_convertS::t2rl;
   } else if(action=="NORM"){
      normC = Nmap;
      comp = &ns_convertS::norm;
   } else if(action=="LOGNORM"){
      normC = log(Nmap);
      comp = &ns_convertS::logNorm;
   } else {
      error("Something went wrong.\n");
      return 1;
   }
   if(!((action=="NORM")||(action=="LOGNORM"))){
      if(trans){
         message("TRANS.\n");
         for(j=0;j<M;j++){
            for(i=0;i<N-1;i++){
               inFile>>sample;
               outFile<<comp(sample,normC,trInfo.effL(j))<<" ";
            }
            inFile>>sample;
            outFile<<comp(sample,normC,trInfo.effL(j))<<endl;
         }
      }else{
         message("NO TRANS.\n");
         for(i=0;i<N;i++){
            for(j=0;j<M-1;j++){
               inFile>>sample;
               outFile<<comp(sample,normC,trInfo.effL(j))<<" ";
            }
            inFile>>sample;
            outFile<<comp(sample,normC,trInfo.effL(j))<<endl;
         }
      }
   }else{
      if(trans){
         for(j=0;j<M;j++){
            for(i=0;i<N-1;i++){
               inFile>>sample;
               outFile<<comp(sample,normC,1)<<" ";
            }
            inFile>>sample;
            outFile<<comp(sample,normC,1)<<endl;
         }
      }else{
         for(i=0;i<N;i++){
            for(j=0;j<M-1;j++){
               inFile>>sample;
               outFile<<comp(sample,normC,1)<<" ";
            }
            inFile>>sample;
            outFile<<comp(sample,normC,1)<<endl;
         }
      }
   }
   inFile.close();
   outFile.close();
   if(args.verbose)message("Done.\n");
   return 0;
}
