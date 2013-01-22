#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

#include "FileHeader.h"
#include "TranscriptInfo.h"
#include "ArgumentParser.h"
#include "common.h"

#define Sof(x) (long)x.size()

int main(int argc,char* argv[]){
   buildTime(argv[0],__DATE__,__TIME__);
   string programDescription=
"Converts or normalizes MCMC expression samples.\n\
   [sampleFile] should contain transposed MCMC samples.";
   // Set options {{{
   ArgumentParser args(programDescription,"[sampleFile]",1);
   args.addOptionS("o","outFile","outFileName",1,"Name of the output file.");
//   args.addOptionB("l","log","log",0,"Use logged values.");
   args.addOptionS("a","action","action",1,"Action to perform options: (T2RL - theta to log-rpkm , C2R - coverage to rpkm, R2C - rpkms 2 coverage, LOGNORM - log+normalize, NORM - normalize.");
   args.addOptionD("","Nmap","Nmap",0,"Total number of aligned reads. Or a normalization constant, when normalizing.");
   args.addOptionS("t","trInfoFile","trInfoFileName",0,"File containing transcript information.");
   if(!args.parse(argc,argv))return 0;
   // }}}
   
   long M=0,i,j,m,N;
   double Nmap=0;
   if(args.isSet("Nmap"))Nmap=args.getD("Nmap");
   bool trans;
   ifstream inFile;
   FileHeader fh;
   string geName,trName,action=args.getS("action");
   TranscriptInfo trInfo;

   // Load TR file if necessary {{{
   if((action=="C2R")||(action=="T2RL")||(action=="R2C")||(action=="RPKM2COVERAGE")){
      if((! args.isSet("trInfoFileName")) || (! trInfo.readInfo(args.getS("trInfoFileName")))){
         error("ERROR: Main: Transcript info file read failed. Please provide valid file with --trInfoFile option.");
         return 1;
      }
      M=trInfo.getM();
   } //}}}

   ofstream outFile(args.getS("outFileName").c_str());
   if(!outFile.is_open()){//{{{
      error("ERROR: Main: Unable to open output file");
      return 1;
   }//}}}
   

   inFile.open(args.args()[0].c_str());
   fh.setFile(&inFile);
   if(!fh.samplesHeader(N,m,trans)){//{{{
      error("ERROR: Main: Unable to open samples file");
      return 1;
/*   }else if((trans)&&(! ((action=="--RPKMtoCOVERAGE")||(action=="-R2C")) )){
      error("File should not be transposed");
      return 0;*/ //}}}
   }else if((m==0)||((M!=0)&&(M!=m))){ //{{{
      error("ERROR: Main: Wrong number of transcripts %ld vs %ld\n",M,m);
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

//   vector<double> samples(M,0);
   double sample;
   if((action=="RPKMtoCOVERAGE")||(action=="R2C")){//{{{
      double normC = 1e-9*Nmap;
      if(trans){
         for(j=0;j<M;j++){
            for(i=0;i<N;i++){
               inFile>>sample;
               outFile<<sample * trInfo.effL(j) * normC<<" ";
            }
            outFile<<endl;
         }
      }else{
         for(i=0;i<N;i++){
            for(j=0;j<M;j++){
               inFile>>sample;
               outFile<<sample * trInfo.effL(j) * normC<<" ";
            }
            outFile<<endl;
         }
      } //}}}
   } else if(action=="C2R"){//{{{
      double normC = 1e-9*Nmap;
      if(trans){
         for(j=0;j<M;j++){
            for(i=0;i<N;i++){
               inFile>>sample;
               outFile<<sample / normC / trInfo.effL(j) <<" ";
            }
            outFile<<endl;
         }
      }else{
         for(i=0;i<N;i++){
            for(j=0;j<M;j++){
               inFile>>sample;
               outFile<<sample / normC / trInfo.effL(j) <<" ";
            }
            outFile<<endl;
         }
      } //}}}
   } else if(action=="T2RL"){//{{{
      double normC = log(1000000000);
      if(trans){
         for(j=0;j<M;j++){
            for(i=0;i<N;i++){
               inFile>>sample;
               outFile<<log(sample/trInfo.effL(j))+normC<<" ";
            }
            outFile<<endl;
         }
      }else{
         for(i=0;i<N;i++){
            for(j=0;j<M;j++){
               inFile>>sample;
               outFile<<log(sample/trInfo.effL(j))+normC<<" ";
            }
            outFile<<endl;
         }
      }//}}}
   } else if(action=="NORM"){//{{{
      if(trans){
         for(j=0;j<M;j++){
            for(i=0;i<N;i++){
               inFile>>sample;
               outFile<<sample*Nmap<<" ";
            }
            outFile<<endl;
         }
      }else{
         for(i=0;i<N;i++){
            for(j=0;j<M;j++){
               inFile>>sample;
               outFile<<sample*Nmap<<" ";
            }
            outFile<<endl;
         }
      }//}}}
   } else if(action=="LOGNORM"){ //{{{
      double normC = log(Nmap);
      if(trans){
         for(j=0;j<M;j++){
            for(i=0;i<N;i++){
               inFile>>sample;
               outFile<<log(sample)-normC<<" ";
            }
            outFile<<endl;
         }
      }else{
         for(i=0;i<N;i++){
            for(j=0;j<M;j++){
               inFile>>sample;
               outFile<<log(sample)-normC<<" ";
            }
            outFile<<endl;
         }
      }
   }//}}}
   inFile.close();
   outFile.close();
   if(args.verbose)message("Done.");
   return 0;
}
