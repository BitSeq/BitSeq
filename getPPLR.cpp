#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <cmath>

using namespace std;

#include "common.h"
#include "PosteriorSamples.h"
#include "ArgumentParser.h"

#define PERCENT 1

long getIndex(long N){ // returns index, without checking for duplicates
   return rand() % N;
}

int main(int argc,char* argv[]){
   buildTime(argv[0],__DATE__,__TIME__);
   string programDescription=
"Computes PPLR from MCMC expression samples.\n\
   (the probability of second condition being up-regulated)\n\
   [sampleFiles] should contain transposed MCMC samples from replicates.\n\
                  (use --log option if they are not logged)";   
   // Set options {{{
   ArgumentParser args(programDescription,"[sampleFiles]",1);
   args.addOptionS("o","outFile","outFileName",1,"Name of the output file.");
   args.addOptionB("l","log","log",0,"Use logged values.");
   args.addOptionB("d","distribution","distribution",0,"Produce whole distribution of differences.");
   args.addOptionS("s","selectFile","selectFileName",0,"File containing list of selected transcript IDs, only these will be reported. Only works with --distribution option.");
   if(!args.parse(argc,argv))return 0;
   // }}}

   long i,m,N,M;
   bool logged = args.flag("log"),getAll=false;
   vector<long> trSelect;
   if(! args.isSet("selectFileName")){
      getAll=true;
   }else{
      ifstream selectF (args.getS("selectFileName").c_str());
      if(! selectF.is_open()){
         cerr<<"ERROR: Main: Failed loading selected transcripts."<<endl;
         return 1;
      }
      selectF>>m;
      while(selectF.good()){
         trSelect.push_back(m);
         selectF>>m;
      }
      selectF.close();
      sort(trSelect.begin(),trSelect.end());
   }

   Conditions cond;
   if(! cond.init(M,N,"NONE",args.args())){
      cerr<<"ERROR: Main: Failed loading conditions."<<endl;
      return 1;
   }
   cout<<"M "<<M<<"   N "<<N<<endl;
   ofstream outFile(args.getS("outFileName").c_str());
   if(! outFile.is_open()){
      cerr<<"ERROR: Main: File write probably failed!"<<endl;
      return 1;
   }
   if(getAll){
      trSelect.resize(M);
      for(i=0;i<M;i++)trSelect[i]=i;
   }
   
   vector<vector<double> > tr(2);
   vector<double> difs;
   // take onlu 3/4 of samples -> they will be shuffled
   long subN = N*3/4;
   double pplr,expr,avfc;
   if(! args.flag("distribution")){
      cout<<"Counting PPLR"<<endl;
      outFile<<"# pplr all"<<endl;
//      outFile<<"# T "<<endl;
      outFile<<"# M "<<M<<endl;
      outFile<<"# pplr average_expression average_foldchange"<<endl;
//      outFile<<"# N "<<N<<endl;
      for(m=0;m<M;m++){
         progressLog(m,M);
         cond.getTranscript(0,m,tr[0],subN);
         cond.getTranscript(1,m,tr[1],subN);
         difs.resize(tr[0].size());
         pplr = 0;
         expr=0;
         avfc=0;
         for(i=0;i<subN;i++){
            if(logged){ 
               difs[i]=log(tr[1][i])-log(tr[0][i]);
               avfc+=tr[1][i]/tr[0][i];
            }else{
               difs[i]=tr[1][i]-tr[0][i];
               avfc+=exp(tr[1][i])/exp(tr[0][i]);
            }
            if(difs[i]>0)pplr+=1;
            if(logged)expr+=log(tr[0][i])+log(tr[1][i]);
            else expr+=tr[0][i]+tr[1][i];
         }
         pplr /= subN;
         expr /= (2.0*subN);
         avfc /= subN;
         outFile<<pplr<<" "<<expr<<" "<<avfc<<" ";
         sort(difs.begin(),difs.end());
         // ADD PERCENTILES
         outFile<<endl;
      }
   }else{
      cout<<"Computing PLR distribution"<<endl;
      long selectM = trSelect.size();
      outFile<<"# plr distribution"<<endl;
      outFile<<"# T "<<endl;
      outFile<<"# M "<<selectM<<endl;
      outFile<<"# N "<<N<<endl;
      outFile<<"# first column - tr number"<<endl;
      for(m=0;m<selectM;m++){
         if((selectM>10)&&(m%((long)(selectM/10))==0)&&(m!=0))cout<<"# "<<m<<" done"<<endl;
         cond.getTranscript(0,trSelect[m],tr[0]);
         cond.getTranscript(1,trSelect[m],tr[1]);
         outFile<<trSelect[m]<<" ";
         for(i=0;i<N;i++){
            if(logged)
               outFile<<log(tr[1][i])-log(tr[0][i])<<" ";
            else
               outFile<<tr[1][i]-tr[0][i]<<" ";
         }
         outFile<<endl;
      }
   }
   outFile.close();
   cond.close();
   return 0;
}
