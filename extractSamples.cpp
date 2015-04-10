/*
 *
 * Extract samples of given transcripts.
 *
 *
 */
#include<iostream>
#include<cstdlib>
#include<algorithm>

using namespace std;

#include "PosteriorSamples.h"
#include "ArgumentParser.h"
#include "common.h"

#define Sof(x) (long)x.size()

vector <long> tokenizeL(const string &input,const string &space = " "){//{{{
   vector <long> ret;
   long pos=0,f=0,n=input.size();
   while((pos<n)&&(f<n)&&(f>=0)){
      f=input.find(space,pos);
      if(f==pos)pos++;
      else{
         if((f <n)&&(f>=0)){
            ret.push_back(atoi(input.substr(pos,f-pos).c_str()));
            pos=f+1;
         }
      }
   }
   if(pos<n)ret.push_back(atoi(input.substr(pos,n-pos).c_str()));
   return ret;
} //}}}

extern "C" int extractSamples(int *argc,char* argv[]){
   srand(time(NULL));
   string programDescription=
"Extracts MCMC samples of selected transcripts.\n\
   [sampleFiles] should contain transposed MCMC samples.";   
   // Set options {{{
   ArgumentParser args(programDescription,"[sampleFiles]",1);
   args.addOptionS("o","outFile","outFileName",1,"Name of the output file.");
   args.addOptionS("L","list","list",0,"Comma delimited list of ZERO-BASED transcript ids (i.e. lines) which should be extracted: 0,17,47,1024,4777");
   args.addOptionL("r","random","randomN",0,"Choose random [randomN] transcripts.");
   if(!args.parse(*argc,argv))return 0;
   if(args.verbose)buildTime(argv[0],__DATE__,__TIME__);
   // }}}
   long i,j,c,C,N,M=0,S;
   vector<long> trList;
   Conditions samples;

   // Initialize samples reader
   if( (!samples.init("NONE", args.args(), &C, &M, &N)) || (C<=0) || (M<=0) || (N<=0)){
      cerr<<"ERROR: Main: Failed loading MCMC samples."<<endl;
      return 1;
   }
   C=samples.getRN();
   if(args.isSet("list")){
      // Process transcripts list:
      trList = tokenizeL(args.getS("list"),",");
      sort(trList.begin(),trList.end());
      // Erase invalid and duplicate IDs
      for(i=0;i<Sof(trList);i++){
         if((trList[i]<0)||(trList[i]>=M)||((i>0)&&(trList[i]==trList[i-1]))){
            trList.erase(trList.begin()+i);
            i--;
         }
      }   
      S=Sof(trList);
      if(S==0){
         cerr<<"ERROR: Main: No valid transcript IDs supplied."<<endl;
         return 1;
      }
   }else if(args.isSet("randomN")){
      // Create list of [randomN] random transcripts
      S = args.getL("randomN");
      if((S<=0)||(S>M)){
         cerr<<"ERROR: Main: Wrong number of transcripts ot output: "<<S<<"."<<endl;
         return 1;
      }
      for(i=0;i<S;i++){
         j = rand()%M;
         while(find(trList.begin(),trList.end(),j)!=trList.end())
            j = rand()%M;
         trList.push_back(j);
      }
      sort(trList.begin(),trList.end());
   }else{
      cerr<<"ERROR: Main: Need to specify at least one of --list or --random."<<endl;
      return 1;
   }
   if(args.verbose)cout<<"C: "<<C<<" samples: "<<N<<"\ntranscripts: "<<M<<"\nselected: "<<S<<endl;
   
   // Open output file and write header
   ofstream outFile(args.getS("outFileName").c_str());
   if(! outFile.is_open()){
      cerr<<"ERROR: Main: File write failed!"<<endl;
      return 1;
   }
   outFile<<"# Selected transcripts from: ";
   for(i=0;i<C;i++)outFile<<args.args()[i]<<",";
   outFile<<endl;
   outFile<<"# transcripts(zero-based): "<<trList[0];
   for(i=1;i<S;i++)outFile<<","<<trList[i];
   outFile<<"\n# T (M rows,N cols)\n";
   outFile<<"# C "<<C<<" (conditions)\n";
   outFile<<"# M "<<S<<" (out of: "<<M<<")\n# N "<<N<<endl;
   outFile.precision(9);
   outFile<<scientific;

   // Copy samples
   vector<double> tr;
   for(j=0;j<S;j++){
      if(args.verbose)cout<<trList[j]<<" ";
      cout.flush();
      for(c=0;c<C;c++){
         samples.getTranscript(c,trList[j], tr);
         for(i=0;i<N;i++)outFile<<tr[i]<<" ";
         outFile<<endl;
      }
   }
   outFile.close();
   if(args.verbose)cout<<"DONE"<<endl;
   return 0;
}
