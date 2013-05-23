#include<algorithm>
#include<cstdlib>
#include<vector>

using namespace std;

#include "common.h"
#include "FileHeader.h"
#include "misc.h"

#include "PosteriorSamples.h"
   
#define Sof(x) (long)x.size()
#define SS second
#define FF first

#define MINUS_INF -47
#define PLUS_INF 1e10

void PosteriorSamples::clear(){//{{{
   N=0;
   M=0;
   norm = 1.0;
   failed=true;
   transposed=true;
   areLogged=false;
}//}}}
bool PosteriorSamples::open(string fileName){//{{{
   if(samplesF.is_open())samplesF.close();
   samplesF.open(fileName.c_str());
   if(!samplesF.is_open()){
      error("PosterioSamples: File open failed: %s\n",(fileName).c_str());
      failed=true;
      return false;
   }
   return true;
}//}}}
bool PosteriorSamples::initSet(long *m,long *n, string fileName){//{{{
   failed=false;
   if(! open(fileName))return false;
   
   FileHeader fh(&samplesF);
   if(!fh.samplesHeader(n,m,&transposed,&areLogged)){
      error("PosteriorSamples: File header reading failed.\n");
      failed=true;
      return false;
   }
   N=*n;
   M=*m;
   return read();
}//}}}
bool PosteriorSamples::read(){//{{{
   if(failed)return false;
   if(transposed){
      lines=vector<long>(M,-1);
      lines[0]=samplesF.tellg();
   }else{
      if(N*M > PS_maxStoredSamples){
         error("PosteriorSamples: Too many samples to store,use trasposed file.\n");
         failed=true;
         return false;
      }
      samples.resize(M,vector<double>(N,0));
      for(long i=0;i<N;i++)
         for(long j=0;j<M;j++)
            samplesF>>samples[j][i];
      if(!samplesF.good()){
         failed=true;
         return false;
      }
   }
   return true;
}//}}}
bool PosteriorSamples::getTranscript(long tr,vector<double> &trSamples){//{{{
   if((tr>=M)||(failed))return false;
   string str;
   bool good=true;
   if(Sof(trSamples)!=N)trSamples.resize(N);
   if(transposed){
      long i;
      if(lines[tr]==-1){
         for(i=0;lines[i+1]!=-1;i++);
         samplesF.seekg(lines[i]);
         while((samplesF.good())&&(i<tr)){
            i++;
            samplesF.ignore(10000000,'\n');
            lines[i]=samplesF.tellg();
         }
      }else{
         samplesF.seekg(lines[tr]);
      }
      for(i=0;(i<N)&&(samplesF.good());i++){
         samplesF>>trSamples[i];
         // apply normalisation.
         trSamples[i] *= norm;
         if(samplesF.eof())break;
         if(samplesF.fail()){
            samplesF.clear();
            samplesF.seekg(-1,ios::cur);
            samplesF>>str;
            if(ns_misc::toLower(str)=="-inf")trSamples[i]=MINUS_INF;
            else if(ns_misc::toLower(str)=="nan")trSamples[i]=PLUS_INF;
            else error("PosteriorSamples: Unknown value: %s in [tr:%ld,pos:%ld]\n",(str).c_str(),tr,i);
            good=false;
         }
      }
      if(i!=N){
         good=false;
         error("PosteriorSamples: Reading failed at position:  [tr:%ld,pos:%ld]\n",tr,i);
      }
   }else{
      trSamples = samples[tr];
      // FIXME(glausp) it is not very efficient to do this every time. 
      // However this part only works for small data files.
      if(norm!=1.0){
         for(long i=0;i<N;i++)trSamples[i] *= norm;
      }
   }
   return good;
}//}}}
void PosteriorSamples::close(){//{{{
   samplesF.close();
   failed=true;
}//}}}


Conditions::Conditions(){//{{{
   mapping=false;
   CN=0;
   C=0;
}//}}}
long Conditions::getIndex(long max){ // {{{returns index, without checking for duplicates
   return rand() % max;
}//}}}
long Conditions::getRC(long c) const { //{{{
   if(c>C)return -1;
   return cIndex[c].SS;
}//}}}
bool Conditions::init(string trFileName, vector<string> filesGot, long *m, long *n){//{{{
   long c;
   return init(trFileName,filesGot,&c,m,n);
}//}}}
bool Conditions::init(string trFileName, vector<string> filesGot, long *c, long *m, long *n){//{{{
   long i,j,x,colN;
   bool sameMs=true;
   vector<string> files;
   cIndex.resize(1,pair<long,long>(0,0));
   for(i=0;i<Sof(filesGot);i++){
      if(filesGot[i]=="C"){
         if((cIndex.end()-1)->SS!=0){
            cIndex.push_back(pair<long,long>(Sof(files),0));
         }
      }else{
         (cIndex.end()-1)->SS++;
         files.push_back(filesGot[i]);
      }
   }
   if((cIndex.end()-1)->SS==0){
      cIndex.pop_back();
   }
   C = Sof(cIndex);
   *c = C;
   //message("File names processed.\n");

   CN = Sof(files);
   samples.resize(CN);
   Ms.resize(CN);
   Ns.resize(CN);
   if(! samples[0].initSet(&Ms[0],&Ns[0],files[0])){
      error("Conditions: file %s failed to open.\n",(files[0]).c_str());
      return false;
   }
   areLogged = samples[0].logged();
   N=Ns[0];
   M=Ms[0];
   for(i=1;i<CN;i++){
      if(! samples[i].initSet(&Ms[i],&Ns[i],files[i])){
         error("Conditions: file %s failed to open.\n",(files[i]).c_str());
         return false;
      }
      if(areLogged != samples[i].logged()){
         error("Conditions: Problem reading %s: some samples are logged and some are not.\n",(files[i]).c_str());
         return false;
      }
      if(M!=Ms[i]){
         sameMs=false;
      }
      if(N>Ns[i])N=Ns[i];
   }
   *n=N;

   ifstream trFile(trFileName.c_str());
   if(! trFile.is_open()){
   // if there is no transcript join file, the we have to make sure that Ms are the same
      if(sameMs){
         M=Ms[0];
         *m=M;
         mapping = false;
         return true;
      }else{
         error("Conditions: Different number of transcripts and missing transcript-join file\n");
         return false;
      }  
   }else{
      FileHeader fh(&trFile);
      if((!fh.transcriptsHeader(&M,&colN))||(M==0)||(colN<CN+1)){
         error("Conditions: Wrong transcript join descriptor file - m: %ld colN: %ld\n",M,colN);
         return false;
      }
      *m=M;
      trMap.resize(M,vector<long>(CN));
      for(i=0;i<M;i++){
         trFile>>x;
         for(j=0;j<colN;j++)
            if(j<CN)trFile >> trMap[i][j];
            else trFile >> x;
      }
      trFile.close();
      sort(trMap.begin(),trMap.end());// sort for faster disc access
      mapping=true;
      return true;
   }
   return false; // we should not get here
}//}}}
bool Conditions::setNorm(vector<double> norms){//{{{
   if((long)norms.size()!=CN){
      error("Conditions: The number of normalization constants does not match number of experiments (files with samples).\n");
      return false;
   }
   for(long i=0;i<CN;i++){
      samples[i].setNorm(norms[i]);
   }
   return true;
}//}}}
bool Conditions::getTranscript(long cond, long rep, long tr, vector<double> &trSamples){//{{{
   if((cond>C)||(rep>cIndex[cond].SS)){
      trSamples.clear();
      return false;
   }
   return getTranscript(rep+cIndex[cond].FF, tr, trSamples);
}//}}}
bool Conditions::getTranscript(long cond, long tr, vector<double> &trSamples){//{{{
   bool status=false;
   static vector<double> tmpSamples;
   if(cond>=CN){
      error("Conditions: Wrong condition request.\n");
      return false;
   }
   if(tr>=M){
      error("Conitions: Wrong transcript request.\n");
      return false;
   }
   if(mapping) tr = trMap[tr][cond];
   if(N != Ns[cond]){
      status = samples[cond].getTranscript(tr, tmpSamples);
      if(Sof(trSamples) != N)trSamples.resize(N);
      for(long i=0;i<N;i++)trSamples[i] = tmpSamples[ getIndex(Ns[cond]) ];
   }else{
      status = samples[cond].getTranscript(tr, trSamples);
   }
   return status;
}//}}}
bool Conditions::getTranscript(long cond, long tr, vector<double> &trSamples, long samplesN){//{{{
   bool status=false;
   static vector<double> tmpSamples;
   if(cond>=CN){
      error("Conditions: Wrong condition request.\n");
      return false;
   }
   if(tr>=M){
      error("Conitions: Wrong transcript request.\n");
      return false;
   }
   if(samplesN > Ns[cond]){
      error("Conitions: Wrong not enough samples.\n");
      return false;
   }
   if(samplesN <1){
      error("Conitions: Wrong number of samples.\n");
      return false;
   }
   if(mapping)tr=trMap[tr][cond];
   if(samplesN != Ns[cond]){
      status = samples[cond].getTranscript(tr, tmpSamples);
      if(Sof(trSamples) != samplesN)trSamples.resize(samplesN);
      for(long i=0;i<samplesN;i++)
         trSamples[i] = tmpSamples[ getIndex(Ns[cond]) ];
   }else{
      status = samples[cond].getTranscript(tr, trSamples);
   }
   return status;
}//}}}
void Conditions::close(){//{{{
   for(long i=0;i<CN;i++){
      samples[i].close();
   }
   cIndex.clear();
}//}}}
