#include<cstdlib>
#include<fstream>
#include<iomanip>
#include<vector>

using namespace std;

#include "common.h"
#include "FileHeader.h"
#include "transposeFiles.h"

bool transposeFiles(vector<string> inFileNames, string outFileName, bool verbose, string message){
   long M=0,fileN=1,i,j,bufMax,bufN,m,n,totalN,maxN=0,f;
   bool trans=false,transposed=false;
   vector<long> N;
   bufMax=BUFFER_DEFAULT;

   ofstream outFile(outFileName.c_str());
   if(!outFile.is_open()){//{{{
      error("TransposeFile: Unable to open output file\n");
      return 0;
   }//}}}
   //{{{ Opening input
   fileN = inFileNames.size();
   ifstream *inFile = new ifstream[fileN];
   totalN=0;
   FileHeader fh;
   for(i=0;i<fileN;i++){
      inFile[i].open(inFileNames[i].c_str());
      fh.setFile(&inFile[i]);
      m = n = 0;
      if((!fh.samplesHeader(&n,&m,&trans)) || (m == 0) || (n == 0)){
         error("TransposeFile: Unable to read header of file: %s\n",(inFileNames[i]).c_str());
         return false;
      }
      if(N.size()==0){
         M=m;
         transposed=trans;
         maxN=n;
      }else if((M!=m)||(transposed!=trans)){
         error("TransposeFile: Different number of transcripts or file %s is in wrong format.\n",(inFileNames[i]).c_str());
         return false;
      }
      outFile<<"# "<<inFileNames[i]<<" "<<n<<endl;
      N.push_back(n);
      if(n>maxN)maxN=n;
      totalN+=n;
   }
   if(bufMax>M)bufMax=M;
   //}}}

   outFile<<message;
   if(!trans)
      outFile<<"# T (M rows,N cols)";
   else 
      outFile<<"# (N rows,M cols)";
   outFile<<"\n# M "<<M<<"\n# N "<<totalN<<endl;
   outFile.precision(9);
   outFile<<scientific;
   if(verbose)message("Transposing files:\n Samples: %ld Transcripts: %ld Buffer size: %ld\n",totalN,M,bufMax);
   if(!trans){ // {{{
      vector< vector<long> > seeks(fileN,vector<long>(maxN,-1));
      vector<vector<string> > valueBuf(bufMax,vector<string>(totalN));
      long lastBuf = 0, done=0;
      bufN=bufMax;
      if(verbose)messageF("(r");
      for(f=0;f<fileN;f++){
         for(i=0;i<N[f];i++){
            for(j=0;j<bufN;j++) inFile[f]>>valueBuf[j][lastBuf];
            lastBuf++;
            seeks[f][i]=inFile[f].tellg();
            inFile[f].ignore(10000000,'\n');
         }
      }
      if(verbose)messageF(">w.");
      for(j=0;j<bufN;j++){
         for(i=0;i < lastBuf - 1;i++)
            outFile<<valueBuf[j][i]<<" ";
         // Write last value without space.
         outFile<<valueBuf[j][i]<<endl;
      }
      lastBuf=0;
      done=bufN;
      while(done<M){
         bufN=bufMax;
         if(M-done<bufMax)bufN=M-done;
         if(verbose)messageF("r");
         for(f=0;f<fileN;f++){
            for(i=0;i<N[f];i++){
               inFile[f].seekg(seeks[f][i]);
               for(j=0;j<bufN;j++) inFile[f]>>valueBuf[j][lastBuf];
               lastBuf++;
               seeks[f][i]=inFile[f].tellg();
            }
         }
         if(verbose)messageF(">w.");
         for(j=0;j<bufN;j++){
            for(i=0;i < lastBuf - 1;i++)
               outFile<<valueBuf[j][i]<<" ";
            // Write last value without space.
            outFile<<valueBuf[j][i]<<endl;
         }
         lastBuf=0;
         done+=bufN;
      }
      for(f=0;f<fileN;f++)inFile[f].close();
      if(verbose)message(")\n");
   } // }}}
   else{ // if(trans) {{{
      vector<long> seeks(M,-1);
      vector<vector<string> > valueBuf(M,vector<string>(bufMax));
      long done;
      if(verbose)message("(");
      for(f=0;f<fileN;f++){
         seeks.assign(M,-1);
         done = 0;
         while(done<N[f]){
            bufN=bufMax;
            if(bufN>N[f]-done)bufN=N[f]-done;
            if(verbose)messageF("r");
            for(j=0;j<M;j++){
               if(seeks[j]!=-1)inFile[f].seekg(seeks[j]);
               for(i=0;i<bufN;i++){
                  inFile[f]>>valueBuf[j][i];
               }
               seeks[j]=inFile[f].tellg();
               if((j+1<M)&&(seeks[j+1]==-1))inFile[f].ignore(100000000,'\n');
            }
            if(verbose)messageF(">w.");
            for(i=0;i<bufN;i++){
               for(j=0;j < M - 1;j++)
                  outFile<<valueBuf[j][i]<<" ";
               // Write last value without space.
               outFile<<valueBuf[j][i]<<endl;
            }
            done+=bufN;
         }
         inFile[f].close();      
      }
      if(verbose)message(")\n");
   } //}}}
   delete[] inFile;
   outFile.close();
   return true;
}
