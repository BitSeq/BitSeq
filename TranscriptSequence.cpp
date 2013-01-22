#include<algorithm>
#include<fstream>
#include<ctime>
#include<sstream>

#include "TranscriptSequence.h"
#include "common.h"

#define WORST_SEARCH 10
#define MAX_LINE 190

TranscriptSequence::TranscriptSequence(){//{{{
   srand(time(NULL));
   M=0;
   cM=0;
   gotGeneNames=false;
}//}}}
TranscriptSequence::TranscriptSequence(string fileName){//{{{
   TranscriptSequence();
   readSequence(fileName);
}//}}}
bool TranscriptSequence::readSequence(string fileName){//{{{
   fastaF.open(fileName.c_str());
   if(!fastaF.is_open()){
      error("TranscriptSequence: problem reading transcript file.\n");
      return false;
   }
   trSeqInfoT newTr;
   newTr.use=0;
   newTr.cache=-1;
   string trDesc,geneName;
   long pos;
   istringstream geneDesc;
   gotGeneNames = true;
   while(fastaF.good()){
      while((fastaF.peek()!='>')&&(fastaF.good()))
         fastaF.ignore(1000,'\n');
      if(! fastaF.good())break;
      // skip description line:
      getline(fastaF, trDesc, '\n');
      // look for gene name:
      pos=trDesc.find("gene:");
      if(pos!=(long)string::npos){
         geneDesc.clear();
         geneDesc.str(trDesc.substr(pos+5));
         geneDesc >> geneName;
         geneNames.push_back(geneName);
      }else{
         gotGeneNames = false;
      }
      // remember position:
      newTr.seek=fastaF.tellg();
      trs.push_back(newTr);
   }
   M = trs.size();
   cache.resize(min(M,(long)TRS_CACHE_MAX));
   cachedTrs.resize(min(M,(long)TRS_CACHE_MAX));
   if(fastaF.bad()){
      error("TranscriptSequence problem reading file.\n");
      return false;
   }
   fastaF.clear();
   return true;
}//}}}
long TranscriptSequence::getM(){//{{{
   return M;
}//}}}
const string* TranscriptSequence::getTr(long tr){//{{{
   if((tr<0)||(tr>=M))return NULL;
   trs[tr].use++;
   return &cache[acquireSequence(tr)];
}//}}}
string TranscriptSequence::getSeq(long tr,long start,long l,bool doReverse){//{{{
   if((tr<0)||(tr>=M))return "";
   trs[tr].use++;
   long trI = acquireSequence(tr);
   
   if(start>=(long)cache[trI].size())return string(l,'N');

   string ret;
   if(start<0){
      ret.assign(-start,'N');
      ret+=cache[trI].substr(0,l+start);
   }else{
      ret = cache[trI].substr(start,l);
      if(((long)ret.size()) < l)ret.append(l-ret.size(), 'N');
   }
   if(!doReverse){
      return ret;
   }else{
      reverse(ret.begin(),ret.end());
      for(long i=0;i<l;i++)
         if((ret[i]=='A')||(ret[i]=='a'))ret[i]='T';
         else if((ret[i]=='T')||(ret[i]=='t'))ret[i]='A';
         else if((ret[i]=='C')||(ret[i]=='c'))ret[i]='G';
         else if((ret[i]=='G')||(ret[i]=='g'))ret[i]='C';
      return ret;
   }
}//}}}
long TranscriptSequence::acquireSequence(long tr){//{{{
   if(trs[tr].cache!=-1)return trs[tr].cache;
   long i,newP,j;
   if(cM<TRS_CACHE_MAX){
      newP=cM;
      cM++;
   }else{
      newP=rand()%cM;
      for(i=0;i<WORST_SEARCH;i++){
         j=rand()%cM;
         if(trs[cachedTrs[newP]].use>trs[cachedTrs[j]].use)newP=j;
      }
      trs[cachedTrs[newP]].cache=-1;
      cache[newP].clear();
   }
   fastaF.seekg(trs[tr].seek);
   char line[MAX_LINE];
   while((fastaF.peek()!='>')&&( fastaF.getline(line,MAX_LINE) )){
      cache[newP]+=line;
   }
   if(fastaF.bad()){
      error("TranscriptSequence: Failed reading transcript %ld\n",tr);
      return 0;
   }
   fastaF.clear();
   cachedTrs[newP]=tr;
   trs[tr].cache=newP;
   return newP;
}//}}}

