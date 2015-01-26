#include<algorithm>
#include<fstream>
#include<set>
#include<sstream>

#include "TranscriptSequence.h"

#include "misc.h"

#include "common.h"

// Number of times we randomly probe for old cache record.
// CR: #define WORST_SEARCH_N 10

TranscriptSequence::TranscriptSequence(){//{{{
   // CR: srand(time(NULL));
   M=0;
   cM=0;
   gotGeneNames=false;
   // CR: useCounter = 0;
}//}}}
TranscriptSequence::TranscriptSequence(string fileName, refFormatT format){//{{{
   TranscriptSequence();
   readSequence(fileName,format);
}//}}}
bool TranscriptSequence::readSequence(string fileName, refFormatT format){//{{{
   fastaF.open(fileName.c_str());
   if(!fastaF.is_open()){
      error("TranscriptSequence: problem reading transcript file.\n");
      return false;
   }
   trSeqInfoT newTr;
   // CR: newTr.lastUse=0;
   // CR: newTr.cache=-1;
   string trDesc,geneName;
   long pos;
   istringstream geneDesc;
   trNames.clear();
   geneNames.clear();
   gotGeneNames = true;
   // Record trNames only from gencode ref.
   gotTrNames = (format == GENCODE);
   while(fastaF.good()){
      while((fastaF.peek()!='>')&&(fastaF.good()))
         fastaF.ignore(1000,'\n');
      if(! fastaF.good())break;
      // Read description line:
      getline(fastaF, trDesc, '\n');
      // look for gene name if previous lines had gene name:
      if(gotGeneNames){
         if(format == GENCODE){
            vector<string> lineTokens = ns_misc::tokenize(trDesc,"|");
            if(lineTokens.size()>1){
               geneNames.push_back(lineTokens[1]);
               trNames.push_back(lineTokens[0].substr(1));
            }else{
               gotGeneNames = false;
               gotTrNames = false;
            }
         }else{ // format == STANDARD
            pos=min(trDesc.find(" gene:"),trDesc.find("gene="));
            if(pos!=(long)string::npos){
               geneDesc.clear();
               geneDesc.str(trDesc.substr(pos+6));
               geneDesc >> geneName;
               geneNames.push_back(geneName);
            }else{
               gotGeneNames = false;
            }
         }
      }
      // remember position:
      newTr.seek=fastaF.tellg();
      trs.push_back(newTr);
   }
   // Exit if there was an error while reading the file.
   if(fastaF.bad()){
      error("TranscriptSequence: problem reading file.\n");
      return false;
   }
   M = trs.size();
   // Allocate cache for all.
   cache.resize(M);
   //cache.resize(min(M,(long)TRS_CACHE_MAX));
   //cachedTrs.resize(min(M,(long)TRS_CACHE_MAX));
   // Clear eof flag from input stream.
   fastaF.clear();
   return loadSequence();
}//}}}
bool TranscriptSequence::loadSequence(){//{{{
   cache.resize(M);
   string seqLine;
   for(long tr=0;tr<M;tr++){
      // Set input stream to transcript's position.
      fastaF.seekg(trs[tr].seek);
      // Read line by line until reaching EOF or next header line '>'.
      while((fastaF.peek()!='>')&&( getline(fastaF,seqLine,'\n').good())){
         cache[tr]+=seqLine;
      }
      if(fastaF.bad()){
         error("TranscriptSequence: Failed reading transcript %ld\n",tr);
         return false;
      }
      // Clear flags (just in case).
      fastaF.clear();
   }
   return true;
}//}}}
long TranscriptSequence::getG() const{//{{{
   if(!gotGeneNames)return 0;
   return (set<string>(geneNames.begin(),geneNames.end())).size();
}//}}}
const string &TranscriptSequence::getTr(long tr) const{//{{{
   if((tr<0)||(tr>=M))return noneTr;
   // Return pointer to the sequence in cache.
   return cache[tr];
   /* Used with cacheing. {{{
   // Update last use info.
   trs[tr].lastUse = useCounter++;
   return cache[acquireSequence(tr)];
   }}} */
}//}}}
string TranscriptSequence::getSeq(long trI,long start,long l,bool doReverse) const{//{{{
   // Return empty string for unknown transcript.
   if((trI<0)||(trI>=M))return "";
   /* Used with cacheing. {{{
   // Update last use info.
   trs[tr].lastUse = useCounter++;
   // Get position within cache.
   long trI = acquireSequence(tr);
   }}} */
   
   // If position is not within the sequence, return Ns.
   if(start>=(long)cache[trI].size())return string(l,'N');

   string ret;
   // Copy appropriate sequence, fill up the rest with Ns.
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
      // For reverse return reversed string with complemented bases.
      reverse(ret.begin(),ret.end());
      for(long i=0;i<l;i++)
         if((ret[i]=='A')||(ret[i]=='a'))ret[i]='T';
         else if((ret[i]=='T')||(ret[i]=='t'))ret[i]='A';
         else if((ret[i]=='C')||(ret[i]=='c'))ret[i]='G';
         else if((ret[i]=='G')||(ret[i]=='g'))ret[i]='C';
      return ret;
   }
}//}}}
/* long TranscriptSequence::acquireSequence(long tr){//{{{
   // If the sequence is stored in cache then just return it's cache index.
   if(trs[tr].cache!=-1)return trs[tr].cache;
   long i,newP,j;
   // See if cache is full.
   if(cM<TRS_CACHE_MAX){
      // If cache limit not reached, just add new sequence.
      newP=cM;
      cM++;
   }else{
      // If cache is full, look at WORST_SEARCH_N positions and choose the one least used.
      newP=rand()%cM;
      for(i=0;i<WORST_SEARCH_N;i++){
         j=rand()%cM;
         if(trs[cachedTrs[newP]].lastUse > trs[cachedTrs[j]].lastUse)newP=j;
      }
      // "remove" the transcript from position newP from cache.
      trs[cachedTrs[newP]].cache=-1;
      cache[newP].clear();
   }
   // Set input stream to transcript's position.
   fastaF.seekg(trs[tr].seek);
   string seqLine;
   // Read line by line until reaching EOF or next header line '>'.
   while((fastaF.peek()!='>')&&( getline(fastaF,seqLine,'\n').good())){
      cache[newP]+=seqLine;
   }
   if(fastaF.bad()){
      error("TranscriptSequence: Failed reading transcript %ld\n",tr);
      return 0;
   }
   // Clear flags.
   fastaF.clear();
   // Update cache information.
   cachedTrs[newP]=tr;
   trs[tr].cache=newP;
   // Return transcripts index within cache.
   return newP;
}//}}} */
