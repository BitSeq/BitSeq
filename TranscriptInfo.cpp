#include<fstream>

#include "common.h"
#include"TranscriptInfo.h"

bool TranscriptInfo::writeInfo(string fileName, bool force){//{{{
   ofstream trF;
   if(! force){
      ifstream testF(fileName.c_str());
      if(testF.is_open()){
         testF.close();
         return false;
      }
      testF.close();
   }
   trF.open(fileName.c_str(),ios::out | ios::trunc);
   if(! trF.is_open() ) return false;
   trF<<"# M "<<M<<endl;
   for(long i=0;i<M;i++)
      trF<<transcripts[i].g<<" "<<transcripts[i].t<<" "<<transcripts[i].l<<" "<<transcripts[i].effL<<endl;
   trF.close();
   return true;
}//}}}
bool TranscriptInfo::setInfo(vector<string> gNames,vector<string> tNames, vector<long> lengths){//{{{
   // The sizes have to be equal.
   if((gNames.size()!=tNames.size())||(tNames.size()!=lengths.size())) return false;
   transcriptT newT;
   M = (long) gNames.size();
   // Create new entry for each transcript.
   for(long i=0;i<M;i++){
      newT.g=gNames[i];
      newT.t=tNames[i];
      newT.l=(int_least32_t)lengths[i];
      newT.effL = lengths[i];
      transcripts.push_back(newT);
   }
   // Initialize gene info based on gene names.
   setGeneInfo();
   isInitialized = true;
   return isInitialized;
}//}}}
void TranscriptInfo::setGeneInfo(){//{{{
   // Map of genes: name -> position within gene vector.
   map<string,long> names;
   geneT tmpG;
   long gi=0,i;
   groupedByGenes = true;
   string previousName = "-noname-";
   for(i=0;i<M;i++){
      // If gene name same as previous, then just add new transcript.
      if(transcripts[i].g == previousName){
         genes[gi].m++;
         genes[gi].trs.push_back(i);
      }else{
         // Check whether the gene name is new or was seen before.
         if(names.count(transcripts[i].g) == 0){
            // Prepare entry for new gene, starting with one (current) transcript.
            tmpG.name = transcripts[i].g;
            tmpG.m = 1;
            tmpG.trs = vector<long>(1,i);
            // Add entry to the gene list.
            genes.push_back(tmpG);
            // Set current gene index.
            gi=genes.size()-1;
            // Map gene name to it's index and update previousName.
            names[transcripts[i].g] = gi;
            previousName=transcripts[i].g;
         }else{
            // If gene name was seen before then transcripts are not grouped by genes.
            groupedByGenes=false;
            //warning("TranscriptInfo: Transcripts of gene %ld are not grouped.\n",transcripts[i].g);
            gi = names[transcripts[i].g];
            genes[gi].m++;
            genes[gi].trs.push_back(i);
         }
      }
   }
   G = genes.size();
}//}}}
TranscriptInfo::TranscriptInfo(){ clearTranscriptInfo(); }
void TranscriptInfo::clearTranscriptInfo(){//{{{
   M=G=0;
   isInitialized=false;
   groupedByGenes=true;
   transcripts.clear();
   genes.clear();
}//}}}
TranscriptInfo::TranscriptInfo(string fileName){//{{{
   // TranscriptInfo();
   readInfo(fileName);
}//}}}
bool TranscriptInfo::readInfo(string fileName){//{{{
   clearTranscriptInfo();
   ifstream trFile(fileName.c_str());
   if(!trFile.is_open()){
      error("TranscriptInfo: problem reading transcript file.\n");
      return false;
   }
   transcriptT newT;
   // Read all lines of file ignoring lines starting with #.
   while(trFile.good()){
      while(trFile.good() && (trFile.peek()=='#'))
         trFile.ignore(100000000,'\n');
      if(!trFile.good()) break;
      // Read gene name, tr name and length.
      trFile>>newT.g>>newT.t>>newT.l;
      // Should not hit EOF or any other error yet.
      if(!trFile.good()) break;
      // Read effective length if present:
      while((trFile.peek() == '\t')||(trFile.peek() == ' ')) trFile.get();
      // If end of line is reached then use length as effective length.
      if((trFile.good()) && (trFile.peek() == '\n')) newT.effL = newT.l;
      else trFile>>newT.effL;
      // If the line was OK, then push new entry (EOF when looking for effective length is allowed).
      if(!trFile.fail())
         transcripts.push_back(newT);
      // Ignore rest of the line.
      trFile.ignore(100000000,'\n');
   }
   trFile.close();
   isInitialized = true;
   M = (long)transcripts.size();
   setGeneInfo();
   return isInitialized;
}//}}}
long TranscriptInfo::getM(){//{{{
   return M;
}//}}}
long TranscriptInfo::getG(){//{{{
   return G;
}//}}}
const vector<long>* TranscriptInfo::getGtrs(long i){//{{{
   if(i>G) return NULL;
   return &genes[i].trs;
}//}}}
double TranscriptInfo::effL(long i){//{{{
   if(isInitialized && (i<M))return transcripts[i].effL;
   return 0;
}//}}}
long TranscriptInfo::L(long i){//{{{
   if(isInitialized && (i<M))return transcripts[i].l;
   return 0;
}//}}}
string TranscriptInfo::trName(long i){//{{{
   if(isInitialized && (i<M))return transcripts[i].t;
   return "";
}//}}}
string TranscriptInfo::geName(long i){//{{{
   if(isInitialized && (i<M))return transcripts[i].g;
   return "";
}//}}}
void TranscriptInfo::setEffectiveLength(vector<double> effL){//{{{
   if((long)effL.size() != M){
      warning("TranscriptInfo: Wrong array size for effective length adjustment.\n");
      return;
   }
   // Adjust effective length to similar scale as normal length
   double sumL = 0,sumN = 0,norm;
   for(long i=0;i<M;i++){
      sumN+=effL[i];
      sumL+=transcripts[i].l;
   }
// don't normalize
//   norm = sumL / sumN;
   norm = 1;
   for(long i=0;i<M;i++){
      transcripts[i].effL = effL[i] * norm;
   }
}//}}}
vector<double> *TranscriptInfo::getShiftedLengths(bool effective){//{{{
   vector<double> *Ls = new vector<double>(M+1);
   for(long i=0;i<M;i++){
      if(effective)(*Ls)[i+1] = transcripts[i].effL;
      else (*Ls)[i+1] = transcripts[i].l;
   }
   return Ls;
}//}}}
bool TranscriptInfo::updateGeneNames(const vector<string> &geneList){//{{{
   if((long)geneList.size() != M){
      warning("TranscriptInfo: Number of items in gene list (%ld) does not match number of transcripts (%ld).",geneList.size(),M);
      return false;
   }
   // Copy gene names in the order they are.
   for(long i=0;i<M;i++){
      transcripts[i].g = geneList[i];
   }
   // Initialize gene info.
   setGeneInfo();
   return true;
}//}}}
bool TranscriptInfo::updateGeneNames(const map<string,string> &trGeneList){//{{{
   if((long)trGeneList.size() != M){
      warning("TranscriptInfo: Number of items in gene list (%ld) does not match number of transcripts (%ld).",trGeneList.size(),M);
      return false;
   }
   // Check all transcripts have associated gene name.
   for(long i=0;i<M;i++){
      if(!trGeneList.count(transcripts[i].t)){
         warning("TranscriptInfo: No gene name for transcript [%s].",transcripts[i].t.c_str());
         return false;
      }
   }
   // Set gene names.
   for(long i=0;i<M;i++){
      transcripts[i].g = trGeneList.find(transcripts[i].t)->second;
   }
   // Initialize gene info.
   setGeneInfo();
   return true;
}//}}}
