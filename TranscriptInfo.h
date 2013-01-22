#ifndef TRANSCRIPTINFO_H
#define TRANSCRIPTINFO_H
#include<string>
#include<vector>
#include<map>
#include<stdint.h>

using namespace std;

struct transcriptT {//{{{
   string g,t;
   int_least32_t l;
   double effL;
   bool operator< (const transcriptT& d2) const{
      if(g==d2.g)return t<d2.t;
      return g<d2.g;
   }
};//}}}

struct geneT {//{{{
   string name;
   int_least32_t m;
   vector<long> trs;
};//}}}

class TranscriptInfo{
   private:
      long M,G;
      bool ok,ordered;
      vector<transcriptT> transcripts;
      vector<geneT> genes;
      void setGeneInfo();
   public:
      TranscriptInfo();
      void clearTranscriptInfo();
      TranscriptInfo(string fileName);
      bool readInfo(string fileName);
      bool writeInfo(string fileName, bool force = false);
      bool setInfo(vector<string> gNames,vector<string> tNames, vector<long> lengths);
      bool isOK();
      long getM();
      long getG();
      const vector<long>* getGtrs(long i);
      long L(long i);
      double effL(long i);
      string trName(long i);
      string geName(long i);
      bool genesOrdered();
      void setEffectiveLength(vector<double> effL);
      vector<double> *getShiftedLengths(bool effective = false);
      bool updateGeneNames(vector<string> geneList);
      bool updateGeneNames(map<string,string> trGeneList);
};

#endif
