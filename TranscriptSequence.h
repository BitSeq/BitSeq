#ifndef TRANSCRIPTSEQUENCE_H
#define TRANSCRIPTSEQUENCE_H
#include<string>
#include<vector>
#include<fstream>

using namespace std;

#define TRS_CACHE_MAX 200000

struct trSeqInfoT{
   streampos seek;
   long use,cache;
};

class TranscriptSequence{
   private:
      long M,cM;
      bool gotGeneNames;
      vector<string> geneNames;
      vector<trSeqInfoT> trs;
      vector<string> cache;
      vector<long> cachedTrs;
      ifstream fastaF;

      long acquireSequence(long tr);

   public:
      TranscriptSequence();
      TranscriptSequence(string fileName);
      bool readSequence(string fileName);
      long getM();
      const string* getTr(long tr);
      string getSeq(long tr, long start, long l,bool doReverse = false); 
      bool hasGeneNames(){ return gotGeneNames; }
      const vector<string>* getGeneNames(){ return &geneNames; }
};

#endif
