#ifndef TRANSCRIPTSEQUENCE_H
#define TRANSCRIPTSEQUENCE_H
#include<fstream>
#include<stdint.h>
#include<string>
#include<vector>

using namespace std;

// Max number f transcripts to be cached at a time.
#define TRS_CACHE_MAX 200000

struct trSeqInfoT{
   streampos seek;
   long cache;
   uint_least64_t lastUse;
};

/*
TranscriptSequence class manages fasta file with transcript sequence.
Only up to TRS_CACHE_MAX transcripts are "cached" at a time.
*/
class TranscriptSequence{
   private:
      // Total number of transcripts and number of chached transcripts.
      long M,cM;
      // Counter for the least recently used entry.
      uint_least64_t useCounter;
      // Flag indicating whether it was possible to obtain gene names from the reference file.
      bool gotGeneNames;
      // Gene names for each transcript. (Note - it's not transcript name)
      vector<string> geneNames;
      // Transcript cache information: seek position, use and cache position.
      vector<trSeqInfoT> trs;
      // Cache of transcript sequences.
      vector<string> cache;
      // IDs of transcripts currently in the cache (same order as cache).
      vector<long> cachedTrs;
      // Input stream for the fasta file.
      ifstream fastaF;

      // Read transcript sequence from the file, save to cache and return it's cache index.
      long acquireSequence(long tr);

   public:
      TranscriptSequence();
      // Initialize class and cass readSequence(fileName).
      TranscriptSequence(string fileName);
      // Process input file fileName and record begining of each transcript.
      // The transcript sequence is not read and cached until it is actually requested.
      bool readSequence(string fileName);
      // Return number of transcripts.
      long getM(){ return M; }
      // Return pointer to the transcript sequence. The reference is not persistent.
      // NULL for unknown transcript.
      const string* getTr(long tr);
      // Return sequence from transcript <tr> starting from <start> of length <l>.
      string getSeq(long tr, long start, long l,bool doReverse = false); 
      // Reports whether gene names were extracted from the sequence file.
      bool hasGeneNames(){ return gotGeneNames; }
      // Return pointer to vector containing the geneNames.
      const vector<string> &getGeneNames(){ return geneNames; }
};

#endif
