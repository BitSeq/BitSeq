#ifndef TRANSCRIPTSEQUENCE_H
#define TRANSCRIPTSEQUENCE_H
#include<fstream>
#include<stdint.h>
#include<string>
#include<vector>

using namespace std;

/*
   Lines commented with CR: -> cache related.
   This was commented out when cacheing was removed.
*/

// Max number f transcripts to be cached at a time.
// CR: #define TRS_CACHE_MAX 200000

struct trSeqInfoT{
   streampos seek;
// CR: long cache;
// CR: uint_least64_t lastUse;
};

enum refFormatT { STANDARD, GENCODE };

/*
TranscriptSequence class manages fasta file with transcript sequence.
// CR: Only up to TRS_CACHE_MAX transcripts are "cached" at a time.
*/
class TranscriptSequence{
   private:
      // Total number of transcripts and number of cached transcripts.
      long M,cM;
      // Flag indicating whether it was possible to obtain gene names from the reference file.
      bool gotGeneNames;
      // Gene names for each transcript. (Note - it's not transcript name)
      vector<string> geneNames;
      // Transcript cache information: seek position, use and cache position.
      vector<trSeqInfoT> trs;
      // Cache of transcript sequences.
      vector<string> cache;
      // Input stream for the fasta file.
      ifstream fastaF;
      // Empty transcript.
      string noneTr;

      // Counter for the least recently used entry.
      // CR: uint_least64_t useCounter;
      // IDs of transcripts currently in the cache (same order as cache).
      // CR: vector<long> cachedTrs;
      // Read transcript sequence from the file, save to cache and return it's cache index.
      // CR: long acquireSequence(long tr);

      bool loadSequence();
   public:
      TranscriptSequence();
      // Initialize class and cass readSequence(fileName).
      TranscriptSequence(string fileName, refFormatT format = STANDARD);
      // Process input file fileName and record beginning of each transcript.
      bool readSequence(string fileName, refFormatT format = STANDARD);
      // Return number of transcripts.
      long getM() const{ return M; }
      // Return number of gene names.
      long getG() const{ return geneNames.size(); }
      // Return pointer to the transcript sequence. The reference is not persistent.
      // NULL for unknown transcript.
      const string &getTr(long tr) const;
      // Return sequence from transcript <tr> starting from <start> of length <l>.
      string getSeq(long trI, long start, long l,bool doReverse = false) const; 
      // Reports whether gene names were extracted from the sequence file.
      bool hasGeneNames() const{ return gotGeneNames; }
      // Return pointer to vector containing the geneNames.
      const vector<string> &getGeneNames() const{ return geneNames; }
};

#endif
