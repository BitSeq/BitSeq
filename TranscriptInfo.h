#ifndef TRANSCRIPTINFO_H
#define TRANSCRIPTINFO_H
#include<string>
#include<vector>
#include<map>
#include<stdint.h>

using namespace std;

struct transcriptT {//{{{
   string g,t;
   int_least32_t l,gI;
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
      // Number of transcripts, genes.
      long M,G;
      // Flags.
      bool isInitialized, groupedByGenes;
      // Transcript information:
      // gene name, transcript name, length, effective length
      vector<transcriptT> transcripts;
      // Gene information:
      // name, number of transcripts, list of transcripts
      // Length is G+1 after initialization (with dummy record at the end).
      vector<geneT> genes;
      // Populate genes variable with gene information based on gene names saved in transcript information.
      void setGeneInfo();
      string noName;

   public:
      TranscriptInfo();
      // Clears all information.
      void clearTranscriptInfo();
      TranscriptInfo(string fileName);
      // Read info from a file name.
      // Header (# M <num>) is ignored. The file is read until EOF.
      bool readInfo(string fileName);
      // Write transcript into into file. Does not overwrite existing file unless force=true.
      bool writeInfo(string fileName, bool force = false) const;
      bool writeGeneInfo(string fileName) const;
      bool setInfo(vector<string> gNames, vector<string> tNames, vector<long> lengths);
      bool isOK() const{ return isInitialized; }
      long getM() const;
      long getG() const;
      const vector<long> &getGtrs(long i) const;
      long L(long i) const;
      double effL(long i) const;
      const string &trName(long i) const;
      const string &geName(long i) const;
      long geId(long i) const;
      bool genesOrdered() const{ return groupedByGenes; }
      void setEffectiveLength(vector<double> effL);
      // Return pointer to a vector of lengths with transcript IDs starting from 1.
      vector<double> *getShiftedLengths(bool effective = false) const;
      // Update transcript names with new names from the list.
      bool updateTrNames(const vector<string> &trList);
      // Sets gene names to transcripts and calls setGeneInfo to initialize gene information.
      bool updateGeneNames(const vector<string> &geneList);
      bool updateGeneNames(const map<string,string> &trGeneList);
};

#endif
