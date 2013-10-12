#ifndef READDISTRIBUTION_H
#define READDISTRIBUTION_H

#include<vector>
#include<map>

using namespace std;

#include "TranscriptInfo.h"
#include "TranscriptSequence.h"
#include "TranscriptExpression.h"

#ifdef BIOC_BUILD

#include "samtoolsHeader.h"
#include <Rinternals.h>

#define bam_init1() ((bam1_t*)S_alloc(1, sizeof(bam1_t)))
// empty destroy, R frees memory itself
#define bam_destroy1(b) 

#else

#include "bam.h"
#include "sam.h"

//#define bam_init1() ((bam1_t*)calloc(1, sizeof(bam1_t)))
/*
#define bam_destroy1(b) do { \
   if (b) { free((b)->data); free(b); }	\
} while (0)
*/

#endif

namespace ns_rD {

// Defaults: {{{
const char LOW_PROB_MISSES  = 6;
const char MAX_NODE_PAR = 2;
const long trSizes [] = { 1334,2104,2977,4389};
const char trSizesN = 4;
const char trNumberOfBins = 20;
const char vlmmNodeDependence [] = { 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, 0};
const char vlmmNodesN = 21;
const char vlmmStartOffset = 8;
const long pows4 [] = {1,4,16,64,256,1024,4096};
//}}}

struct fragmentT{//{{{
   bam1_t *first,*second;
   bool paired;
   fragmentT(){
      first = bam_init1();
      second = bam_init1();
      paired = true;
   }
   ~fragmentT(){
      bam_destroy1(first);
      bam_destroy1(second);
   }
};

typedef fragmentT* fragmentP;
//}}}

class VlmmNode{//{{{
   private:
      long parentsN;
      vector<double> probs;

   public:
      VlmmNode(){parentsN = 0;}
      VlmmNode(long p);
      void setParentsN(long p);
      void update(double Iexp, char b, char bp, char bpp);
      void normalize();
      double getP(char b, char bp, char bpp) const;
      double getPsum(char b) const;
};//}}}

enum biasT { readM_5, readM_3, uniformM_5, uniformM_3, weight_5, weight_3};
enum readT { mate_5, mate_3, FullPair };

} // namespace ns_rD

class ReadDistribution{
   private:
      long procN,M,fragSeen,singleReadLength,minFragLen;
      double lMu,lSigma,logLengthSum,logLengthSqSum;
      long lowProbMismatches;
      bool verbose,uniform,unstranded,lengthSet,gotExpression,normalized;
      bool validLength;
      long warnPos, warnTIDmismatch, warnUnknownTID, noteFirstMateDown;
      TranscriptInfo* trInf;
      TranscriptSequence* trSeq;
      TranscriptExpression* trExp;
      // for each transcript, remember seen fragments in map: length->(sum of probs)
      vector<map<long,double> > trFragSeen5,trFragSeen3;
      // cache for already computed weight norms for:
      //    (single reads 5',3', Pair) x Transcript x Length
      vector<vector<map<long, double> > > weightNorms;
      // position probability arrays (RE-FACTOR to array of 4 vectors)
      vector<vector<vector<double> > > posProb;
      vector<vector<ns_rD::VlmmNode> > seqProb;
      // Cache probabilities for Phred score.
      vector<double> lProbMis;
      vector<double> lProbHit;
      // Mismatch likelihods along read.
      vector<double> lFreqHit;
      vector<double> lFreqMis;
      // Cache length probabilities.
      vector<double> lLengthP,lLengthNorm;
   
      double getLengthLP(long len) const;
      double computeLengthLP(double len) const;
      double getLengthLNorm(long trLen) const;
      void computeLengthProb();
      void updateMismatchFreq(bam1_t *samA);
      void updatePosBias(long pos, ns_rD::biasT bias, long tid, double Iexp);
      void updateSeqBias(long pos, ns_rD::biasT bias, long tid, double Iexp);
      double getPosBias(long start, long end, ns_rD::readT read,
                        long trLen) const;
      double getSeqBias(long pos, ns_rD::readT read, long tid) const;
      inline char getBase(long pos, const string &fSeq) const;
      double getSeqBias(long start, long end, ns_rD::readT read,
                        const string &fSeq) const;
      //inline char complementBase(char base) const;
      double getWeightNorm(long len, ns_rD::readT read, long tid);
      pair<double, double> getSequenceLProb(bam1_t *samA) const;
   public:
      ReadDistribution();
      void setProcN(long procN);
      void writeWarnings();
      bool init(long m, TranscriptInfo* trI, TranscriptSequence* trS, TranscriptExpression* trE, bool unstranded, bool verb = true);
      bool initUniform(long m, TranscriptInfo* trI, TranscriptSequence* trS, bool verb = true);
      void setLowProbMismatches(long m);
      void setLength(double mu, double sigma);
      bool observed(ns_rD::fragmentP frag);
      void normalize();
      void logProfiles(string logFileName = "");
      bool getP(ns_rD::fragmentP frag,double &prob,double &probNoise);
      long getWeightNormCount() const;
      vector<double> getEffectiveLengths();
}; 

#endif
