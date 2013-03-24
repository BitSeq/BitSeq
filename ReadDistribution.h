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


// Defaults: {{{
#define LOW_PROB_MISSES 6
#define MAX_NODE_PAR 2
const long trSizes [] = { 1334,2104,2977,4389};
const long trSizesN = 4;
const long trNumberOfBins = 20;
const long vlmmNodeDependence [] = { 0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0, 0};
const long vlmmNodesN = 21;
const long vlmmStartOffset = 8;
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
      double getP(char b, char bp, char bpp);
      double getPsum(char b);
};//}}}

enum biasT { readM_5, readM_3, uniformM_5, uniformM_3, weight_5, weight_3};
enum readT { mate_5, mate_3, FullPair };

class ReadDistribution{
   private:
      long M,fragSeen,singleReadLength,minFragLen;
      double lMu,lSigma,logLengthSum,logLengthSqSum;
      long lowProbMismatches;
      bool verbose,uniform,lengthSet,gotExpression,normalized;
      bool warnPos, warnTIDmismatch, validLength;
      TranscriptInfo* trInf;
      TranscriptSequence* trSeq;
      TranscriptExpression* trExp;
      // for each transcript, remember seen fragments in map: length->(sum of probs)
      vector<map<long,double> > trFragSeen5,trFragSeen3;
      // cache for already computed weight norms for single reads 4',3', Pair x Transcript x Length
      vector<vector<map<long, double> > > weightNorms;
      // position probability arrays (RE-FACTOR to array of 4 vectors)
      vector<vector<vector<double> > > posProb;
      vector<vector<VlmmNode> > seqProb;
   
      double getLengthLP(double len);
      double getLengthLNorm(double trLen);
      void updatePosBias(long pos, biasT bias, long tid, double Iexp);
      void updateSeqBias(long pos, biasT bias, long tid, double Iexp);
      double getPosBias(long pos, readT read, long tid);
      double getSeqBias(long pos, readT read, long tid);
      double getWeightNorm(long len, readT, long tid);
      pair<double, double> getSequenceLProb(bam1_t *samA);
   public:
      ReadDistribution(long m);
      bool init(TranscriptInfo* trI, TranscriptSequence* trS, TranscriptExpression* trE, bool verb = true);
      bool initUniform(TranscriptInfo* trI, TranscriptSequence* trS, bool verb = true);
      void setLowProbMismatches(long m);
      void setLength(double mu, double sigma);
      void observed(fragmentP frag);
      void normalize();
      void logProfiles(string logFileName = "");
      void getP(fragmentP frag,double &prob,double &probNoise);
      vector<double> getEffectiveLengths();
}; 

#endif
