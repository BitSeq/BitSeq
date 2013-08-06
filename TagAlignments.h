#ifndef TAGALIGNMENTS_H
#define TAGALIGNMENTS_H

#include<stdint.h>
#include<vector>

using namespace std;

// Probabilities are stored in log scale.

class TagAlignments{
   private:
      vector<int_least32_t> trIds;
      vector<double> probs;
      vector<int_least32_t> readIndex;
      vector<int_least32_t> readsInIsoform;

      bool storeLog,knowNtotal,knowNreads;
      long M,Ntotal,Nreads,currentRead,reservedN;
   public:
      // Constructor, can specify whether the probabilities should be stored in log space.
      TagAlignments(bool storeL = true);
      // Initialize reader. For non-zero arguments, also reserves some memory.
      void init(long Nreads = 0,long Ntotal = 0,long M = 0);
      // Add alignment for currently processed read.
      void pushAlignment(long trId, double prob);
      // Add alignment for currently processed read, with probability in log scale.
      void pushAlignmentL(long trId, double lProb);
      // Finish processing current read and move onto new read. 
      void pushRead();
      // Finalizes reading reads and sets N, Nreads, Ntotal.
      void finalizeRead(long *M, long *Nreads, long *Ntotal);
      // Return TrID of i-th alignment.
      int_least32_t getTrId(long i) const;
      // Return alignment probability of i-th alignment as it is stored.
      // (if it is stored in log space, return log-probability)
      double getProb(long i) const;
      // Get index for i-th read's alignments.
      int_least32_t getReadsI(long i) const;
      // Get number of reads.
      long getNreads() const { return Nreads;}
      // Get number of transcripts.
      long getM() const { return M; }
      // Get (total) number of alignments.
      long getNhits() const { return Ntotal; }
}; 

#endif
