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

      bool knowNtotal,knowNreads;
      long M,Ntotal,Nreads,currentRead,reservedN;
   public:
      TagAlignments();
      void init(long Nreads = 0,long Ntotal = 0,long M = 0);
      void pushAlignment(long trId, double prob);
      void pushRead();
      void finalizeRead(long *M, long *Nreads, long *Ntotal);
      int_least32_t getTrId(long i) const;
      double getProb(long i) const;
      int_least32_t getReadsI(long i) const;
      long getNreads() const { return Nreads;}
      //void setProb(long i,double p);
}; 

#endif
