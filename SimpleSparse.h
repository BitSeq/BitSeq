#ifndef SIMPLESPARSE_H
#define SIMPLESPARSE_H

#include<stdint.h>

//#define setVal(x,i,y) {for(i=0;i<x->T;i++)x->val[i]=y;}

class SimpleSparse {
   private:
   bool base;
   public:
   long N,M,T; // reads, transcripts, total
   int_least32_t *rowStart,*colStart,*col;
   double *val;

   SimpleSparse(long n,long m, long t);
   SimpleSparse(SimpleSparse *m0);
   ~SimpleSparse();
   void softmax(SimpleSparse *res) const;
   void softmaxInplace(SimpleSparse *res);
   long countAboveDelta(double delta = 0.99) const;
   void sumCols(double res[]) const;
   void sumRows(double res[]) const;
   double logSumExpVal(long st, long en) const;
};

#endif
