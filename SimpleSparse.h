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
   void softmax(SimpleSparse *res);
   void softmaxInplace(SimpleSparse *res);
   void sumCols(double res[]);
   void sumRows(double res[]);
};

#endif
