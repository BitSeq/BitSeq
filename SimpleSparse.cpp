#include<omp.h>
#include<cstring>
#include<cmath>

#include "SimpleSparse.h"

void SimpleSparse::sumRows(double res[]){//{{{
   long i,r;
   for(r=0;r<N;r++){
      res[r]=0;
      for(i=rowStart[r];i<rowStart[r+1];i++){
         res[r]+=val[i];
      }
   }
}//}}}
void SimpleSparse::sumCols(double res[]){//{{{
   memset(res,0,M*sizeof(double));
   for(long i=0;i<T;i++)res[col[i]]+=val[i];
}//}}}

void SimpleSparse::softmaxInplace(SimpleSparse *res){//{{{
   double rowSum = 0;
   long i,r;
   #pragma omp parallel for private(i,rowSum)
   for(r=0;r<N;r++){
      rowSum = 0;
      for(i=rowStart[r];i<rowStart[r+1];i++){
         res->val[i] = exp(val[i]);
         rowSum+=res->val[i];
      }
      for(i=rowStart[r];i<rowStart[r+1];i++){
         res->val[i] /= rowSum;
         val[i] = log( res->val[i] );
      }
   }
}//}}}
void SimpleSparse::softmax(SimpleSparse *res){//{{{
   double rowSum = 0;
   long i,r;
   for(r=0;r<N;r++){
      rowSum = 0;
      for(i=rowStart[r];i<rowStart[r+1];i++){
         res->val[i] = exp(val[i]);
         rowSum+=res->val[i];
      }
      for(i=rowStart[r];i<rowStart[r+1];i++)res->val[i]/=rowSum;
   }
}//}}}

SimpleSparse::SimpleSparse(long n,long m, long t){//{{{
   N=n;
   M=m;
   T=t;
   val = new double[T];
   base = true; // base matrix with it's own col & rowStart information
   col = new int_least32_t[T];
   rowStart = new int_least32_t[N+1];
   //colStart = new long[M+1];
}//}}}
SimpleSparse::SimpleSparse(SimpleSparse *m0){//{{{
   N=m0->N;
   M=m0->M;
   T=m0->T;
   val = new double[T];
   base = false; // use col & rowStart information from the base matrix m0
   col = m0->col;
   rowStart = m0->rowStart;
   /*col = new long[T];
   rowStart = new long[N+1];
   memcpy(col, m0->col, T*sizeof(long));
   memcpy(rowStart, m0->rowStart, (N+1)*sizeof(long));
   */
}//}}}
SimpleSparse::~SimpleSparse(){//{{{
   delete[] val;
   if(base){
      // BEWARE there could be other matrices using this data 
      delete[] col;
      delete[] rowStart;
   }
}//}}}
