#include<omp.h>
#include<cstring>
#include<cmath>

#include "SimpleSparse.h"

double SimpleSparse::logSumExpVal(long st, long en) const{//{{{
   if(st<0)st = 0;
   if((en == -1) || (en > T)) en = T;
   if(st >= en) return 0;
   long i;
   double sumE = 0, m = val[st];
   for(i = st; i < en; i++)if(val[i] > m)m = val[i];
   for(i = st; i < en; i++)
      sumE += exp(val[i] - m);
   return  m + log(sumE);
}//}}}
void SimpleSparse::sumRows(double res[]) const{//{{{
   long i,r;
   for(r=0;r<N;r++){
      res[r]=0;
      for(i=rowStart[r];i<rowStart[r+1];i++){
         res[r]+=val[i];
      }
   }
}//}}}
void SimpleSparse::sumCols(double res[]) const{//{{{
   memset(res,0,M*sizeof(double));
   for(long i=0;i<T;i++)res[col[i]]+=val[i];
}//}}}
long SimpleSparse::countAboveDelta(double delta) const{//{{{
   long i,count=0;
   #pragma omp parallel for reduction(+:count)
   for(i=0;i<T;i++){
      if(val[i]>delta)count++;
   }
   return count;
}//}}}

void SimpleSparse::getFixed(double delta, vector<long> *counts, vector<long> *unfixed) const{//{{{
   counts->assign(M,0);
   unfixed->clear();
   long i,r;
   for(r=0;r<N;r++){
      for(i=rowStart[r];i<rowStart[r+1];i++){
         if(val[i]>delta)break;
      }
      if(i == rowStart[r+1]){
         unfixed->push_back(r);
      }else{
         (*counts)[col[i]]++;
      }
   }
}//}}}

void SimpleSparse::softmaxInplace(SimpleSparse *res){//{{{
   double logRowSum = 0;
   long i,r;
   #pragma omp parallel for private(i,logRowSum)
   for(r=0;r<N;r++){
      logRowSum = logSumExpVal(rowStart[r],rowStart[r+1]);
      for(i=rowStart[r];i<rowStart[r+1];i++){
         val[i] = val[i] - logRowSum;
         res->val[i] = exp( val[i] );
      }
   }
}//}}}
void SimpleSparse::softmax(SimpleSparse *res) const{//{{{
   double logRowSum = 0;
   long i,r;
   #pragma omp parallel for private(i,logRowSum)
   for(r=0;r<N;r++){
      logRowSum = logSumExpVal(rowStart[r],rowStart[r+1]);
      for(i=rowStart[r];i<rowStart[r+1];i++){
         res->val[i] = exp(val[i] - logRowSum);
      }
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
