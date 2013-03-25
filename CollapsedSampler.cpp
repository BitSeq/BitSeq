#include<sys/time.h>

#include "CollapsedSampler.h"
#include "common.h"
#include "misc.h"

#define DEBUG(x) 

CollapsedSampler::CollapsedSampler(){ //{{{
//   message("COLLAPSED\n");
}//}}}
CollapsedSampler::~CollapsedSampler(){ //{{{
//   message("COLLAPSED DIE\n");
//   Sampler::~Sampler();
}//}}}
/*void CollapsedSampler::init(long n, long m, long samplesTotal, long samplesOut, long Nmap, long Nunmap, const vector<long> &alignI, const vector<TagAlignment> &alignments, const distributionParameters &betaPar, const distributionParameters &dirPar,long seed){//{{{

   Sampler::init(n,m,samplesTotal,samplesOut,Nmap,Nunmap,alignI,alignments,betaPar,dirPar,long seed);

   Z.assign(n,0);
}//}}}*/
void CollapsedSampler::sampleZ(){//{{{
   DEBUG(message("sampleZ\n");)
   int_least32_t i,j,k;
   if(Sof(Z)!=Nmap){
      Z.assign(Nmap,0);
      // init Z&C
      for(i=0;i<Nmap;i++){
         //choose random transcript;
         k = (int_least32_t) (m * uniformDistribution(rng_mt));
         Z[i]=k;
         C[k]++;
      }
   }
#ifdef DoSTATS
   nZ++;
   struct timeval start, end;
   gettimeofday(&start, NULL);
#endif
   vector<double> lphi(m,0); 
   // phi of size M should be enough 
   // because of summing the probabilities for each isoform when reading the data
   double lProbNorm,r,sum,const1;
   int_least32_t readsAlignmentsN;

   const1=m * dir->alpha + Nmap - 1;
   // randomize order: ???
   for(i=0;i<Nmap;i++){  // XXX Nmap-1 ?
      C[Z[i]]--; // use counts without the current one 
      readsAlignmentsN = alignments->getReadsI(i+1) - alignments->getReadsI(i);
      for(j=0, k=alignments->getReadsI(i); j<readsAlignmentsN; j++, k++){
         //message("%ld %lf ",(*alignments)[k].getTrId(),(*alignments)[k].getProb());
         if(alignments->getTrId(k) == 0){
            lphi[j] = alignments->getProb(k) +
               log(beta->beta + Nunmap + C[0]) + 
               log(const1 - C[0]); // this comes from division in "false part"
         }else{
            lphi[j] = alignments->getProb(k) + 
               log(beta->alpha + Nmap - 1 - C[0]) + 
               log(dir->alpha + C[ alignments->getTrId(k) ]); 
               /* 
               / (m * dir->alpha + Nmap - 1 - C[0]) ;
               this term was replaced by (const1 - C[0]) 
               and moved into "true part" as multiplication 
               */
         }
         //message("%lf\n",phi[j]);
      }
      // Normalization constant:
      lProbNorm = ns_math::logSumExp(lphi);
      r = uniformDistribution(rng_mt);
      for(j = 0, sum = 0 ; (sum<r) && (j<=readsAlignmentsN); j++){
         sum += exp(lphi[j] - lProbNorm);
      }
      if(j==0){
         // e.g. if probNorm == 0
         // assign to noise.
         Z[i] = 0;
      } else {
         Z[i] = alignments->getTrId(alignments->getReadsI(i) + j -1);
      }
      C[ Z[i] ]++;
   }
#ifdef DoSTATS
   gettimeofday(&end, NULL);
   tZ += (end.tv_sec-start.tv_sec)*1000*1000+(end.tv_usec-start.tv_usec);
#endif
   DEBUG(message("Z sampled\n");)
}//}}}

void CollapsedSampler::update(){//{{{
   Sampler::update();

   sampleTheta();

   updateSums();
   if((doLog)&&(save))appendFile();
}//}}}
void CollapsedSampler::sample(){//{{{
   Sampler::sample();

   sampleZ();
}//}}}
