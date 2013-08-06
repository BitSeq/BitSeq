#ifdef DoSTATS
#include<sys/time.h>
#endif

#include "CollapsedSampler.h"
#include "common.h"

void CollapsedSampler::sampleZ(){//{{{
   int_least32_t i,j,k,unfI;
   // Resize Z and initialize if not big enough. {{{
   if((long)Z.size() != Nmap){
      Z.assign(Nmap,0);
      // init Z&C
      unfI = 0;
      for(i=0;i<Nmap;i++){
         if(someFixed){
         // If some are fixed, only used indexes in unfixed.
            if(unfI>=(long)unfixed.size())break;
            if(i<unfixed[unfI])continue;
            if(i>unfixed[unfI]){unfI++;continue;}
            unfI++;
         }
         //choose random transcript;
         k = (int_least32_t) (m * uniformDistribution(rng_mt));
         Z[i]=k;
         C[k]++;
      }
   }//}}}
   // TimeStats {{{
#ifdef DoSTATS
   nZ++;
   struct timeval start, end;
   gettimeofday(&start, NULL);
#endif
   // }}}
   vector<double> phi(m,0); 
   // phi of size M should be enough 
   // because of summing the probabilities for each isoform when reading the data
   double probNorm,r,sum,const1a,const1b,const2a;
   int_least32_t readsAlignmentsN;

   const1a = beta->beta + Nunmap;
   const1b = m * dir->alpha + Nmap - 1;
   const2a = beta->alpha + Nmap - 1;
   // randomize order: ???
   unfI = 0;
   for(i=0;i<Nmap;i++){
      if(someFixed){
      // If some are fixed, only used indexes in unfixed.
         if(unfI>=(long)unfixed.size())break;
         if(i<unfixed[unfI])continue;
         if(i>unfixed[unfI]){unfI++;continue;}
         unfI++;
      }
      probNorm=0;
      C[Z[i]]--; // use counts without the current one 
      readsAlignmentsN = alignments->getReadsI(i+1) - alignments->getReadsI(i);
      for(j=0, k=alignments->getReadsI(i); j<readsAlignmentsN; j++, k++){
         //message("%ld %lf ",(*alignments)[k].getTrId(),(*alignments)[k].getProb());
         if(alignments->getTrId(k) == 0){
            phi[j] = alignments->getProb(k) *
               (const1a + C[0]) *
               (const1b - C[0]); // this comes from division in "false part"
         }else{
            phi[j] = alignments->getProb(k) *
               (const2a - C[0]) *
               (dir->alpha + C[ alignments->getTrId(k) ]); 
               /* 
               /(m * dir->alpha + Nmap - 1 - C[0]) ;
               this term was replaced by *(const1b - C[0]) 
               and moved into "true part" as multiplication 
               */
         }
         probNorm += phi[j];
      }
      r = uniformDistribution(rng_mt);
      // Apply Normalization constant:
      r *= probNorm;
      for(j = 0, sum = 0 ; (sum<r) && (j<readsAlignmentsN); j++){
         sum += phi[j];
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
   // TimeStats {{{
#ifdef DoSTATS
   gettimeofday(&end, NULL);
   tZ += (end.tv_sec-start.tv_sec)*1000*1000+(end.tv_usec-start.tv_usec);
#endif
   // }}}
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

void CollapsedSampler::fixReads(const vector<long> &counts, const vector<long> &unfixed){
   someFixed = true;
   this->C.assign(counts.begin(),counts.end());
   this->unfixed.assign(unfixed.begin(),unfixed.end());
}
