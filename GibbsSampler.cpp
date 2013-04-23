#ifdef DoSTATS
#include<sys/time.h>
#endif

#include "GibbsSampler.h"
#include "common.h"

GibbsSampler::GibbsSampler(){ //{{{
   thetaAct=0;
}//}}}
void GibbsSampler::sampleZ(){//{{{
   // TimeStats {{{
#ifdef DoSTATS
   nZ++;
   struct timeval start, end;
   gettimeofday(&start, NULL);
#endif
   // }}}
   long i,j,k;
   vector<double> phi(m,0); 
   // phi of size M should be enough 
   // because of summing the probabilities for each isoform when reading the data
   double probNorm,r,sum;
   int_least32_t readsAlignmentsN;

   // Reset C to zeros.
   C.assign(C.size(),0);
   // Assign reads.
   for(i=0;i<Nmap;i++){
      probNorm=0;
      readsAlignmentsN = alignments->getReadsI(i+1) - alignments->getReadsI(i);
      for(j=0, k=alignments->getReadsI(i); j < readsAlignmentsN; j++, k++){
         if(alignments->getTrId(k) == 0){
            phi[j] = alignments->getProb(k) * (1 - thetaAct);
         }else{
            phi[j] = alignments->getProb(k) *
               thetaAct * theta[alignments->getTrId(k)];
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
         C[0]++;
      }else{
         // Assign to the chosen transcript.
         C[ alignments->getTrId( alignments->getReadsI(i)+j-1 ) ]++;
      }
   }
   // TimeStats {{{
#ifdef DoSTATS
   gettimeofday(&end, NULL);
   tZ += (end.tv_sec-start.tv_sec)*1000*1000+(end.tv_usec-start.tv_usec);
#endif
   // }}}
}//}}}
void GibbsSampler::sampleThetaAct(){//{{{
#ifdef DoSTATS
   nTa++;
   struct timeval start, end;
   gettimeofday(&start, NULL);
#endif
   double C0=C[0]+Nunmap,X,Y; 
   // counting C_0 from all reads
   // generate thetaAct~Beta(a,b) as thetaAct = X/(X+Y) ; X~Gamma(a,1), Y~Gamma(b,1)
   gammaDistribution.param(gDP(beta->alpha + Nmap+Nunmap - C0, 1));
   X = gammaDistribution(rng_mt);
   gammaDistribution.param(gDP(beta->beta + C0, 1));
   Y = gammaDistribution(rng_mt);
   
   thetaAct = X / (X+Y);
#ifdef DoSTATS
   gettimeofday(&end, NULL);
   tTa += (end.tv_sec-start.tv_sec)*1000*1000+(end.tv_usec-start.tv_usec);
#endif
}//}}}
void GibbsSampler::update(){//{{{
   Sampler::update();

   theta[0]=thetaAct; // save thetaAct as theta_0

   updateSums();
   if((doLog)&&(save))appendFile();
}//}}}
void GibbsSampler::sample(){//{{{
   Sampler::sample();

   sampleTheta();
   sampleThetaAct();
   sampleZ();
}//}}}
