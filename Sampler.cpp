#ifdef DoSTATS
#include<sys/time.h>
#endif

#include "Sampler.h"
#include "common.h"

#define Sof(x) (long)x.size()

Sampler::Sampler(){ //{{{
   m=samplesN=samplesLogged=samplesTotal=samplesOut=Nmap=Nunmap=0;
   isoformLengths = NULL;
#ifdef DoSTATS
   tT=tTa=tZ=0;
   nT=nTa=nZ=0;
#endif
}//}}}
Sampler::~Sampler(){ //{{{
#ifdef DoSTATS
   message("---------------------------\nSTATISTICS:\n");
   message("Theta: %lg   %lgm   av:%lgs\n",nT,tT/60000.0,tT/1000.0/nT);
   message("Z: %lg   %lgm   av:%lgs\n",nZ,tZ/60000.0,tZ/1000.0/nZ);
   if(nTa>0)message("Theta Act: %ld   %lgm   av:%lgs\n",nTa,tTa/60000.0,tTa/1000.0/nTa);
   message("Total time: %lgm\n",(tT+tZ)/60000.0);
#endif
}//}}}
void Sampler::init(long m, long samplesTotal, long samplesOut, long Nunmap,const TagAlignments *alignments, const distributionParameters &betaPar, const distributionParameters &dirPar, long &seed){//{{{
//   this->n=n;
   this->m=m;
   this->samplesOut=samplesOut;
   this->Nmap=alignments->getNreads();
   this->Nunmap=Nunmap;
   this->alignments=alignments;
   beta=&betaPar;
   dir=&dirPar;
   //dir=new distributionParameters;
   //dir->alpha=1.0/m;
   //dir->beta=dirPar.beta;
   rng_mt.seed(seed);
   seed = (long) (1717171717.17*uniformDistribution(rng_mt));

   resetSampler(samplesTotal);

   theta.assign(m,0);
   C.assign(m,0);
}//}}}
void Sampler::resetSampler(long samplesTotal){//{{{
   this->samplesTotal=samplesTotal;
   samplesN = 0;
   samplesLogged = 0;
   logRate=(double)samplesOut/samplesTotal;
   sumC0 = 0;
   sumNorm.first = sumNorm.second = 0;
   thetaSum.assign(m,pairD(0,0));
   thetaSqSum.assign(m,pairD(0,0));
}//}}}
long Sampler::getAverageC0(){//{{{
   return (long) (sumC0 / sumNorm.first);
}//}}}
void Sampler::getAverage(vector<pairD> &av){//{{{
   long i;
   if(Sof(av)<m)
      av.assign(m,pairD(0,0));
   for(i=0;i<m;i++){
      if(sumNorm.first != 0)
         av[i].first=thetaSum[i].first/sumNorm.first;
      if(sumNorm.second != 0)
         av[i].second=thetaSum[i].second/sumNorm.second;
   }
}//}}}
pairD Sampler::getAverage(long i){//{{{
   double av1,av2;
   av1=(sumNorm.first==0)?0:thetaSum[i].first/sumNorm.first;
   av2=(sumNorm.second==0)?0:thetaSum[i].second/sumNorm.second;
   return pairD(av1,av2);
}//}}}
pairD Sampler::getWithinVariance(long i){//{{{
   double va1,va2;
   if(sumNorm.first==0)
      va1=0;
   else
      va1=thetaSqSum[i].first/(sumNorm.first-1.0) - 
           (thetaSum[i].first/(sumNorm.first-1.0))*
           (thetaSum[i].first/sumNorm.first);
   if(sumNorm.second==0)
      va2=0;
   else 
      va2=thetaSqSum[i].second/(sumNorm.second-1.0) - 
           (thetaSum[i].second/(sumNorm.second-1.0))*
           (thetaSum[i].second/sumNorm.second);
   if(va1<0)message("minus %lg %lg %lg\n",thetaSqSum[i].first,thetaSum[i].first,sumNorm.first);
   return pairD(va1,va2);
}//}}}
void Sampler::getThetaSums(long i, double *thSqSum, double *thSum, double *sumN){//{{{
   if(i >= m){
      (*thSqSum) = (*thSum) = (*sumN) = 0;
      return;
   }
   *thSqSum = thetaSqSum[i].first;
   *thSum = thetaSum[i].first;
   *sumN = sumNorm.first;
}//}}}
void Sampler::getTau(vector<double> &tau, double norm){//{{{
   long i;
   double tauSum=0;

   if ((Sof(theta) > Sof(tau)) || (Sof((*isoformLengths)) != Sof(tau)))
     error("Sampler failed");

   tau.assign(Sof(tau),0);

   tau[0]=theta[0]; // set thetaAct
   // divide by length:
   for(i=1;i<Sof(theta);i++){
      tau[ i ] = theta[i] / (*isoformLengths)[ i ] * norm;
      tauSum += tau[i];
   }
   // DO normalize:
   for(i=1;i<Sof(tau);i++)
      if(tau[i]>0) tau[i] /= tauSum;
}//}}}
void Sampler::appendFile(){//{{{
   long i;
   double norm=saveNorm;
   if((!save) || (outFile == NULL))return;
   thetaActLog.push_back(theta[0]);
   outFile->precision(9);
   (*outFile)<<scientific;
   if(saveType == "counts"){
      if(norm == 0)norm = Nmap;
      for(i=1;i<m;i++)
         (*outFile)<<theta[i]*norm<<" ";
   }else if(saveType == "rpkm"){
      if(norm == 0)norm = 1000000000.0;
      for(i=1;i<m;i++)
         if((*isoformLengths)[i]>0)
            (*outFile)<<theta[i]*norm/(*isoformLengths)[i]<<" ";
         else
            (*outFile)<<theta[i]*norm<<" ";
   }else if(saveType == "theta"){
      if(norm == 0)norm=1.0;
      for(i=1;i<m;i++)
         (*outFile)<<theta[i]*norm<<" ";
   }else if(saveType == "tau"){
      if(norm == 0)norm=1.0;
      vector<double> tau(m);
      getTau(tau,norm);
      for(i=1;i<m;i++)
         (*outFile)<<tau[i]<<" ";
   }
   (*outFile)<<endl;
}//}}}
void Sampler::updateSums(){//{{{
   long i;
   double s;
   for(i=0;i<m;i++){
      thetaSum[i].first+=theta[i];
      thetaSqSum[i].first+=theta[i]*theta[i];
   }
   sumC0+=C[0];
   sumNorm.first++;
   //if(doLog){
   for(i=0;i<m;i++){
      s = log(theta[i]) - log(1-theta[i]);//LOGIT
      thetaSum[i].second += s;
      thetaSqSum[i].second += s * s;
   }
   sumNorm.second++;
   //}
}//}}}
void Sampler::saveSamples(ofstream *outFile, const vector<double> *isoformLengths, const string &saveType, double norm){//{{{
   this->outFile = outFile;
   this->isoformLengths = isoformLengths;
   this->saveType = saveType;
   saveNorm = norm;
   save = true;
   thetaActLog.clear();
}//}}}
void Sampler::noSave(){//{{{
   save = false;
   outFile = NULL;
   if(isoformLengths){
      delete isoformLengths;
      isoformLengths = NULL;
   }
}//}}}

void Sampler::sampleTheta(){//{{{
#ifdef DoSTATS
   nT++;
   struct timeval start, end;
   gettimeofday(&start, NULL);
#endif
   vector<double> gamma(m,0);
   double gammaSum=0;
   long i;
   for(i=1;i<m;i++){
      gammaDistribution.param(gDP(dir->alpha + C[i], dir->beta));
      gamma[i]=gammaDistribution(rng_mt);
      gammaSum+=gamma[i];
   }
   if (gammaSum<=0) // at least something should be more than zero
     error("Sampler failed");

   for(i=1;i<m;i++){
      theta[i]=gamma[i]/gammaSum;
   }
#ifdef DoSTATS
   gettimeofday(&end, NULL);
   tT += (end.tv_sec-start.tv_sec)*1000*1000+(end.tv_usec-start.tv_usec);
#endif
}//}}}
void Sampler::sample(){//{{{
   samplesN++;
}//}}}
void Sampler::update(){//{{{
   doLog = false;
   if(samplesOut-samplesLogged>0){
      if(samplesTotal-samplesN<=samplesOut-samplesLogged)doLog=true;
      else if((long)(logRate * samplesN) > samplesLogged)doLog=true;
   }
   if(doLog) samplesLogged ++;
}//}}}
