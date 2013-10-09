#include<cmath>
#include<cstring>
#include<fstream>
#include<iomanip>
#include<omp.h>
#include "boost/random/normal_distribution.hpp"
#include "boost/random/gamma_distribution.hpp"

#include "asa103/asa103.hpp"
#include "common.h"

#include "VariationalBayes.h"

#define SWAPD(x,y) {tmpD=x;x=y;y=tmpD;}
#define ZERO_LIMIT 1e-12

typedef boost::random::gamma_distribution<double>::param_type gDP;

void VariationalBayes::setLog(string logFileName,MyTimer *timer){//{{{
   this->logFileName=logFileName;
   this->logTimer=timer;
}//}}}
VariationalBayes::VariationalBayes(SimpleSparse *_beta,double *_alpha,long seed,long procN){//{{{
/*
 As bitseq_vb::__init__(self, alpha, beta) in python
 Python difference:
  - python version excludes beta.data <= 1e-40
*/
   quiet = false;
   logFileName = "tmp.convLog";
   logTimer = NULL;
#ifdef SUPPORT_OPENMP
   omp_set_num_threads(procN);
#endif
   long i;
   beta=_beta;
   N=beta->N;
   M=beta->M;
   T=beta->T;
   
   //logBeta= new SimpleSparse(beta);
   //beta already contains log probabilities.

   if(_alpha){
      alpha = _alpha;
   }else{
      alpha = new double[M];
      for(i=0;i<M;i++)alpha[i]=1.;
   }
   phiHat = new double[M];
   digA_pH = new double[M];
   
   rng_mt.seed(seed);
   boost::random::normal_distribution<long double> normalD;
   //typedef boost::random::normal_distribution<long double>::param_type nDP;
   //normalD.param(nDP(0,1));

   phi_sm = new SimpleSparse(beta);
   for(i=0;i<T;i++)phi_sm->val[i] = normalD(rng_mt);
   phi = new SimpleSparse(beta);
   // PyDif make phi a copy of phi_sm <- not important because of unpack() coming up next
   
   unpack(phi_sm->val); //unpack(pack()); 

   double alphaS=0,gAlphaS=0;
   for(i=0;i<M;i++){
      alphaS+=alpha[i];
      gAlphaS+=lgamma(alpha[i]);
   }
   boundConstant = lgamma(alphaS) - gAlphaS - lgamma(alphaS+N);
}//}}}
VariationalBayes::~VariationalBayes(){//{{{
   delete[] alpha;
   delete[] phiHat;
   delete[] digA_pH;
   delete phi_sm;
   delete phi;
}//}}}
void VariationalBayes::unpack(double vals[],double adds[]){//{{{
   if(adds==NULL){
      if(vals!=phi_sm->val)memcpy(phi_sm->val,vals,T*sizeof(double));
   }else{
      long i;
      #pragma omp parallel for
      for(i=0;i<T;i++)phi_sm->val[i] = vals[i]+adds[i];
   }
   phi_sm->softmaxInplace(phi); //softmax  phi_sm into phi; and set phi_sm to log(phi)
   phi->sumCols(phiHat); // sumCols of phi into phiHat
}//}}}

void VariationalBayes::negGradient(double res[]){//{{{
   long i;
   int err=0,totalError=0;
   #pragma omp parallel for private(err) reduction(+:totalError)
   for(i=0;i<M;i++){
      digA_pH[i]=digama(alpha[i]+phiHat[i], &err);
      totalError += err;
   }
   if(totalError){error("VariationalBayes: Digamma error (%d).\n",totalError); }
   // beta is logged now
   #pragma omp parallel for
   for(i=0;i<T;i++)res[i]= - (beta->val[i] - phi_sm->val[i] - 1.0 + digA_pH[beta->col[i]]);
}//}}}
double VariationalBayes::getBound(){//{{{
   // the lower bound on the model likelihood
   double A=0,B=0,C=0;
   long i;
   #pragma omp parallel for reduction(+:A,B)
   for(i=0;i<T;i++){
      // beta is logged now.
      A += phi->val[i] * beta->val[i];
      // PyDif use nansum instead of ZERO_LIMIT (nansum sums all elements treating NaN as zero
      if(phi->val[i]>ZERO_LIMIT){
         B += phi->val[i] * phi_sm->val[i];
      }
   }
   #pragma omp parallel for reduction(+:C)
   for(i=0;i<M;i++){
      C += lgamma(alpha[i]+phiHat[i]);
   }
   return A+B+C+boundConstant;
}//}}}

void VariationalBayes::optimize(bool verbose,OPT_TYPE method,long maxIter,double ftol, double gtol){//{{{
   bool usedSteepest;
   long iteration=0,i,r;
   double boundOld,bound,squareNorm,squareNormOld=1,valBeta=0,valBetaDiv,natGrad_i,gradGamma_i,phiGradPhiSum_r;
   double *gradPhi,*natGrad,*gradGamma,*searchDir,*tmpD,*phiOld;
   gradPhi=natGrad=gradGamma=searchDir=tmpD=phiOld=NULL;
   // allocate stuff {{{
   //SimpleSparse *phiGradPhi=new SimpleSparse(beta);
   gradPhi = new double[T];
   // phiOld = new double[T]; will use gradPhi memory for this
   phiOld = NULL;
   natGrad = new double[T];
   if(method == OPTT_HS)
      gradGamma = new double[T];
   searchDir = new double[T];
   //searchDirOld = new double[T];
   //phiGradPhi_sum = new double[N];
   // }}}
#ifdef LOG_CONV
   ofstream logF(logFileName.c_str());
   logF.precision(15);
   logF<<"# iter bound squareNorm time(m) [M*means M*vars]"<<endl;
   if(logTimer)logTimer->setQuiet();
   #ifdef LONG_LOG
   vector<double> dirAlpha(M);
   #endif
#endif
   boundOld=getBound();
   while(true){
      negGradient(gradPhi);
      // "yuck"
      //setVal(phiGradPhi,i,phi->val[i]*gradPhi[i]);
      //phiGradPhi->sumRows(phiGradPhi_sum);
      // removed need for phiGradPhi matrix:
      // removed need for phiGradPhi_sum
      /*for(r=0;r<N;r++){
         phiGradPhi_sum[r] = 0;
         for(i=phi->rowStart[r];i<phi->rowStart[r+1];i++) phiGradPhi_sum[r] += phi->val[i] * gradPhi[i];
      }*/

      // set natGrad & gradGamma
      squareNorm=0;
      valBeta = 0;
      valBetaDiv = 0;
      #pragma omp parallel for private(i,phiGradPhiSum_r,natGrad_i,gradGamma_i) reduction(+:squareNorm,valBeta,valBetaDiv)
      for(r=0;r<N;r++){
         phiGradPhiSum_r = 0;
         for(i = phi->rowStart[r]; i < phi->rowStart[r+1]; i++) 
            phiGradPhiSum_r += phi->val[i] * gradPhi[i];
         
         for(i = phi->rowStart[r]; i < phi->rowStart[r+1]; i++){
            natGrad_i = gradPhi[i] - phiGradPhiSum_r;
            gradGamma_i = natGrad_i * phi->val[i];
            squareNorm += natGrad_i * gradGamma_i;
            
            if(method==OPTT_PR){
               valBeta += (natGrad_i - natGrad[i])*gradGamma_i;
            }
            if(method==OPTT_HS){
               valBeta += (natGrad_i-natGrad[i])*gradGamma_i;
               valBetaDiv += (natGrad_i-natGrad[i])*gradGamma[i];
               gradGamma[i] = gradGamma_i;
            }
            natGrad[i] = natGrad_i;
         }
      }
      
      if((method==OPTT_STEEPEST) || (iteration % (N*M)==0)){
         valBeta=0;
      }else if(method==OPTT_PR ){
         // already computed:
         // valBeta=0;
         // for(i=0;i<T;i++)valBeta+= (natGrad[i]-natGradOld[i])*gradGamma[i];
         valBeta /= squareNormOld;
      }else if(method==OPTT_FR ){
         valBeta = squareNorm / squareNormOld;
      }else if(method==OPTT_HS ){
         // already computed:
         //valBeta=div=0;
         //for(i=0;i<T;i++){
         //   valBeta += (natGrad[i]-natGradOld[i])*gradGamma[i];
         //   div += (natGrad[i]-natGradOld[i])*gradGammaOld[i];
         //}
         if(valBetaDiv!=0)valBeta /= valBetaDiv;
         else valBeta = 0;
      }

      if(valBeta>0){
         usedSteepest = false;
         //for(i=0;i<T;i++)searchDir[i]= -natGrad[i] + valBeta*searchDirOld[i];
         // removed need for searchDirOld:
         #pragma omp parallel for
         for(i=0;i<T;i++)
            searchDir[i]= -natGrad[i] + valBeta*searchDir[i];
      }else{
         usedSteepest = true;
         #pragma omp parallel for
         for(i=0;i<T;i++)
            searchDir[i]= -natGrad[i];
      }

      //try conjugate step
      SWAPD(gradPhi,phiOld);
      memcpy(phiOld,phi_sm->val,T*sizeof(double)); // memcpy(phiOld,pack(),T*sizeof(double));
      unpack(phiOld,searchDir);
      bound = getBound();
      iteration++;
      // make sure there is an increase in L, else revert to steepest
      if((bound<boundOld) && (valBeta>0)){
         usedSteepest = true;
         #pragma omp parallel for
         for(i=0;i<T;i++)
            searchDir[i]= -natGrad[i];
         unpack(phiOld,searchDir);
         bound = getBound();
         // this should not be increased: iteration++;
      }
      if(bound<boundOld) { // If bound decreased even after using steepest, step back and quit.
         unpack(phiOld);
      }
      SWAPD(gradPhi,phiOld);
      if(verbose){
         #ifdef SHOW_FIXED
            messageF("iter(%c): %5.ld  bound: %.3lf grad: %.7lf  beta: %.7lf  fixed: %ld\n",(usedSteepest?'s':'o'),iteration,bound,squareNorm,valBeta,phi->countAboveDelta(0.999));
         #else
            messageF("iter(%c): %5.ld  bound: %.3lf grad: %.7lf  beta: %.7lf\n",(usedSteepest?'s':'o'),iteration,bound,squareNorm,valBeta);
         #endif
      }else if(!quiet){
         messageF("\riter(%c): %5.ld  bound: %.3lf grad: %.7lf  beta: %.7lf      ",(usedSteepest?'s':'o'),iteration,bound,squareNorm,valBeta);
      }
#ifdef LOG_CONV
   if((iteration%100==0) ||
      ((iteration<500) && (iteration%50==0)) ||
      ((iteration<150) && (iteration%10==0)) ||
      ((iteration<50) && (iteration%5==0))){
      logF<<iteration<<" "<<bound<<" "<<squareNorm;
      if(logTimer)logF<<" "<<logTimer->current(0,'m');
      #ifdef LONG_LOG
      double alphaSum = 0, alphaVarNorm;
      // True 'alpha' - Dirichlet parameter is alpha+phiHat.
      for(i=1;i<M;i++){
         dirAlpha[i] = alpha[i] + phiHat[i];
         alphaSum += dirAlpha[i];
      }
      for(i=1;i<M;i++)logF<< " " << dirAlpha[i] / alphaSum;
      alphaVarNorm = alphaSum*alphaSum*(alphaSum+1);
      for(i=1;i<M;i++)logF<<" "<<dirAlpha[i]*(alphaSum-dirAlpha[i])/alphaVarNorm;
      #endif
      logF<<endl;
   }
#endif

      // convergence check {{{
      if(bound<boundOld){
         message("\nEnd: bound decrease\n");
         break;
      }
      if(abs(bound-boundOld)<=ftol){
         message("\nEnd: converged (ftol)\n");
         break;
      }
      if(squareNorm<=gtol){
         message("\nEnd: converged (gtol)\n");
         break;
      }
      if(iteration>=maxIter){
         message("\nEnd: maxIter exceeded\n");
         break;
      }
      // }}}
      // store essentials {{{
      squareNormOld=squareNorm;
      boundOld=bound;
      // }}}
      R_INTERUPT;
   }
   if(quiet){
      messageF("iter(%c): %5.ld  bound: %.3lf grad: %.7lf  beta: %.7lf\n",(usedSteepest?'s':'o'),iteration,bound,squareNorm,valBeta);
   }
#ifdef LOG_CONV
   logF<<iteration<<" "<<bound<<" "<<squareNorm;
   if(logTimer)logF<<" "<<logTimer->current(0,'m');
   #ifdef LONG_LOG
   double alphaSum = 0, alphaVarNorm;
   // True 'alpha' - Dirichlet parameter is alpha+phiHat.
   for(i=1;i<M;i++){
      dirAlpha[i] = alpha[i] + phiHat[i];
      alphaSum += dirAlpha[i];
   }
   for(i=1;i<M;i++)logF<< " " << dirAlpha[i] / alphaSum;
   alphaVarNorm = alphaSum*alphaSum*(alphaSum+1);
   for(i=1;i<M;i++)logF<<" "<<dirAlpha[i]*(alphaSum-dirAlpha[i])/alphaVarNorm;
   #endif
   logF<<endl;
   if(logTimer)logTimer->setVerbose();
   logF.close();
#endif
   // free memory {{{
   //delete phiGradPhi;
   delete[] gradPhi;
   delete[] natGrad;
   if(method == OPTT_HS)
      delete[] gradGamma;
   delete[] searchDir;
   //delete[] searchDirOld;
   //delete[] phiGradPhi_sum;
   // }}}
}//}}}

double *VariationalBayes::getAlphas(){//{{{
   double *alphas = new double[M];
   for(long i=0;i<M;i++)alphas[i] = alpha[i] + phiHat[i];
   return alphas;
}//}}}

void VariationalBayes::generateSamples(long samplesN, const string &outTypeS, const vector<double> *isoformLengths, ofstream *outF) {//{{{
   vector<double> gamma(M,0);
   vector<gDP> alphaParam;
   boost::random::gamma_distribution<double> gammaDistribution;
   long n,m;
   double gammaSum, norm, normC = 1.0;
   // Set normalisation.
   if(outTypeS == "counts") normC = N; // N is Nmap.
   if(outTypeS == "rpkm") normC = 1e9;
   // Pre-compute Dirichlet's alpha and save them as parameters for Gamma.
   for(m=0;m<M;m++)alphaParam.push_back(gDP(alpha[m] + phiHat[m], 1.0));
   // Sample.
   outF->precision(9);
   (*outF)<<scientific;
   for(n=0;n<samplesN;n++){
      // Compute M gammas and sum. Ignore 0 - noise transcript.
      gammaSum = 0;
      for(m=1;m<M;m++){
         gammaDistribution.param(alphaParam[m]);
         gamma[m] = gammaDistribution(rng_mt);
         gammaSum += gamma[m];
      }
      // For rpkm normalize by length.
      if(outTypeS == "rpkm"){
         if((long)isoformLengths->size() < M){
            error("VariationalBayes: Too few isoform lengths for RPKM computation.");
            return;
         }
         for(m=1;m<M;m++)
            if((*isoformLengths)[m]>0)
               gamma[m] /= (*isoformLengths)[m];
      }
      norm = normC / gammaSum;
      for(m=1;m < M;m++){
         (*outF)<<gamma[m] * norm<<" ";
      }
      (*outF)<<endl;
      R_INTERUPT;
   }
   // Delete lengths.
   delete isoformLengths;
}//}}}
