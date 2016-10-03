#ifndef VARIATIONALBAYES_H
#define VARIATIONALBAYES_H

#include "boost/random/mersenne_twister.hpp"

#include "MyTimer.h"
#include "SimpleSparse.h"

//#define LOG_CONV
//#define LONG_LOG
//#define SHOW_FIXED

enum OPT_TYPE { OPTT_STEEPEST, OPTT_PR, OPTT_FR, OPTT_HS};

class VariationalBayes {
   private:
      long N,M,T; // N - read num K- read number?
      double * alpha; // prior over expression
      double * phiHat;
      double * digA_pH;
      double boundConstant;
      SimpleSparse *beta,*phi_sm,*phi;
      // logBeta replaced by logging beta itself
      string logFileName;
      MyTimer *logTimer;
      // mersen twister random number generator
      boost::random::mt11213b rng_mt;
      bool quiet;

   public:
      VariationalBayes(SimpleSparse *_beta,double *_alpha=NULL,long seed = 0,long procN = 1);
      ~VariationalBayes();
      //double *pack(){return phi_sm->val;} 
      void unpack(double vals[], double adds[] = NULL); // set phi_m, phi=softmax(phi_m), phi_hat=sumOverCols(phi)
      void negGradient(double res[]);
      double getBound();
      void optimize(bool verbose=false, OPT_TYPE method=OPTT_STEEPEST,long maxIter=10000,double ftol=1e-5, double gtol=1e-5);
      double *getAlphas();
      SimpleSparse *getPhi();
      void setLog(string logFileName,MyTimer *timer);
      // Generates samples from the distribution. The 0 (noise) transcript is left out.
      void generateSamples(long samplesN, const string &outTypeS, const vector<double> *isoformLengths, ofstream *outF);
      void beQuiet(){ quiet = true; }
};

#endif
