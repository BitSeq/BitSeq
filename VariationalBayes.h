#ifndef VARIATIONALBAYES_H
#define VARIATIONALBAYES_H

#include "SimpleSparse.h"
#include "MyTimer.h"

#define LOG_CONV

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

   public:
      VariationalBayes(SimpleSparse *_beta,double *_alpha=NULL,long seed = 0,long procN = 1);
      ~VariationalBayes();
      //double *pack(){return phi_sm->val;} 
      void unpack(double vals[], double adds[] = NULL); // set phi_m, phi=softmax(phi_m), phi_hat=sumOverCols(phi)
      void negGradient(double res[]);
      double getBound();
      void optimize(bool verbose=false, OPT_TYPE method=OPTT_STEEPEST,long maxIter=50000,double ftol=1e-6, double gtol=1e-6);
      double *getAlphas();
      void setLog(string logFileName,MyTimer *timer);
};

#endif
