#ifndef GIBBSPARAMETERS_H
#define GIBBSPARAMETERS_H

#include<string>

using namespace std;

#include "ArgumentParser.h"

struct distributionParameters{//{{{
   double alpha,beta;
};//}}}

class gibbsParameters{
   private:
      long gs_burnIn, gs_samplesN, gs_chainsN, gs_samplesNmax, gs_samplesSave;
      double gs_targetScaleReduction;
      bool verbose;
      distributionParameters dirP, betaP;
      string gs_samplesFile,gs_meansFile,paramFileName;
      void parameter(string name, bool &variable, double value);
      void parameter(string name, long &variable, double value);
      void parameter(string name, double &variable, double value);
   public:
      gibbsParameters(bool verbose = true);
      bool setParameters(string paramFileName);
      bool setParameters(ArgumentParser &args);
      bool readParameters();
      void getAllParameters();
      long burnIn() const {return gs_burnIn;}
      long samplesN() const {return gs_samplesN;}
      long samplesSave() const {return gs_samplesSave;}
      long samplesNmax() const {return gs_samplesNmax;}
      long chainsN() const {return gs_chainsN;}
      const distributionParameters& dir() const {return dirP;}
      const distributionParameters& beta()const {return betaP;}
      double targetScaleReduction() const {return gs_targetScaleReduction;}
//      string samplesFile() const {return gs_samplesFile;}
//      string meansFile() const {return gs_meansFile;}
//      void setLogFiles(string tau,string tauMeans);
//      outputType output() const {return (outputType)gs_output;}
};


#endif
