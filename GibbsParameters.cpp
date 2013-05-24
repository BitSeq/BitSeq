#include<fstream>

using namespace std;

#include "GibbsParameters.h"
#include "common.h"

#define DEBUGGP(x) 
#define Sof(x) (long)x.size()


/*void gibbsParameters::setLogFiles(string tau,string tauMeans){//{{{
   gs_samplesFile=tau;
   gs_meansFile=tauMeans; 
}//}}}*/
void gibbsParameters::getAllParameters(){//{{{
   message("Parameters:\n   burnIn: %ld\
\n   samplesN: %ld\n   samplesSave: %ld\
\n   samplesNmax: %ld\n   chainsN: %ld\
\n   targetScaleReduction: %lf\n   dirAlpha: %lf\
\n   dirBeta: %lf\n   betaAlpha: %lf\n   betaBeta: %lf\n",
gs_burnIn,gs_samplesN,gs_samplesSave,gs_samplesNmax,gs_chainsN,gs_targetScaleReduction,dirP.alpha,dirP.beta,betaP.alpha,betaP.beta);
}//}}}
bool gibbsParameters::setParameters(string paramFileName){//{{{
   this->paramFileName = paramFileName;
   return readParameters();
}//}}}
bool gibbsParameters::setParameters(ArgumentParser &args){//{{{
   if(args.isSet("MCMC_burnIn"))gs_burnIn=args.getL("MCMC_burnIn");
   if(args.isSet("MCMC_samplesN"))gs_samplesN=args.getL("MCMC_samplesN");
   if(args.isSet("MCMC_samplesSave"))gs_samplesSave=args.getL("MCMC_samplesSave");
   if(args.isSet("MCMC_samplesNmax"))gs_samplesNmax=args.getL("MCMC_samplesNmax");
   if(args.isSet("MCMC_chainsN"))gs_chainsN=args.getL("MCMC_chainsN");
   if(args.isSet("MCMC_scaleReduction"))gs_targetScaleReduction=args.getD("MCMC_scaleReduction");
   if(args.isSet("MCMC_dirAlpha"))dirP.alpha=args.getD("MCMC_dirAlpha");
   return true;
}//}}}
bool gibbsParameters::readParameters(){//{{{
   ifstream pFile;
   string param;
   double val;
   char tmp[256];
   pFile.open(paramFileName.c_str());
   while((pFile.is_open())&&(! pFile.eof())){
      if((! (pFile>>param)) || (Sof(param)==0) || (param[0]=='#')){
         pFile.getline(tmp,256);
         continue;
      }
      pFile>>val;
      if(pFile.good()){
         DEBUGGP(message("# DEBUG gPar ||%s==%lf||\n",(param).c_str(),val);)
         if(param=="burnIn")parameter("burnIn",gs_burnIn,val);
         if(param=="samplesN")parameter("samplesN",gs_samplesN,val);
         if(param=="samplesSave")parameter("samplesSave",gs_samplesSave,val);
         if(param=="samplesNmax")parameter("samplesNmax",gs_samplesNmax,val);
         if(param=="chainsN")parameter("chainsN",gs_chainsN,val);
         if(param=="targetScaleReduction")parameter("targetScaleReduction",gs_targetScaleReduction,val);
         if(param=="dirAlpha")parameter("dirAlpha",dirP.alpha,val);
         if(param=="dirBeta")parameter("dirBeta",dirP.beta,val);
         if(param=="betaAlpha")parameter("betaAlpha",betaP.alpha,val);
         if(param=="betaBeta")parameter("betaBeta",betaP.beta,val);
         //if(param=="output")parameter("output",gs_output,val);
      }
      pFile.getline(tmp,256);
   }
   //if(gs_samplesN>gs_samplesNmax)gs_samplesNmax=gs_samplesN;
   pFile.close();
   return true;
}//}}}
void gibbsParameters::parameter(string name, double &variable, double value){//{{{
   bool output=false;
   if(verbose && (variable != value))output = true;
   variable = value;
   if(output){
      message("### %s: %lf\n",(name).c_str(),variable);
   }
}//}}}
void gibbsParameters::parameter(string name, long &variable, double value){//{{{
   bool output=false;
   if(verbose && (variable != (long) value))output = true;
   variable = (long) value;
   if(output){
      message("### %s: %ld\n",(name).c_str(),variable);
   }
}//}}}
void gibbsParameters::parameter(string name, bool &variable, double value){//{{{
   bool output=false;
   if(verbose && (variable !=(bool)((long) value)))output = true;
   variable = (bool)((long)value);
   if(output){
      message("### %s: %d\n",(name).c_str(),variable);
   }
}//}}}

gibbsParameters::gibbsParameters(bool verbose){//{{{
   this->verbose = verbose;
   gs_burnIn=1000;
   gs_samplesN=1000;
   gs_samplesNmax=50000;
   gs_samplesSave=500;
   gs_chainsN=4;
   gs_targetScaleReduction=1.2;
   dirP.alpha=1;
   dirP.beta=1;
   betaP.alpha=10;
   betaP.beta=2;
   gs_samplesFile="gibbs_log.rpkmS";
   gs_meansFile="gibbs_log.thetaMeans";
   //gs_output=RPKM;
}//}}}
