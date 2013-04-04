#ifndef SAMPLER_H
#define SAMPLER_H

#include<vector>
#include<fstream>
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/gamma_distribution.hpp"
#include "boost/random/uniform_01.hpp"

using namespace std;

#include "GibbsParameters.h"
#include "TagAlignments.h"

#define Sof(x) (long)x.size()
// compute statistics
//#define DoSTATS
//#define DoDebug


typedef pair<double,double> pairD;

class Sampler{
   protected:
   long m, samplesN, samplesLogged, samplesTotal, samplesOut, Nmap, Nunmap;
   const distributionParameters *beta,*dir;
//   distributionParameters *dir;
   const TagAlignments *alignments;
   const vector<double> *isoformLengths;
   boost::random::mt11213b rng_mt;
   boost::random::gamma_distribution<double> gammaDistribution;
   typedef boost::random::gamma_distribution<double>::param_type gDP;
   // for kids only:
   boost::random::uniform_01<double> uniformDistribution;
   
   bool doLog,save;
   outputType saveType;
   ofstream *outFile;
   double saveNorm,logRate;
#ifdef DoSTATS   
   long long nT,nZ,nTa;
   double tT,tZ,tTa;
#endif

   vector<long> C;
   double sumC0;
   vector<double> theta;
   vector<double> thetaActLog;
   vector<pairD> thetaSum;
   vector<pairD> thetaSqSum;
   pairD sumNorm;

   void sampleTheta();
   void getTau(vector <double> &tau, double norm);
   void appendFile();
   void updateSums();

   public:

   Sampler();
   virtual ~Sampler();
   virtual void init(long m, long samplesTotal, long samplesOut, long Nunmap,const TagAlignments *alignments, const distributionParameters &betaPar, const distributionParameters &dirPar, long &seed);
   
   void resetSampler(long samplesTotal);
   long getAverageC0();
   void getAverage(vector<pairD> &av);
   pairD getAverage(long i);
   void getWithinVariance(vector<pairD> &va);
   pairD getWithinVariance(long i);
   void saveSamples(ofstream *outFile, const vector<double> *isoformLengths, outputType saveType, double norm = 0);
   void noSave();
   const vector<double>& getThetaActLog(){return thetaActLog;}

   
   virtual void sample(); // produces new McMc samples
   virtual void update(); // samples theta if necessary and updates sums
};


#endif
