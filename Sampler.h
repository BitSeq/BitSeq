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

// compute statistics
//#define DoSTATS
//#define DoDebug

typedef pair<double,double> pairD;

class Sampler{
   protected:
   long m, samplesN, samplesLogged, samplesTotal, samplesOut, Nmap, Nunmap;
   const distributionParameters *beta,*dir;
   const TagAlignments *alignments;
   const vector<double> *isoformLengths;
   boost::random::mt11213b rng_mt;
   boost::random::gamma_distribution<double> gammaDistribution;
   typedef boost::random::gamma_distribution<double>::param_type gDP;
   // Need by children:
   boost::random::uniform_01<double> uniformDistribution;
   
   bool doLog,save;
   string saveType;
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

   // Sample theta.
   void sampleTheta();
   // Compute tau.
   void getTau(vector <double> &tau, double norm);
   // Append current expression samples into file opened for saving samples.
   void appendFile();
   // Update sums of theta and theta^2.
   void updateSums();

   public:

   Sampler();
   virtual ~Sampler();
   // Initialize sampler, set seed and use it to generate new seed.
   void init(long m, long samplesTotal, long samplesOut, long Nunmap,
             const TagAlignments *alignments,
             const distributionParameters &betaPar, 
             const distributionParameters &dirPar, 
             long &seed);
   // Reset sampler's stats before new iteration
   void resetSampler(long samplesTotal);
   // Return mean C[0].
   long getAverageC0();
   // Get vector of mean theta expression. Has "two columns" first is calculated
   // from all samples, the second only from the "thinned" samples.
   void getAverage(vector<pairD> &av);
   // Get mean for transcript i.
   pairD getAverage(long i);
   // Get within variance for transcript i.
   pairD getWithinVariance(long i);
   // Get sum of theta^2, sum of theta, and their norm for transcript i.
   void getThetaSums(long i, double *thSqSum, double *thSum, double *sumN);
   // Return norms for theta sums.
   pairD getSumNorms() const { return sumNorm; }
   // Set sampler into state where samples are saved into the outFile.
   void saveSamples(ofstream *outFile, const vector<double> *isoformLengths,
                    const string &saveType, double norm = 0);
   // Stop saving samples into the file.
   void noSave();
   // Get theta act logged values.
   const vector<double>& getThetaActLog(){return thetaActLog;}

   // Produce new McMc samples.
   virtual void sample();
   // If necessary ("thinned sample") sample theta; update sums.
   virtual void update();
};


#endif
