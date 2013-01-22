
#include "Sampler.h"

class GibbsSampler : public Sampler{
   private:
   double thetaAct;

   void sampleThetaAct();
   void sampleZ();
      
   // USING uniformDistribution from class Sampler
   // USING gammaDistribution from class Sampler
   
   public:

   GibbsSampler();
   virtual ~GibbsSampler();
   
//   virtual void init(long n, long m, long samplesTotal, long samplesOut, long Nmap, long Nunmap, const vector<long> &alignI, const vector<TagAlignment> &alignments, const distributionParameters &betaPar, const distributionParameters &dirPar, long seed);
   virtual void update();
   virtual void sample();
};
