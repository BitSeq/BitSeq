#include<stdint.h>

#include "Sampler.h"

class CollapsedSampler : public Sampler{
   private:
   vector<int_least32_t> Z;
   vector<double> par1, par2, par3;

   void sampleZ();

   // USING unifromDistribution from class Sampler.
   
   public:

   CollapsedSampler();
   virtual ~CollapsedSampler();
   
   virtual void init(long m, long samplesTotal, long samplesOut, long Nunmap,const TagAlignments *alignments, const distributionParameters &betaPar, const distributionParameters &dirPar, long &seed);
   virtual void update();
   virtual void sample();
   
};
