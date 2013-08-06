#include<stdint.h>

#include "Sampler.h"

class CollapsedSampler : public Sampler{
   private:
   vector<int_least32_t> Z;
   bool someFixed;
   vector<long> unfixed;

   void sampleZ();

   public:
   CollapsedSampler() : Sampler() { someFixed=false; }
   void fixReads(const vector<long> &counts, const vector<long> &unfixed);

   virtual void update();
   virtual void sample();
   
};
