#include<stdint.h>

#include "Sampler.h"

class CollapsedSampler : public Sampler{
   private:
   vector<int_least32_t> Z;

   void sampleZ();

   public:

   virtual void update();
   virtual void sample();
   
};
