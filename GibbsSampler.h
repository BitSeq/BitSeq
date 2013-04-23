
#include "Sampler.h"

class GibbsSampler : public Sampler{
   private:
   double thetaAct;

   void sampleThetaAct();
   void sampleZ();
      
   public:

   GibbsSampler();
   
   virtual void update();
   virtual void sample();
};
