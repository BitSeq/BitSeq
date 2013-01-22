#ifndef TAGALIGNMENT_H
#define TAGALIGNMENT_H



class TagAlignment{
   protected:
      long trId;
//      bool strand; // true = forward; false = reverse
      long double prob;
   public:
      TagAlignment(long t=0,long double p = 0){
         trId=t;
//         strand=s;
         prob=p;
      }
      long getTrId()const {return trId;}
      double getProb()const {return prob;}
      void setProb(double p){prob=p;}
}; 

class TagAlignment2: public TagAlignment {
   private:
      long double lowProb;
   public:
      //TagAlignment(long t=0,bool s=true,long double p = 0,long double lp = 0){
      TagAlignment2(long t=0,long double p = 0,long double lp = 0){
         trId=t;
//         strand=s;
         prob=p;
         lowProb = lp;
      }
      double getLowProb()const {return lowProb;}
}; 


#endif
