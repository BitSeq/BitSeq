#ifndef MYTIMER_H
#define MYTIMER_H

#include<vector>

using namespace std;

class MyTimer{
 private:
   vector<time_t> times;
   long N;
   bool quiet;
   // Adjust time to format m or h.
   void adjust(double &time,char f);
   // Write time in format.
   void write(double time,char f);
 public:
   MyTimer();
   void setQuiet(){quiet=true;}
   void setVerbose(){quiet=false;}
   void start(long timer=0);
   double split(long timer=0, char f='s');
   double getTime(long timer=0, char f='s');
   double current(long timer=0, char f='s');
   double stop(long timer=0, char f='s');
};

#endif
