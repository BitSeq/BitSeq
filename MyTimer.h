#ifndef MYTIMER_H
#define MYTIMER_H

#include<vector>
#include<cstdlib>
#include<ctime>

using namespace std;

#include "common.h"

class MyTimer{
   private:
   vector<time_t> times;
   long N;
   bool quiet;
      void adjust(double &time,char f){//{{{
         if(f=='m')time/=60.0;
         if(f=='h')time/=3600.0;
      }//}}}
      void write(double time,char f){//{{{
         if(!quiet)message("[time: +%lf %c]\n",time,f);
      }//}}}
   public:
   MyTimer(){//{{{
      N=1;
      quiet=false;
      times.resize(N);
      times[0]=time(NULL);
   }//}}}
   void setQuiet(){quiet=true;}
   void setVerbose(){quiet=false;}
   void start(long timer=0){//{{{
      if(timer>=N){
         N=timer+1;
         times.resize(N);
      }
      times[timer]=time(NULL);
   }//}}}
   double split(long timer=0, char f='s'){//{{{
      if(timer>=N)return 0;
      double ret;
      ret=time(NULL)-times[timer];
      adjust(ret,f);
      write(ret,f);
      times[timer]=time(NULL);
      return ret;
   }//}}}
   double current(long timer=0, char f='s'){//{{{
      if(timer>=N)return 0;
      double ret;
      ret=time(NULL)-times[timer];
      adjust(ret,f);
      write(ret,f);
      return ret;
   }//}}}
   double stop(long timer=0, char f='s'){//{{{
      if(timer>=N)return 0;
      double ret;
      ret=time(NULL)-times[timer];
      adjust(ret,f);
      write(ret,f);
      times[timer]=time(NULL);
      return ret;
   }//}}}
};

#endif
