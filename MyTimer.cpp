#include<ctime>

#include "MyTimer.h"

#include "common.h"

void MyTimer::adjust(double &time,char f){//{{{
   if(f=='m')time/=60.0;
   if(f=='h')time/=3600.0;
}//}}}
void MyTimer::write(double time,char f){//{{{
   if(!quiet)messageF("[time: +%.2lf %c]\n",time,f);
}//}}}
MyTimer::MyTimer(){//{{{
   N=1;
   quiet=false;
   times.resize(N);
   times[0]=time(NULL);
}//}}}
void MyTimer::start(long timer){//{{{
   if(timer>=N){
      N=timer+1;
      times.resize(N);
   }
   times[timer]=time(NULL);
}//}}}
double MyTimer::split(long timer, char f){//{{{
   if(timer>=N)return 0;
   double ret;
   ret=time(NULL)-times[timer];
   adjust(ret,f);
   write(ret,f);
   times[timer]=time(NULL);
   return ret;
}//}}}
double MyTimer::getTime(long timer, char f){//{{{
   if(timer>=N)return 0;
   double ret;
   ret=time(NULL)-times[timer];
   adjust(ret,f);
   return ret;
}//}}}
double MyTimer::current(long timer, char f){//{{{
   if(timer>=N)return 0;
   double ret;
   ret=time(NULL)-times[timer];
   adjust(ret,f);
   write(ret,f);
   return ret;
}//}}}
double MyTimer::stop(long timer, char f){//{{{
   if(timer>=N)return 0;
   double ret;
   ret=time(NULL)-times[timer];
   adjust(ret,f);
   write(ret,f);
   times[timer]=time(NULL);
   return ret;
}//}}}

