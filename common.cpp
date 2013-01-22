#include <cstdlib>
#include <string>

#include "common.h"

using namespace std;

void buildTime(char *argv0, string compileDate, string compileTime){
#ifdef BIOC_BUILD
   return ; // dont want to print compile information
#endif
   message("### %s build: %s %s\n",argv0,compileDate.c_str(),compileTime.c_str());
}

bool progressLog(long cur,long outOf, long steps) {
   // output progress status every 10%
   if((outOf>steps)&&(cur%((long)(outOf/steps))==0)&&(cur!=0)){
      message("# %ld done.\n",cur);
      return true;
   }
   return false;
}
