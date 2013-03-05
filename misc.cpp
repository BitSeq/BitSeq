#include "misc.h"
#include <algorithm>
#include <ctime>

#include "common.h"
#include "FileHeader.h"


namespace ns_misc {

long getSeed(const ArgumentParser &args){//{{{
   long seed;
   if(args.isSet("seed"))seed=args.getL("seed");
   else seed = time(NULL);
   if(args.verbose)message("seed: %ld\n",seed);
   return seed;
}//}}}
bool openOutput(const ArgumentParser &args, ofstream *outF){//{{{
   outF->open(args.getS("outFileName").c_str());
   if(!outF->is_open()){
      error("Main: Output file open failed.\n");
      return false;
   }
   return true;
}//}}}
bool openOutput(const string &name, ofstream *outF) {//{{{
   outF->open(name.c_str());
   if(!outF->is_open()){
      error("Main: File '%s' open failed.\n",name.c_str());
      return false;
   }
   return true;
}//}}}

bool readConditions(const ArgumentParser &args, long *C, long *M, long *N, Conditions *cond){//{{{
   if(! cond->init("NONE", args.args(), C, M, N)){
      error("Main: Failed loading MCMC samples.\n");
      return false;
   }
   if(args.isSet("normalization")){
      if(! cond->setNorm(args.getTokenizedS2D("normalization"))){
         error("Main: Applying normalization constants failed.\n");
         return false;
      }
   }
   if(!cond->logged() && args.verb()){
      message("Samples are not logged. (will log for you)\n");
      message("Using %lg as minimum instead of log(0).\n",LOG_ZERO);
   }
   if(args.verb())message("Files with samples loaded.\n");
   return true;
}//}}}

} // namespace ns_misc

namespace ns_params {
bool readParams(const string &name, vector<paramT> *params, ofstream *outF){//{{{
   long parN;
   ifstream parFile(name.c_str());
   FileHeader fh(&parFile);
   if(!fh.paramsHeader(&parN, outF)){
      error("Main: Problem loading parameters file %s\n",name.c_str());
      return false;
   }
   // Vector of parameters: (mean expression, (alpha, beta) )
   paramT param;
   while(parFile.good()){
      while((parFile.good())&&(parFile.peek()=='#')){
         parFile.ignore(10000000,'\n');
      }
      parFile>>param.alpha>>param.beta>>param.expr;
      if(parFile.good())
         params->push_back(param);
      parFile.ignore(10000000,'\n');
   }
   if((parN>0)&&(parN != (long)params->size())){
      warning("Main: declared number of parameters does not match number of lines read (%ld %ld).\n", parN, (long)params->size());
   }
   fh.close();
   sort(params->begin(),params->end());
   return true;
}//}}}

} // namespace ns_params
