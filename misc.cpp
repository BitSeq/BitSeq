#include "misc.h"
#include <algorithm>
#include <cmath>
#include <ctime>

#include "common.h"
#include "FileHeader.h"

namespace ns_math {
double logAddExp(double a, double b){ //{{{
   if(a>b){
      return a+log1p(exp(b-a));
   }else {
      return b+log1p(exp(a-b));
   }
} //}}}
double logSumExp(const vector<double> &vals, long st, long en){ //{{{
   if(st<0)st = 0;
   if((en == -1) || (en > (long)vals.size())) en = vals.size();
   if(st >= en)return 0;
   double sumE = 0, m = *max_element(vals.begin() + st,vals.begin() + en);
   for(long i = st; i < en; i++)
      sumE += exp(vals[i] - m);
   return  m + log(sumE);
} //}}}
} // namespace ns_math

namespace ns_expression {

string getOutputType(const ArgumentParser &args, const string &defaultType){ //{{{
   string type = ns_misc::toLower(args.getS("outputType"));
   if((type!="theta") && (type!="rpkm") && (type!="counts") && (type!="tau")){
      type = defaultType;
      warning("Using output type %s.",type.c_str());
   }
   return type;
} //}}}
} // namespace ns_expression

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

void computeCI(double cf, vector<double> *difs, double *ciLow, double *ciHigh){//{{{
   cf = (100 - cf) / 2.0;
   double N = difs->size();
   sort(difs->begin(),difs->end());
   *ciLow = (*difs)[(long)(N/100.*cf)];
   *ciHigh = (*difs)[(long)(N-N/100.*cf)];
}//}}}

string toLower(string str){//{{{
   for(size_t i=0;i<str.size();i++)
      if((str[i]>='A')&&(str[i]<='Z'))str[i]=str[i]-'A'+'a';
   return str;
}//}}}
} // namespace ns_misc

namespace ns_genes {
bool getLog(const ArgumentParser &args){// {{{
   if(args.flag("log")){
      if(args.verb())message("Using logged values.\n");
      return true;
   }
   if(args.verb())message("NOT using logged values.\n");
   return false;
}// }}}

bool prepareInput(const ArgumentParser &args, TranscriptInfo *trInfo, PosteriorSamples *samples, long *M, long *N, long *G){// {{{
   if(! trInfo->readInfo(args.getS("trInfoFileName"))) return false;
   *G = trInfo->getG();
   if((! samples->initSet(M,N,args.args()[0]))||(*M<=0)||(*N<=0)){
      error("Main: Failed loading MCMC samples.\n");
      return false;
   }
   if(*M!=trInfo->getM()){
      error("Main: Number of transcripts in the info file and samples file are different: %ld vs %ld\n",trInfo->getM(),*M);
      return false;
   }
   if(args.verb())messageF("Transcripts: %ld\n",*M);
   return true;
}// }}}

bool updateGenes(const ArgumentParser &args, TranscriptInfo *trInfo, long *G){//{{{
   if(!(args.isSet("trMapFile") || args.isSet("geneListFile")))return true;
   if(args.isSet("trMapFile") && args.isSet("geneListFile")){
      error("Main: Please provide only one of trMapFile and geneListFile, both serve the same function.\n");
      return false;
   }
   bool isMap;
   ifstream mapFile;
   if(args.isSet("trMapFile")){
      isMap = true;
      mapFile.open(args.getS("trMapFile").c_str());
   }else {
      isMap = false;
      mapFile.open(args.getS("geneListFile").c_str());
   }
   if(!mapFile.is_open()){
      if(isMap){
         error("Main: Failed reading file with transcript to gene mapping.\n");
      }else{
         error("Main: Failed reading file with gene names.\n");
      }
      return false;
   }
   map<string,string> trMap;
   vector<string> geneList;
   string trName,geName;
   while(mapFile.good()){
      while(mapFile.good() && (mapFile.peek()=='#'))
         mapFile.ignore(100000000,'\n');
      if(!mapFile.good()) break;
      mapFile>>geName;
      if(isMap){
         mapFile>>trName;
      }
      if(!mapFile.fail()){
         if(isMap){
            trMap[trName]=geName;
         }else{
            geneList.push_back(geName);
         }
      }
      mapFile.ignore(100000000,'\n');
   }
   mapFile.close();
   bool succ;
   if(isMap)succ = trInfo->updateGeneNames(trMap);
   else succ = trInfo->updateGeneNames(geneList);
   if(!succ){
      error("Main: Filed setting gene information.\n");
      return false;
   }
   *G = trInfo->getG();
   return true;
}//}}}

bool checkGeneCount(long G, long M){//{{{
   if((G != 1) && (G != M)) return true;
   if(G==1){
      error("Main: All transcripts share just one gene.\n");
   }else{
      error("Main: There are no transcripts sharing one gene.\n");
   }
   message("Please provide valid transcript to gene mapping (trMapFile or geneListFile).\n"
           "   (trMap file should contain rows in format: <geneName> <transcriptName>.)\n"
           "   (geneList file should contain rows with gene names, one per transcript.)\n");
   return false;
}//}}}
} // namespace ns_genes

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
