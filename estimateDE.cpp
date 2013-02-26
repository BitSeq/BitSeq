/*
 * Original model applying the DE model to individual sets of samples independently.
 * One set of samples == 1 sample from each replicate of each condition.
 */
#include<algorithm>
#include<cmath>
#include<fstream>
#include "boost/random/gamma_distribution.hpp"
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/normal_distribution.hpp"

using namespace std;

#include "ArgumentParser.h"
#include "common.h"
#include "FileHeader.h"
#include "MyTimer.h"
#include "PosteriorSamples.h"

#define FF first
#define SS second

//#define PERCENT 0.9

#define LAMBDA_0 2.0
const double LOG_ZERO=-1000;

namespace ns_estimateDE {

struct paramT {//{{{
   double expr, alpha, beta;
   bool operator< (const paramT &p2) const{
      return expr<p2.expr;
   }
};//}}}

// Open and write headers into appropriate output files.
// The size of outFiles[] should be C+1.
// Retruns true if everything went OK.
bool initializeOutputFile(long M, long N, long C, const ArgumentParser &args, ofstream *outF, ofstream outFiles[]);
// For a given mean expression expr finds alpha and beta for which were estimated for a closes expression.
void getParams(double expr,const vector<paramT> &params, double *alpha, double *beta);
bool readParams(const ArgumentParser &args, long *parN, vector<paramT> *params);

}

extern "C" int estimateDE(int *argc,char* argv[]){
string programDescription =
"Estimate differential expression from the dataset.\n\
   [sampleFiles] should contain transposed MCMC samples from replicates.\n\
   To distinguish conditions use C between them e.g.:\n\
      samplesC1-R1.rpkm samplesC1-R2.rpkm C samplesC2-R1.rpkm samplesC2-R2.rpkm";
   // Intro: {{{
   buildTime(argv[0],__DATE__,__TIME__);
   ArgumentParser args(programDescription,"[sampleFiles]",1);
   args.addOptionS("o","outPrefix","outFilePrefix",1,"Prefix for the output files.");
   args.addOptionS("p","parameters","parFileName",1,"File containing estimated hyperparameters.");
   args.addOptionB("s","samples","samples",0,"Produce samples of condition mean expression apart from PPLR and confidence.");
   args.addOptionD("l","lambda0","lambda0",0,"Parameter lambda_0.",LAMBDA_0);
   args.addOptionD("c","confidencePerc","cf",0,"Percentage for confidence intervals.", 5);
   args.addOptionS("","norm","normalization",0,"Normalization constants for each input file provided as comma separated list of doubles (e.g. 1.0017,1.0,0.9999 ).");
   if(!args.parse(*argc,argv))return 0;
   //}}}
   long C,M,N,m,n,c,r,parN;
   // Open file with hyper parameters and read those in.
   vector<ns_estimateDE::paramT> params(parN);
   if(!readParams(args, &parN, &params)) return 1;
   // samples files: {{{
   // Initialize sample files handled by object cond.
   Conditions cond;
   if(! cond.init(C,M,N,"NONE",args.args())){
      error("Main: Failed loading conditions.\n");
      return 0;
   }
   if(args.isSet("normalization")){
      if(! cond.setNorm(args.getTokenizedS2D("normalization"))){
         error("Main: Appying normalization constants failed.\n");
         return 1;
      }
   }
   bool logged = cond.logged();
   if( (!logged) && args.verbose){
      message("Samples are not logged. (will log for you)\n");
      message("Using %lg as minimum instead of log(0).\n",LOG_ZERO);
   }
   // CR = cond.getRN();
   if(args.verbose)message("Sample files loaded.\n");
   // }}}
   // Initialize output files.
   ofstream outF;
   // Use standard array as we don't want to bother with vector of pointers.
   ofstream outFiles[C+1];
   if(!ns_estimateDE::initializeOutputFile(M, N, C, args, &outF, outFiles)) return 1;

   // variables {{{
   vector<vector<vector<double> > > tr(C);
   vector<vector<double> > samples(C,vector<double>(N));
   vector<double> vars(N);
   vector<double> normMu(C);
   vector<double> mu_0(C);
//   vector<vector<double> > mus(C,vector<double>(N,0));
//   vector<double> vars(N);
   double prec,sum,sumSq,alpha,beta,betaPar,mu,al0,be0,mu_00,divC,divT;
   double lambda0 = args.getD("lambda0");
   long RC;
   MyTimer timer;
   boost::random::mt11213b rng_mt(time(NULL));
   boost::random::gamma_distribution<long double> gammaDistribution;
   typedef boost::random::gamma_distribution<long double>::param_type gDP;
   boost::random::normal_distribution<long double> normalDistribution;
   typedef boost::random::normal_distribution<long double>::param_type nDP;
   // }}}

   double logFC, pplr, cfLow, cfHigh;
   vector<long double> difs(N);
   /*
    * N - number of samples in one replicate (the smallest number for replicates with different N_r)
    * M - number of transcripts
    * C - number of conditions
    * not used: CR - total number of replicates in all conditions
    *
    */
   if(args.verbose){ //{{{
      timer.split();
      message("Sampling condition mean expression.\n");
   }//}}}
   for(m=0;m<M;m++){
      if(progressLog(m,M))timer.split();
      // Read and prepare {{{
      mu_00 = divT = 0;
      for(c=0;c<C;c++){
         mu_0[c]=0;
         divC=0;
         RC = cond.getRC(c);
         if((long)tr[c].size() < RC){
            tr[c].resize( RC );
         }
         for(r=0;r< RC;r++){
            if(cond.getTranscript(c, r , m, tr[c][r]), N){
               for(n=0;n<N;n++){
                  if(!logged)tr[c][r][n] = (tr[c][r][n] == 0)? LOG_ZERO : log (tr[c][r][n] ); // NO LOGGING
                  mu_0[c]+=tr[c][r][n];
               }
               divC+=1;
            }else{
               warning("Main: Condition %ld replicate %ld does not seem to have transcript %ld.\n",c,r,m);
            }
         }
#ifdef BIOC_BUILD
	 R_CheckUserInterrupt();
#endif
         mu_0[c] /= (divC * N); // take mean over all replicates
         mu_00+=mu_0[c];
         if(mu_0[c]!=0)divT++;
      }
      mu_00/=divT; 
      //}}}
      // Sample condition mean expressions {{{
      for(n=0;n<N;n++){
         for(c=0;c<C;c++){
            RC = cond.getRC(c);
            ns_estimateDE::getParams(mu_0[c],params,&al0,&be0);
            alpha = al0 + RC / 2.0;
            betaPar = lambda0*mu_00*mu_00;

            sum=0;
            sumSq=0;
            for(r=0;r< RC;r++){
               sum += tr[c][r][n];
               sumSq += tr[c][r][n]*tr[c][r][n];
            }
            betaPar += sumSq - (lambda0*mu_00 + sum)*(lambda0*mu_00 + sum) /
               (lambda0 + RC);
            normMu[c]= (lambda0*mu_00 + sum) / (lambda0 + RC);
            beta = be0 + betaPar / 2 ;
            gammaDistribution.param(gDP(alpha, 1.0/beta));
            prec=gammaDistribution(rng_mt);

            normalDistribution.param(nDP(normMu[c], 1/sqrt(prec *(lambda0 + RC))));
            mu = normalDistribution(rng_mt);
            // save sample
            samples[c][n] = mu;
            vars[n] = 1/(prec *(lambda0 + RC));
         }
#ifdef BIOC_BUILD
	 R_CheckUserInterrupt();
#endif
      }
      // }}}
      for(c=0;c<C;c++){
         mu_0[c] = 0;
         for(n=0;n<N;n++)mu_0[c] +=samples[c][n];
         mu_0[c] /= N;
      }
      pplr = 0;
      logFC = 0;
      for(n=0;n<N;n++){
         if(samples[1][n] > samples[0][n])pplr+=1;
         logFC += samples[1][n]-samples[0][n];
         difs[n] = samples[1][n]-samples[0][n];
      }
      // Use log2FC
      logFC /= log(2);
      logFC /= N;
      pplr /= N;
      sort(difs.begin(),difs.end());
      cfLow = difs[(long)(N/100.*args.getD("cf"))];
      cfHigh = difs[(long)(N-N/100.*args.getD("cf"))];
      outF<<pplr<<" "<<logFC<<" "<<cfLow<<" "<<cfHigh;
      for(c=0;c<C;c++)outF<<" "<<mu_0[c];
      outF<<endl;
      if(args.flag("samples")){//{{{
         for(c=0;c<C;c++){
            for(n=0;n<N;n++)outFiles[c]<<samples[c][n]<<" ";
            outFiles[c]<<endl;
         }
         for(n=0;n<N;n++){
            outFiles[C]<<vars[n]<<" ";
         }
         outFiles[C]<<endl;
      }//}}}
   }
   // Close and exit {{{
   if(args.flag("samples")){
      for(c=0;c<C+1;c++)outFiles[c].close();
   }
   outF.close();
   if(args.verbose)message("DONE\n");
   // }}}
   return 0;
}

#ifndef BIOC_BUILD
int main(int argc,char* argv[]){
   return estimateDE(&argc,argv);
}
#endif

namespace ns_estimateDE {

bool initializeOutputFile(long M, long N, long C, const ArgumentParser &args, ofstream *outF, ofstream outFiles[]){//{{{
   if(args.flag("samples")){
      // If samples flag is set, then write condition mean expression samples into -C?.est files.
      // Also write variance samples into samples file.
      stringstream fnStream;
      string fileName;
      // Initialize samples files.
      for(long c=0;c<C;c++){
         fnStream.str("");
         fnStream<<args.getS("outFilePrefix")<<"-C"<<c<<".est";
         fileName = fnStream.str();
         outFiles[c].open(fileName.c_str());
         if(! outFiles[c].is_open()){
            error("Unable to open output file %s\n",fileName.c_str());
            return false;
         }
         // Write header for samples file.
         outFiles[c]<<"# Inferred means\n";
         outFiles[c]<<"# condition "<<c<<endl;
         outFiles[c]<<"# ";
         for(long i=0;i<(long)args.args().size();i++){
            outFiles[c]<<args.args()[i]<<" ";
         }
         outFiles[c]<<"\n# lambda_0 "<<args.getD("lambda0")<<"\n# T \n# M "<<M<<"\n# N "<<N<<endl;
      }
      // Initialize file for variances.
      string varFileName = args.getS("outFilePrefix")+".estVar";
      outFiles[C].open(varFileName.c_str());
      if(! outFiles[C].is_open()){
         error("Unable to open output file %s\n",varFileName.c_str());
         return false;
      }
      // Write header for variance file.
      outFiles[C]<<"# infered variances\n";
      outFiles[C]<<"\n# lambda_0 "<<args.getD("lambda0")<<"\n# T \n# M "<<M<<"\n# N "<<N<<endl;
   }
   // Initialize PPLR file.
   string outFileName = args.getS("outFilePrefix")+".pplr";
   outF->open(outFileName.c_str());
   if(! outF->is_open()){
      error("Unable to open output file %s\n",outFileName.c_str());
      return false;
   }
   // Write header for PPLR file.
   *outF<<"# ";
   for(long i=0;i<(long)args.args().size();i++){
      *outF<<args.args()[i]<<" ";
   }
   *outF<<"\n# lambda_0 "<<args.getD("lambda0")<<"\n# T \n# M "<<M<<"\n# N "<<N<<"\n# Columns:\n";
   *outF<<"# PPLR log2FC ConfidenceLow ConfidenceHigh [mean condition mean expressions]"<<endl;
   return true;
}//}}}

void getParams(double expr,const vector<paramT> &params, double *alpha, double *beta){//{{{
   long i=0,j=params.size()-1,k;
   if(expr<=params[0].expr){
      *alpha=params[0].alpha;
      *beta=params[0].beta;
      return;
   }
   if(expr>=params[j].expr){
      *alpha=params[j].alpha;
      *beta=params[j].beta;
      return;
   }
   while(j-i>1){
      k=(i+j)/2;
      if(params[k].expr<=expr)i=k;
      else j=k;
   }
   if(expr-params[i].expr<params[j].expr-expr)k=i;
   else k=j;
   
   *alpha=params[k].alpha;
   *beta=params[k].beta;
}//}}}

bool readParams(const ArgumentParser &args, long *parN, vector<paramT> *params){//{{{
   ifstream paramF(args.getS("parFileName").c_str());
   FileHeader fh(&paramF);
   if((!fh.paramsHeader(parN))||(*parN==0)){
      error("Main: Problem loading parameters file %s\n",args.getS("parFileName").c_str());
      return false;
   }
   // Vector of parameters: (mean expression, (alpha, beta) )
   params->resize(*parN);
   for(long i=0;i<*parN;i++){
      paramF>>(*params)[i].alpha>>(*params)[i].beta>>(*params)[i].expr;
   }
   fh.close();
   sort(params->begin(),params->end());
   if(args.verb())message("Parameters loaded.\n");
   return true;
}//}}}
}
