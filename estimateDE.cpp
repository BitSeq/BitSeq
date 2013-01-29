/*
 * Original model applying the DE model to individual sets (from all replicates) of samples independently  
 *
 */
#include<cmath>
#include<algorithm>
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/gamma_distribution.hpp"
#include "boost/random/normal_distribution.hpp"

using namespace std;

#include "PosteriorSamples.h"
#include "FileHeader.h"
#include "MyTimer.h"
#include "ArgumentParser.h"
#include "common.h"

#define FF first
#define SS second

//#define PERCENT 0.9

#define LAMBDA_0 2.0
//#define MU_0 (10/M)
// CHECK: ------------------------------
//#define ALPHA 1.2
//#define BETA 3.2
// CHECK: ------------------------------
const double LOG_ZERO=-1000;

void getParams(double &alpha, double &beta, double expr, vector<pair<double,pair<double,double> > > &params){//{{{
   long i=0,j=params.size()-1,k;
   if(expr<=params[0].FF){
      alpha=params[0].SS.FF;
      beta=params[0].SS.SS;
      return;
   }
   if(expr>=params[j].FF){
      alpha=params[j].SS.FF;
      beta=params[j].SS.SS;
      return;
   }
   while(j-i>1){
      k=(i+j)/2;
      if(params[k].FF<=expr)i=k;
      else j=k;
   }
   if(expr-params[i].FF<params[j].FF-expr)k=i;
   else k=j;
   
   alpha=params[k].SS.FF;
   beta=params[k].SS.SS;
}//}}}

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
   long i,C,M,N,m,n,c,r,parN;
   // param file:{{{
   ifstream paramF(args.getS("parFileName").c_str());
   FileHeader fh(&paramF);
   if((!fh.paramsHeader(parN))||(parN==0)){
      error("Main: Problem loading parameters file %s\n",args.getS("parFileName").c_str());
      return 1;
   }
   vector<pair<double, pair<double,double> > > params(parN);
   for(i=0;i<parN;i++){
      paramF>>params[i].SS.FF>>params[i].SS.SS>>params[i].FF;
   }
   fh.close();
   sort(params.begin(),params.end());
   if(args.verbose)message("Parameters loaded.\n");
   //}}}
   // samples files: {{{
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
   // output files:{{{
   ofstream outF;
   ofstream outFiles[C+1];
   if(args.flag("samples")){
      stringstream fName;
      for(c=0;c<C;c++){
         fName.str("");
         fName<<args.getS("outFilePrefix")<<"-C"<<c<<".est";
         outFiles[c].open(fName.str().c_str());
         if(! outFiles[c].is_open()){
            error("Unable to open output file %s\n",(fName.str()).c_str());
            return 1;
         }
         outFiles[c]<<"# Inferred means\n";
         outFiles[c]<<"# condition "<<c<<endl;
         outFiles[c]<<"# ";
         for(i=0;i<(long)args.args().size();i++){
            outFiles[c]<<args.args()[i]<<" ";
         }
         outFiles[c]<<"\n# lambda_0 "<<args.getD("lambda0")<<"\n# T \n# M "<<M<<"\n# N "<<N<<endl;
      }
      outFiles[C].open((args.getS("outFilePrefix")+".estVar").c_str());
      if(! outFiles[C].is_open()){
         error("Unable to open output file %s\n",((args.getS("outFilePrefix")+".estVar")).c_str());
         return 1;
      }
      outFiles[C]<<"# infered variances\n";
      outFiles[C]<<"\n# lambda_0 "<<args.getD("lambda0")<<"\n# T \n# M "<<M<<"\n# N "<<N<<endl;
   }
   outF.open((args.getS("outFilePrefix")+".pplr").c_str());
   if(! outF.is_open()){
      error("Unable to open output file %s\n",((args.getS("outFilePrefix")+".pplr")).c_str());
      return 1;
   }
   outF<<"# ";
   for(i=0;i<(long)args.args().size();i++){
      outF<<args.args()[i]<<" ";
   }
   outF<<"\n# lambda_0 "<<args.getD("lambda0")<<"\n# T \n# M "<<M<<"\n# N "<<N<<"\n# Columns:\n";
   outF<<"# PPLR log2FC ConfidenceLow ConfidenceHigh [mean condition mean expressions]"<<endl;
   // }}}

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
            getParams(al0,be0,mu_0[c],params);
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
