/*
 * Original model applying the DE model to individual sets of samples independently.
 * One set of samples == 1 sample from each replicate of each condition.
 */
#include<algorithm>
#include<cmath>
#include<fstream>
#include<sstream>
#include "boost/random/gamma_distribution.hpp"
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/normal_distribution.hpp"

using namespace std;

#include "ArgumentParser.h"
#include "common.h"
#include "misc.h"
#include "MyTimer.h"
#include "PosteriorSamples.h"

#define FF first
#define SS second

//#define PERCENT 0.9

#define LAMBDA_0 2.0

using ns_params::paramT;

namespace ns_estimateDE {

// Open and write headers into appropriate output files.
// The size of outFiles[] should be C+1.
// Returns true if everything went OK.
bool initializeOutputFile(long C, long M, long N, const ArgumentParser &args, ofstream *outF, ofstream outFiles[]);
// For a given mean expression expr finds alpha and beta for which were estimated for a closes expression.
void getParams(double expr,const vector<paramT> &params, paramT *par);
// Read transcript m into tr and prepare mu_0 and mu_00, cond does not really change.
void readNextTranscript(long m, long C, long N, Conditions *cond, const vector<paramT> &params, vector<vector<vector<double> > > *tr, vector<paramT> *curParams, double *mu_00);

}

extern "C" int estimateDE(int *argc,char* argv[]){
string programDescription =
"Estimate differential expression from the dataset.\n\
   [sampleFiles] should contain transposed MCMC samples from replicates.\n\
   To distinguish conditions use C between them e.g.:\n\
      samplesC1-R1.rpkm samplesC1-R2.rpkm C samplesC2-R1.rpkm samplesC2-R2.rpkm";
   // Intro: {{{
   ArgumentParser args(programDescription,"[sampleFiles]",1);
   args.addOptionS("o","outPrefix","outFilePrefix",1,"Prefix for the output files.");
   args.addOptionS("p","parameters","parFileName",1,"File containing estimated hyperparameters.");
   args.addOptionB("s","samples","samples",0,"Produce samples of condition mean expression apart from PPLR and confidence.");
   args.addOptionD("l","lambda0","lambda0",0,"Parameter lambda_0.",LAMBDA_0);
   args.addOptionD("","confidenceInterval","cf",0,"Percentage for confidence intervals.", 95);
   args.addOptionS("","norm","normalization",0,"Normalization constants for each input file provided as comma separated list of doubles (e.g. 1.0017,1.0,0.9999 ).");
   args.addOptionL("","seed","seed",0,"Random initialization seed.");
   if(!args.parse(*argc,argv))return 0;
   if(args.verbose)buildTime(argv[0],__DATE__,__TIME__);
   //}}}
   /*
    * N - number of samples in one replicate (the smallest number for replicates with different N_r)
    * M - number of transcripts
    * C - number of conditions
    */
   long C,M,N;
   vector<paramT> params;
   Conditions cond;
   // Open file with hyper parameters and read those in.
   if(!ns_params::readParams(args.getS("parFileName"), &params)) return 1;
   if(args.verb())message("Parameters loaded.\n");
   // Initialize sample files handled by object cond.
   if(!ns_misc::readConditions(args, &C, &M, &N, &cond)) return 1;
   // Initialize output files.
   ofstream outF;
   ofstream outFiles[C+1];
   // Use standard array as we don't want to bother with vector of pointers.
   if(!ns_estimateDE::initializeOutputFile(C, M, N, args, &outF, outFiles)) return 1;

   // variables {{{
   vector<vector<vector<double> > > tr(C);
   vector<paramT> curParams(C);
   vector<vector<double> > samples(C,vector<double>(N));
   vector<double> vars(N);
   vector<double> mu_c(C);
//   vector<vector<double> > mus(C,vector<double>(N,0));
//   vector<double> vars(N);
   long c,c2,m,n,r;
   double prec,var,sum,sumSq,alpha,beta,betaPar,mu_00,normMu;
   double lambda0 = args.getD("lambda0");
   long RC;
   MyTimer timer;
   boost::random::mt11213b rng_mt(ns_misc::getSeed(args));
   boost::random::gamma_distribution<long double> gammaDistribution;
   typedef boost::random::gamma_distribution<long double>::param_type gDP;
   boost::random::normal_distribution<long double> normalDistribution;
   typedef boost::random::normal_distribution<long double>::param_type nDP;
   double log2FC, pplr, ciLow, ciHigh;
   vector<double> difs(N);
   // }}}

   if(args.verbose){ //{{{
      timer.split();
      message("Sampling condition mean expression.\n");
   }//}}}
   for(m=0;m<M;m++){
      if(progressLog(m,M))timer.split();
      // Read into tr and assign hyperparameters into curParams, initialize mu_00.
      // cond does not really change, just reads more data from file.
      ns_estimateDE::readNextTranscript(m, C, N, &cond, params, &tr, &curParams, &mu_00);
      // Zero "mean condition mean expression".
      mu_c.assign(C,0);
      // Sample condition mean expressions {{{
      for(n=0;n<N;n++){
         for(c=0;c<C;c++){
            RC = cond.getRC(c);
            alpha = curParams[c].alpha + RC / 2.0;
            betaPar = lambda0*mu_00*mu_00;

            sum=0;
            sumSq=0;
            for(r=0;r< RC;r++){
               sum += tr[c][r][n];
               sumSq += tr[c][r][n]*tr[c][r][n];
            }
            betaPar += sumSq - (lambda0*mu_00 + sum)*(lambda0*mu_00 + sum) /
               (lambda0 + RC);
            normMu= (lambda0*mu_00 + sum) / (lambda0 + RC);
            beta = curParams[c].beta + betaPar / 2 ;
            // Set parameters of gamma distribution.
            gammaDistribution.param(gDP(alpha, 1.0/beta));
            // Sample precision.
            prec = gammaDistribution(rng_mt);
            // Variance, the precision is scaled by (lambda0+RC).
            var = 1/(prec *(lambda0 + RC));
            vars[n] = var;

            // Set parameter for normal distribution.
            normalDistribution.param(nDP(normMu, sqrt(var)));
            // Sample condition mean.
            samples[c][n] = normalDistribution(rng_mt);
            mu_c[c] += samples[c][n];
         }
         R_INTERUPT;
      }
      // }}}
      // Compute condition mean for each condition.
      for(c=0;c<C;c++) mu_c[c] /= N;
      // Calculate and write pplr for each pair of conditions. {{{
      for(c=0;c<C;c++){
         for(c2=c+1;c2<C;c2++){
            pplr = 0;
            for(n=0;n<N;n++)
               if(samples[c2][n] > samples[c][n])pplr+=1;
            pplr/=N;
            outF<<pplr<<" ";
         }
      }
      // }}}
      // Calculate log2FC; write log2FC and CIs for each pair of conditions. {{{
      for(c=0;c<C;c++){
         for(c2=c+1;c2<C;c2++){
            for(n=0;n<N;n++)
               difs[n] = samples[c2][n]-samples[c][n];
            ns_misc::computeCI(args.getD("cf"), &difs, &ciLow, &ciHigh);
            ciLow /= log(2);
            ciHigh /= log(2);
            log2FC = (mu_c[c2] - mu_c[c])/log(2);
            outF<<log2FC<<" "<<ciLow<<" "<<ciHigh<<" ";
         }
      }
      // }}}
      // Write logged condition mean for each condition. No space before EOL. {{{
      for(c = 0; c < C-1; c++)outF<<mu_c[c]<<" ";
      outF<<mu_c[C-1]<<endl;
      // }}}
      // Write samples if necessary. {{{ 
      if(args.flag("samples")){
         for(c=0;c<C;c++){
            for(n=0;n<N;n++)outFiles[c]<<samples[c][n]<<" ";
            outFiles[c]<<endl;
         }
         // Save sampled variance as well.
         for(n=0;n<N;n++) outFiles[C]<<vars[n]<<" ";
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

bool initializeOutputFile(long C, long M, long N, const ArgumentParser &args, ofstream *outF, ofstream outFiles[]){//{{{
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
         outFiles[c]<<"# Inferred condition mean log expression.\n"
                      "# condition "<<c+1
                    <<"\n# ";
         for(long i=0;i<(long)args.args().size();i++){
            outFiles[c]<<args.args()[i]<<" ";
         }
         outFiles[c]<<"\n# lambda_0 "<<args.getD("lambda0")<<"\n# T (Mrows_Ncols) L (logged)\n# M "<<M<<"\n# N "<<N<<endl;
      }
      // Initialize file for variances.
      string varFileName = args.getS("outFilePrefix")+".estVar";
      outFiles[C].open(varFileName.c_str());
      if(! outFiles[C].is_open()){
         error("Unable to open output file %s\n",varFileName.c_str());
         return false;
      }
      // Write header for variance file.
      outFiles[C]<<"# Inferred variances in last condition.\n"
                   "# lambda_0 "<<args.getD("lambda0")
                 <<"\n# T \n# M "<<M<<"\n# N "<<N
                 <<endl;
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
   *outF<<"\n# lambda_0 "<<args.getD("lambda0")<<"\n# T \n# M "<<M<<"\n# N "<<N<<"\n"
        <<"# Conditions: C "<<C<<" Condition pairs("<<C*(C-1)/2<<"): ";
   for(long c=0;c<C;c++)
      for(long c2=c+1;c2<C;c2++)
      *outF<<c+1<<"~"<<c2+1<<" ";
   *outF<<"\n# Columns contain PPLR for each pair of conditions, "
          "log2 fold change with confidence intervals for each pair of conditions and "
          "log mean condition mean expression for each condition.\n"
          "# CPxPPLR CPx(log2FC ConfidenceLow ConfidenceHigh) "
          "Cx(log mean condition mean expressions)"
        <<endl;
   return true;
}//}}}

void getParams(double expr,const vector<paramT> &params, paramT *par){//{{{
   long i=0,j=params.size()-1,k;
   if(expr<=params[0].expr){
      par->alpha=params[0].alpha;
      par->beta=params[0].beta;
      return;
   }
   if(expr>=params[j].expr){
      par->alpha=params[j].alpha;
      par->beta=params[j].beta;
      return;
   }
   while(j-i>1){
      k=(i+j)/2;
      if(params[k].expr<=expr)i=k;
      else j=k;
   }
   if(expr-params[i].expr<params[j].expr-expr)k=i;
   else k=j;
   
   par->alpha=params[k].alpha;
   par->beta=params[k].beta;
}//}}}

void readNextTranscript(long m, long C, long N, Conditions *cond, const vector<paramT> &params, vector<vector<vector<double> > > *tr, vector<paramT> *curParams, double *mu_00){//{{{
   double divT = 0, divC, mu_0;
   long c,r,n,RC;
   *mu_00 = 0;
   for(c=0;c<C;c++){
      mu_0=0;
      divC=0;
      RC = cond->getRC(c);
      if((long)(*tr)[c].size() < RC){
         (*tr)[c].resize( RC );
      }
      for(r=0;r<RC;r++){
         if(cond->getTranscript(c, r , m, (*tr)[c][r]), N){
            for(n=0;n<N;n++){
               // Log the expression samples if the files don't have logged flag set.
               if(!cond->logged())(*tr)[c][r][n] = ((*tr)[c][r][n] == 0)? ns_misc::LOG_ZERO : log ((*tr)[c][r][n] );
               mu_0+=(*tr)[c][r][n];
            }
            divC+=1;
         }else{
            warning("Main: Condition %ld replicate %ld does not seem to have transcript %ld.\n",c,r,m);
         }
      }
      R_INTERUPT;
      if(divC>0){
         mu_0 /= (divC * N); // take mean over all replicates
         *mu_00+=mu_0;
         divT++;
      }
      getParams(mu_0, params, &(*curParams)[c]);
   }
   *mu_00/=divT; 
}//}}}

}
