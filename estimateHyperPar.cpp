/*
 * Hyperparameter model in estimate[*]HyperPar.cpp always depends on the model used in 
 *  relevant estimate[*]DE.cpp
 */
// DECLARATIONS: {{{
#include <cmath>
#include <algorithm>
#include "boost/random/normal_distribution.hpp"
#include "boost/random/uniform_01.hpp"
#include "boost/random/mersenne_twister.hpp"
using namespace std;

#include "PosteriorSamples.h"
#include "MyTimer.h"
#include "ArgumentParser.h"
#include "lowess.h"
#include "TranscriptExpression.h"
#include "common.h"
//}}}
// Defaults: {{{
#define ALPHA_PROP 0.1
#define BETA_PROP 0.08
#define subM_MIN 10
#define subM_MAX 5000
#define SAMPLES_N 2
#define MAX_ITER 1000
#define MAX_RETRIES 10
#define MAX_PARAM 5000
#define LOG_ZERO -1000
//}}}

struct paramT {//{{{
   double e,a,b;
   bool operator< (const paramT& d2) const{
      return e<d2.e;
   }
};//}}}


extern "C" int estimateHyperPar(int *argc,char* argv[]){
string programDescription =
"Estimate expression dependent hyperparameters from the dataset.\n\
   [sample Files] should contain transposed MCMC samples from replicates.\n\
   To distinguish conditions use C between them e.g.:\n\
      samplesC1-R1.rpkm samplesC1-R2.rpkm C samplesC2-R1.rpkm samplesC2-R2.rpkm";
   // Intro: {{{
   buildTime(argv[0],__DATE__,__TIME__);
   // Set options {{{
   ArgumentParser args(programDescription,"[sampleFiles]",1);
   args.addOptionB("V","veryVerbose","veryVerbose",0,"More verbose output.");
   args.addOptionS("o","outFile","outFileName",1,"Name of the output file.");
   args.addOptionS("p","paramsAllFile","paramsAllFileName",0,"Name of the file to which to store all parameter values generated prior to lowess smoothing.");
   args.addOptionS("","meanFile","meanFileName",0,"Name of the file containing joint mean and variance.");
   args.addOptionL("g","groupsNumber","groupsN",0,"Number of groups of transcript of similar size.",200);
   args.addOptionL("s","samplesNumber","samplesN",0,"Number of samples generated for each group.",SAMPLES_N);
   args.addOptionD("l","lambda0","lambda0",0,"Precision scaling parameter lambda0.",2.0);
   args.addOptionD("","exThreshold","exT",0,"Threshold of lowest expression for which the estimation is done.",-5);
   args.addOptionB("S","smoothOnly","smoothOnly",0,"Input file contains previously sampled hyperparameters which should smoothed only.");
   args.addOptionD("","lowess-f","lowess-f",0,"Parameter F for lowess smoothing specifying amount of smoothing.",0.2);
   args.addOptionL("","lowess-steps","lowess-steps",0,"Parameter Nsteps for lowess smoothing specifying number of iterations.",5);
   args.addOptionB("","noforce","noforce",0,"Do not force smoothing of the parameters.",true);
   args.addOptionS("","norm","normalization",0,"Normalization constants for each input file provided as comma separated list of doubles (e.g. 1.0017,1.0,0.9999 ).");
   if(!args.parse(*argc,argv))return 0;
   // }}}

   MyTimer timer;
   timer.start(1);
   long i,M=0,N,RTN,C;
   bool logged=false,storeAll=args.isSet("paramsAllFileName");
   vector<paramT> params;
   paramT param;
   TranscriptExpression trExp;

   if(! args.flag("smoothOnly")){
      if(! args.isSet("meanFileName")){
         error("Main: Please provide mean file name.\n");
         return 1;
      }
      trExp.readExpression(args.getS("meanFileName"), MEAN_VARIANCE);
      M = trExp.getM();
      if(args.verbose)message("M: %ld\n",M);
      logged = trExp.isLogged();
      if(args.verbose){
         if(logged)message("Using logged values.\n");
         else message("NOT using logged values.\n");
      }
      trExp.doSort(true);
   } 
   
   ofstream outF(args.getS("outFileName").c_str());
   if(!outF.is_open()){
      error("Main: Out file open failed.\n");
      return 1;
   }
   ///}}}

   if(args.flag("smoothOnly")){ 
      // Reading previously sampled parameters. {{{
      ifstream paramsF(args.args()[0].c_str());
      if(!paramsF.is_open()){
         error("Main: Input file open failed.\n");
         return 1;
      }
      // Copy header lines.
      string strBuffer;
      while((paramsF.good())&&(paramsF.peek()=='#')){
         getline(paramsF,strBuffer,'\n');
         outF<<strBuffer<<endl;
      }
      // Read parameters.
      while(paramsF.good()){
         while((paramsF.good())&&(paramsF.peek()=='#')){
            paramsF.ignore(10000000,'\n');
         }
         paramsF>>param.a>>param.b>>param.e;
         if(paramsF.good())
            params.push_back(param);
         paramsF.ignore(10000000,'\n');
      }
      paramsF.close();
      // }}}
   }else{ 
      // Sampling parameters based on data
      // Read conditions {{{   
      Conditions cond;
      if(! cond.init("NONE", args.args(), &C, &M, &N)){
         error("Main: Failed loading MCMC samples.\n");
         return 0;
      }
      if(args.isSet("normalization")){
         if(! cond.setNorm(args.getTokenizedS2D("normalization"))){
            error("Main: Appying normalization constants failed.\n");
            return 1;
         }
      }
      RTN = cond.getRN();
      if(args.verbose)message("Total replicates: %ld\n",RTN);

      ofstream paramsF;
      if(storeAll){
         paramsF.open(args.getS("paramsAllFileName").c_str());
         if(!paramsF.is_open()){
            error("Main: Failed opening %s.\n",(args.getS("paramsAllFileName")).c_str());
            return 0;
         }
         paramsF<<"# lambda0 "<<args.getD("lambda0")<<endl;
      }
      // }}}
      // Declarations {{{
      vector<long double> mu0(subM_MAX,0);
      vector<vector<vector<double> > > tr(subM_MAX,vector<vector<double> >(RTN));
      vector<vector<long double> > bAdd(subM_MAX,vector<long double> (C,0));
      boost::random::mt11213b rng_mt(time(NULL));
      boost::random::uniform_01<long double> uniformDistribution;
      boost::random::normal_distribution<long double> normalDistributionA,normalDistributionB;
      typedef boost::random::normal_distribution<long double>::param_type nDP;
      
      long double alpha,beta,alphaP,betaP,prob,probAll,probC,mean,old_mult,proposalMultiplier,acceptR,sum,sumS,lambda0,exDelta,exLast;
      long samp,samplesN,samplesREDO,maxIter,r,c,m,curM,Rc,subM;
      bool breaked=false,good=false;
      //}}}
      // Initial values {{{ 
      alpha=uniformDistribution(rng_mt)*10.0;
      beta=uniformDistribution(rng_mt)*5.0;
      old_mult=0;
      proposalMultiplier=2.0;
      prob = 0;
      lambda0 = args.getD("lambda0");
      samplesN = args.getL("samplesN");
      curM=0;
      exDelta = (trExp.exp(0)-trExp.exp(M-1))/args.getL("groupsN");
      exLast = trExp.exp(0);
      if(args.verbose)message("Expression step: %Lg\n",exDelta);
      // }}}
      timer.split();
      if(args.verbose)message("Running sampler.\n");
      while(curM<M){
         // Reading next group of transcripts {{{
         mean=0;
         m = 0;
         while((curM<M)&&(m<subM_MAX)){
            if(trExp.exp(curM)<args.getD("exT")){
               if(args.verbose)message("skipping expression: %lg\n",trExp.exp(curM));
               break;
            }
            for(r=0;r<RTN;r++){
               good = cond.getTranscript(r, trExp.id(curM), tr[m][r],samplesN+MAX_RETRIES);
               if(!good)break;
               if(logged) // should check whether samples are logged as well
                  for(samp=0;samp<samplesN+MAX_RETRIES;samp++){
                     tr[m][r][samp] = (tr[m][r][samp] == 0)? LOG_ZERO:log(tr[m][r][samp]);
                  }
            }
            if(good){
               mu0[m]=trExp.exp(curM);
               mean+=mu0[m];
               m++;
            }
            curM++;
            if(args.flag("veryVerbose"))if(progressLog(curM,M))timer.split(0,'m');
            if((m>=subM_MIN)&&(exDelta<exLast-trExp.exp(curM-1)))break;
         }
         exLast = trExp.exp(curM-1);
         if(m<subM_MIN)break;
         subM = m;
         mean/=subM;
         if(args.flag("veryVerbose"))message("# mean: %Lg  subM: %ld\n",mean,subM);
         if(storeAll)paramsF<<"# mean: "<<mean<<"  subM: "<<subM<<endl;
         samplesREDO = 0;
         //}}}
         for(samp=0;samp<samplesN+samplesREDO;samp++){
            // Computing Badd_gc and initializing {{{
            for(m=0;m<subM;m++){
               i=0; // counter over all replicates;
               for(c=0;c<C;c++){
                  sum = 0;
                  sumS = 0;
                  Rc=cond.getRC(c);
                  for(r=0;r<Rc;r++){
                     sum += tr[m][i][samp];
                     sumS += tr[m][i][samp]*tr[m][i][samp];
                     i++;
                  }
                  bAdd[m][c]=0.5*(sumS + mu0[m]*mu0[m]*lambda0 - 
                        (sum+mu0[m]*lambda0)*(sum+mu0[m]*lambda0)/(lambda0+Rc));
               }
            } 
            acceptR=0;
            old_mult=0;
            proposalMultiplier=proposalMultiplier*2.0;
            normalDistributionA.param(nDP(0,ALPHA_PROP*proposalMultiplier));
            normalDistributionB.param(nDP(0,BETA_PROP*proposalMultiplier));
            maxIter=0;
            breaked = false;
#ifdef BIOC_BUILD
	    R_CheckUserInterrupt();
#endif
            //}}}
            while((acceptR<0.25)||(acceptR>0.5)||(old_mult!=proposalMultiplier)){
               // Convergence control based on acceptance ratio. {{{
               maxIter++;
               if(maxIter>MAX_ITER){
                  if(args.flag("veryVerbose"))
                     message("(BREAKED acceptR %Lg mult %Lg)\n",acceptR,proposalMultiplier);
                  if(storeAll)
                     paramsF<<"#(BREAKED acceptR "<<acceptR<<" mult "<<proposalMultiplier<<")"<<endl;
                  breaked=true;
                  break;
               }
               if((alpha>MAX_PARAM)||(beta>MAX_PARAM)){
                  if(args.flag("veryVerbose"))
                     message("(OVERFLOW acceptR %Lg mult %Lg)\n",acceptR,proposalMultiplier);
                  if(storeAll)
                     paramsF<<"#(OVERFLOW acceptR "<<acceptR<<" mult "<<proposalMultiplier<<")"<<endl;
                  breaked=true;
                  break;
               }
               old_mult=proposalMultiplier;
               if(acceptR<0.25)proposalMultiplier/=1.02;
               if(acceptR>0.5)proposalMultiplier*=1.02;
               if(old_mult!=proposalMultiplier){
                  normalDistributionA.param(nDP(0,ALPHA_PROP*proposalMultiplier));
                  normalDistributionB.param(nDP(0,BETA_PROP*proposalMultiplier));
               }
               //}}}
               acceptR=0;
#ifdef BIOC_BUILD
	       R_CheckUserInterrupt();
#endif
               for(i=0;i<1000;i++){ // Sampling 1000 samples {{{
                  alphaP = alpha + normalDistributionA(rng_mt);
                  if(alphaP<0)alphaP = -alphaP;
                  betaP= beta + normalDistributionB(rng_mt);
                  if(betaP<0)betaP = -betaP;
                  if((alphaP==0)||(betaP==0)){
                     prob=0;
                  }else{
                     prob = 1.0;
                     probAll = pow(betaP,alphaP) / pow(beta,alpha);
                     for(c=0;c<C;c++){
                        probC = lgamma(alphaP + cond.getRC(c)/2.0)+
                            lgamma(alpha) -
                            lgamma(alpha + cond.getRC(c)/2.0) -
                            lgamma(alphaP);
                        probC = probAll * exp(probC);
                        for(m=0;m<subM;m++){
      //                  message(" (var_g %lg) (pow %lg %lg %lg) ",bAdd[g]/2.0,pow(beta+bAdd[g]/2, alpha),pow(betaP+bAdd[g]/2, alphaP),pow((beta+bAdd[g]/2)/(betaP+bAdd[g]/2),SUB_N/2));
                           prob *= probC;
                           prob *= pow(beta+bAdd[m][c], alpha) / 
                                   pow(betaP+bAdd[m][c], alphaP); 
                           prob *= pow( (beta+bAdd[m][c])/(betaP+bAdd[m][c]), (long double)(cond.getRC(c)/2.0));
                        }
                     }
                     if((prob>1.0)||(uniformDistribution(rng_mt)< prob)){
                        alpha=alphaP;
                        beta=betaP;
                        acceptR++;
                     }
                  }
               } //}}}
               acceptR/=i;
            }
            // Save generated parameters {{{
            if(storeAll)
               paramsF<<"#(acceptR "<<acceptR<<" mult "<<proposalMultiplier<<" iter "<<maxIter<<")"<<endl;
            if(!breaked){
               if(args.flag("veryVerbose")) message("%Lg  %Lg\n",alpha,beta);
               if(storeAll) paramsF<<alpha<<" "<<beta<<" "<<mean<<endl;
               param.e=mean;
               param.a=alpha;
               param.b=beta;
               params.push_back(param);
            }else{
               if(args.flag("veryVerbose")) message("# %Lg %Lg %Lg\n",alpha,beta,mean);
               if(storeAll) paramsF<<"# "<<alpha<<"  "<<beta<<endl;
               proposalMultiplier=2;
               normalDistributionA.param(nDP(0,ALPHA_PROP*proposalMultiplier));
               normalDistributionB.param(nDP(0,BETA_PROP*proposalMultiplier));
               alpha=uniformDistribution(rng_mt)*10.0;
               beta=uniformDistribution(rng_mt)*5.0;
               if(samplesREDO<MAX_RETRIES){
                  samplesREDO++;
               }
            }
            //}}}
         }
         if((args.verbose)&&(!args.flag("veryVerbose"))){
            message(".");
#ifndef BIOC_BUILD
            fflush(stdout);
#endif
         }
      }
      cond.close();
      if(storeAll)paramsF.close();
      outF<<"# lambda0 "<<args.getD("lambda0")<<endl;
      if(args.verbose)message("\nSampling done.\n");
   }
   sort(params.begin(),params.end());
   long pAll=(long)params.size(), pDistinct;
   if(args.verbose)message("Have %ld parameters to smooth.\n",pAll);
   vector<double> exp(pAll),alp(pAll),bet(pAll),alpS,betS;
   for(i=0;i<pAll;i++){
      exp[i]=params[i].e;
      alp[i]=params[i].a;
      bet[i]=params[i].b;
   }
   double f = args.getD("lowess-f");
   long iter = args.getL("lowess-steps"),iterAdd;
   bool redoSmooth;
   for(iterAdd=0;iterAdd<6;iterAdd++){ // Increase iteration if anything is <=0
      redoSmooth = false;
      lowess(exp,alp,f,iter+iterAdd,alpS);
      for(i=0;i<pAll;i++)
         if(alpS[i]<=0){
            redoSmooth = true;
            if(args.flag("veryVerbose"))message(" negative alpha: %lg exp: %lg\n",alpS[i],exp[i]);
         }
      if(!redoSmooth)break;
      if(args.verbose)message("Re-Smoothing alpha.\n");
   }
   outF<<"# alphaSmooth f: "<<f<<" nSteps: "<<iter+iterAdd<<endl;
   if(args.verbose)message("# alphaSmooth f: %lg nSteps: %ld\n",f,iter+iterAdd);
   if((iterAdd==6)&&(args.flag("noforce"))){
      error("Main: Unable to produce smooth alpha >0.\nTry adjusting the parameter lowess-f.\n");
      outF.close();
      remove(args.getS("outFileName").c_str());
      return 0;
   }
   for(iterAdd=0;iterAdd<6;iterAdd++){ // Increase iteration if anything is <=0
      redoSmooth = false;
      lowess(exp,bet,f,iter+iterAdd,betS);
      for(i=0;i<pAll;i++)
         if(betS[i]<=0){
            redoSmooth = true;
            if(args.flag("veryVerbose"))message(" negative beta: %lg exp: %lg\n",betS[i],exp[i]);
         }
      if(!redoSmooth)break;
      if(args.verbose)message("Re-Smoothing beta.\n");
   }
   outF<<"# betaSmooth f: "<<f<<" nSteps: "<<iter+iterAdd<<endl;
   if(args.verbose)message("# betaSmooth f: %lg nSteps: %ld\n",f,iter+iterAdd);
   if((iterAdd==6)&&(args.flag("noforce"))){
      error("Main: Unable to produce smooth beta >0.\nTry adjusting the parameter lowess-f.\n");
      outF.close();
      remove(args.getS("outFileName").c_str());
      return 0;
   }
   if(!args.flag("noforce")){
      for(i=0;i<pAll;i++)
         while((i<pAll)&&((alpS[i]<=0)||(betS[i]<=0))){
            message("Removing: %lg %lg %lg\n",alpS[i],betS[i],exp[i]);
            alpS.erase(alpS.begin()+i); betS.erase(betS.begin()+i); exp.erase(exp.begin()+i);
            pAll = alpS.size();
         }
   }
   pDistinct = 1;
   for(i=1;i<pAll;i++)if(exp[i]!=exp[i-1])pDistinct++;
   outF<<"# PN "<<pDistinct<<" hyperparameters"<<endl;
   outF<<"# columns: alpha beta expression "<<endl;
   outF<<alpS[0]<<" "<<betS[0]<<" "<<exp[0]<<endl;
   for(i=1;i<pAll;i++)
      if(exp[i]!=exp[i-1])outF<<alpS[i]<<" "<<betS[i]<<" "<<exp[i]<<endl;
   outF.close();
   if(args.verbose){message("DONE.\n");timer.stop(1,'m');}
   return 0;
}

#ifndef BIOC_BUILD
int main(int argc,char* argv[]){
   return estimateHyperPar(&argc,argv);
}
#endif
