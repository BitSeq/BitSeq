#include<omp.h>
#include<cmath>
#include<algorithm>

#include "CollapsedSampler.h"
#include "GibbsSampler.h"
#include "Sampler.h"
#include "FileHeader.h"
#include "MyTimer.h"
#include "ArgumentParser.h"
#include "TranscriptInfo.h"
#include "transposeFiles.h"
#include "TagAlignments.h"
#include "common.h"


#define DEBUG(x)
#define FF first
#define SS second

//#define LOG_NEED

TranscriptInfo trInfo;

long  M;//, mAll; // M : number of transcripts (include transcript 0 ~ Noise)
//long N, 
long Nunmap; // N: number of read, unmappable read, mappable reads

string outTypeS;
outputType outTypeI;
vector<string> samplesFileNames;
string failedMessage;

void clearDataEE(){
   samplesFileNames.clear();
}

TagAlignments* readData(ArgumentParser &args) {//{{{
   long i,j,num,tid;
   double prb;
   long Ntotal=0,Nmap=0;
   string readId,strand,blank;
   ifstream inFile;
   MyTimer timer;
   TagAlignments *alignments = new TagAlignments();

   // Read alignment probabilities {{{
   inFile.open(args.args()[0].c_str());
   FileHeader fh(&inFile);
   bool newformat=true;
   if((!fh.probHeader(Nmap,Ntotal,newformat)) || (Nmap ==0)){//{{{
      error("Prob file header read failed.\n");
      return NULL;
   }//}}}
   message("Mappings: %ld\n",Nmap);
   message("Ntotal: %ld\n",Ntotal);
   if(Ntotal>Nmap)Nunmap=Ntotal-Nmap;
   else Nunmap=1; //no valid count file assume only one not aligned properly
   alignments->init(Nmap,0,M);
   long mod=10000;
   long bad = 0;
   timer.start();
   for(i = 0; i < Nmap; i++) {
      inFile>>readId>>num;
      if(!newformat)inFile>>blank;
      if(!inFile.good())break;
     //    message("%s %ld\n",(readId).c_str(),num);
      for(j = 0; j < num; j++) {
         if(newformat) inFile>>tid>>prb;
         else inFile>>tid>>strand>>prb;
         if(inFile.fail()){
            inFile.clear();
            // ignore rest of line
            //inFile.ignore(10000000,'\n');
            // ignore other read's alignments
            j=num;
            // this read goes to noise
            tid=0;
            prb=100;
            bad++;
         }
         if((!newformat) && (tid!=0)){
            // these lengths are not shifted
            prb /= trInfo.L(tid-1);
         }
         alignments->pushAlignment(tid, prb);         
      }
      // ignore rest of line
      inFile.ignore(10000000,'\n');

      alignments->pushRead();

#ifdef BIOC_BUILD
      R_CheckUserInterrupt();
#endif
      if((i % mod == 0)&&(i>0)){
         message("  %ld ",i);
         timer.split();
         mod*=10;
      }
   }
   //message("Bad: %ld\n",bad);
   inFile.close();
   long Nhits,NreadsReal;
   alignments->finalizeRead(&M, &NreadsReal, &Nhits);
   // If the transcript info is initialized, check that the number of transcripts has not changed.
   // The number can't be smaller as it starts off with trInfo->M
   if((trInfo.isOK())&&(M > trInfo.getM() + 1)){
      if(outTypeI == RPKM){
         error("Main: Number of transcripts in .prob file is higher than in the .tr file (%ld %ld)!\n",M,trInfo.getM() + 1);
         delete alignments;
         return NULL;
      }else{
         warning("Main: Number of transcripts in .prob file is higher than in the .tr file (%ld %ld)!\n   This can cause problems later on!\n",M,trInfo.getM() + 1);
      }
   }
   //}}}
   if(i<Nmap)message("Read only %ld reads.\n",NreadsReal);
   message("Finished Reading!\nTotal hits = %ld\n",Nhits);
   message("Isoforms: %ld\n",M);
   Nmap = NreadsReal;
   return alignments;
   /* {{{ remapping isoforms to ignore those without any hits
   M = mAll;
   M = isoformsHit;
   isoformMap.assign(M);
   for(i=0,j=0;i<mAll;i++)
      if(readInIsoform[i]!=-1){
         readInIsoform[i]=j;
         isoformMap[j]=i;
         j++;
      }
   for(i=0;i<Sof(alignments);i++){
      alignments[i].setTrId( readInIsoform[ alignments[i].getTrId() ] );
   }
   }}}*/
}//}}}

void MCMC(TagAlignments *alignments,gibbsParameters &gPar,ArgumentParser &args){//{{{
   // Declarations: {{{
   DEBUG(message("Declaratinos:\n"));
   long i,j,samplesHave=0,totalSamples=0,samplesN,chainsN,samplesSave,seed;
   pairD rMean,tmpA,tmpV;
   double rH1,rH2;
   ofstream meansFile,samplesFile[gPar.chainsN()];
   MyTimer timer;
   bool quitNext = false;
   vector<pairD> betwVar(M),withVar(M),s2j(M),totAverage(M),av,var;
   vector<pair<pairD,long> > rHat2(M);
   // }}}
   // Init: {{{
   DEBUG(message("Initialization:\n"));
   samplesN=gPar.samplesN();
   chainsN=gPar.chainsN();
   samplesSave=(gPar.samplesSave()-1)/chainsN+1;

   vector<Sampler*> samplers(chainsN);
   if( ! args.flag("gibbs")){
      for(j=0;j<chainsN;j++)
         samplers[j] = new CollapsedSampler;
   }else{
      for(j=0;j<chainsN;j++)
         samplers[j] = new GibbsSampler;
   }

   timer.start();;
   // parallel block: 
   // make sure that all functions used are CONST and variables are beaing READ or private
   // private: samplesHave
   if(args.isSet("seed"))seed=args.getL("seed");
   else seed = time(NULL);
   if(args.verbose)message("seed: %ld\n",seed);
   for(i=0;i<chainsN;i++){
      // Init samplers
      DEBUG(message("Sampler %ld init.\n",i);)
      samplers[i]->noSave();
      DEBUG(message("init\n");)
      samplers[i]->init(M, samplesN, samplesSave, Nunmap, alignments, gPar.beta(), gPar.dir(), seed);
      DEBUG(message("   seed: %ld\n",seed);)
      // sampler is initialized with 'seed' and then sets 'seed' to new random seed for the next sampler
   }
   #pragma omp parallel for private(samplesHave)
   for(i=0;i<chainsN;i++){
      DEBUG(message(" burn in\n");) 
      for(samplesHave=0;samplesHave<gPar.burnIn();samplesHave++){
        samplers[i]->sample();
#ifdef BIOC_BUILD
         R_CheckUserInterrupt();
#endif
      }
   }
   message("Burn in: %ld DONE. ",gPar.burnIn());
   DEBUG(message(" reseting samplers after burnin\n"));
   for(i=0;i<chainsN;i++){
      samplers[i]->resetSampler(samplesN);
   }
   timer.split(0,'m');
   //}}}
   // Main sampling loop:
   while(1){
      timer.start();
      // Sample: {{{
      // parallel block:
      // make sure that all functions used are CONST and variables are being READ or private
      // private: samplesHave, samplesSkipped
      #pragma omp parallel for private(samplesHave)
      for(i=0;i<chainsN;i++){
         for(samplesHave = 0;samplesHave<samplesN;samplesHave++){
            samplers[i]->sample();
            samplers[i]->update();
#ifdef BIOC_BUILD
            R_CheckUserInterrupt();
#endif
         }
      }
      totalSamples+=samplesN*chainsN;
      message("\nSampling DONE. ");
      timer.split(0,'m');
      //}}}
      // Check for change of parameters: {{{
      gPar.readParameters();
      // }}}
      // Compute convergence statistics {{{
      totAverage.assign(M,pairD(0,0));
      betwVar.assign(M,pairD(0,0));
      withVar.assign(M,pairD(0,0));
      samplesHave = samplesN;
      for(i=0;i<M;i++){
         for(j=0;j<chainsN;j++){
            tmpA = samplers[j]->getAverage(i);
            tmpV = samplers[j]->getWithinVariance(i);
            totAverage[i].FF += tmpA.FF;
            totAverage[i].SS += tmpA.SS;
            withVar[i].FF += tmpV.FF;
            withVar[i].SS += tmpV.SS;
         }
         totAverage[i].FF /= chainsN;
         totAverage[i].SS /= chainsN;
         withVar[i].FF /= chainsN;
         withVar[i].SS /= chainsN;
         for(j=0;j<chainsN;j++){
            tmpA = samplers[j]->getAverage(i);
            betwVar[i].FF += (totAverage[i].FF - tmpA.FF)*(totAverage[i].FF - tmpA.FF);
            betwVar[i].SS += (totAverage[i].SS - tmpA.SS)*(totAverage[i].SS - tmpA.SS);
         }
         betwVar[i].FF /= (chainsN-1.0);
         betwVar[i].SS /= (chainsN-1.0);
      }
      for(i=0;i<M;i++){
         // betwVar[i] *= samplesHave / (chainsN - 1.0);
         rHat2[i].SS=i;
         if(withVar[i].FF == 0 ){
            rHat2[i].FF.FF = 0;
            rHat2[i].FF.SS = 0;
         } else { 
            rHat2[i].FF.FF = (samplesHave - 1.0) / samplesHave + betwVar[i].FF / withVar[i].FF ;
            rHat2[i].FF.SS = (samplesHave - 1.0) / samplesHave + betwVar[i].SS / withVar[i].SS ;
               //betwVar[i] / ( samplesHave * withVar[i] );
         }
      }
      sort(rHat2.rbegin(),rHat2.rend());
      message("rHat (for %ld samples) \n",samplesN);
      rMean.FF=0;
      rMean.SS=0;
      message("   rHat (rHat from subset |    tid | mean theta)\n");
      for(i=0;(i<10) && (i<M);i++){
         rH1 = sqrt(rHat2[i].FF.FF);
         rH2 = sqrt(rHat2[i].FF.SS);
         rMean.FF+=rH1;
         rMean.SS+=rH2;
//         message("   %lf (%lf | %ld | %lf|%lf|%lf)",rHat2[i].FF.FF,rHat2[i].FF.SS,rHat2[i].SS,totAverage[rHat2[i].SS].FF,withVar[rHat2[i].SS].FF,betwVar[rHat2[i].SS].FF/samplesHave);
         if((i<3) || args.verbose){
            message("   %7.4lf (%7.4lf | %6ld | %7.4lf)",rH1,rH2,rHat2[i].SS-1,totAverage[rHat2[i].SS].FF);
            message("\n");
         }
//                  message("   %lf",sqrt(rHat2[i].FF));
      }                  
      rMean.FF /= 10.0;
      rMean.SS /= 10.0;
      message("  Mean rHat of worst 10 transcripts: %lf\n",rMean.FF);
      if(args.flag("scaleReduction"))message("   (target: %.3lf)\n",gPar.targetScaleReduction());
      message("  Mean C0: (");
      for(j=0;j<chainsN;j++)message("%ld ",samplers[j]->getAverageC0());
      message("). Nunmap: %ld\n",Nunmap);
      if(args.flag("gibbs"))message("  Mean thetaAct (noise parameter)\n   %lf\n",totAverage[0].FF);
      message("\n");
      //}}}
      // Increase sample size and start over: {{{
      if(quitNext){// Sampling iterations end {{{
         if(sqrt(rHat2[0].FF.FF) > gPar.targetScaleReduction()){
            message("WARNING: Following transcripts failed to converge entirely\n   (however the estimates might still be usable):\n");
            stringstream sstr;
            sstr.str();
            sstr<<"# unconverged_transcripts: ";
            for(i=0;(i<M) && (sqrt(rHat2[i].FF.FF) > gPar.targetScaleReduction());i++){
               sstr<<rHat2[i].SS<<" ("<<sqrt(rHat2[i].FF.FF)<<") ";
               message("   %s( %ld , %lf )\n",(trInfo.trName(rHat2[i].SS-1)).c_str(),rHat2[i].SS-1,sqrt(rHat2[i].FF.FF));
            }
            sstr<<"\n";
            failedMessage=sstr.str();
         }
         for(j=0;j<chainsN;j++){
            samplers[j]->noSave();
            samplesFile[j].close();
         }
         break;
      }//}}}
      if(! args.flag("scaleReduction")){
         vector<double> needS(M,0);
         for(i=1;i<M;i++){
            // between variance was not multiplied by samplesHave===n
            // there is no chainsN in the denominator because samplesSave was already divided by chainsN
            needS[i] = samplesSave * samplesHave/
                     ((samplesHave-1.0)/samplesHave*withVar[i].FF/betwVar[i].FF+1.0);
         } 
         // log the number of effective samples, only when testing... //{{{
         #ifdef LOG_NEED
            stringstream sstr;
            sstr<<args.getS("outFilePrefix")<<".effLog";
            ofstream effLog(sstr.str().c_str());
            for(i=1;i<M;i++){
               effLog<<needS[rHat2[i].SS]<<" "<<sqrt(rHat2[i].FF.FF)<<" "<<samplesHave*betwVar[rHat2[i].SS].FF<<" "<<withVar[rHat2[i].SS].FF<<" "<<rHat2[i].SS<<endl;
            }
            effLog.close();
         #endif
         //}}}
         sort(needS.begin(),needS.end());
         i = (long)(M*0.95)+1; // make at least 95% transcripts converged 
         samplesN = max((long)needS[i],samplesSave);
         quitNext = true;
      }else{
            // Prepare for producing samples if Rhat^2<target scale reduction
            // OR reached samplesNmax
            // OR produced too many samples (>500 000)
         if((totalSamples < 5000000) && (rMean.FF > gPar.targetScaleReduction())){
            samplesN *= 2;
         }else{
            quitNext = true;
         }
         if((samplesN >= gPar.samplesNmax()) || args.flag("MCMC_samplesDOmax")){
            samplesN=gPar.samplesNmax();
            quitNext = true;
         }
      }
      // if next iteration is the last one, prepare the files and make samples write samples
      if(quitNext){ 
         message("Producing %ld final samples.\n",samplesN);
         // if samplesN<samplesSave, only samplesN samples will be saved
         if(samplesN<samplesSave)samplesSave = samplesN;
         stringstream sstr;
         for(j=0;j<chainsN;j++){
            sstr.str("");
            sstr<<args.getS("outFilePrefix")<<"."<<outTypeS<<"S-"<<j;
            samplesFileNames.push_back(sstr.str());
            samplesFile[j].open(samplesFileNames[j].c_str());
            if(! samplesFile[j].is_open()){
               error("Main: Unable to open output file '%s'.\n",(sstr.str()).c_str());
            }else{
               samplesFile[j]<<"#\n# M "<<M-1<<"\n# N "<<samplesSave<<endl;
               samplers[j]->saveSamples(&samplesFile[j],trInfo.getShiftedLengths(true),outTypeI);
            }
         }
      }
      for(j=0;j<chainsN;j++){
         samplers[j]->resetSampler(samplesN);
      }
      samplesHave=0;
      //}}}
   }
   // Write means: {{{
   meansFile.open((args.getS("outFilePrefix")+".thetaMeans").c_str());
   if(meansFile.is_open()){
      meansFile<<"# T => Mrows \n# M "<<M-1<<endl;
      meansFile<<"# file containing the mean value of theta - realtive abundace of fragments and counts\n"
                 "# (overall mean, overall counts, mean of saved samples, and mean from every chain are reported)\n"
                 "# columns:\n"
                 "# <transcriptID> <meanThetaOverall> <meanReadCountOverall> <meanThetaSaved>";
      for(j=0;j<chainsN;j++)meansFile<<" <chain"<<j+1<<"mean>";
      meansFile<<endl;
      meansFile<<scientific;
      meansFile.precision(9);
      double sum,sum2;
      for(i=0;i<M;i++){
         sum=sum2=0;
         for(j=0;j<chainsN;j++){
            sum+=samplers[j]->getAverage(i).FF; 
            sum2+=samplers[j]->getAverage(i).SS; 
         }
         if(i==0){
            meansFile<<"#thetaAct:";
         }else{
            meansFile<<i;
         }
         meansFile<<" "<<sum/chainsN<<" "<<(long)floor(sum/chainsN*alignments->getNreads()+0.5)<<" "<<sum2/chainsN;
         for(j=0;j<chainsN;j++)
            meansFile<<" "<<samplers[j]->getAverage(i).FF;
         meansFile<<endl;
      }
      meansFile.close();
   }else{
      warning("Main: Unable to write thetaMeans into: %s\n",(args.getS("outFilePrefix")+".thetaMeans").c_str());
   }
   //}}}
   // Write thetaAct: {{{
   if(args.isSet("thetaActFileName")){
      ofstream actFile(args.getS("thetaActFileName").c_str());
      if(actFile.is_open()){
         actFile<<"# samples of thetaAct parameter (only generated when using gibbs sampling)\n";
         actFile<<"# N "<<chainsN*samplesSave<<endl;
         for(j=0;j<chainsN;j++){
            for(i=0;i<(long)samplers[j]->getThetaActLog().size();i++)
               actFile<<samplers[j]->getThetaActLog()[i]<<" ";
         }
         actFile<<endl;
         actFile.close();
      }else{
         warning("Main: Unable to write thetaAct log: %s.\n",(args.getS("thetaActFileName")).c_str());
      }
   }
   // }}}
   // Free memory: {{{
   for(j=0;j<chainsN;j++){
      delete samplers[j];
   }
//   delete [] samplers;
   //}}}
   message("Total samples: %ld\n",totalSamples);
}//}}}

extern "C" int estimateExpression(int *argc, char* argv[]) {//{{{
clearDataEE();
string programDescription =
"Estimates expression given precomputed probabilities of (observed) reads' alignments.\n\
   Uses MCMC sampling algorithm to produce relative abundance or RPKM.\n";
   buildTime(argv[0],__DATE__,__TIME__);
   // Set options {{{
   ArgumentParser args;
   args.init(programDescription,"[prob file]",1);
   args.addOptionS("o","outPrefix","outFilePrefix",1,"Prefix for the output files.");
   args.addOptionS("O","outType","outputType",0,"Output type (theta, RPKM, counts, tau).","theta");
   args.addOptionB("G","gibbs","gibbs",0,"Use gibbs sampling instead of collapsed gibbs sampling.");
   args.addOptionS("p","parFile","parFileName",0,"File containing parameters for the sampler, which can be otherwise specified by --MCMC* options. As the file is checked after every MCMC iteration, the parameters can be adjusted while running.");
   args.addOptionS("t","trInfoFile","trInfoFileName",0,"File containing transcript information. (Necessary for RPKM)");
   args.addOptionL("P","procN","procN",0,"Limit the maximum number of threads to be used. (Default is the number of MCMC chains.)");
   args.addOptionS("","thetaActFile","thetaActFileName",0,"File for logging noise parameter theta^{act}.");
   args.addOptionL("","MCMC_burnIn","MCMC_burnIn",0,"Length of sampler's burn in period.",1000);
   args.addOptionL("","MCMC_samplesN","MCMC_samplesN",0,"Initial number of samples produced. Doubles after every iteration.",1000);
   args.addOptionL("","MCMC_samplesSave","MCMC_samplesSave",0,"Number of samples recorder in total.",1000);
   args.addOptionL("","MCMC_samplesNmax","MCMC_samplesNmax",0,"Maximum number of samples produced in one iteration. After producing samplesNmax samples sampler finishes.",50000);
   args.addOptionB("","MCMC_samplesDOmax","MCMC_samplesDOmax",0,"Produce maximum number of samples (samplesNmax) in second iteration and quit.");
   args.addOptionL("","MCMC_chainsN","MCMC_chainsN",0,"Number of parallel chains used. At least two chains will be used.",4);
   args.addOptionD("","MCMC_scaleReduction","MCMC_scaleReduction",0,"Target scale reduction, sampler finishes after this value is met.",1.2);
   args.addOptionD("","MCMC_dirAlpha","MCMC_dirAlpha",0,"Alpha parameter for the Dirichlet distribution.",1.0);
   args.addOptionB("","scaleReduction","scaleReduction",0,"Use scale reduction as stopping criterion, instead of computing effective sample size.");
   args.addOptionL("s","seed","seed",0,"Random initialization seed.");
   if(!args.parse(*argc,argv))return 0;
   // }}}
   MyTimer timer;
#ifdef SUPPORT_OPENMP
   if(args.isSet("procN"))
      omp_set_num_threads(args.getL("procN"));
   else
      omp_set_num_threads(args.getL("MCMC_chainsN"));
#endif

   gibbsParameters gPar;
   TagAlignments *alignments=NULL;
//{{{ Initialization:

   gPar.setParameters(args);
   if(args.isSet("parFileName")){
      gPar.setParameters(args.getS("parFileName"));
   }
   
   if((args.getS("outputType") == "theta")||(args.getS("outputType") == "THETA")){
      outTypeI=THETA;
      outTypeS="theta";
   }else if((args.getS("outputType") == "RPKM")||(args.getS("outputType") == "rpkm")){
      outTypeI=RPKM;
      outTypeS="rpkm";
   }else if(args.getS("outputType") == "tau"){
      outTypeI=TAU;
      outTypeS="tau";
   }else{
      outTypeI=COVERAGE;
      outTypeS="counts";
      if(args.getS("outputType") != "counts")
         message("Using output type \"counts\"\n");
   }

   if(args.verbose)gPar.getAllParameters();

   //}}}
   // {{{ Read transcriptInfo and .prob file 
   if(args.verbose)message("Reading data.\n");
   if((!args.isSet("trInfoFileName"))||(!trInfo.readInfo(args.getS("trInfoFileName")))){
      if(outTypeI==RPKM){
         error("Main: Missing transcript info file. This will cause problems if producing RPKM.");
         return 1;
      }
   }else{
      M = trInfo.getM()+1;
   }
   alignments = readData(args);
   if(! alignments){
      error("Main: Reading alignments failed.\n");
      return 1;
   }
   if(M<=0){
      error("Main: Invalid number of transcripts in .prob file.\n");
      return 1;
   }
   // }}}

   if(args.verbose)timer.split();
   if(args.verbose)message("Starting the sampler.\n");
   MCMC(alignments,gPar,args);
   // {{{ Transpose and merge sample file 
   if(transposeFiles(samplesFileNames,args.getS("outFilePrefix")+"."+outTypeS,args.verbose,failedMessage)){
      if(args.verbose)message("Sample files transposed. Deleting.\n");
      for(long i=0;i<(long)samplesFileNames.size();i++){
         remove(samplesFileNames[i].c_str());
      }
   }else{
      message("Transposing files failed. Please check the files and try using trasposeLargeFile to transpose & merge the files into single file.\n");
   }
   //}}}
   delete alignments;
   if(args.verbose){message("DONE. "); timer.split(0,'m');}
   return 0;
}//}}}

#ifndef BIOC_BUILD
int main(int argc, char* argv[]) {
   return estimateExpression(&argc,argv);
}
#endif
