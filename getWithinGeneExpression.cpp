/*
 *
 * Produce relative expression within gene
 *
 */
#include<cmath>

using namespace std;

#include "ArgumentParser.h"
#include "common.h"
#include "misc.h"
#include "PosteriorSamples.h"
#include "TranscriptInfo.h"

extern "C" int getWithinGeneExpression(int *argc,char* argv[]){
   string programDescription=
"Computes relative expression of transcripts within genes.\n\
   [samplesFile] should contain transposed MCMC expression samples.\n\
   program can produce means and variance and write them into [sumFile]\n\
   or individual MCMC samples which are written into [outFile].";   
   // Set options {{{
   ArgumentParser args(programDescription,"[samplesFile]",1);
   args.addOptionS("t","trInfoFile","trInfoFileName",1,"Name of the transcript file.");
   args.addOptionB("a","adjustByLength","adjust",0,"Adjust expression by transcripts length.");
   args.addOptionS("o","outFile","outFileName",0,"Name of the output file.");
   args.addOptionS("s","sumFile","sumFileName",0,"Name of summarization file where true mean, true variance and relative mean and relative variance are saved.");
   args.addOptionB("l","log","log",0,"Use logged values.");
//   args.addOptionS("t","type","type",0,"Type of variance, possible values: [sample,sqDif] for sample variance or sqared difference.","sample");
   if(!args.parse(*argc,argv))return 0;
   if(args.verbose)buildTime(argv[0],__DATE__,__TIME__);
   // }}}
   bool doLog,doOut=args.isSet("outFileName"),doSum=args.isSet("sumFileName");
   if(! (doOut || doSum)){
      error("Main: Have to secify at least one of --outFile/--sumFile.\n");
      return 1;
   }
   doLog = ns_genes::getLog(args);
   
   long N=0,M=0,G;
   TranscriptInfo trInfo;
   PosteriorSamples  samples;
   if(!ns_genes::prepareInput(args, &trInfo, &samples, &M, &N, &G))return 1;

   ofstream outFile,sumFile;
   if(doOut){
      if(!ns_misc::openOutput(args, &outFile))return 1;;
      // Write output header {{{
      outFile<<"# from: "<<args.args()[0]<<"\n# samples of within gene expression\n";
      if(! trInfo.genesOrdered()){
         warning("Main: Transcripts in output file will be reordered and grouped by genes.\n");
         outFile<<"# WARNING: transcripts in output file are reordered and grouped by genes.\n";
      }
      if(doLog)outFile<<"# L \n";
      outFile<<"# T (M rows,N cols)\n";
      outFile<<"# M "<<M<<"\n# N "<<N<<endl;
      // }}}
   }
   if(doSum){
      if(!ns_misc::openOutput(args.getS("sumFileName"), &sumFile))return 1;
      sumFile<<"# from: "<<args.args()[0]<<"\n# <mean> <variance> <mean of within gene expression>  <variance of within gene expression>\n# M "<<M<<endl;
   }
   vector<long double> mean(M),mean2(M),sqSum(M),sqSum2(M);
   vector< vector<double> > trs;
   vector<long double> normals(N,0);
   long double x,sum,var,var2,l;
   long i,j,g,gM,m;
   if(args.flag("adjust")&&(doSum)){
      vector<double> tr(M);
      if(args.verbose)message("Computing normalization constants, because of length adjustment.\n");
      for(j=0;j<M;j++){
         if(args.verbose)progressLog(j,M);
         samples.getTranscript(j,tr);
         for(i=0;i<N;i++)
            normals[i] += tr[i]/trInfo.L(j);
      }
   }
   if(args.verbose)message("Computing within gene relative expression.\n");
   for(g=0;g<G;g++){
      if(args.verbose)progressLog(g,G);
      gM = trInfo.getGtrs(g)->size();
      if((long)trs.size()<gM)trs.resize(gM);
      //message("%ld\n",gM);
      for(j=0;j<gM;j++){
         m = (*trInfo.getGtrs(g))[j];
         samples.getTranscript( m , trs[j]);
         mean[m] = mean2[m] = sqSum[m] = sqSum2[m] = 0;
      }
      for(i=0;i<N;i++){
         sum = 0;
         for(j=0;j<gM;j++){
            if(args.flag("adjust")){
               m = (*trInfo.getGtrs(g))[j];
               l = trInfo.L(m);
               trs[j][i]/=l;
               if(normals[i]!=0)trs[j][i]/=normals[i];
            }
            x=trs[j][i];
            sum += x;
            if(doSum){
               m = (*trInfo.getGtrs(g))[j];
               if(doLog)x=log(trs[j][i]);
               mean[m] += x;
               sqSum[m] += x*x;
            }
         }
         for(j=0;j<gM;j++){
            trs[j][i] /= sum;
            if(doLog)trs[j][i] = log(trs[j][i]);
            if(doSum){
               m = (*trInfo.getGtrs(g))[j];
               mean2[m] += trs[j][i];
               sqSum2[m] += trs[j][i]*trs[j][i];
            }
         }
      }
      if(doOut){
         for(j=0;j<gM;j++){
            for(i=0;i<N;i++)
               outFile<<trs[j][i]<<" ";
            outFile<<endl;
         }
      }
   }
   if(doOut)outFile.close();
   if(doSum){
      for(j=0;j<M;j++){
         mean[j] /= N;
         var = sqSum[j]/N - mean[j]*mean[j];
         mean2[j] /= N;
         var2 = sqSum2[j]/N - mean2[j]*mean2[j];
         sumFile<<mean[j]<<" "<<var<<" "<<mean2[j]<<" "<<var2<<endl;
      }
      sumFile.close();
   }
   if(args.verbose)message("DONE\n");
   return 0;
}


#ifndef BIOC_BUILD
int main(int argc,char* argv[]){
   return getWithinGeneExpression(&argc,argv);
}
#endif
