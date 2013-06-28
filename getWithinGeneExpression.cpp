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

namespace ns_withinGene {

// Read transcripts of gene g.
void readTranscripts(long g, const TranscriptInfo &trInfo, PosteriorSamples *samples, long *gM, vector< vector<double> > *trs);

// Adjust expression samples by transcript length.
void adjustExpression(long g, const TranscriptInfo &trInfo, vector< vector<double> > *trs);

// Compute sum of samples of transcripts from one gene.
void getSum(long gM, long N, const vector< vector<double> > &trs, vector<double> *sum);

// Update 'mean' and squareSum with new value.
void updateSummaries(double x, long double *mean, long double *sqSum, double norm = 1);

// Append samples of a transcript into output file.
void writeTr(long N, const vector<double> &tr, ofstream *outFile);

} // namespace ns_withinGene

extern "C" int getWithinGeneExpression(int *argc,char* argv[]){
   string programDescription=
"Computes relative expression of transcripts within genes.\n\
   [samplesFile] should contain transposed MCMC expression samples.\n\
   program can produce means and variance and write them into [sumFile]\n\
   or individual MCMC samples which are written into [outFile].";   
   // Set options {{{
   ArgumentParser args(programDescription,"[samplesFile] -t [trInfoFileName]",1);
   args.addOptionS("t","trInfoFile","trInfoFileName",1,"Name of the transcript file.");
   args.addOptionB("a","adjustByLength","adjust",0,"Adjust expression by transcripts length.");
   args.addOptionS("o","outFile","outFileName",0,"Name of the output file.");
   args.addOptionS("s","sumFile","sumFileName",0,"Name of summarization file where true mean, true variance and relative mean and relative variance are saved.");
   args.addOptionB("l","log","log",0,"Use logged values.");
   args.addOptionS("T","trMap","trMapFile",0,"Name of the file containing transcript to gene mapping.");
   args.addOptionB("","groupByGene","groupByGene",0,"Group transcripts by genes (this can change the default order of output.");
   if(!args.parse(*argc,argv))return 0;
   if(args.verbose)buildTime(argv[0],__DATE__,__TIME__);
   // }}}
   bool doLog,doOut=args.isSet("outFileName"),doSummaries=args.isSet("sumFileName");
   if(! (doOut || doSummaries)){
      error("Main: Have to specify at least one of --outFile/--sumFile.\n");
      return 1;
   }
   doLog = ns_genes::getLog(args);
   
   long N=0,M=0,G;
   TranscriptInfo trInfo;
   PosteriorSamples  samples;
   if(!ns_genes::prepareInput(args, &trInfo, &samples, &M, &N, &G))return 1;
   if(args.isSet("trMapFile") && (!ns_genes::updateGenes(args.getS("trMapFile"), &trInfo, &G)))return 1;
   if(args.verb())messageF("Genes: %ld\n",G);
   if(!ns_genes::checkGeneCount(G,M))return 1;
   if(args.isSet("trMapFile")){
      if(args.verb())message("Updating transcript info file with gene names provided in trMapFile.\n");
      if(!trInfo.writeInfo(args.getS("trInfoFileName"), true)){
         if(args.verb())warning("Main: Updating trInfoFile failed.\n");
      }
   }

   ofstream outFile,sumFile;
   if(doOut){
      if(!ns_misc::openOutput(args, &outFile))return 1;;
      // Write output header {{{
      outFile<<"# from: "<<args.args()[0]<<"\n# samples of within gene expression\n";
      if(! trInfo.genesOrdered()){
         if(args.flag("groupByGene")){
            warning("Main: Transcripts in output file will be reordered and grouped by genes.\n");
            outFile<<"# WARNING: transcripts in output file are reordered and grouped by genes.\n";
         }else{
            warning("Main: Transcripts are not grouped by genes.\n"
                    "   The transcript order will be kept the same but computation will be slower.\n");
         }
      }
      if(doLog)outFile<<"# L \n";
      outFile<<"# T (M rows,N cols)\n";
      outFile<<"# M "<<M<<"\n# N "<<N<<endl;
      // }}}
   }
   if(doSummaries){
      if(!ns_misc::openOutput(args.getS("sumFileName"), &sumFile))return 1;
      sumFile<<"# from: "<<args.args()[0]<<"\n# <mean> <variance> <mean of within gene expression>  <variance of within gene expression>\n# M "<<M<<endl;
   }
   vector<long double> mean(M,0),mean2(M,0),sqSum(M,0),sqSum2(M,0);
   vector< vector<double> > trs;
   vector<double> sum;
   vector<double> normals(N,1);
   long i,j,g,gM,m;
   if(args.flag("adjust")&&(doSummaries)){
      // 'normals' are only precomputed so that non-relative mean and variance are computed from
      // length adjusted and normalised expression.
      vector<double> tr(M);
      if(args.verbose)message("Computing normalization constants, because of length adjustment.\n");
      normals.assign(N,0);
      for(j=0;j<M;j++){
         if(args.verbose)progressLog(j,M);
         samples.getTranscript(j,tr);
         for(i=0;i<N;i++)
            normals[i] += tr[i]/trInfo.L(j);
      }
   }
   if(args.verbose)message("Computing within gene relative expression.\n");
   g = -2;
   if(!args.flag("groupByGene")){
      long curJ=0;
      // Here we iterate over transcripts:
      //  For each transript: load all transcripts of a gene of current transcripts
      //  If gene is same as for previous, then just reuse information
      for(m=0;m<M;m++){
         if(args.verbose)progressLog(m,M);
         if(trInfo.geId(m) == g){
            for(j=0;j<gM;j++)if(trInfo.getGtrs(g)[j] == m){curJ = j; break;}
         }else{
            g = trInfo.geId(m);
            ns_withinGene::readTranscripts(g, trInfo, &samples, &gM, &trs);
            curJ = 0;
            for(j=0;j<gM;j++)if(trInfo.getGtrs(g)[j] == m){curJ = j; break;}
            if(args.flag("adjust"))ns_withinGene::adjustExpression(g, trInfo, &trs);
            ns_withinGene::getSum(gM, N, trs, &sum);
         }
         for(i=0;i<N;i++){
            if(doLog)trs[curJ][i] = log(trs[curJ][i]);
            if(doSummaries) ns_withinGene::updateSummaries(trs[curJ][i], &mean[m], &sqSum[m], normals[i]);
            if(doLog)trs[curJ][i] -= log(sum[i]);
            else trs[curJ][i] /= sum[i];
            if(doSummaries) ns_withinGene::updateSummaries(trs[curJ][i], &mean2[m], &sqSum2[m]);
         }
         if(doOut){
            ns_withinGene::writeTr(N, trs[curJ], &outFile);
         }
      }
   }else{
      // Here we iterate over genes:
      //  Calculate values for all their transcripts
      //  Write all transcripts of current gene
      for(g=0;g<G;g++){
         if(args.verbose)progressLog(g,G);
         ns_withinGene::readTranscripts(g, trInfo, &samples, &gM, &trs);
         if(args.flag("adjust"))ns_withinGene::adjustExpression(g, trInfo, &trs);
         ns_withinGene::getSum(gM, N, trs, &sum);
         for(i=0;i<N;i++){
            for(j=0;j<gM;j++){
               m = trInfo.getGtrs(g)[j];
               if(doLog)trs[j][i] = log(trs[j][i]);
               if(doSummaries) ns_withinGene::updateSummaries(trs[j][i], &mean[m], &sqSum[m], normals[i]);
               if(doLog)trs[j][i] -= log(sum[i]);
               else trs[j][i] /= sum[i];
               if(doSummaries) ns_withinGene::updateSummaries(trs[j][i], &mean2[m], &sqSum2[m]);
            }
         }
         if(doOut){
            for(j=0;j<gM;j++){
               ns_withinGene::writeTr(N, trs[j], &outFile);
            }
         }
      }
   }
   if(doOut)outFile.close();
   if(doSummaries){
      long double var,var2;
      for(m=0;m<M;m++){
         mean[m] /= N;
         var = sqSum[m]/N - mean[m]*mean[m];
         mean2[m] /= N;
         var2 = sqSum2[m]/N - mean2[m]*mean2[m];
         sumFile<<mean[m]<<" "<<var<<" "<<mean2[m]<<" "<<var2<<endl;
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


namespace ns_withinGene {

void readTranscripts(long g, const TranscriptInfo &trInfo, PosteriorSamples *samples, long *gM, vector< vector<double> > *trs){//{{{
   *gM = trInfo.getGtrs(g).size();
   if((long)trs->size() < *gM)trs->resize(*gM);
   for(long j = 0; j < *gM; j++){
      samples->getTranscript( trInfo.getGtrs(g)[j] , (*trs)[j]);
   }
}// }}}

void adjustExpression(long g, const TranscriptInfo &trInfo, vector< vector<double> > *trs){//{{{
   long N,gM = trInfo.getGtrs(g).size();
   double l;
   for(long j=0; j<gM; j++){
      l = trInfo.L(trInfo.getGtrs(g)[j]);
      N = (*trs)[j].size();
      for(long n=0; n<N; n++){
         (*trs)[j][n] /= l;
      }
   }
}// }}}

void getSum(long gM, long N, const vector< vector<double> > &trs, vector<double> *sum){//{{{
   sum->assign(N,0);
   for(long j=0; j<gM; j++)
      for(long n=0; n<N; n++)(*sum)[n] += trs[j][n];
}// }}}

void updateSummaries(double x, long double *mean, long double *sqSum, double norm){//{{{
   if(x>0){// division is safe
      x /= norm;
      *mean += x;
      *sqSum += x*x;
   }
}// }}}

void writeTr(long N, const vector<double> &tr, ofstream *outFile){//{{{
   for(long n=0; n<N-1; n++)
      (*outFile)<<tr[n]<<" ";
   (*outFile)<<tr[N-1]<<endl;
}// }}}

} // namespace ns_withinGene

