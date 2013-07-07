/*
 *
 * Produce overall gene expression
 *
 */
#include<cmath>

using namespace std;

#include "ArgumentParser.h"
#include "common.h"
#include "misc.h"
#include "PosteriorSamples.h"
#include "TranscriptInfo.h"

extern "C" int getGeneExpression(int *argc,char* argv[]){
   string programDescription=
"Computes expression of whole genes.\n\
   [samplesFile] should contain transposed MCMC samples which will be transformed into gene expression samples.";   
   // Set options {{{
   ArgumentParser args(programDescription,"[samplesFile]",1);
   args.addOptionS("t","trInfoFile","trInfoFileName",1,"Name of the transcript file.");
   args.addOptionB("a","adjustByLength","adjust",0,"Adjust expression by transcripts length.");
   args.addOptionB("","theta2rpkm","rpkm",0,"Transform transcript expression in theta to gene expression in RPKM.");
   args.addOptionS("o","outFile","outFileName",1,"Name of the output file.");
   args.addOptionB("l","log","log",0,"Output logged values.");
   args.addOptionS("T","trMap","trMapFile",0,"Name of the file containing transcript to gene mapping.");
   args.addOptionS("G","geneList","geneListFile",0,"Name of the file containing list of gene names (one for each transcript).");
   args.addOptionB("","updateTrFile","updateTrFile",0,"Update trInfoFile if new gene names were provided (with trMapFile or geneListFile).");
   args.addOptionS("g","geneInfoFile","geneInfoFile",0,"Name of while to which gene information will be saved.");
   if(!args.parse(*argc,argv))return 0;
   if(args.verbose)buildTime(argv[0],__DATE__,__TIME__);
   // }}}
   bool doLog,doAdjust=args.flag("adjust")||args.flag("rpkm"),doRPKM=args.flag("rpkm");
   doLog = ns_genes::getLog(args);
   
   long N=0,M=0,G;
   TranscriptInfo trInfo;
   PosteriorSamples  samples;
   if(!ns_genes::prepareInput(args, &trInfo, &samples, &M, &N, &G))return 1;
   if(!ns_genes::updateGenes(args, &trInfo, &G))return 1;
   if(args.verb())messageF("Genes: %ld\n",G);
   if(!ns_genes::checkGeneCount(G,M))return 1;
   if(args.flag("updateTrFile") && (args.isSet("trMapFile") || args.isSet("geneListFile"))){
      if(args.verb())message("Updating transcript info file with new gene names.\n");
      if(!trInfo.writeInfo(args.getS("trInfoFileName"), true)){
         if(args.verb())warning("Main: Updating trInfoFile failed.\n");
      }
   }
   if(args.isSet("geneInfoFile")){
      if(args.verb())message("Saving gene information into: %s.\n",args.getS("geneInfoFile").c_str());
      if(!trInfo.writeGeneInfo(args.getS("geneInfoFile"))){
         warning("Main: Writing gene information failed.\n");
      }
   }

   ofstream outFile;
   if(!ns_misc::openOutput(args, &outFile))return 1;;
   // Write ouput header {{{
   outFile<<"# from: "<<args.args()[0]<<"\n# samples of gene expression\n";
   if(args.verbose)message("Genes will be ordered as they first appear in %s.\n",(args.getS("trInfoFileName")).c_str());
   outFile<<"# Genes will be ordered as they first appear in "<<args.getS("trInfoFileName")<<"\n";
   if(doRPKM)outFile<<"# data in RPKM\n";
   if(doLog)outFile<<"# L \n";
   outFile<<"# T (M rows,N cols)\n";
   outFile<<"# G = M "<<G<<"\n# N "<<N<<endl;
   // Set precision.
   outFile.precision(9);
   outFile<<scientific;
   // }}}
   vector< vector<double> > trs;
   vector<long double> normals(N,0);
   long double sum;
   long i,j,g,gM,m;
   if(doAdjust){
      vector<double> tr(M);
      if(args.verbose)message("Computing normalization constants, because of length adjustment.\n");
      for(j=0;j<M;j++){
         if(args.verbose)progressLog(j,M);
         samples.getTranscript(j,tr);
         for(i=0;i<N;i++)
            normals[i] += tr[i]/trInfo.L(j);
      }
   }
   if(args.verbose)message("Computing gene expression.\n");
   for(g=0;g<G;g++){
      if(args.verbose)progressLog(g,G);
      gM = trInfo.getGtrs(g).size();
      if((long)trs.size()<gM)trs.resize(gM);
      for(j=0;j<gM;j++){
         m = trInfo.getGtrs(g)[j];
         samples.getTranscript( m , trs[j]);
      }
      for(i=0;i<N;i++){
         sum = 0;
         for(j=0;j<gM;j++){
            if(doAdjust&&(normals[i]>0)){
               m = trInfo.getGtrs(g)[j];
               sum+=(trs[j][i] / trInfo.L(m)) / normals[i];
            }else{
               sum+=trs[j][i];
            }
         }
         if(doRPKM)sum=sum*10e9;
         if(doLog)sum=log(sum);
         outFile<<sum<<" ";
      }
      outFile<<endl;
   }
   outFile.close();
   if(args.verbose)message("DONE\n");
   return 0;
}

#ifndef BIOC_BUILD
int main(int argc,char* argv[]){
   return getGeneExpression(&argc,argv);
}
#endif
