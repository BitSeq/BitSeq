#include<cmath>
#include<omp.h>

#include "common.h"
#include "misc.h"
#include "MyTimer.h"
#include "ReadDistribution.h"

void inline progressLogRD(long cur,long outOf) {//{{{
   // output progress status every 10%
   if((outOf>10)&&(cur%((long)(outOf/10))==0)&&(cur!=0))message("# %ld done.\n",cur);
}//}}}

#define DEBUG(x) 

inline char int2base(int B){//{{{
   switch(B){
      case 0: return 'a';
      case 1: return 'c';
      case 2: return 'g';
      case 3: return 't';
      default: return 'n';
   }
}//}}}
inline long base2int(char B){//{{{
   switch(B){
      case 'A': case 'a': return 0; 
      case 'C': case 'c': return 1; 
      case 'G': case 'g': return 2; 
      case 'T': case 't': return 3;
      default: return -1;
   }
}//}}}
inline long bamBase2int(int B){//{{{
   switch(B){
      case 1: return 0;
      case 2: return 1;
      case 4: return 2;
      case 8: return 3;
      default: return -1;
   }
}//}}}
inline void mapAdd(map<long,double > &m, long key, double val){//{{{
   if(m.count(key)==0)
      m[key] = val;
   else
      m[key] += val;
}//}}}

ReadDistribution::ReadDistribution(){ //{{{
   M=0;
   uniform = lengthSet = gotExpression = normalized = validLength = false;
   warnPos = warnTIDmismatch = warnUnknownTID = noteFirstMateDown = 0;
   lMu=100;
   lSigma=10;
   verbose = true;
   singleReadLength = 0;
   minFragLen=10000;
   lowProbMismatches = LOW_PROB_MISSES;
}//}}}
void ReadDistribution::writeWarnings() {//{{{
   if(warnPos>0){
      warning("ReadDistribution: %ld upstream reads from a pair did not align to the sense strand of transcript.\n", warnPos);
   }
   if(warnTIDmismatch>0){
      warning("ReadDistribution: %ld pair reads were aligned to different transcripts.\n", warnTIDmismatch);
   }
   if(warnUnknownTID>0){
      warning("ReadDistribution: %ld fragments were aligned to unknown transcripts.\n", warnUnknownTID);
   }
   if(noteFirstMateDown){
      message("NOTE: ReadDistribution: First mate from a pair was downstream (%ld times).\n", noteFirstMateDown);
   }
   warnPos = warnTIDmismatch = warnUnknownTID = noteFirstMateDown = 0;
}//}}}
bool ReadDistribution::init(long m, TranscriptInfo* trI, TranscriptSequence* trS, TranscriptExpression* trE, bool verb){ //{{{
   M = m;
   verbose = verb;
   if(trI==NULL){
      error("ReadDistribution: Missing TranscriptInfo.\n");
      return false;
   }
   if(trS==NULL){
      error("ReadDistribution: Missing TranscriptSequence.\n");
      return false;
   }
   uniform = false;
   trInf=trI;
   trSeq=trS;
   trExp=trE;
   if(trExp) gotExpression = true;
   else gotExpression = false;
   lengthSet = false;
   logLengthSum = logLengthSqSum = 0;
   fragSeen = 0;
   // Initialize tr - frag_length - expression maps:
   trFragSeen5.resize(M);
   trFragSeen3.resize(M);
   weightNorms.resize(3,vector<map<long, double> >(M));
   // Initialize position bias matrices:
   posProb.resize( 6, vector<vector<double> >(trSizesN + 1, vector<double>(trNumberOfBins,0.01/trNumberOfBins)));
   // Initialize sequence bias VLMMs: 
   seqProb.resize(4);
   for(long i=0;i<vlmmNodesN;i++){
      for(long j=0;j<4;j++)
         seqProb[j].push_back(VlmmNode(vlmmNodeDependence[i]));
   }
   return true;
}//}}}
bool ReadDistribution::initUniform(long m, TranscriptInfo* trI, TranscriptSequence* trS, bool verb){ //{{{
   M = m;
   verbose = verb;
   if(trI==NULL){
      error("ReadDistribution: Missing TranscriptInfo.\n");
      return false;
   }
   trInf = trI;
   trSeq = trS;
   trExp = NULL;
   uniform = true;
   lengthSet = false;
   gotExpression = false;
   logLengthSum = logLengthSqSum = 0;
   fragSeen = 0;
   return true;
}//}}}
void ReadDistribution::setLowProbMismatches(long m){//{{{
   lowProbMismatches = m>1 ? m:1;
}//}}}
void ReadDistribution::setLength(double mu, double sigma){ //{{{
   lMu=mu;
   lSigma=sigma;
   lengthSet=true;
}//}}}
void ReadDistribution::observed(fragmentP frag){ //{{{
   DEBUG(message("%s===%s\n",bam1_qname(frag->first),bam1_qname(frag->second));)
   long tid = frag->first->core.tid;
   if((frag->paired)&&(tid!=frag->second->core.tid)){
      warnTIDmismatch++;
      return;
   }
   if((tid < 0)||(tid>=M)){
      warnUnknownTID++;
      return;
   }
   // Set inverse expression
   double Iexp = (gotExpression)? 1.0/trExp->exp(tid) : 1.0;
   // Calculate reads' true end position:
   long frag_first_endPos, frag_second_endPos=0;
   frag_first_endPos = bam_calend(&frag->first->core, bam1_cigar(frag->first));
   if(frag->paired){
      frag_second_endPos = bam_calend(&frag->second->core, bam1_cigar(frag->second));
   }
   // update lengths: //{{{
   DEBUG(message("   length update\n");)
   double len,logLen;
   if(frag->paired){
      fragSeen ++;
      if(frag->second->core.pos>frag->first->core.pos)
         len = frag_second_endPos  - frag->first->core.pos;
      else{
         len = frag_first_endPos - frag->second->core.pos;
      }
      if(minFragLen>(long)len)minFragLen = (long) len;
      logLen = log(len);
      logLengthSum += logLen;
      logLengthSqSum += logLen*logLen;
   }else{
      len = frag_first_endPos - frag->first->core.pos;
      singleReadLength = (long)len;
      if(singleReadLength<minFragLen)minFragLen = singleReadLength;
   } //}}}
   // for uniform distribution ignore other estimation:
   if(uniform) return;

   // check mates relative position: {{{
   if((frag->paired) && (frag->first->core.pos > frag->second->core.pos)){
      noteFirstMateDown ++;
      bam1_t *tmp = frag->second;
      frag->second = frag->first;
      frag->first = tmp;
   }
   if((frag->paired) && (frag->first->core.flag & BAM_FREVERSE)){
      warnPos ++;
      return;
   }//}}}
   // positional bias:
   // sequence bias:
   DEBUG(message("   positional & sequence bias\n");)
   if(! frag->paired){
      if(frag->first->core.flag & BAM_FREVERSE){
         // Antisense strand of transcript is 5'end of fragment
         updatePosBias(frag_first_endPos, readM_5, tid, Iexp);
         // readM_5 and uniformM_5 are always "second mates" 
         // this is assumed also in getP(...);
         updateSeqBias(frag_first_endPos, readM_5, tid, Iexp);
         // update sum of expression of  fragments of given length
         mapAdd(trFragSeen5[tid], (long)len, Iexp);
      }else{
         // Sense strand of transcript is 3'end of fragment
         updatePosBias( frag->first->core.pos, readM_3, tid, Iexp);
         updateSeqBias( frag->first->core.pos, readM_3, tid, Iexp);
         mapAdd(trFragSeen3[tid], (long)len, Iexp);
      }
   }else{
      updatePosBias( frag->first->core.pos, readM_3, tid, Iexp);
      updateSeqBias( frag->first->core.pos, readM_3, tid, Iexp);
      mapAdd(trFragSeen3[tid], (long)len, Iexp);
         
      updatePosBias( frag_second_endPos, readM_5, tid, Iexp);
      updateSeqBias( frag_second_endPos, readM_5, tid, Iexp);
      mapAdd(trFragSeen5[tid], (long)len, Iexp);
   }
}//}}}
void ReadDistribution::normalize(){ //{{{
   // length distribution: {{{
   double newMu=0, newSigma=0;
  
   if(fragSeen>10){
      newMu = logLengthSum / fragSeen;
      newSigma = sqrt(logLengthSqSum / fragSeen - newMu*newMu);
      if(verbose)message("ReadDistribution: fragment length mu: %lg sigma: %lg\n",newMu,newSigma);
      validLength = true;
   }else{
      validLength = false;
   }
   if(lengthSet){
      // check difference between estimated mean and provided mean
      if(abs(newMu-lMu)>lSigma){
         warning("ReadDistribution: Estimated length mean (%lg) differs too much from the one provided (%lg).\n",newMu,lMu);
      }
   }else{
      lMu = newMu;
      lSigma = newSigma;
   }
   // }}}
   if(uniform) return;
   map<long,double>::iterator mIt;
   long i,j,m,group,trLen,fragLen;
   double Iexp,norm;
   // set Uniform position position bias: //{{{
   if(verbose)message("ReadDistribution: Computing uniform positional bias.\n");
   for(m=0;m<M;m++){
      if(verbose)progressLogRD(m,M);
      trLen = trInf->L(m);
      if(trLen<trNumberOfBins)continue;
      //message(" %ld %ld %ld\n",m,trLen,trFragSeen[m].size());
      for(group=0;group<trSizesN;group++)
         if(trLen<trSizes[group])break;
      // update 3' positional bias
      for( mIt=trFragSeen3[m].begin(); mIt != trFragSeen3[m].end(); mIt++){
         fragLen = mIt->first;
         Iexp = mIt->second / (trLen - fragLen + 1);
         for(i=0;i<trNumberOfBins;i++){
            // update probability of each bin by Iexp*"effective length of current bin"
            if((i+1)*trLen/trNumberOfBins < fragLen)continue;
            if(i*trLen/trNumberOfBins < fragLen){
               posProb[uniformM_3][group][trNumberOfBins -1 -i]+= Iexp * ((double)(i+1)*trLen / trNumberOfBins - fragLen + 1);
            }else{
               posProb[uniformM_3][group][trNumberOfBins -1 -i]+= Iexp * ((double)trLen / trNumberOfBins);
            }
         }
      }  
      // update 5' positional bias
      for( mIt=trFragSeen5[m].begin(); mIt != trFragSeen5[m].end(); mIt++){
         fragLen = mIt->first;
         Iexp = mIt->second / (trLen - fragLen + 1);
         for(i=0;i<trNumberOfBins;i++){
            // update probability of each bin by Iexp*"effective length of current bin"
            if((i+1)*trLen/trNumberOfBins < fragLen)continue;
            if(i*trLen/trNumberOfBins < fragLen){
               posProb[uniformM_5][group][i]+= Iexp * ((double)(i+1)*trLen / trNumberOfBins - fragLen + 1);
            }else{
               posProb[uniformM_5][group][i]+= Iexp * ((double)trLen / trNumberOfBins);
            }
         }
      }  
   }// }}}
   // pre-compute position bias weights: {{{
   for(j=0;j<4;j++)
      for(group=0;group<=trSizesN;group++){
         norm = 0;
         for(i=0;i<trNumberOfBins;i++)norm += posProb[j][group][i];
         for(i=0;i<trNumberOfBins;i++)posProb[j][group][i] /= norm;
      }
   for(group=0;group <= trSizesN;group++){
      for(i=0;i<trNumberOfBins;i++){
         // FIX HERE
         posProb[weight_5][group][i] = posProb[readM_5][group][i]/posProb[uniformM_5][group][i];
         // FIX HERE
         posProb[weight_3][group][i] = posProb[readM_3][group][i]/posProb[uniformM_3][group][i];
      }
   }//}}}
   //set Uniform sequence bias: {{{
   if(verbose)message("ReadDistribution: Computing uniform sequence bias.\n");
   double IexpSum5,IexpSum3;
   map<long,double>::reverse_iterator mItR;
   long p;
   for(m=0;m<M;m++){
      if(verbose)progressLogRD(m,M);
      trLen = trInf->L(m);
      IexpSum3=0;
      for(mIt=trFragSeen3[m].begin();mIt!= trFragSeen3[m].end();mIt++)
         IexpSum3+=mIt->second / (trLen - mIt->first + 1);
      IexpSum5=0;
      mItR=trFragSeen3[m].rbegin();
      mIt=trFragSeen5[m].begin();
      // STL map iterator IS sorted by key <=> length
      for(p=0;p<trLen;p++){
         while((mIt!=trFragSeen5[m].end())&&(mIt->first <= p+1)){IexpSum5+=mIt->second/ (trLen - mIt->first + 1); mIt++;}
         while((mItR!=trFragSeen3[m].rend())&&(trLen-p < mItR->first)){IexpSum3-= mItR->second / (trLen - mItR->first + 1) ; mItR++;}
         // 5' end is expected to be "after"
         updateSeqBias(p+1, uniformM_5, m, IexpSum5);
         updateSeqBias(p, uniformM_3, m, IexpSum3);
      }
   }//}}}
   // normalize VLMM nodes: {{{
   for(i=0;i<vlmmNodesN;i++){
      for(long j=0;j<4;j++)
         seqProb[j][i].normalize();
   }//}}} 
}//}}}
void ReadDistribution::logProfiles(string logFileName){//{{{
   ofstream outF;
   outF.open(logFileName.c_str());
   if(!outF.is_open()){
      error("ReadDistribution: Unable to open profile file: %s\n",(logFileName).c_str());
      return;
   }
   long i,j,g;
   outF<<"# BASES: (readM_5, readM_3, uniformM_5, uniformM_3)"<<endl;
   for(j=0;j<4;j++){
      outF<<"# "<<endl;
      for(i=0;i<vlmmNodesN;i++){
         outF<<seqProb[j][i].getPsum('A')<<" "<<seqProb[j][i].getPsum('C')<<" "<<seqProb[j][i].getPsum('G')<<" "<<seqProb[j][i].getPsum('T')<<endl;
      }
   }

   outF<<"#\n# Position: (readM_5, readM_3, uniformM_5, uniformM_3, weight_5, weight_3)"<<endl;
   for(j=0;j<6;j++){
      outF<<"# "<<endl;
      for(g=0;g<=trSizesN;g++){
         for(i=0;i<trNumberOfBins;i++)
            outF<<posProb[j][g][i]<<" ";
         outF<<endl;
      }
   }
   outF.close();
}//}}}
pair<double,double> ReadDistribution::getSequenceLProb(bam1_t *samA){//{{{
   if(! samA) return pair<double, double>(0,0);
   double lProb=0,lowLProb=0,lProbMis;
   bam1_core_t *samC = &samA->core;
   uint8_t *qualP=bam1_qual(samA);
   long i,j,misses,len=samC->l_qseq;
   long deletionN=0;
   // Count number of deletions-insertions
   for(i=0;i<samC->n_cigar;i++){
      switch(bam1_cigar(samA)[i]&BAM_CIGAR_MASK){
         case BAM_CDEL:
            deletionN += (long)(bam1_cigar(samA)[i]>>BAM_CIGAR_SHIFT);
            break;
         case BAM_CINS:
            deletionN -= (long)(bam1_cigar(samA)[i]>>BAM_CIGAR_SHIFT);
            break;
      }
   }
   string seq = trSeq->getSeq(samC->tid, samC->pos, len+deletionN, false);
   misses=lowProbMismatches;
   long cigarOp,cigarI,cigarOpCount;
   cigarOp=cigarI=cigarOpCount=0;
   // i - iterates within reference sequence, j - iterates within read
   i=j=0;
   while((i<len+deletionN) && (j<len)){
      if(cigarOpCount == 0){
         if(cigarI >= samC->n_cigar) break;
         cigarOp = bam1_cigar(samA)[cigarI]&BAM_CIGAR_MASK;
         cigarOpCount = (long)(bam1_cigar(samA)[cigarI]>>BAM_CIGAR_SHIFT);
         cigarI++;
      }
      cigarOpCount --;
      switch(cigarOp){
         case BAM_CDEL:
            i++; break;
         case BAM_CINS:
            j++; break;
         case BAM_CMATCH:
         case BAM_CEQUAL:
         case BAM_CDIFF:
            if((base2int(seq[i]) == -1)||
               (base2int(seq[i]) != bamBase2int(bam1_seqi(bam1_seq(samA),j))))misses--;
            i++;
            j++;
      }
   }
   if(misses<=0)misses=1;
   // start from the end so the "mismatched" bases for lowProb are at the end
   i=len+deletionN-1;
   j=len-1;
   cigarI = samC->n_cigar-1;
   while((i>=0) && (j>=0)){
      if(cigarOpCount == 0){
         if(cigarI < 0) break;
         cigarOp = bam1_cigar(samA)[cigarI]&BAM_CIGAR_MASK;
         cigarOpCount = (long)(bam1_cigar(samA)[cigarI]>>BAM_CIGAR_SHIFT);
         cigarI--;
      }
      cigarOpCount --;
      switch(cigarOp){
         case BAM_CDEL:
            i--; continue;
         case BAM_CINS:
            j--; continue;
         /*case BAM_CMATCH:
         case BAM_CEQUAL:
         case BAM_CDIFF:*/
      }
      lProbMis = (((double) qualP[j])/-10.0) * log(10.0);
         
      if((base2int(seq[i]) == -1)||(base2int(seq[i]) != bamBase2int(bam1_seqi(bam1_seq(samA),j)))){
         // If bases don't match, multiply probability by probability of error.
         lProb += lProbMis;
         lowLProb += lProbMis;
      }else{
         // If bases do match, multiple probability by inverse of probability of error.
         double lProbHit = log1p(-exp(lProbMis));
         lProb += lProbHit;
         if(misses>0){
            // If there are some misses left add a 'miss' to the 'low probability'.
            lowLProb += lProbMis;
            misses--;
         }else{
            lowLProb += lProbHit;
         }
      }
      i--;
      j--;
   }
   return pair<double, double>(lProb,lowLProb);
}//}}}
bool ReadDistribution::getP(fragmentP frag,double &lProb,double &lProbNoise){ //{{{
   double lP = 0;
   lProb = ns_misc::LOG_ZERO;
   lProbNoise = ns_misc::LOG_ZERO;
   pair<double, double> lpSeq1,lpSeq2;
   // Get probability based on base mismatches: {{{
   lpSeq1 = getSequenceLProb(frag->first);
   lpSeq2 = getSequenceLProb(frag->second);
   long tid = frag->first->core.tid;
   if((frag->paired)&&(tid!=frag->second->core.tid)){
      warnTIDmismatch++;
      return false;
   }
   if((tid < 0)||(tid>=M)){
      warnUnknownTID++;
      return false;
   }
   // Calculate reads' true end position:
   long frag_first_endPos, frag_second_endPos=0;
   frag_first_endPos = bam_calend(&frag->first->core, bam1_cigar(frag->first));
   if(frag->paired){
      frag_second_endPos = bam_calend(&frag->second->core, bam1_cigar(frag->second));
   }
   double trLen = trInf->L(tid),len;
   // }}}
   if(frag->paired){
   // Get probability of length {{{
      if(frag->second->core.pos>frag->first->core.pos)
         len = frag_second_endPos - frag->first->core.pos;
      else{
         len = frag_first_endPos - frag->second->core.pos;
      }
      // compute length probability and normalize by probability of all possible lengths (cdf):
      // P*=lengthP/lengthNorm
      if(validLength) lP += getLengthLP(len) - getLengthLNorm(trLen);
      // }}}
   }else{
      len = frag_first_endPos - frag->first->core.pos;
   }
   if(uniform){
      // Get probability of position for uniform distribution
      // P*=1/(trLen-len+1)
      lP -= log(trLen-len+1);
   }else{ // Positional & Sequence bias {{{
      // Get probability of position given read bias model
      // check mates' relative position:
      if( frag->paired && (frag->first->core.pos > frag->second->core.pos)){
         noteFirstMateDown ++;
         bam1_t *tmp = frag->second;
         frag->second = frag->first;
         frag->first = tmp;
      }
      // check strand of the first read:
      if(frag->paired && (frag->first->core.flag & BAM_FREVERSE)){
         warnPos++;
         return false;
      }
      if(!frag->paired){
         if(frag->first->core.flag & BAM_FREVERSE){
            // P*=posBias5'*seqBias5'/weightNorm5'
            lP += log(getPosBias(frag_first_endPos, mate_5, tid)) +
               log(getSeqBias(frag_first_endPos, mate_5, tid )) -
               log(getWeightNorm( (long) len, mate_5, tid));
         }else{
            // P*=posBias3'*seqBias3'/weightNorm3'
            lP += log(getPosBias(frag->first->core.pos , mate_3, tid)) +
               log(getSeqBias(frag->first->core.pos , mate_3, tid )) -
               log(getWeightNorm( (long) len, mate_3, tid));
         }
      }else{
//#pragma omp parallel sections num_threads (2) reduction(*:P)
//{
//   #pragma omp section
         // P*=1/weightNormFull
         lP -= log(getWeightNorm( (long) len, FullPair, tid));
//   #pragma omp section
//   {
         // P*=posBias5'*posBias3'*seqBias5'*seqBias3'
         lP += log(getPosBias(frag_second_endPos, mate_5, tid))
          + log(getPosBias(frag->first->core.pos , mate_3, tid))
          + log(getSeqBias(frag_second_endPos, mate_5, tid ))
          + log(getSeqBias(frag->first->core.pos , mate_3, tid )); 
//   }
//}
      }
   } //}}}
   lProb = lP + lpSeq1.first+lpSeq2.first;
   lProbNoise = lP + lpSeq1.second+lpSeq2.second;
   return true;
}//}}}
void ReadDistribution::updatePosBias(long pos, biasT bias, long tid, double Iexp){ //{{{
   if((bias==readM_5)||(bias==uniformM_5))pos--;
   long group, rel, trLen;
   trLen = trInf->L(tid);
   // transcript too short:
   if(trLen < trNumberOfBins) return;
   // choose group:
   for(group = 0;group < trSizesN;group++)
      if(trLen<trSizes[group])break;
   // find relative position:
   rel = (pos * trNumberOfBins) / trLen;
   if(rel>=trNumberOfBins)rel=trNumberOfBins-1;
   //add inverse expression:
   posProb[bias][ group ][ rel ] += Iexp;
}//}}}
void ReadDistribution::updateSeqBias(long pos, biasT bias, long tid, double Iexp){ //{{{
   if(Iexp<=0)return;
   if(bias>3)return; //this should not happen
   long start ;
   string seq;

   if((bias == uniformM_3)||(bias == readM_3)){
      start = pos - vlmmStartOffset - MAX_NODE_PAR;
      seq = trSeq->getSeq(tid, start, vlmmNodesN + MAX_NODE_PAR);
   }else{ // else get reverse complement
      start = pos + vlmmStartOffset - vlmmNodesN ;
      seq = trSeq->getSeq(tid, start, vlmmNodesN + MAX_NODE_PAR, true);
   }
   for(long i=0;i<vlmmNodesN;i++){
      seqProb[bias][i].update( Iexp, seq[i+2], seq[i+1], seq[i]);
   }
}//}}}
double ReadDistribution::getPosBias(long pos, readT read, long tid){ //{{{
   long group, rel, trLen;
   if(read == mate_5)pos--;
   trLen = trInf->L(tid);
   // transcript too short:
   if(trLen < trNumberOfBins) return 1;
   // choose group:
   for(group = 0;group < trSizesN;group++)
      if(trLen<trSizes[group])break;
   // find relative position:
   rel = (pos * trNumberOfBins) / trLen;
   if(rel>=trNumberOfBins)rel=trNumberOfBins-1;
   // return bias weight
   if(read == mate_5)
      return posProb[ weight_5 ][ group ][ rel ];
   if(read == mate_3)
      return posProb[ weight_3 ][ group ][ rel ];
   // shouldn't happen
   return 0;
}//}}}
double ReadDistribution::getSeqBias(long pos, readT read, long tid){ //{{{
   if(read==FullPair)return 0; // this should never happen
   string seq;
   long start;
   double B = 1;
   // Assuming the 5' mate is the second mate (which determines the start of interesting sequence)
   if(read == mate_3) start = pos - vlmmStartOffset - MAX_NODE_PAR;
   else start = pos + vlmmStartOffset - vlmmNodesN;

   if(read == mate_3){
#pragma omp critical
{
      seq = trSeq->getSeq(tid, start, vlmmNodesN + MAX_NODE_PAR);
}
      for(long i=0;i<vlmmNodesN;i++)
         // FIX HERE
         B *= seqProb[readM_3][i].getP( seq[i+2], seq[i+1], seq[i]) /
              seqProb[uniformM_3][i].getP( seq[i+2], seq[i+1], seq[i]);
   }else{ // else get reverse complement
#pragma omp critical
{
      seq = trSeq->getSeq(tid, start, vlmmNodesN + MAX_NODE_PAR, true);
}
      for(long i=0;i<vlmmNodesN;i++)
         // FIX HERE
         B *= seqProb[readM_5][i].getP( seq[i+2], seq[i+1], seq[i]) /
              seqProb[uniformM_5][i].getP( seq[i+2], seq[i+1], seq[i]);
   }
   return B;
}//}}}
double ReadDistribution::getWeightNorm(long len, readT read, long tid){ //{{{
   if(len == 0)return 1;
   if(weightNorms[read][tid].count(len) == 0){
      long trLen = trInf->L(tid),pos;
      double norm = 0,w;
#pragma omp parallel for private(w) reduction(+:norm)
      for(pos = 0;pos <= trLen-len;pos++){
         w=1.0;
         if((read == FullPair)||(read == mate_3)){
            w*=getPosBias(pos, mate_3, tid)*getSeqBias(pos,mate_3,tid);
         }
         if((read == FullPair)||(read == mate_5)){
            w*=getPosBias(pos + len, mate_5, tid)*getSeqBias(pos + len, mate_5, tid);
         }
         norm+=w;
      }
      weightNorms[read][tid][len] = norm;
//      message("w: %ld %ld %ld  %ld%lf\n",read,tid,len,trLen<"   ",norm);
      return norm;
   }
   return weightNorms[read][tid][len];
}//}}}
double ReadDistribution::getLengthLP(double len){//{{{
   //return 1./(len*lSigma*sqrt_2_pi)*exp(-pow(log(len) - lMu, (double)2.0)/(2 * pow(lSigma, (double)2)));
   const double log_sqrt_2_pi = .918938533192; // log(sqrt(2*pi))
   const double lLen = log(len);
   return - (lLen + 
             log(lSigma) + 
             log_sqrt_2_pi + 
             pow( (lLen - lMu) / lSigma, 2.0) / 2.0 );
}//}}}
double ReadDistribution::getLengthLNorm(double trLen){//{{{
   // erfc needs compiler with C99 standard 
   // other option might be to use boost/math/special_functions/erf.hpp
   return log(0.5)+log(erfc(-(log(trLen)-lMu)/(lSigma*1.41421356237309)));
}//}}}
vector<double> ReadDistribution::getEffectiveLengths(){ //{{{
   vector<double> effL(M,0);
   long m,len,trLen,pos;
   double eL, lCdfNorm,lenP, wNorm;
   MyTimer timer;
   timer.start();
#pragma omp parallel for private(len,trLen,pos,eL,lenP,wNorm,lCdfNorm)
   for(m=0;m<M;m++){
      if(verbose && (m!=0) && (m%(M/10)==0)){
#pragma omp critical
         {
            message("# %ld done. ",m);
            timer.current();
         }
      }
      trLen = trInf->L(m);
      if(!validLength){
         if(trLen>singleReadLength*2) effL[m] = trLen - singleReadLength; 
         else if(trLen>singleReadLength) effL[m] = singleReadLength;
         else effL[m] = trLen;
         continue;
      }
      lCdfNorm = getLengthLNorm(trLen);
// always computing the effective length using fragLen only
      if(uniform){
         eL = 0;
         for(len=1;len<=trLen;len++){
            eL += exp(getLengthLP(len)-lCdfNorm) * (trLen-len);
         }
         // dont go below minimal fragment length
         effL[m] = eL>minFragLen?eL:trLen;
      }else{
         vector<double> posBias5(trLen),posBias3(trLen);
         for(pos = 0;pos<trLen;pos++){
            posBias5[pos] = getPosBias(pos+1, mate_5, m)*getSeqBias(pos+1, mate_5, m);
            posBias3[pos] = getPosBias(pos, mate_3, m)*getSeqBias(pos, mate_3, m);
            //if(m==0)message(" %ld %lf %lf\n",pos,posBias5[pos],posBias3[pos]);
         }
         eL=0;
         for(len=1;len<=trLen;len++){
            wNorm = 0;
            for(pos=0;pos <= trLen - len;pos++){
               wNorm += posBias3[pos] * posBias5[pos+len-1];
            }
            lenP = exp(getLengthLP( len ) - lCdfNorm);
            //if(m==0)message("   %ld  %lf   %lf\n",len,lenP,wNorm);
            eL += lenP * wNorm;
         }
         // dont go below minimal fragment length
         effL[m] = eL>minFragLen?eL:trLen;
      }
   }
   if(! uniform){
      // normalize effective length to same sum as original length
      double effSum=0,lSum=0;
      for(m=0;m<M;m++){
         lSum+=trInf->L(m);
         effSum+=effL[m];
      }
      for(m=0;m<M;m++)effL[m] *= lSum/effSum;
   }
   for(m=0;m<M;m++)if(effL[m]<=0) effL[m]=1;
   return effL;
}//}}}

double VlmmNode::getPsum(char b) {//{{{
   if(base2int(b) == -1) return 1/4;
   if(parentsN == 2)return getP(b,'N','N')*16;
   if(parentsN == 1)return getP(b,'N','N')*4;
   return probs[base2int(b)];
}//}}}
VlmmNode::VlmmNode(long p) {//{{{
   setParentsN(p);
}//}}}
void VlmmNode::setParentsN(long p) {//{{{
   parentsN = p;
   if(parentsN>2){
      warning("VlmmNode: Code not read for using more than 2 parents.\n");
      parentsN = 2;
   }
   // initialize probability matrix, set pseudocount:
   probs.assign(pows4[parentsN+1], 0.01/pows4[parentsN+1]);
}//}}}
void VlmmNode::update(double Iexp, char b, char bp, char bpp) {//{{{
   double expDiv = 1.0;
   if(base2int(b) == -1)expDiv *=4.0;
   if((parentsN>0)&&(base2int(bp) == -1))expDiv *=4.0;
   if((parentsN>1)&&(base2int(bpp) == -1))expDiv *=4.0;
   if(expDiv == 1){
      // All bases are known:
      long i=0;
      switch(parentsN){
         case 2: 
            i += pows4[2]*base2int(bpp);
         case 1: 
            i += pows4[1]*base2int(bp);
         default: 
            i += base2int(b);
      }
      probs[ i ] += Iexp;
   }else{
      long i=0,j=0,k=0;
      Iexp /= expDiv;
      if(parentsN==2){
         for(i=0;i<4;i++)
            if((base2int(bpp) == i) || (base2int(bpp) == -1))
               for(j=0;j<4;j++)
                  if((base2int(bp) == j) || (base2int(bp) == -1))
                     for(k=0;k<4;k++)
                        if((base2int(b) == k) || (base2int(b) == -1))
                           probs[pows4[2]*i + pows4[1]*j+ k]+=Iexp;
      }else if(parentsN==1){
         for(j=0;j<4;j++)
            if((base2int(bp) == j) || (base2int(bp) == -1))
               for(k=0;k<4;k++)
                  if((base2int(b) == k) || (base2int(b) == -1))
                     probs[pows4[1]*j+ k]+=Iexp;
      }else{
         for(k=0;k<4;k++)
            // if((base2int(b) == k) || (base2int(b) == -1)); WE KNOW THAT b == 'N'
               probs[k]+=Iexp;
      }
   }
}//}}}
void VlmmNode::normalize() {//{{{
   double sum=0;
   long i,j,k,index;
   if(parentsN == 2){
      for(k=0;k<4;k++)
         for(j=0;j<4;j++){
            index = pows4[2]*k + pows4[1]*j;
            sum = 0;
            for(i=0;i<4;i++)sum += probs[i + index];
            for(i=0;i<4;i++)probs[i + index] /= sum;
         }
   }else if(parentsN == 1){
      for(j=0;j<4;j++){
         index = pows4[1]*j;
         sum = 0;
         for(i=0;i<4;i++)sum += probs[i + index];
         for(i=0;i<4;i++)probs[i + index] /= sum;
      }
   }else{
      sum = 0;
      for(i=0;i<pows4[parentsN+1];i++)sum += probs[i];
      for(i=0;i<pows4[parentsN+1];i++)probs[i] /= sum;
   }
}//}}}
double VlmmNode::getP(char b, char bp, char bpp) {//{{{
   if(base2int(b) == -1)return 1.0/4.0;
   double probDiv = 1.0;
   if((parentsN>0)&&(base2int(bp) == -1))probDiv *=4.0;
   if((parentsN>1)&&(base2int(bpp) == -1))probDiv *=4.0;
   if(probDiv == 1.0){
      // All bases are known:
      long i=0;
      switch(parentsN){
         case 2: 
            i += pows4[2]*base2int(bpp);
         case 1: 
            i += pows4[1]*base2int(bp);
         default: 
            i += base2int(b);
      }
      return probs[ i ];
   }else{
      long i=0,j=0,k=0;
      double prob = 0;
      // either one ore both parents are unknown==undefined
      if(parentsN==2){
         k = base2int(b);
         for(i=0;i<4;i++)
            if((base2int(bpp) == i) || (base2int(bpp) == -1))
               for(j=0;j<4;j++)
                  if((base2int(bp) == j) || (base2int(bp) == -1))
                     prob += probs[pows4[2]*i + pows4[1]*j+ k];
      }else if(parentsN==1){
         // there was an unknown => we know that parent is unknown
         k = base2int(b);
         for(j=0;j<4;j++)
            prob += probs[pows4[1]*j+ k];
      }else ;// Covered by all bases unknown;
      return prob / probDiv;
   }
}//}}}

