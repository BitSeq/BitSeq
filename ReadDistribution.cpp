#include<algorithm>
#include<cmath>
#ifdef _OPENMP
#include<omp.h>
#endif

#include "ReadDistribution.h"

#include "misc.h"
#include "MyTimer.h"

#include "common.h"

#define DEBUG(x) 

namespace ns_rD {
// Base 2 Int mapping.
vector<char> tableB2I;
vector<int> tableB2BI;

/*void inline progressLogRD(long cur,long outOf) {//{{{
   // output progress status every 10%
   if((outOf>10)&&(cur%((long)(outOf/10))==0)&&(cur!=0))message("# %ld done.\n",cur);
}//}}} */
void fillTable() {//{{{
   tableB2I.assign(256,-1);
   tableB2I['A'] = tableB2I['a'] = 0;
   tableB2I['C'] = tableB2I['c'] = 1;
   tableB2I['G'] = tableB2I['g'] = 2;
   tableB2I['T'] = tableB2I['t'] = 3;
   tableB2BI.assign(256,15);
   tableB2BI['A'] = tableB2BI['a'] = 1;
   tableB2BI['C'] = tableB2BI['c'] = 2;
   tableB2BI['G'] = tableB2BI['g'] = 4;
   tableB2BI['T'] = tableB2BI['t'] = 8;
}//}}}
inline char base2int(char B){//{{{
   /* switch(B){
      case 'A': case 'a': return 0; 
      case 'C': case 'c': return 1; 
      case 'G': case 'g': return 2; 
      case 'T': case 't': return 3;
      default: return -1;
   } */
   return tableB2I[B];
}//}}}
inline int base2BAMint(char B){//{{{
   return tableB2BI[B];
}//}}}
inline void mapAdd(map<long,double > &m, long key, double val){//{{{
   if(m.count(key)==0)
      m[key] = val;
   else
      m[key] += val;
}//}}}
inline bool readHasPhred(const bam1_t *samA){//{{{
   if(samA->core.l_qseq < 1) return false;
   return bam1_qual(samA)[0] != 0xff;
}//}}}
// Count (number of deleteions) - (number of insertions).
long countDeletions(const bam1_t *samA){//{{{
   long deletionN = 0;
   for(long i=0;i<samA->core.n_cigar;i++){
      switch(bam1_cigar(samA)[i]&BAM_CIGAR_MASK){
         case BAM_CDEL:
            deletionN += (long)(bam1_cigar(samA)[i]>>BAM_CIGAR_SHIFT);
            break;
         case BAM_CINS:
            deletionN -= (long)(bam1_cigar(samA)[i]>>BAM_CIGAR_SHIFT);
            break;
      }
   }
   return deletionN;
}//}}}
inline bool getCigarOp(const bam1_t *samA, long cigarI, long *cigarOp,
                       long *cigarOpCount){//{{{
   if((cigarI<0) || (cigarI >= samA->core.n_cigar)) return false;
   *cigarOp = bam1_cigar(samA)[cigarI]&BAM_CIGAR_MASK;
   *cigarOpCount = (long)(bam1_cigar(samA)[cigarI]>>BAM_CIGAR_SHIFT);
   return true;
}//}}}
} // namespace ns_rD

using namespace ns_rD;

ReadDistribution::ReadDistribution(){ //{{{
   M=0;
   uniform = lengthSet = gotExpression = normalized = validLength = false;
   warnPos = warnTIDmismatch = warnUnknownTID = noteFirstMateDown = 0;
   procN = 1; 
#ifdef _OPENMP
   omp_set_num_threads(procN);
#endif
   lMu=100;
   lSigma=10;
   verbose = true;
   singleReadLength = 0;
   minFragLen=10000;
   lowProbMismatches = LOW_PROB_MISSES;
   lProbMis.resize(256,0);
   lProbHit.resize(256,0);
   for(long i=0; i<256; i++){
      lProbMis[i] = - i / 10.0 * log(10.0);
      lProbHit[i] = log1p(-exp(lProbMis[i]));
   }
   fillTable();
}//}}}
void ReadDistribution::writeWarnings() {//{{{
   if(warnPos>0){
      warning("ReadDistribution: %ld reads from a pair did not align to the expected strand of a transcript.\n   Use --unstranded option in case the 5' and 3' mate are not expected to be from sense and anti-sense strands respectively.\n", warnPos);
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
void ReadDistribution::setProcN(long procN){//{{{
   if(procN<0)procN=1;
   if(procN>32)procN=4;
#ifdef _OPENMP
   this->procN = procN;
   omp_set_num_threads(procN);
#else
   this->procN = 1;
#endif
}//}}}
bool ReadDistribution::init(long m, TranscriptInfo* trI, TranscriptSequence* trS, TranscriptExpression* trE, bool unstranded, bool verb){ //{{{
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
   this->unstranded = unstranded;
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
   validLength=true;
   computeLengthProb();
}//}}}
bool ReadDistribution::observed(fragmentP frag){ //{{{
   DEBUG(message("%s===%s\n",bam1_qname(frag->first),bam1_qname(frag->second));)
   long tid = frag->first->core.tid;
   if((frag->paired)&&(tid!=frag->second->core.tid)){
      warnTIDmismatch++;
      return false;
   }
   if((tid < 0)||(tid>=M)){
      warnUnknownTID++;
      return false;
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
   // Update Mismatch frequencies if no Phred. //{{{
   if((!readHasPhred(frag->first)) || (frag->paired && !readHasPhred(frag->second))){
      updateMismatchFreq(frag->first);
      if(frag->paired)updateMismatchFreq(frag->second);
   }
   // }}}
   // for uniform distribution ignore other estimation:
   if(uniform) return true;

   // check mates relative position: {{{
   if((frag->paired) && (frag->first->core.pos > frag->second->core.pos)){
      noteFirstMateDown ++;
      bam1_t *tmp = frag->second;
      frag->second = frag->first;
      frag->first = tmp;
   }
   if((frag->paired) && (!unstranded) && 
      ((frag->first->core.flag & BAM_FREVERSE) ||
       (! frag->second->core.flag & BAM_FREVERSE))){
      warnPos ++;
      return false;
   }//}}}
   // positional bias:
   // sequence bias:
   DEBUG(message("   positional & sequence bias\n");)
   if(! frag->paired){
      if(frag->first->core.flag & BAM_FREVERSE){
         // Antisense strand of transcript is 3'end of fragment
         updatePosBias(frag_first_endPos, readM_3, tid, Iexp);
         // readM_5 and uniformM_5 are always "second mates" 
         // this is assumed also in getP(...);
         updateSeqBias(frag_first_endPos, readM_3, tid, Iexp);
         // update sum of expression of  fragments of given length
         mapAdd(trFragSeen3[tid], (long)len, Iexp);
      }else{
         // Sense strand of transcript is 5'end of fragment
         updatePosBias( frag->first->core.pos, readM_5, tid, Iexp);
         updateSeqBias( frag->first->core.pos, readM_5, tid, Iexp);
         mapAdd(trFragSeen5[tid], (long)len, Iexp);
      }
   }else{
      updatePosBias( frag->first->core.pos, readM_5, tid, Iexp);
      updateSeqBias( frag->first->core.pos, readM_5, tid, Iexp);
      mapAdd(trFragSeen5[tid], (long)len, Iexp);
         
      updatePosBias( frag_second_endPos, readM_3, tid, Iexp);
      updateSeqBias( frag_second_endPos, readM_3, tid, Iexp);
      mapAdd(trFragSeen3[tid], (long)len, Iexp);
   }
   return true;
}//}}}
void ReadDistribution::normalize(){ //{{{
   // length distribution: {{{
   double newMu=0, newSigma=0;
  
   if(fragSeen>10){
      // Estimate mean and sigma for length distribution.
      newMu = logLengthSum / fragSeen;
      newSigma = sqrt(logLengthSqSum / fragSeen - newMu*newMu);
      if(verbose)message("ReadDistribution: fragment length mu: %lg sigma: %lg\n",newMu,newSigma);
      validLength = true;
   }
   if(lengthSet){
      // check difference between estimated mean and provided mean
      if(abs(newMu-lMu)>lSigma){
         warning("ReadDistribution: Estimated length mean (%lg) differs too much from the one provided (%lg).\n",newMu,lMu);
      }
   }else{
      // Use estimated mean and sigma;
      lMu = newMu;
      lSigma = newSigma;
      if(validLength)computeLengthProb();
   }
   // }}}
   // mismatch frequencies: {{{
   double lFreqSum;
   for(size_t i=0;i<lFreqHit.size();i++){
      lFreqSum = log(lFreqHit[i]+lFreqMis[i]);
      lFreqHit[i] = log(lFreqHit[i]) - lFreqSum;
      lFreqMis[i] = log(lFreqMis[i]) - lFreqSum;
   }
   // }}}
   if(uniform) return;
   map<long,double>::iterator mIt;
   long i,j,m,group,trLen,fragLen;
   double Iexp,norm;
   double binSize;
   // set Uniform position position bias: //{{{
   if(verbose)message("ReadDistribution: Computing uniform positional bias.\n");
   for(m=0;m<M;m++){
      //if(verbose)progressLogRD(m,M);
      trLen = trInf->L(m);
      if(trLen<trNumberOfBins)continue;
      binSize = (double)trLen / trNumberOfBins;
      //message(" %ld %ld %ld\n",m,trLen,trFragSeen[m].size());
      for(group=0;group<trSizesN;group++)
         if(trLen<trSizes[group])break;
      // update 5' positional bias
      for( mIt=trFragSeen5[m].begin(); mIt != trFragSeen5[m].end(); mIt++){
         fragLen = mIt->first;
         Iexp = mIt->second / (trLen - fragLen + 1);
         for(i=0;i<trNumberOfBins;i++){
            // update probability of each bin by Iexp*"effective length of current bin"
            if((i+1) * binSize <= fragLen)continue;
            if(i * binSize < fragLen){
               posProb[uniformM_5][group][trNumberOfBins -1 -i] +=
                  Iexp * ((i+1) * binSize - fragLen + 1);
            }else{
               posProb[uniformM_5][group][trNumberOfBins -1 -i] +=
                  Iexp * binSize;
            }
         }
      }  
      // update 3' positional bias
      for( mIt=trFragSeen3[m].begin(); mIt != trFragSeen3[m].end(); mIt++){
         fragLen = mIt->first;
         Iexp = mIt->second / (trLen - fragLen + 1);
         for(i=0;i<trNumberOfBins;i++){
            // update probability of each bin by Iexp*"effective length of current bin"
            if((i+1) * binSize <= fragLen)continue;
            if(i * binSize < fragLen){
               posProb[uniformM_3][group][i] +=
                  Iexp * ((i+1) * binSize - fragLen + 1);
            }else{
               posProb[uniformM_3][group][i] +=
                  Iexp * binSize;
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
      //if(verbose)progressLogRD(m,M);
      trLen = trInf->L(m);
      IexpSum5=0;
      for(mIt=trFragSeen5[m].begin();mIt!= trFragSeen5[m].end();mIt++)
         IexpSum5+=mIt->second / (trLen - mIt->first + 1);
      IexpSum3=0;
      mItR=trFragSeen5[m].rbegin();
      mIt=trFragSeen3[m].begin();
      // STL map iterator IS sorted by key <=> length
      for(p=0;p<trLen;p++){
         while((mIt!=trFragSeen3[m].end())&&(mIt->first <= p+1)){IexpSum3+=mIt->second/ (trLen - mIt->first + 1); mIt++;}
         while((mItR!=trFragSeen5[m].rend())&&(trLen-p < mItR->first)){IexpSum5-= mItR->second / (trLen - mItR->first + 1) ; mItR++;}
         updateSeqBias(p, uniformM_5, m, IexpSum5);
         // 3' end is expected to be "after"
         updateSeqBias(p+1, uniformM_3, m, IexpSum3);
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
   outF.precision(6);
   outF<<scientific;
   if(!outF.is_open()){
      error("ReadDistribution: Unable to open profile file: %s\n",(logFileName).c_str());
      return;
   }
   long i,j,g;
   outF<<"# BASES: (readM_5, readM_3, uniformM_5, uniformM_3)"<<endl;
   if(!uniform){
      for(j=0;j<4;j++){
         outF<<"# "<<endl;
         for(i=0;i<vlmmNodesN;i++){
            outF<<seqProb[j][i].getPsum('A')<<" "<<seqProb[j][i].getPsum('C')<<" "<<seqProb[j][i].getPsum('G')<<" "<<seqProb[j][i].getPsum('T')<<endl;
         }
      }
   }

   outF<<"#\n# Position: (readM_5, readM_3, uniformM_5, uniformM_3, weight_5, weight_3)"<<endl;
   if(!uniform){
      for(j=0;j<6;j++){
         outF<<"# "<<endl;
         for(g=0;g<=trSizesN;g++){
            for(i=0;i<trNumberOfBins;i++)
               outF<<posProb[j][g][i]<<" ";
            outF<<endl;
         }
      }
   }
   outF<<"#\n# Mismatch likelihood: (probHit, probMis)"<<endl;
   for(i=0;i<(long)lFreqHit.size();i++)outF<<exp(lFreqHit[i])<<" ";
   outF<<endl;
   for(i=0;i<(long)lFreqMis.size();i++)outF<<exp(lFreqMis[i])<<" ";
   outF<<endl;
   outF.close();
}//}}}
void ReadDistribution::updateMismatchFreq(bam1_t *samA) {//{{{
   if(! samA) return;
   bam1_core_t *samC = &samA->core;
   long i,j,k,kStart,kDir,len=samC->l_qseq;
   // Make sure we have place for storing data.
   if(len>(long)lFreqHit.size()){
      lFreqHit.resize(len,1.0);
      lFreqMis.resize(len,1.0);
   }
   // Set direction for storing mismatches depending on read orientation.
   if(samC->flag & BAM_FREVERSE){
      kStart = len - 1;
      kDir = -1;
   }else{
      kStart = 0;
      kDir = +1;
   }
   long deletionN = countDeletions(samA);
   string seq = trSeq->getSeq(samC->tid, samC->pos, len+deletionN, false);
   long cigarOp,cigarI,cigarOpCount;
   cigarOp=cigarI=cigarOpCount=0;
   // i - iterates within reference sequence
   // j - iterates within read
   // k - iterates within frequency arrays, can be reversed
   for(i=j=0,k=kStart;(i<len+deletionN) && (j<len);){
      if(cigarOpCount == 0){
         if(! getCigarOp(samA, cigarI, &cigarOp, &cigarOpCount))break;
         cigarI++;
      }
      switch(cigarOp){
         case BAM_CDEL: i+=cigarOpCount; cigarOpCount=0; continue;
         case BAM_CINS:
            j+= cigarOpCount; 
            k+= kDir * cigarOpCount;
            cigarOpCount=0; 
            continue;
      }
      if(base2int(seq[i]) > -1){
         if(base2BAMint(seq[i]) != bam1_seqi(bam1_seq(samA),j))lFreqMis[k]+=1;
         else lFreqHit[k]+=1;
      }
      i++;
      j++;
      k+=kDir;
      cigarOpCount --;
   }
}//}}}
pair<double,double> ReadDistribution::getSequenceLProb(bam1_t *samA) const{//{{{
   if(! samA) return pair<double, double>(0,0);
   double lProb=0,lowLProb=0, lPHit, lPMis;
   bam1_core_t *samC = &samA->core;
   uint8_t *qualP=bam1_qual(samA);
   bool hasPhred = readHasPhred(samA);
   long i,j,k,len=samC->l_qseq;
   long deletionN = countDeletions(samA);
   string seq = trSeq->getSeq(samC->tid, samC->pos, len+deletionN, false);
   long hitC, misC, addMisC;
   long cigarOp,cigarI,cigarOpCount;
   bool reversed = (samC->flag & BAM_FREVERSE);

   // First count the number fo misses to add for low probability. {{{
   cigarOp = cigarI = cigarOpCount = 0;
   hitC = misC = 0;
   // i - iterates within reference sequence
   // j - iterates within read
   for(i=j=0;(i<len+deletionN) && (j<len);){
      if(cigarOpCount == 0){
         if(! getCigarOp(samA, cigarI, &cigarOp, &cigarOpCount))break;
         cigarI++;
      }
      switch(cigarOp){
         case BAM_CDEL: i+=cigarOpCount; cigarOpCount=0; continue;
         case BAM_CINS: j+=cigarOpCount; cigarOpCount=0; continue;
      }
      if((base2int(seq[i]) == -1)||
         (base2BAMint(seq[i]) != bam1_seqi(bam1_seq(samA),j)))misC++;
      else hitC++;
      i++;
      j++;
      cigarOpCount --;
   }
   addMisC = max((long)1, lowProbMismatches - misC);
   // }}}
   
   cigarOp = cigarI = cigarOpCount = 0;
   for(i=j=0;(i<len+deletionN) && (j<len);){
      if(cigarOpCount == 0){
         if(! getCigarOp(samA, cigarI, &cigarOp, &cigarOpCount))break;
         cigarI++;
      }
      switch(cigarOp){
         case BAM_CDEL: i+=cigarOpCount; cigarOpCount=0; continue;
         case BAM_CINS: j+=cigarOpCount; cigarOpCount=0; continue;
         /*case BAM_CMATCH:
         case BAM_CEQUAL:
         case BAM_CDIFF:*/
      }
      if(hasPhred){
         lPHit = lProbHit[qualP[j]];
         lPMis = lProbMis[qualP[j]];
      }else{
         if(!reversed)k = j;
         else k = len-j-1;
         if((k>=0)&&(k<(long)lFreqHit.size())){
            lPHit = lFreqHit[k];
            lPMis = lFreqMis[k];
         }else{
            lPHit = lPMis = 0.5;
         }
      }
      if((base2int(seq[i]) == -1) ||
         (base2BAMint(seq[i]) != bam1_seqi(bam1_seq(samA),j))){
         // If bases don't match, multiply probability by probability of error.
         lProb += lPMis;
         lowLProb += lPMis;
      }else{
         lProb += lPHit;
         hitC --;
         if((addMisC>0) && (reversed || (addMisC>hitC))){
            // If there are some misses left add a 'miss' to the 'low probability'.
            lowLProb += lPMis;
            addMisC--;
         }else{
            lowLProb += lPHit;
         }
      }
      i++;
      j++;
      cigarOpCount --;
   }
   return pair<double, double>(lProb,lowLProb);
}//}}}
bool ReadDistribution::getP(fragmentP frag,double &lProb,double &lProbNoise){ //{{{
   lProb = ns_misc::LOG_ZERO;
   lProbNoise = ns_misc::LOG_ZERO;
   long tid = frag->first->core.tid;
   long trLen = trInf->L(tid),len;
   // Check transcript IDs {{{
   if((frag->paired)&&(tid!=frag->second->core.tid)){
      warnTIDmismatch++;
      return false;
   }
   if((tid < 0)||(tid>=M)){
      warnUnknownTID++;
      return false;
   }
   //}}}
   double lP = 0;
   // Get probability based on base mismatches: {{{
   pair<double, double> lpSeq1(0,0),lpSeq2(0,0);
   lpSeq1 = getSequenceLProb(frag->first);
   if(frag->paired)lpSeq2 = getSequenceLProb(frag->second);
   // }}}
   // Calculate reads' true end position: {{{
   long frag_first_endPos, frag_second_endPos=0;
   frag_first_endPos = bam_calend(&frag->first->core, bam1_cigar(frag->first));
   if(frag->paired){
      frag_second_endPos = bam_calend(&frag->second->core, bam1_cigar(frag->second));
   }
   // }}}
   if(frag->paired){
   // Get probability of length {{{
      if(frag->second->core.pos > frag->first->core.pos)
         len = frag_second_endPos - frag->first->core.pos;
      else{
         len = frag_first_endPos - frag->second->core.pos;
      }
      // compute length probability and normalize by probability of all possible lengths (cdf):
      // P*=lengthP/lengthNorm
      // }}}
      if(validLength) lP += getLengthLP(len) - getLengthLNorm(trLen);
   }else{
      len = frag_first_endPos - frag->first->core.pos;
   }
   if(uniform){
      // Get probability of position for uniform distribution
      // P*=1/(trLen-len+1)
      lP -= log(trLen - len + 1.0);
   }else{ // Positional & Sequence bias {{{
      // Get probability of position given read bias model
      // check mates' relative position:
      if( frag->paired && (frag->first->core.pos > frag->second->core.pos)){
         noteFirstMateDown ++;
         bam1_t *tmp = frag->second;
         frag->second = frag->first;
         frag->first = tmp;
      }
      if(!frag->paired){
         if(frag->first->core.flag & BAM_FREVERSE){
            // If read was reverse complement, then it's 3' mate.
            // P*=posBias3'*seqBias3'/weightNorm3'
            lP += log(getPosBias(frag->first->core.pos, frag_first_endPos, 
                                 mate_3, trLen)) +
               log(getSeqBias(frag_first_endPos , mate_3, tid )) -
               log(getWeightNorm( (long) len, mate_3, tid));
         }else{
            // P*=posBias5'*seqBias5'/weightNorm5'
            lP += log(getPosBias(frag->first->core.pos, frag_first_endPos,
                                 mate_5, trLen)) +
               log(getSeqBias(frag->first->core.pos, mate_5, tid )) -
               log(getWeightNorm( (long) len, mate_5, tid));
         }
      }else{
         // check strand of the reads:
         if((!unstranded) && 
            ((frag->first->core.flag & BAM_FREVERSE) ||
            (! frag->second->core.flag & BAM_FREVERSE))){
               warnPos ++;
               return false;
         }
//#pragma omp parallel sections num_threads (2) reduction(*:P)
//{
//   #pragma omp section
         // P*=1/weightNormFull
         lP -= log(getWeightNorm( (long) len, FullPair, tid));
//   #pragma omp section
//   {
         // P*=posBias5'*posBias3'*seqBias5'*seqBias3'
         lP += log(getPosBias(frag->first->core.pos, frag_second_endPos,
                              FullPair, trLen))
          + log(getSeqBias(frag->first->core.pos, mate_5, tid ))
          + log(getSeqBias(frag_second_endPos , mate_3, tid )); 
//   }
//}
      }
   } //}}}
   lProb = lP + lpSeq1.first+lpSeq2.first;
   lProbNoise = lP + lpSeq1.second+lpSeq2.second;
   return true;
}//}}}
void ReadDistribution::updatePosBias(long pos, biasT bias, long tid, double Iexp){ //{{{
   if(bias == readM_3)pos--;
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
   // Set correct start based on orientation.
   if((bias == readM_5)||(bias == uniformM_5)){
      start = pos - vlmmStartOffset - MAX_NODE_PAR;
      seq = trSeq->getSeq(tid, start, vlmmNodesN + MAX_NODE_PAR);
   }else{
      start = pos + vlmmStartOffset - vlmmNodesN ;
      // Get don't need complementing as it is always complement.
      seq = trSeq->getSeq(tid, start, vlmmNodesN + MAX_NODE_PAR);
      // Only reverse the sequence.
      reverse(seq.begin(),seq.end());
   }
   // Update bias weights.
   for(long i=0;i<vlmmNodesN;i++){
      seqProb[bias][i].update( Iexp, seq[i+2], seq[i+1], seq[i]);
   }
}//}}}
double ReadDistribution::getPosBias(long start, long end, readT read, long trLen) const { //{{{
   end --;
   // transcript too short:
   if(trLen < trNumberOfBins) return 1;
   long group, relS, relE;
   // choose group:
   for(group = 0;group < trSizesN;group++)
      if(trLen<trSizes[group])break;
   // find relative positions:
   relS = (start * trNumberOfBins) / trLen;
   if(relS>=trNumberOfBins)relS=trNumberOfBins-1;
   relE = (end * trNumberOfBins) / trLen;
   if(relE>=trNumberOfBins)relE=trNumberOfBins-1;
   double posBias = 1;
   // return bias weight
   if((read == FullPair) || (read == mate_5))
      posBias *= posProb[ weight_5 ][ group ][ relS ];
   if((read == FullPair) || (read == mate_3))
      posBias *= posProb[ weight_3 ][ group ][ relE ];
   return posBias;
}//}}}
double ReadDistribution::getSeqBias(long pos, readT read, long tid) const{ //{{{
   if(read==FullPair)return 0; // this should never happen
   long start;
   biasT bias,biasNorm;
   // Get sequence based on which fragment end we are dealing with.
   if(read == mate_5){
      start = pos - vlmmStartOffset - MAX_NODE_PAR;
   }else{
      start = pos + vlmmStartOffset - vlmmNodesN;
   }
   string seq = trSeq->getSeq(tid, start, vlmmNodesN + MAX_NODE_PAR);
   if(read == mate_5){
      bias = readM_5;
      biasNorm = uniformM_5;
   }else{
      bias = readM_3;
      biasNorm = uniformM_3;
      // Reverse the sequence for 3' end.
      reverse(seq.begin(),seq.end());
   }
   double B = 1;
   for(long i=0;i<vlmmNodesN;i++)
      // FIX HERE (probably that we are always doing 'same' division)
      B *= seqProb[bias][i].getP( seq[i+2], seq[i+1], seq[i]) /
           seqProb[biasNorm][i].getP( seq[i+2], seq[i+1], seq[i]);
   return B;
}//}}}
inline char ReadDistribution::getBase(long pos, const string &fSeq) const{ //{{{
   if((pos<0)||(pos>=(long)fSeq.size()))return 'N';
   return fSeq[pos];
}//}}}
double ReadDistribution::getSeqBias(long start, long end, readT read, const string &fSeq) const{ //{{{
   start = start - vlmmStartOffset - MAX_NODE_PAR;
   end = end + vlmmStartOffset + MAX_NODE_PAR - 1;
   
   double B = 1;
   long i,j;
   if((read==FullPair) || (read == mate_5)){
      for(i=0,j=start; i<vlmmNodesN; i++, j++)
         // FIX HERE (probably that we are always doing 'same' division)
         B *= seqProb[readM_5][i].getP( getBase(j+2,fSeq), getBase(j+1,fSeq), getBase(j,fSeq)) /
              seqProb[uniformM_5][i].getP( getBase(j+2,fSeq), getBase(j+1,fSeq), getBase(j,fSeq));
   }
   if((read==FullPair) || (read == mate_3)){
      // For 3' bias we go from 'end' position backwards.
      for(i=0,j=end; i<vlmmNodesN; i++, j--)
         // FIX HERE (probably that we are always doing 'same' division)
         B *= seqProb[readM_3][i].getP( getBase(j-2,fSeq), getBase(j-1,fSeq), getBase(j,fSeq)) /
              seqProb[uniformM_3][i].getP( getBase(j-2,fSeq), getBase(j-1,fSeq), getBase(j,fSeq));
   }
   return B;
}//}}}
/* inline char ReadDistribution::complementBase(char base) const{ //{{{
   if((base=='A')||(base=='a')) return'T';
   if((base=='T')||(base=='t')) return 'A';
   if((base=='C')||(base=='c')) return 'G';
   if((base=='G')||(base=='g')) return 'C';
   return 'N';
}//}}} */
double ReadDistribution::getWeightNorm(long len, readT read, long tid){ //{{{
   if(len == 0)return 1;
   if(weightNorms[read][tid].count(len) == 0){
      const string &trS = trSeq->getTr(tid);
      // We are not complementing.
      //for(size_t i=0;i<trRS.size();i++)trRS[i] = complementBase(trRS[i]);
      long trLen = trInf->L(tid), pos;
      double norm = 0,w;
      #pragma omp parallel for \
         private(w) \
         reduction(+:norm)
      for(pos = 0;pos <= trLen-len;pos++){
         w = getPosBias(pos, pos + len, read, trLen) *
             getSeqBias(pos, pos + len, read, trS);
         norm+=w;
      }
      weightNorms[read][tid][len] = norm;
//      message("w: %ld %ld %ld  %ld%lf\n",read,tid,len,trLen<"   ",norm);
      return norm;
   }
   return weightNorms[read][tid][len];
}//}}}
long ReadDistribution::getWeightNormCount() const{//{{{
   long length_sum=0;
   for(size_t i=0;i<weightNorms.size();i++)
      for(size_t j=0;j<weightNorms[i].size();j++)
         length_sum+=weightNorms[i][j].size();
   return length_sum;
}//}}}
double ReadDistribution::getLengthLP(long len) const{//{{{
   if(len>=(double)lLengthP.size())return computeLengthLP(len);
   return lLengthP[len];
}//}}}
double ReadDistribution::computeLengthLP(double len) const{//{{{
   //return 1./(len*lSigma*sqrt_2_pi)*exp(-pow(log(len) - lMu, (double)2.0)/(2 * pow(lSigma, (double)2)));
   if(len == 0)return ns_misc::LOG_ZERO;
   const double log_sqrt_2_pi = .918938533192; // log(sqrt(2*pi))
   const double lLen = log(len);
   return - (lLen +
             log(lSigma) + 
             log_sqrt_2_pi + 
             pow( (lLen - lMu) / lSigma, 2.0) / 2.0 );
}//}}}
double ReadDistribution::getLengthLNorm(long trLen) const{//{{{
   if(trLen<(double)lLengthNorm.size())return lLengthNorm[trLen];

   // erfc needs compiler with C99 standard 
   // other option might be to use boost/math/special_functions/erf.hpp
   const long double sqrt_2 = 1.41421356237309;
   long double CDF2 = erfcl((lMu-log((long double)trLen)) / (lSigma * sqrt_2));
   if(CDF2 == 0)return log(0.5)+ns_misc::LOG_ZERO;
   return (double)(log(0.5)+log(CDF2));
}//}}}
void ReadDistribution::computeLengthProb() {//{{{
   MyTimer timer;
   if(verbose){
      message("Pre-computing length probabilities. ");
      timer.start();
   }
   long max=0;
   if(trInf){
      for(long i=0;i<M;i++)if(trInf->L(i)>max)max=trInf->L(i);
      max = min(max,(long)150000);
   }else{
      max = 100000;
   }
   lLengthP.assign(max+1,ns_misc::LOG_ZERO);
   lLengthNorm.assign(max+1,ns_misc::LOG_ZERO);
   bool normIsOne = false;
   for(long i=1;i<=max;i++){
      if(normIsOne){
         // lP is LOG_ZERO already, set norm to log(1).
         lLengthNorm[i] = 0;
         continue;
      }
      lLengthP[i] = computeLengthLP(i);
      lLengthNorm[i] = ns_math::logAddExp(lLengthNorm[i-1],lLengthP[i]);
      if(lLengthNorm[i] > -1e-15){
         normIsOne=true;
      }
   }
   if(verbose)timer.current();
}//}}}
vector<double> ReadDistribution::getEffectiveLengths(){ //{{{
   vector<double> effL(M,0);
   long m,len,trLen,pos;
   double eL, lCdfNorm,lenP, wNorm;
   string trRS;
   // Make one caching array for each process.
   vector<vector<double> > posBias5All(procN),posBias3All(procN);
   MyTimer timer;
   timer.start();
   DEBUG(message("Eff length: validLength %d ; minFragLen: %ld.\n",(int)validLength,minFragLen));
   #pragma omp parallel for \
      schedule (dynamic,5) \
      private (len,trLen,pos,eL,lenP,wNorm,lCdfNorm,trRS)
   for(m=0;m<M;m++){
      if(verbose && (m!=0) && (M>20) && (m%(M/10)==0)){
         #pragma omp critical
         {
            message("# %ld done. ",m);
            timer.current();
         }
      }
      long threadID = 0;
#ifdef _OPENMP
      threadID = omp_get_thread_num();
#endif
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
         DEBUG(message("Copy sequence.\n"));
         const string &trS = trSeq->getTr(m);
         vector<double> &posBias5 = posBias5All[threadID];
         vector<double> &posBias3 = posBias3All[threadID];
         posBias5.resize(trLen);
         posBias3.resize(trLen);
         DEBUG(message("Precomputing posBias.\n"));
         for(pos = 0;pos<trLen;pos++){
            // Don't care about end position.
            posBias5[pos] = getPosBias(pos, trLen, mate_5, trLen)*
                            getSeqBias(pos, trLen, mate_5, trS);
            // Don't care about start position.
            posBias3[pos] = getPosBias(0, pos+1, mate_3, trLen)*
                            getSeqBias(0, pos+1, mate_3, trS);
         }
         eL=0;
         DEBUG(message("Computing norms.\n"));
         for(len=1;len<=trLen;len++){
            wNorm = 0;
            for(pos=0;pos <= trLen - len;pos++){
               wNorm += posBias5[pos] * posBias3[pos+len-1];
            }
            lenP = exp(getLengthLP( len ) - lCdfNorm);
            eL += lenP * wNorm;
         }
         // Check for weirdness and don't go below 0 (some transcripts already had 5 bases).
         // Function isnormal assumes C99 or C++11.
         if((!isnormal(eL)) || (eL <= 1)){
            effL[m] = trLen;
            DEBUG(message("weird: %lf %ld %ld\n",eL,trLen,m));
         }else{
            effL[m] = eL;
         }
      }
   }
   DEBUG(long same = 0);
   if(! uniform){
      // normalize effective length to same sum as original length
      double effSum=0,lSum=0;
      for(m=0;m<M;m++){
         DEBUG(if(effL[m] == trInf->L(m))same++);
         lSum+=trInf->L(m);
         effSum+=effL[m];
      }
      for(m=0;m<M;m++)effL[m] *= lSum/effSum;
   }
   DEBUG(message(" same: %ld.\n",same));
   for(m=0;m<M;m++)if(effL[m]<=0) effL[m]=trInf->L(m);
   return effL;
}//}}}

double VlmmNode::getPsum(char b) const{//{{{
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
double VlmmNode::getP(char b, char bp, char bpp) const{//{{{
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

