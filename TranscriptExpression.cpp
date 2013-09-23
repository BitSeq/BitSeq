#include<algorithm>

#include "TranscriptExpression.h"
#include "FileHeader.h"
#include "common.h"

TE_FileType TranscriptExpression::guessFileType(const string &fileName){//{{{
   string extension = fileName.substr(fileName.rfind(".")+1);
   if(extension == "thetaMeans") return SAMPLER_MEANS;
   if(extension == "m_alphas") return M_ALPHAS;
   // Ends with 'mean' or 'variance' or 'var'
   if((extension.rfind("mean") == extension.size() - 4) ||
      (extension.rfind("variance") == extension.size() - 8) ||
      (extension.rfind("var") == extension.size() - 3)) return MEAN_VARIANCE;
   // Default is SAMPLER_MEANS.
   return SAMPLER_MEANS;
}//}}}
TranscriptExpression::TranscriptExpression(){//{{{
   M=0;
   logged=false;
}//}}}
TranscriptExpression::TranscriptExpression(const string &fileName, TE_FileType fileType){//{{{
   TranscriptExpression();
   readExpression(fileName,fileType);
}//}}}
bool TranscriptExpression::readExpression(const string &fileName, TE_FileType fileType){//{{{
   long i;
   if(fileType == GUESS)fileType = guessFileType(fileName);
   ifstream varFile(fileName.c_str());
   FileHeader fh(&varFile);
   if((!fh.varianceHeader(&M,&logged))||(M==0)){
      error("TranscriptExpression: Problem loading variance file %s\n",(fileName).c_str());
      return false;
   }
   // M_ALPHAS file contains nosie transcript.
   if(fileType == M_ALPHAS) M--;
   trs.resize(M);
   if(fileType == SAMPLER_MEANS){
      double count,mean2;
      for(i=0;i<M;i++){
         varFile>>trs[i].id>>trs[i].exp>>count>>mean2>>trs[i].var;
         // IDs in SAMPLER_MEANS file are shifted by 1
         trs[i].id--;
         varFile.ignore(1000,'\n');
         if(varFile.bad()){
            error("TranscriptExpression: Problem reading transcript %ld.\n",i);
            return false;
         }
      }
   }else if(fileType == MEAN_VARIANCE){
      for(i=0;i<M;i++){
         trs[i].id=i;
         varFile>>trs[i].exp>>trs[i].var;
         varFile.ignore(1000,'\n');
         if(varFile.bad()){
            error("TranscriptExpression: Problem reading transcript %ld.\n",i);
            return false;
         }
      }
   }else if(fileType == M_ALPHAS){
      double alpha, beta, beta0;
      // Skip first entry - noise transcript.
      varFile>>trs[0].exp>>alpha>>beta0;
      varFile.ignore(1000,'\n');
      for(i=0;i<M;i++){
         trs[i].id=i;
         varFile>>trs[i].exp>>alpha>>beta;
         // Beta0 is the sum of all except noise.
         trs[i].exp = alpha / beta0;
         trs[i].var = alpha * (beta0-alpha) / (beta0 * beta0 * (beta0 + 1));
         varFile.ignore(1000,'\n');
         if(varFile.bad()){
            error("TranscriptExpression: Problem reading transcript %ld.\n",i);
            return false;
         }
      }
   }
   fh.close();
   return true;
}//}}}
void TranscriptExpression::doSort(bool reverse){//{{{
   if(! reverse)
      sort(trs.begin(),trs.end());
   else
      sort(trs.rbegin(),trs.rend());
}//}}}

