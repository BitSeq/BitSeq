#include<algorithm>

#include "TranscriptExpression.h"
#include "FileHeader.h"
#include "common.h"

TranscriptExpression::TranscriptExpression(){//{{{
   M=0;
   logged=false;
}//}}}
TranscriptExpression::TranscriptExpression(string fileName, TE_FileType fileType){//{{{
   TranscriptExpression();
   readExpression(fileName,fileType);
}//}}}
bool TranscriptExpression::readExpression(string fileName, TE_FileType fileType){//{{{
   long i;
   ifstream varFile(fileName.c_str());
   FileHeader fh(&varFile);
   if((!fh.varianceHeader(M,logged))||(M==0)){
      error("TranscriptExpression: Problem loading variance file %s\n",(fileName).c_str());
      return false;
   }
   trs.resize(M);
   if(fileType == SAMPLER_MEANS){
      for(i=0;i<M;i++){
         trs[i].var=0;
         varFile>>trs[i].id>>trs[i].exp;
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

