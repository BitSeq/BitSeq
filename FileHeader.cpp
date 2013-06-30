#include<cstdlib>

#include "FileHeader.h"
#include "common.h"

using namespace ns_fileHeader;

void FileHeader::skipEmptyLines() {//{{{
   if(!file) return;
   while(file->good() &&
         ((file->peek() == ' ') ||
          (file->peek() == '\n')))
      file->get();
}//}}}

// static
vector<string> FileHeader::tokenizer(const string &input,const string &space){//{{{
   vector<string> ret;
   long pos=0,f=0,n=input.size();
   while((pos<n)&&(f<n)&&(f>=0)){
      f=input.find(space,pos);
      if(f==pos)pos++;
      else{
         if((f<n)&&(f>=0)){
            ret.push_back(input.substr(pos,f-pos));
            pos=f+1;
         }
      }
   }
   if(pos<n)ret.push_back(input.substr(pos,n-pos));
   return ret;
} //}}}

bool FileHeader::readValues(ofstream *outF){//{{{
   if((file==NULL)||(!file->is_open())){
      error("FileHeader: Input file not opened for reading.\n");
      return false;
   }
   string line;
   vector<string> words;
   long value;
   char *chP;
   skipEmptyLines();
   while(file->good() && (file->peek() == '#')){
      // Read line.
      getline(*file, line);
      // If outF is defined, copy the header there.
      if(outF!=NULL)(*outF)<<line<<endl;
      skipEmptyLines();
      // Tokenize line into words.
      words = tokenizer(line);
      // Store words as flags. Start with 1st word as the 0th one are hashes.
      // If word is followed by a numeric value, use it as a value for the flag.
      for(long i=1;i<(long)words.size();i++){
         // Only add new entry if it wasn't there already.
         if(values.count(words[i])==0)
            values[words[i]] = no_value;
         // See if next word is numeric and if so, then use it as a value.
         if(i+1<(long)words.size()){
            value = strtol(words[i+1].c_str(), &chP, 10);
            // Conversion was succesful the value is non-zero OR the pointer should point to end of string (null character).
            if((value!=0)||(*chP=='\0')) {
               // Save value and skip the number.
               values[words[i]] = value;
               i++;
            }
         }
      }
   }
   return true;
}//}}}

bool FileHeader::samplesHeader(long *n, long *m, bool *transposed, bool *logged){//{{{
   if(!readValues()){
      *n=0;
      *m=0;
      return false;
   }
   if(logged!=NULL)if(values.count("L"))*logged = true;
   if(values.count("T"))*transposed = true;
   if(values.count("M") && (values["M"]!=no_value))*m = values["M"];
   if(values.count("N") && (values["N"]!=no_value))*n = values["N"];
   return true;
}//}}}

bool FileHeader::transcriptsHeader(long *m, long *colN){//{{{
   if(!readValues()){
      *m=0;
      return false;
   }
   if(values.count("M") && (values["M"]!=no_value))*m = values["M"];
   if(colN!=NULL)
      if(values.count("colN") && (values["colN"]!=no_value))*colN = values["colN"];
   return true;
}//}}}

bool FileHeader::probHeader(long *Nmap,long *Ntotal, AlignmentFileType *format){//{{{
   if(!readValues()){
      *Nmap=0;
      return false;
   }
   if(values.count("LOGFORMAT")){*format = LOG_FORMAT;}
   else if(values.count("NEWFORMAT")){*format = NEW_FORMAT;}
   else *format = OLD_FORMAT;
   if(values.count("Ntotal") && (values["Ntotal"]!=no_value))*Ntotal = values["Ntotal"];
   if(values.count("Nmap") && (values["Nmap"]!=no_value))*Nmap = values["Nmap"];
   return true;
}//}}}

bool FileHeader::varianceHeader(long *m,bool *logged){//{{{
   if(!readValues()){
      *m=0;
      return false;
   }
   if(logged!=NULL)if(values.count("L"))*logged = true;
   if(values.count("M") && (values["M"]!=no_value))*m = values["M"];
   return true;
}//}}}

bool FileHeader::paramsHeader(long *parN, ofstream *outF){//{{{
   if(!readValues(outF)){
      *parN=0;
      return false;
   }
   *parN = 0;
   if(values.count("PN") && (values["PN"]!=no_value))*parN = values["PN"];
   return true;
}//}}}
   


