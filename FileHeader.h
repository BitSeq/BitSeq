#ifndef FILEHEADER_H
#define FILEHEADER_H

#include<fstream>
#include<map>
#include<vector>

using namespace std;

const long no_value = -4747;

// FileHeader class parses file headers (lines starting with # at the beginning of the file).
// Every word (space separated string) is considered a possible FLAG.
// If a FLAG is followed by a numeric value, than the value is stored as the FLAG's value.
// The individual functions then just look whether FLAG was present, and in case of integers, whether it had some value assigned to it.
class FileHeader {
 private:
   ifstream *file;
   map<string,long> values;
   bool readValues(ofstream *outF = NULL);

   void skipEmptyLines();
   static vector<string> tokenizer(const string &input,const string &space = " ");
 public:
   FileHeader(ifstream *f = NULL) {
      file = f;
   }
   void setFile(ifstream *f){
      file = f;
   }
   void close(){
      file->close();
      file=NULL;
   }
   bool samplesHeader(long *n, long *m, bool *transposed, bool *logged = NULL);
   bool transcriptsHeader(long *m, long *colN);
   bool probHeader(long *Nmap,long *Ntotal,bool *newformat);
   bool varianceHeader(long *m,bool *logged);
   bool paramsHeader(long *parN, ofstream *outF);
};

#endif
