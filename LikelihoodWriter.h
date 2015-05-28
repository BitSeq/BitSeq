#ifndef PROBSFILE_H
#define PROBSFILE_H

#include <string>
#include <vector>
#include <fstream>
#include "parseAlignment.h"

using namespace std;

namespace ns_parseAlignment {
class TagAlignment{//{{{
   protected:
      int_least32_t trId;
//      bool strand; // true = forward; false = reverse
      double prob,lowProb;
   public:
      TagAlignment(long t=0,double p = 0,double lp = 0){
         trId=(int_least32_t)t;
//         strand=s;
         prob=p;
         lowProb=lp;
      }
      long getTrId()const {return trId;}
      double getProb()const {return prob;}
      double getLowProb()const {return lowProb;}
      void setProb(double p){prob=p;}
}; //}}}
} // namespace ns_parseAlignment

class LikelihoodWriter {
public:
  virtual void writeRead(string name, vector<ns_parseAlignment::TagAlignment> alignments) = 0;
  virtual void writeDummy(string name) = 0;
  virtual void finalize() = 0;
};

class ProbWriter: public LikelihoodWriter {
private:
  ofstream outF;
public:
  ProbWriter(string fname, long Ntotal, long Nmap, long M);
  ~ProbWriter() {
    if (outF.is_open())
      outF.close();
  }
  void writeRead(string name, vector<ns_parseAlignment::TagAlignment> alignments);
  void writeDummy(string name);
  void finalize();
};

class MemoryWriter: public LikelihoodWriter {
private:
  AlignmentLikelihoods *likelihoods;
public:
  MemoryWriter(AlignmentLikelihoods *outalign, long Ntotal, long Nmap, long M);
  ~MemoryWriter() { }
  void writeRead(string name, vector<ns_parseAlignment::TagAlignment> alignments);
  void writeDummy(string name);
  void finalize();
};

#endif
