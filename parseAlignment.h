#ifndef PARSEALIGNMENT_H
#define PARSEALIGNMENT_H

#include "TagAlignments.h"

class AlignmentLikelihoods {
public:
  TagAlignments *alignments;
  long Ntotal, Nmap, M;
  AlignmentLikelihoods() {
    alignments=new TagAlignments();
    Ntotal = 0;
    Nmap = 0;
    M = 0;
  }
  ~AlignmentLikelihoods() {
    delete alignments;
  }
};

int parseAlignmentReal(int *argc,char* argv[],AlignmentLikelihoods *outalignments=0);
extern "C" int parseAlignment(int *argc,char* argv[]);

#endif
