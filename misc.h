#ifndef MISC_H
#define MISC_H

#include<fstream>

#include "ArgumentParser.h"
#include "PosteriorSamples.h"

namespace ns_misc {

const double LOG_ZERO=-1000;

long getSeed(const ArgumentParser &args);

bool openOutput(const ArgumentParser &args, ofstream *outF);
bool openOutput(const string &name, ofstream *outF);

// Reads and initializes files containing samples fro each condition and each replicate.
bool readConditions(const ArgumentParser &args, long *C, long *M, long *N, Conditions *cond);

// Compute confidence intervals.
void computeCI(double cf, vector<double> *difs, double *ciLow, double *ciHigh);

}

namespace ns_params{

struct paramT {//{{{
   double expr, alpha, beta;
   bool operator< (const paramT &p2) const{
      return expr<p2.expr;
   }
};//}}}

// Read hyperparameters from a file specified by file name.
// If outF is not NULL, it copies header from input file to outF.
// The vector is sorted by expression at the end.
bool readParams(const string &name, vector<paramT> *params, ofstream *outF = NULL);

}
#endif
