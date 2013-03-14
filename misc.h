#ifndef MISC_H
#define MISC_H

#include<fstream>

#include "ArgumentParser.h"
#include "PosteriorSamples.h"

namespace ns_math {

// For a=log(x), b=log(y); compute log(x+y).
double logAddExp(double a, double b);

// For vals_i = log(x_i); compute log(sum(x_i)).
double logSumExp(const vector<double> &vals);

}

namespace ns_misc {

// Value to use instead of log(0).
const double LOG_ZERO=-1000;

// Return seed; either using seed set in args, or by using time(NULL) as seed.
long getSeed(const ArgumentParser &args);

// Open output file based on standard argument --outFile=<outFileName>.
bool openOutput(const ArgumentParser &args, ofstream *outF);
// Open output file of a give name.
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
