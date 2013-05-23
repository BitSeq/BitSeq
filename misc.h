#ifndef MISC_H
#define MISC_H

#include<fstream>

#include "ArgumentParser.h"
#include "PosteriorSamples.h"
#include "TranscriptInfo.h"

namespace ns_math {

// For a=log(x), b=log(y); compute log(x+y).
double logAddExp(double a, double b);

// For vals_i = log(x_i); compute log(sum(x_i)) for st<=i<en.
double logSumExp(const vector<double> &vals, long st = 0, long en = -1);

}

namespace ns_expression {

// Return output type based on the command line argument (one of theta/rpkm/counts/tau).
string getOutputType(const ArgumentParser &args, const string &defaultType = "rpkm");
}

namespace ns_misc {

// Value to use instead of log(0).
const double LOG_ZERO=-100;

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

// Convert string into lower case.
string toLower(string str);
}

namespace ns_genes {
bool getLog(const ArgumentParser &args);

bool prepareInput(const ArgumentParser &args, TranscriptInfo *trInfo, PosteriorSamples *samples, long *M, long *N, long *G);
} // namespace ns_genes

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
