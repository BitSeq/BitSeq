#ifndef POSTERIORSAMPLES_H
#define POSTERIORSAMPLES_H

#include<vector>
#include<fstream>
#include<string>

using namespace std;

const long PS_maxStoredSamples = 100000000;

class PosteriorSamples{//{{{
   private:
      long N,M;
      double norm;
      bool transposed,failed,areLogged;
      ifstream samplesF;
      vector<long> lines;
      vector<vector<double> > samples;

      bool open(string fileName);
      bool read();
   public:
   PosteriorSamples() { clear(); }
   ~PosteriorSamples() { close(); }
   // Copy constructor and assginment. Both just create new class. For vectors only.
   PosteriorSamples(const PosteriorSamples &other) { clear(); }
   PosteriorSamples& operator=(const PosteriorSamples & other) { //{{{
      close();
      clear();
      return *this;
   } //}}}
   void clear();
   bool initSet(long *m, long *n, string fileName);
   bool getTranscript(long tr, vector<double> &trSamples);
   void close();
   bool logged(){return areLogged;}
   void setNorm(double norm){this->norm = norm;}
};//}}}

class Conditions{//{{{
   private:
      long M,N,CN,C;
      bool mapping,areLogged;
      vector<long> Ms,Ns;
      vector<vector <long> > trMap;
      vector<PosteriorSamples> samples;
      vector<pair<long,long> > cIndex;
      
      long getIndex(long max); // return index without checking for duplicats
   public:
      Conditions();
      void close();
      long getRC(long c) const;
      long getRN() const { return CN;}
      long getC() const { return C;}
      bool init(string trFileName, vector<string> filesGot, long *c, long *m, long *n);
      bool init(string trFileName, vector<string> filesGot, long *m, long *n);
      bool setNorm(vector<double> norms);
      bool getTranscript(long cond, long rep, long tr, vector<double> &trSamples);
      bool getTranscript(long cond, long tr, vector<double> &trSamples);
      bool getTranscript(long cond, long tr, vector<double> &trSamples, long samplesN);
      double probFC(vector<double> x, vector<double> y, double logThreshold);
      bool transcriptStat(vector<double> sam, long tr, vector< vector<double> > &stat);
      bool logged() const { return areLogged; }
};//}}}

#endif
