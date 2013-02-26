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
   PosteriorSamples();
//   bool init(long n, long m, bool t, string fileName);
   bool initSet(long &m, long &n, string fileName);
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
      PosteriorSamples *samples;
      vector<pair<long,long> > cIndex;
      
      long getIndex(long max); // return index without checking for duplicats
   public:
      Conditions();
      void close();
      long getRC(long c);
      long getRN(){ return CN;}
      long getC(){ return C;}
      bool init(string trFileName, vector<string> filesGot, long *c, long *m, long *n);
      bool init(string trFileName, vector<string> filesGot, long *m, long *n);
      bool setNorm(vector<double> norms);
      bool getTranscript(long cond, long rep, long tr, vector<double> &trSamples);
      bool getTranscript(long cond, long tr, vector<double> &trSamples);
      bool getTranscript(long cond, long tr, vector<double> &trSamples, long samplesN);
      bool logged(){return areLogged;}
};//}}}
