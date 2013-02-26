#ifndef ARGUMENTPARSER_H
#define ARGUMENTPARSER_H

#include<map>
#include<vector>
#include<string>

using namespace std;

enum OptionType {OTSTRING, OTLONG, OTBOOL, OTDOUBLE};
struct Option{//{{{
   OptionType type;
   string shortName,longName,description;
};//}}}

class ArgumentParser{
   private:
      map<string,string> mapS;
      map<string,long> mapL;
      map<string,bool> mapB;
      map<string,double> mapD;
      map<string,string> names;
      map<string,Option> validOptions;
      string programName, argumentDesc, programDesc;
      vector<string> arguments;
      vector<string> compulsory;
      long minimumArguments;

      bool existsOption(string name) const;
      bool existsName(string name) const;
      template <typename valueType>
      void appendDescription(string &desc,valueType defValue);
   public:
      bool verbose;

      ArgumentParser(const string &pD="",const string &aD="[FILES]", long minArgs = 1){//{{{
         verbose = false;
         init(pD,aD,minArgs);
      }//}}}
      void init(const string &pD="",const string &aD="[FILES]", long minArgs = 1){//{{{
         programDesc=pD; 
         argumentDesc=aD; 
         minimumArguments = minArgs;
      }//}}}
      bool parse(int n,char * argv[]); 
      void addOptionS(string shortName, string longName, string name, bool comp, string description="", string defValue="noDefault");
      void addOptionL(string shortName, string longName, string name, bool comp, string description="", long defValue=-47);
      void addOptionD(string shortName, string longName, string name, bool comp, string description="", double defValue=-47.47);
      void addOptionB(string shortName, string longName, string name, bool comp, string description="", bool defValue=false);
      bool verb() const { return verbose; }
      string getS(string name) const;
      long getL(string name) const;
      double getD(string name) const;
      const vector<string>& args() const { return arguments; }
      vector<double> getTokenizedS2D(string name) const;
      bool flag(string name) const;
      bool isSet(string name) const;
      void usage();
      void writeAll();
};

#endif
