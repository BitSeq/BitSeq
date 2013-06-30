#ifndef ARGUMENTPARSER_H
#define ARGUMENTPARSER_H

#include<map>
#include<string>
#include<vector>

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

      bool existsOption(const string &name) const;
      bool existsName(const string &name) const;
      template <typename valueType>
      void appendDescription(string *desc,valueType defValue);
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
      void addOptionS(const string &shortName, const string &longName, const string &name, bool comp, const string &description="", const string &defValue="noDefault");
      void addOptionL(const string &shortName, const string &longName, const string &name, bool comp, const string &description="", long defValue=-47);
      void addOptionD(const string &shortName, const string &longName, const string &name, bool comp, const string &description="", double defValue=-47.47);
      void addOptionB(const string &shortName, const string &longName, const string &name, bool comp, const string &description="", bool defValue=false);
      bool verb() const { return verbose; }
      string getS(const string &name) const;
      long getL(const string &name) const;
      double getD(const string &name) const;
      const vector<string>& args() const { return arguments; }
      vector<double> getTokenizedS2D(const string &name) const;
      bool flag(const string &name) const;
      bool isSet(const string &name) const;
      void usage();
      void writeAll();
      // Update value of existing option.
      void updateS(const string &name, const string &value);
};

#endif
