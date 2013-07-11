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

      bool existsOption(const string &name, bool warn = false) const;
      bool existsName(const string &name) const;
      template <typename valueType>
      void appendDescription(string *desc,valueType defValue);
   public:
      // The value of verbose option for direct access.
      bool verbose;

      // Constructor for the class sets: programDescription, additional string
      // and minimum number of required arguments.
      ArgumentParser(const string &pD="",const string &aD="[FILES]", long minArgs = 1){//{{{
         verbose = false;
         init(pD,aD,minArgs);
      }//}}}
      // Init function for initialization, sets the same values as constructor.
      void init(const string &pD="",const string &aD="[FILES]", long minArgs = 1){//{{{
         programDesc=pD; 
         argumentDesc=aD; 
         minimumArguments = minArgs;
      }//}}}
      // Parse function given number of arguments and array of arguments 
      // it processes the arguments and makes options available through 
      // get[S/L/D] functions and args() function.
      bool parse(int n,char * argv[]); 
      /*
       * SETTERS:
       */
      // Add option (string) adds new option, name is the name used for referring 
      // to it.
      void addOptionS(const string &shortName,
                      const string &longName,
                      const string &name,
                      bool comp,
                      const string &description="", 
                      const string &defValue="noDefault");
      // Add option (long).
      void addOptionL(const string &shortName, const string &longName,
                      const string &name, bool comp, const string &description="",
                      long defValue=-47);
      // Add option (double).
      void addOptionD(const string &shortName, const string &longName,
                      const string &name, bool comp, const string &description="",
                      double defValue=-47.47);
      // Add option (boolean or 'flag').
      void addOptionB(const string &shortName, const string &longName, 
                      const string &name, bool comp, const string &description="",
                      bool defValue=false);
      /*
       * GETTERS:
       */
      // Return reference to vector of arguments
      // (i.e. the strings provided with no -/-- modifier).
      const vector<string>& args() const { return arguments; }
      // Return true if option <name> was set.
      bool isSet(const string &name) const;
      // Return value of string option <name>.
      string getS(const string &name) const;
      // Return value of string option <name> in lower case.
      string getLowerS(const string &name) const;
      // Return value of integer option <name>.
      long getL(const string &name) const;
      // Return value of double option <name>.
      double getD(const string &name) const;
      // Return value of bool option <name>.
      bool flag(const string &name) const;
      // Return value of verbose.
      bool verb() const { return verbose; }

      /*
       * OTHER:
       */
      // (Advanced get) Return tokenized (comma separated) string as vector of doubles.
      vector<double> getTokenizedS2D(const string &name) const;
      // Write usage string.
      void usage();
      // Write all options.
      void writeAll();
      // Update value of existing string option.
      void updateS(const string &name, const string &value);
};

#endif
