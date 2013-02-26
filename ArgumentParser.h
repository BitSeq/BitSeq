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
/* Unused code {{{
class caseFreeCompare{
   private:
      char lowerC(char c)inline {
         if((c>='A')&&(c<='Z'))return c-'A'+'a';
         return c;
      }
   public:
      bool operator()(string a,string b){
         for(long i=0;i<(long)a.size();i++){
            if(i>=(long)b.size())return false;
            if(a[i]==b[i])continue;
            if(lowerC(a[i])<lowerC(b[i]))return true;
            if(lowerC(a[i])>lowerC(b[i]))return false;
            return ((a[i]>='a')&&(a[i]<='z'));
         }
         return true;
      }
} }}}*/


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

      ArgumentParser(string pD="", string aD="[FILES]", long minArgs = 1){//{{{
         verbose = false;
         programDesc=pD; 
         argumentDesc=aD; 
         minimumArguments = minArgs;
      }//}}}
      void init(string pD="", string aD="[FILES]", long minArgs = 1){//{{{
         programDesc=pD; 
         argumentDesc=aD; 
         minimumArguments = minArgs;
      }//}}}
      bool parse(int n,char * argv[]); 
      void addOptionS(string shortName, string longName, string name, bool comp, string description="", string defValue="noDefault");
      void addOptionL(string shortName, string longName, string name, bool comp, string description="", long defValue=-47);
      void addOptionD(string shortName, string longName, string name, bool comp, string description="", double defValue=-47.47);
      void addOptionB(string shortName, string longName, string name, bool comp, string description="", bool defValue=false);
      inline bool verb() const { return verbose; }
      string getS(string name) const;
      long getL(string name) const;
      double getD(string name) const;
      const vector<string>& args() const;
      vector<double> getTokenizedS2D(string name);
      bool flag(string name) const;
      bool isSet(string name) const;
      void usage();
      void writeAll();
      /*
      ArgumentParser(){//{{{
         init();
      }//}}}
      void init(string pD="", string aD="[FILES]", long minArgs = 1){//{{{
         programDesc=pD; 
         argumentDesc=aD; 
         minimumArguments = minArgs;
      }//}}}*/
};

#endif
