#include <string>
#include <utility>
#include <map>
#include <vector>
#include <cstring>
#include <iostream>

#include "common.h"

extern "C" int parseAlignment(int *argc, char* argv[]);
extern "C" int estimateExpression(int *argc, char* argv[]);
extern "C" int estimateVBExpression(int *argc, char* argv[]);
extern "C" int estimateHyperPar(int *argc, char* argv[]);
extern "C" int estimateDE(int *argc, char* argv[]);
extern "C" int estimateFCProb(int *argc, char* argv[]);
extern "C" int getGeneExpression(int *argc,char* argv[]);
extern "C" int getVariance(int *argc,char* argv[]);
extern "C" int getWithinGeneExpression(int *argc,char* argv[]);

typedef int (*mainfun)(int *, char **);
typedef std::map<std::string, mainfun> mainfun_map;
typedef std::vector<std::pair<std::string, std::string> > doc_vector;

class Commands {
private:
  mainfun_map commandmap;
  doc_vector docs;
public:
  int call(std::string command, int *argc, char* argv[]) {
    mainfun_map::const_iterator it = commandmap.find(command);
    if (it==commandmap.end())
      return -1;
    return (*it->second)(argc, argv);
  }
  void write_docs() {
    for(doc_vector::const_iterator it = docs.begin(); it != docs.end(); ++it) {
      message("%s\n\t%s\n", it->first.c_str(), it->second.c_str());
    }
  }
  
  Commands() {
    commandmap.insert(std::make_pair("parseAlignment", &parseAlignment));
    docs.push_back(std::make_pair("parseAlignment", "Preprocess alignment results for expression estimation"));
    commandmap.insert(std::make_pair("estimateExpression", &estimateExpression));
    docs.push_back(std::make_pair("estimateExpression", "Accurate expression estimation using MCMC sampling"));
    commandmap.insert(std::make_pair("estimateVBExpression", &estimateVBExpression));
    docs.push_back(std::make_pair("estimateVBExpression", "Fast expression estimation using variational inference"));
    commandmap.insert(std::make_pair("estimateHyperPar", &estimateHyperPar));
    docs.push_back(std::make_pair("estimateHyperPar", "Estimate hyperparameters"));
    commandmap.insert(std::make_pair("estimateDE", &estimateDE));
    docs.push_back(std::make_pair("estimateDE", "Estimate differential expression from replicated samples"));
    commandmap.insert(std::make_pair("estimateFCProb", &estimateFCProb));
    docs.push_back(std::make_pair("estimateFCProb", "Estimate fold change probability from two samples"));
    commandmap.insert(std::make_pair("getGeneExpression", &getGeneExpression));
    docs.push_back(std::make_pair("getGeneExpression", "Compute gene expression from transcript expression"));
    commandmap.insert(std::make_pair("getVariance", &getVariance));
    docs.push_back(std::make_pair("getVariance", "Estimate variance from MCMC samples"));
    commandmap.insert(std::make_pair("getWithinGeneExpression", &getWithinGeneExpression));
    docs.push_back(std::make_pair("getWithinGeneExpression", "Compute transcript relative expression within gene"));
  }
};


void writeHelp(Commands commands)
{
  message("Usage: bitseq command options\n\n");
  message("For help on a command: bitseq command\n\n");
  message("Commands:\n\n");
  commands.write_docs();
}


int main(int argc, char* argv[]) {
  Commands commands;
  int newargc = argc - 1;
  char **newargv = new char*[argc];

  for (int i=0; i<argc-1; i++) {
    newargv[i] = new char[strlen(argv[i+1])+1];
    strcpy(newargv[i], argv[i+1]);
  }
  if (argc < 2) {
    writeHelp(commands);
    return -1;
  }
  return commands.call(argv[1], &newargc, newargv);
}
