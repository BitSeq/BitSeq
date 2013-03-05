#include "ArgumentParser.h"
#include "transposeFiles.h"
#include "common.h"

int main(int argc,char* argv[]){
   string programDescription = 
"Transposes [input files] into [outFileName] so that there are M lines with N columns each.";
   ArgumentParser args(programDescription,"[input files]",1);
   args.addOptionS("o","outFile","outFileName",1,"Name of the output file.");
   if(!args.parse(argc,argv))return 0;
   if(args.verbose)buildTime(argv[0],__DATE__,__TIME__);

   if(transposeFiles(args.args(),args.getS("outFileName"),args.verbose)){
      if(args.verbose)message("DONE.");
      return 0;
   }else{
      error("Failed.");
      return 1;
   }
}


