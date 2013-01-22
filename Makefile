CXX = g++
#CFLAGS = -Wall 
HOSTNAME = $(shell hostname)
ARCH = -mtune=generic



ifeq ($(HOSTNAME), valiant)
	ARCH = -march=core2
endif
ifeq ($(HOSTNAME), master1.templar.internal)
	ARCH = -march=core2
endif
ifeq ($(HOSTNAME), gauss)
	ARCH = -march=native
endif
ifeq ($(HOSTNAME), rpc465.cs.man.ac.uk)
	ARCH = -march=native
endif

COFLAGS = $(ARCH) -O3 -pipe
# -ffast-math segfaults with old gcc
#COFLAGS = $(ARCH) -O3 -ffast-math -pipe
DBGFLAGS = -ggdb
CXXFLAGS = -Wall $(COFLAGS)
#COFLAGS = -Wall -U_FORTIFY_SOURCE -O2 -ffast-math 
PROGRAMS = estimateVBExpression extractSamples getGeneExpression getWithinGeneExpression estimateExpression transposeLargeFile estimateDE getVariance estimateHyperPar getPPLR convertSamples parseAlignment getFoldChange
# getCount estimate0_4HyperPar estimate0_4DE 
BOOSTFLAGS = -I .
OPENMP = -fopenmp -DSUPPORT_OPENMP

#host:
#	echo $(HOSTNAME)

all: $(PROGRAMS)

#%.cpp%.o:  
#	$(CXX) $(COFLAGS) -c $<
#.PHONY: test

#estimate0_4HyperPar: estimate0_4HyperPar.cpp common.o PosteriorSamples.o ArgumentParser.o lowess.o TranscriptExpression.o 
#	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) estimate0_4HyperPar.cpp common.o PosteriorSamples.o ArgumentParser.o TranscriptExpression.o lowess.o -o estimate0_4HyperPar

#estimate0_4DE:  estimate0_4DE.cpp common.o PosteriorSamples.o ArgumentParser.o
#	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) estimate0_4DE.cpp common.o PosteriorSamples.o ArgumentParser.o -o estimate0_4DE

estimateVBExpression:  estimateVBExpression.cpp common.o ArgumentParser.o VariationalBayes.o SimpleSparse.o TagAlignments.o asa103/asa103.o
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) $(OPENMP) estimateVBExpression.cpp common.o ArgumentParser.o VariationalBayes.o SimpleSparse.o TagAlignments.o asa103/asa103.o -o estimateVBExpression

VariationalBayes.o: VariationalBayes.cpp VariationalBayes.h SimpleSparse.o
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) $(OPENMP) -c VariationalBayes.cpp 

SimpleSparse.o: SimpleSparse.cpp SimpleSparse.h
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) $(OPENMP) -c SimpleSparse.cpp

extractSamples: extractSamples.cpp common.o ArgumentParser.o PosteriorSamples.o common.h
	$(CXX) $(CXXFLAGS) extractSamples.cpp common.o ArgumentParser.o PosteriorSamples.o -o extractSamples

getFoldChange: getFoldChange.cpp common.o ArgumentParser.o PosteriorSamples.o 
	$(CXX) $(CXXFLAGS) getFoldChange.cpp common.o ArgumentParser.o PosteriorSamples.o -o getFoldChange

getGeneExpression: getGeneExpression.cpp common.o ArgumentParser.o TranscriptInfo.o PosteriorSamples.o common.h
	$(CXX) $(CXXFLAGS) getGeneExpression.cpp common.o ArgumentParser.o TranscriptInfo.o PosteriorSamples.o -o getGeneExpression

getWithinGeneExpression: getWithinGeneExpression.cpp common.o ArgumentParser.o TranscriptInfo.o PosteriorSamples.o common.h
	$(CXX) $(CXXFLAGS) getWithinGeneExpression.cpp common.o ArgumentParser.o TranscriptInfo.o PosteriorSamples.o -o getWithinGeneExpression

parseAlignment: parseAlignment.cpp common.o ArgumentParser.o TranscriptInfo.o TranscriptSequence.o TranscriptExpression.o ReadDistribution.o samtools/sam.o
	$(CXX) $(CXXFLAGS) $(OPENMP) -Isamtools samtools/*.o common.o ArgumentParser.o TranscriptInfo.o TranscriptSequence.o TranscriptExpression.o ReadDistribution.o parseAlignment.cpp -lz -o parseAlignment

samtools/sam.o:
	make --directory samtools

convertSamples: convertSamples.cpp common.o ArgumentParser.o TranscriptInfo.o
	$(CXX) $(CXXFLAGS) convertSamples.cpp common.o ArgumentParser.o  TranscriptInfo.o -o convertSamples

getPPLR:  getPPLR.cpp common.o ArgumentParser.o PosteriorSamples.o
	$(CXX) $(CXXFLAGS) getPPLR.cpp common.o ArgumentParser.o PosteriorSamples.o -o getPPLR

estimateHyperPar: estimateHyperPar.cpp common.o PosteriorSamples.o ArgumentParser.o lowess.o TranscriptExpression.o 
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) estimateHyperPar.cpp common.o PosteriorSamples.o ArgumentParser.o TranscriptExpression.o lowess.o -o estimateHyperPar

estimateDE:  estimateDE.cpp common.o PosteriorSamples.o ArgumentParser.o
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) estimateDE.cpp common.o PosteriorSamples.o ArgumentParser.o -o estimateDE

estimateExpression:  estimateExpression.cpp common.o ArgumentParser.o Sampler.o GibbsSampler.o CollapsedSampler.o GibbsParameters.o TranscriptInfo.o transposeFiles.o TagAlignments.o
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) $(OPENMP) estimateExpression.cpp common.o ArgumentParser.o Sampler.o GibbsSampler.o CollapsedSampler.o GibbsParameters.o TranscriptInfo.o transposeFiles.o TagAlignments.o -o estimateExpression

getVariance:  getVariance.cpp common.o PosteriorSamples.o ArgumentParser.o
	$(CXX) $(CXXFLAGS) getVariance.cpp common.o PosteriorSamples.o ArgumentParser.o -o getVariance

transposeLargeFile: common.o transposeFiles.o ArgumentParser.o
	$(CXX) $(CXXFLAGS) transposeLargeFile.cpp common.o transposeFiles.o ArgumentParser.o -o transposeLargeFile

GibbsSampler.o: GibbsSampler.h GibbsSampler.cpp Sampler.o GibbsParameters.o
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) -c GibbsSampler.cpp

CollapsedSampler.o: CollapsedSampler.h CollapsedSampler.cpp Sampler.o GibbsParameters.o
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) -c CollapsedSampler.cpp

Sampler.o: Sampler.cpp Sampler.h GibbsParameters.o
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) -c Sampler.cpp

ReadDistribution.o: ReadDistribution.cpp ReadDistribution.h TranscriptInfo.o TranscriptSequence.o TranscriptExpression.o 
	$(CXX) $(CXXFLAGS) $(OPENMP) -Isamtools -c ReadDistribution.cpp

asa103/asa103.o:
	make --directory asa103
TagAlignments.o: TagAlignments.cpp TagAlignments.h
transposeFiles.o: transposeFiles.cpp transposeFiles.h FileHeader.h
GibbsParameters.o: GibbsParameters.cpp GibbsParameters.h
ArgumentParser.o: ArgumentParser.cpp ArgumentParser.h
lowess.o: lowess.cpp lowess.h
TranscriptInfo.o: TranscriptInfo.cpp TranscriptInfo.h
TranscriptSequence.o: TranscriptSequence.cpp TranscriptSequence.h
TranscriptExpression.o: TranscriptExpression.cpp TranscriptExpression.h
PosteriorSamples.o: PosteriorSamples.cpp PosteriorSamples.h FileHeader.h
common.o: common.cpp common.h

clean:
	rm samtools/*.o asa103/*.o *.o $(PROGRAMS)

