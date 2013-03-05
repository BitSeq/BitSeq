CXX = g++
HOSTNAME = $(shell hostname)
ARCH = -mtune=generic

ifeq ($(HOSTNAME), valiant)
	ARCH = -march=core2
endif
ifeq ($(HOSTNAME), gauss)
	ARCH = -march=native
endif
ifeq ($(HOSTNAME), rpc465.cs.man.ac.uk)
	ARCH = -march=native
endif
ifeq ($(HOSTNAME), chopok)
	ARCH = -march=native
endif

DBGFLAGS = -ggdb -U_FORTIFY_SOURCE
COFLAGS = $(ARCH) -O3 -pipe
# -ffast-math segfaults with old gcc
#COFLAGS = $(ARCH) -O3 -ffast-math -pipe
CXXFLAGS = -Wall $(COFLAGS)
BOOSTFLAGS = -I .
OPENMP = -fopenmp -DSUPPORT_OPENMP

PROGRAMS = \
   convertSamples \
   estimateDE \
   estimateExpression \
   estimateHyperPar \
   estimateVBExpression \
   extractSamples \
   getFoldChange \
   getGeneExpression \
   getPPLR \
   getVariance \
   getWithinGeneExpression \
   parseAlignment \
   transposeLargeFile

all: $(PROGRAMS)

# PROGRAMS:
convertSamples: convertSamples.cpp ArgumentParser.o common.o TranscriptInfo.o
	$(CXX) $(CXXFLAGS) convertSamples.cpp ArgumentParser.o common.o TranscriptInfo.o -o convertSamples

estimateDE: estimateDE.cpp ArgumentParser.o common.o misc.o PosteriorSamples.o
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) estimateDE.cpp ArgumentParser.o common.o misc.o PosteriorSamples.o -o estimateDE

estimateExpression: estimateExpression.cpp ArgumentParser.o CollapsedSampler.o common.o GibbsParameters.o GibbsSampler.o Sampler.o TagAlignments.o TranscriptInfo.o transposeFiles.o
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) $(OPENMP) estimateExpression.cpp ArgumentParser.o CollapsedSampler.o common.o GibbsParameters.o GibbsSampler.o Sampler.o TagAlignments.o TranscriptInfo.o transposeFiles.o -o estimateExpression

estimateHyperPar: estimateHyperPar.cpp ArgumentParser.o common.o lowess.o misc.o PosteriorSamples.o TranscriptExpression.o 
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) estimateHyperPar.cpp ArgumentParser.o common.o lowess.o misc.o PosteriorSamples.o TranscriptExpression.o -o estimateHyperPar

estimateVBExpression: estimateVBExpression.cpp asa103/asa103.o ArgumentParser.o common.o SimpleSparse.o TagAlignments.o VariationalBayes.o
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) $(OPENMP) estimateVBExpression.cpp asa103/asa103.o ArgumentParser.o common.o SimpleSparse.o TagAlignments.o VariationalBayes.o -o estimateVBExpression

extractSamples: extractSamples.cpp ArgumentParser.o common.o PosteriorSamples.o
	$(CXX) $(CXXFLAGS) extractSamples.cpp ArgumentParser.o common.o PosteriorSamples.o -o extractSamples

getFoldChange: getFoldChange.cpp ArgumentParser.o common.o PosteriorSamples.o 
	$(CXX) $(CXXFLAGS) getFoldChange.cpp ArgumentParser.o common.o PosteriorSamples.o -o getFoldChange

getGeneExpression: getGeneExpression.cpp ArgumentParser.o common.o PosteriorSamples.o TranscriptInfo.o
	$(CXX) $(CXXFLAGS) getGeneExpression.cpp ArgumentParser.o common.o PosteriorSamples.o TranscriptInfo.o -o getGeneExpression

getPPLR: getPPLR.cpp ArgumentParser.o common.o PosteriorSamples.o
	$(CXX) $(CXXFLAGS) getPPLR.cpp ArgumentParser.o common.o PosteriorSamples.o -o getPPLR

getVariance: getVariance.cpp ArgumentParser.o common.o PosteriorSamples.o
	$(CXX) $(CXXFLAGS) getVariance.cpp ArgumentParser.o common.o PosteriorSamples.o -o getVariance

getWithinGeneExpression: getWithinGeneExpression.cpp ArgumentParser.o common.o PosteriorSamples.o TranscriptInfo.o
	$(CXX) $(CXXFLAGS) getWithinGeneExpression.cpp ArgumentParser.o common.o PosteriorSamples.o TranscriptInfo.o -o getWithinGeneExpression

parseAlignment: parseAlignment.cpp ArgumentParser.o common.o ReadDistribution.o samtools/sam.o TranscriptExpression.o TranscriptInfo.o TranscriptSequence.o
	$(CXX) $(CXXFLAGS) $(OPENMP) -Isamtools parseAlignment.cpp ArgumentParser.o common.o ReadDistribution.o samtools/*.o TranscriptExpression.o TranscriptInfo.o TranscriptSequence.o -lz -o parseAlignment

transposeLargeFile: transposeLargeFile.cpp ArgumentParser.o common.o transposeFiles.o
	$(CXX) $(CXXFLAGS) transposeLargeFile.cpp ArgumentParser.o common.o transposeFiles.o -o transposeLargeFile

# LIBRARIES:
CollapsedSampler.o: CollapsedSampler.cpp CollapsedSampler.h GibbsParameters.h Sampler.h
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) -c CollapsedSampler.cpp

GibbsSampler.o: GibbsSampler.cpp GibbsSampler.h GibbsParameters.h Sampler.h
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) -c GibbsSampler.cpp

ReadDistribution.o: ReadDistribution.cpp ReadDistribution.h TranscriptExpression.h TranscriptInfo.h TranscriptSequence.h 
	$(CXX) $(CXXFLAGS) $(OPENMP) -Isamtools -c ReadDistribution.cpp

Sampler.o: Sampler.cpp Sampler.h GibbsParameters.h
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) -c Sampler.cpp

SimpleSparse.o: SimpleSparse.cpp SimpleSparse.h
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) $(OPENMP) -c SimpleSparse.cpp

VariationalBayes.o: VariationalBayes.cpp VariationalBayes.h SimpleSparse.h
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) $(OPENMP) -c VariationalBayes.cpp 

ArgumentParser.o: ArgumentParser.cpp ArgumentParser.h
common.o: common.cpp common.h
GibbsParameters.o: ArgumentParser.h GibbsParameters.cpp GibbsParameters.h
lowess.o: lowess.cpp lowess.h
misc.o: ArgumentParser.h PosteriorSamples.h misc.cpp misc.h
PosteriorSamples.o: PosteriorSamples.cpp PosteriorSamples.h FileHeader.h
TagAlignments.o: TagAlignments.cpp TagAlignments.h
TranscriptExpression.o: TranscriptExpression.cpp TranscriptExpression.h
TranscriptInfo.o: TranscriptInfo.cpp TranscriptInfo.h
TranscriptSequence.o: TranscriptSequence.cpp TranscriptSequence.h
transposeFiles.o: transposeFiles.cpp transposeFiles.h FileHeader.h

# EXTERNAL LIBRARIES:
asa103/asa103.o:
	make --directory asa103

samtools/sam.o:
	make --directory samtools

# CLEAN:
clean:
	rm samtools/*.o asa103/*.o *.o $(PROGRAMS)

