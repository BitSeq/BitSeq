CXX = g++
ARCH = -mtune=generic
VERSION = 0.7.5
#	ARCH = -march=core2
#	ARCH = -march=native


# Use O1 for debuiggging so it's not totally slow.
DBGFLAGS = -O1 -ggdb -U_FORTIFY_SOURCE
COFLAGS = $(ARCH) -O2 -pipe
CXXFLAGS = -DBS_VERSION=\"$(VERSION)\" -Wall $(COFLAGS)
# -Wvla does not work with old gcc
# -ffast-math segfaults with old gcc, don't use.
LDFLAGS = -Wl,-gc-sections
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
   transposeLargeFile \
   gtftool

all: $(PROGRAMS)

COMMON_DEPS = ArgumentParser.o common.o FileHeader.o misc.o MyTimer.o
# PROGRAMS:
convertSamples: convertSamples.cpp $(COMMON_DEPS) TranscriptInfo.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) convertSamples.cpp $(COMMON_DEPS) TranscriptInfo.o -o convertSamples

estimateDE: estimateDE.cpp $(COMMON_DEPS) PosteriorSamples.o
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) $(LDFLAGS) estimateDE.cpp $(COMMON_DEPS) PosteriorSamples.o -o estimateDE

estimateExpression: estimateExpression.cpp $(COMMON_DEPS) CollapsedSampler.o GibbsParameters.o GibbsSampler.o Sampler.o TagAlignments.o TranscriptInfo.o transposeFiles.o
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) $(OPENMP) $(LDFLAGS) estimateExpression.cpp $(COMMON_DEPS) CollapsedSampler.o GibbsParameters.o GibbsSampler.o Sampler.o TagAlignments.o TranscriptInfo.o transposeFiles.o -o estimateExpression

estimateHyperPar: estimateHyperPar.cpp $(COMMON_DEPS) lowess.o PosteriorSamples.o TranscriptExpression.o 
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) $(LDFLAGS) estimateHyperPar.cpp $(COMMON_DEPS) lowess.o PosteriorSamples.o TranscriptExpression.o -o estimateHyperPar

estimateVBExpression: estimateVBExpression.cpp $(COMMON_DEPS) SimpleSparse.o TagAlignments.o TranscriptInfo.o transposeFiles.o VariationalBayes.o
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) $(OPENMP) $(LDFLAGS) estimateVBExpression.cpp $(COMMON_DEPS) SimpleSparse.o TagAlignments.o TranscriptInfo.o transposeFiles.o VariationalBayes.o -o estimateVBExpression

extractSamples: extractSamples.cpp $(COMMON_DEPS) PosteriorSamples.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) extractSamples.cpp $(COMMON_DEPS) PosteriorSamples.o -o extractSamples

getFoldChange: getFoldChange.cpp $(COMMON_DEPS) PosteriorSamples.o 
	$(CXX) $(CXXFLAGS) $(LDFLAGS) getFoldChange.cpp $(COMMON_DEPS) PosteriorSamples.o -o getFoldChange

getGeneExpression: getGeneExpression.cpp $(COMMON_DEPS) PosteriorSamples.o TranscriptInfo.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) getGeneExpression.cpp $(COMMON_DEPS) PosteriorSamples.o TranscriptInfo.o -o getGeneExpression

getPPLR: getPPLR.cpp $(COMMON_DEPS) PosteriorSamples.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) getPPLR.cpp $(COMMON_DEPS) PosteriorSamples.o -o getPPLR

getVariance: getVariance.cpp $(COMMON_DEPS) PosteriorSamples.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) getVariance.cpp $(COMMON_DEPS) PosteriorSamples.o -o getVariance

getWithinGeneExpression: getWithinGeneExpression.cpp $(COMMON_DEPS) PosteriorSamples.o TranscriptInfo.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) getWithinGeneExpression.cpp $(COMMON_DEPS) PosteriorSamples.o TranscriptInfo.o -o getWithinGeneExpression

parseAlignment: parseAlignment.cpp $(COMMON_DEPS) ReadDistribution.o samtools/sam.o TranscriptExpression.o TranscriptInfo.o TranscriptSequence.o
	$(CXX) $(CXXFLAGS) $(OPENMP) $(LDFLAGS) -pthread parseAlignment.cpp $(COMMON_DEPS) ReadDistribution.o samtools/*.o TranscriptExpression.o TranscriptInfo.o TranscriptSequence.o -lz -o parseAlignment

transposeLargeFile: transposeLargeFile.cpp $(COMMON_DEPS) transposeFiles.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) transposeLargeFile.cpp $(COMMON_DEPS) transposeFiles.o -o transposeLargeFile

gtftool: gtftool.cpp $(COMMON_DEPS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) gtftool.cpp $(COMMON_DEPS) -o gtftool

# LIBRARIES:
ArgumentParser.o: ArgumentParser.cpp ArgumentParser.h
	$(CXX) $(CXXFLAGS) -ffunction-sections -fdata-sections -c ArgumentParser.cpp

CollapsedSampler.o: CollapsedSampler.cpp CollapsedSampler.h GibbsParameters.h Sampler.h
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) -c CollapsedSampler.cpp

FileHeader.o: common.h misc.h FileHeader.cpp FileHeader.h
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) -ffunction-sections -fdata-sections -c FileHeader.cpp

GibbsSampler.o: GibbsSampler.cpp GibbsSampler.h GibbsParameters.h Sampler.h
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) -c GibbsSampler.cpp

misc.o: ArgumentParser.h PosteriorSamples.h misc.cpp misc.h
	$(CXX) $(CXXFLAGS) -ffunction-sections -fdata-sections -c misc.cpp

MyTimer.o: MyTimer.h MyTimer.cpp
	$(CXX) $(CXXFLAGS) -ffunction-sections -fdata-sections -c MyTimer.cpp

PosteriorSamples.o: PosteriorSamples.cpp PosteriorSamples.h FileHeader.h
	$(CXX) $(CXXFLAGS) -ffunction-sections -fdata-sections -c PosteriorSamples.cpp

ReadDistribution.o: ReadDistribution.cpp ReadDistribution.h TranscriptExpression.h TranscriptInfo.h TranscriptSequence.h 
	$(CXX) $(CXXFLAGS) $(OPENMP) -c ReadDistribution.cpp

Sampler.o: Sampler.cpp Sampler.h GibbsParameters.h
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) -c Sampler.cpp

SimpleSparse.o: SimpleSparse.cpp SimpleSparse.h
	$(CXX) $(CXXFLAGS) $(OPENMP) -c SimpleSparse.cpp

VariationalBayes.o: VariationalBayes.cpp VariationalBayes.h SimpleSparse.h
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) $(OPENMP) -c VariationalBayes.cpp 

common.o: common.cpp common.h
GibbsParameters.o: ArgumentParser.h GibbsParameters.cpp GibbsParameters.h
lowess.o: lowess.cpp lowess.h
TagAlignments.o: TagAlignments.cpp TagAlignments.h
TranscriptExpression.o: TranscriptExpression.cpp TranscriptExpression.h
TranscriptInfo.o: TranscriptInfo.cpp TranscriptInfo.h
TranscriptSequence.o: TranscriptSequence.cpp TranscriptSequence.h
transposeFiles.o: transposeFiles.cpp transposeFiles.h FileHeader.h

# EXTERNAL LIBRARIES:
samtools/sam.o:
	make --directory samtools

# CLEAN:
clean:
	rm *.o $(PROGRAMS)

clean-all:
	rm samtools/*.o *.o $(PROGRAMS)

