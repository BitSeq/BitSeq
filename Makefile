CXX = g++
ARCH = -mtune=generic
VERSION = 0.7.9
#	ARCH = -march=core2
#	ARCH = -march=native

UNAME := $(shell sh -c 'uname -s 2>/dev/null || echo not')

# Use O1 for debuiggging so it's not totally slow.
DBGFLAGS = -O1 -ggdb -U_FORTIFY_SOURCE
COFLAGS = $(ARCH) -O2 -pipe
CXXFLAGS = -DBS_VERSION=\"$(VERSION)\" -Wall $(COFLAGS)
# -Wvla does not work with old gcc
# -ffast-math segfaults with old gcc, don't use.

BOOSTFLAGS = -I .

# Modern Mac OS X and Clang do not support OpenMP
ifeq ($(UNAME), Darwin)
OPENMP =
else
LDFLAGS = -Wl,-gc-sections
OPENMP = -fopenmp -DSUPPORT_OPENMP
endif

PROGRAMS = bitseq

all: $(PROGRAMS)

COMMON_DEPS = ArgumentParser.o common.o FileHeader.o misc.o MyTimer.o TranscriptInfo.o PosteriorSamples.o

HELPER_OBJECTS = SimpleSparse.o CollapsedSampler.o GibbsParameters.o GibbsSampler.o lowess.o LikelihoodWriter.o ReadDistribution.o Sampler.o TagAlignments.o TranscriptExpression.o TranscriptSequence.o transposeFiles.o VariationalBayes.o

COMMAND_OBJECTS = estimateExpression.o estimateVBExpression.o parseAlignment.o estimateDE.o estimateFCProb.o estimateHyperPar.o getGeneExpression.o getVariance.o getWithinGeneExpression.o getFoldChange.o getPPLR.o convertSamples.o extractSamples.o transposeLargeFile.o gtftool.o

# PROGRAMS:
bitseq: bitseq.cpp $(COMMON_DEPS) $(COMMAND_OBJECTS) $(HELPER_OBJECTS) samtools/sam.o
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) $(OPENMP) $(LDFLAGS) -pthread bitseq.cpp $(COMMON_DEPS) $(COMMAND_OBJECTS) $(HELPER_OBJECTS) samtools/*.o -lz -o bitseq

# default compilation
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) -c $<

# COMMANDS:
estimateExpression.o: estimateExpression.cpp
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) $(OPENMP) -c estimateExpression.cpp

estimateVBExpression.o: estimateVBExpression.cpp
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) $(OPENMP) -c estimateVBExpression.cpp

parseAlignment.o: parseAlignment.cpp
	$(CXX) $(CXXFLAGS) $(BOOSTFLAGS) $(OPENMP) -c parseAlignment.cpp

# LIBRARIES:
ArgumentParser.o: ArgumentParser.cpp ArgumentParser.h
	$(CXX) $(CXXFLAGS) -ffunction-sections -fdata-sections -c ArgumentParser.cpp

CollapsedSampler.o: CollapsedSampler.cpp CollapsedSampler.h GibbsParameters.h Sampler.h

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
LikelihoodWriter.o: LikelihoodWriter.cpp LikelihoodWriter.h
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

