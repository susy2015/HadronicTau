IDIR       = .
ODIR       = obj
SDIR       = $(CMSSW_BASE)/src/SusyAnaTools/Tools


OPENCV_DIRECTORY=$(CMSSW_BASE)/src/opencv

CXX        = g++

CXXFLAGS  += -I. -I$(CMSSW_BASE)/src -I$(OPENCV_DIRECTORY)/include/ -I$(OPENCV_DIRECTORY)/modules/core/include -I$(OPENCV_DIRECTORY)/modules/video/include -I$(OPENCV_DIRECTORY)/modules/objdetect/include -I$(OPENCV_DIRECTORY)/modules/ml/include/ -I$(OPENCV_DIRECTORY)/modules/photo/include/ -I$(OPENCV_DIRECTORY)/modules/imgproc/include/ -std=c++0x

## Optimization flag
CXXFLAGS += -g
## Enable the maximun warning
#CXXFLAGS += -Wall -Wextra -Weffc++ -g

## Include ROOT
CXXFLAGS  += $(shell root-config --cflags)

CXXDEPFLAGS = -MMD -MP

LD         = g++
LDFLAGS    =

ROOTLIBS   = $(shell root-config --glibs)
MT2LIB     = -L$(CMSSW_BASE)/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2
TOPLIB     = -L$(CMSSW_BASE)/lib/${SCRAM_ARCH}/ -lTopTaggerTopTagger -lTopTaggerCfgParser
OPENVCLIBS = -L$(OPENCV_DIRECTORY)/lib/ -lopencv_ml -lopencv_core

#OBJS       = $(patsubst %, $(ODIR)/%, $(OBJ))

PROGRAMS = CS HadTauLL

all: mkobj sampPyWrap $(PROGRAMS)

mkobj:
	@mkdir -p obj

#code to compile shared library to link samples to python                                                                                                                               
sampPyWrap: $(ODIR)/samplesModule.so

$(ODIR)/samplesModule.so: $(ODIR)/samplesPyWrap.o $(ODIR)/samplesModulePyWrap.o
	$(CXX) -shared -o $@ $^

$(ODIR)/samplesPyWrap.o: $(SDIR)/samples.cc $(SDIR)/samples.h
	$(CXX) --std=c++11 -c -fPIC -o $@ $<

$(ODIR)/samplesModulePyWrap.o: $(SDIR)/samplesModule.cc
	$(CXX) --std=c++11 -c -fPIC -o $@ $<


$(ODIR)/%.o : $(IDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(IDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<


CS: $(ODIR)/CS.o $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/samples.o $(ODIR)/searchBins.o $(ODIR)/customize.o $(ODIR)/baselineDef.o $(ODIR)/BTagCorrector.o $(ODIR)/BTagCalibrationStandalone.o $(ODIR)/ISRCorrector.o
	$(LD) $^ $(TOPLIB) $(OPENVCLIBS) $(MT2LIB) $(ROOTLIBS) -o $@

HadTauLL: $(ODIR)/HadTauLL.o $(ODIR)/NTupleReader.o $(ODIR)/SATException.o $(ODIR)/samples.o $(ODIR)/searchBins.o $(ODIR)/customize.o $(ODIR)/baselineDef.o $(ODIR)/BTagCorrector.o $(ODIR)/BTagCalibrationStandalone.o
	$(LD) $^ $(TOPLIB) $(OPENVCLIBS) $(MT2LIB) $(ROOTLIBS) -o $@


clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.so $(ODIR)/*.d $(PROGRAMS) core 

-include $(ODIR)/*.d
