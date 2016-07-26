IDIR = .
ODIR = obj
SDIR = ../../SusyAnaTools/Tools/

CXX = g++
CXXFLAGS += -I. -std=c++0x -I${CMSSW_BASE}/src/ -I$(LHAPDF_DATA_PATH)/../../include
## Optimization flag
CXXFLAGS += -g
## Enable the maximun warning
#CXXFLAGS += -Wall -Wextra -Weffc++ -g

## Include ROOT
CXXFLAGS += $(shell root-config --cflags)
CXXDEPFLAGS = -MMD -MP
LD = g++
LDFLAGS =
LIBS = $(shell root-config --glibs)

LHAPDFLIB = -L$(LHAPDF_DATA_PATH)/../../lib -lLHAPDF
#OBJS = $(patsubst %, $(ODIR)/%, $(OBJ))

PROGRAMS = Template Closure Acc FakeRate IsoTrack DataMC Data Systematic SignalCont PrintSearchBins

all: sampPyWrap $(PROGRAMS)

#mkobj:
#	@mkdir -p obj

#code to compile shared library to link samples to python
sampPyWrap: $(ODIR)/samplesModule.so

$(ODIR)/samplesModule.so: $(ODIR)/samplesPyWrap.o $(ODIR)/samplesModulePyWrap.o
	$(CXX) -shared -o $@ $^

$(ODIR)/samplesPyWrap.o: $(SDIR)/samples.cc $(SDIR)/samples.h
	$(CXX) --std=c++11 -c -fPIC -o $@ $<

$(ODIR)/samplesModulePyWrap.o: $(SDIR)/samplesModule.cc
	$(CXX) --std=c++11 -c -fPIC -o $@ $<


$(ODIR)/%.o : $(SDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(SDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(IDIR)/%.C
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

$(ODIR)/%.o : $(IDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(CXXDEPFLAGS) -o $@ -c $<

Template: $(ODIR)/Template.o $(ODIR)/NTupleReader.o $(ODIR)/samples.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

Closure: $(ODIR)/Closure.o $(ODIR)/NTupleReader.o $(ODIR)/samples.o $(ODIR)/searchBins.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

Acc: $(ODIR)/Acc.o $(ODIR)/NTupleReader.o $(ODIR)/PDFUncertainty.o $(ODIR)/samples.o $(ODIR)/searchBins.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 $(LHAPDFLIB) -o $@

FakeRate: $(ODIR)/FakeRate.o $(ODIR)/NTupleReader.o $(ODIR)/samples.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

IsoTrack: $(ODIR)/IsoTrack.o $(ODIR)/NTupleReader.o $(ODIR)/samples.o $(ODIR)/searchBins.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

DataMC: $(ODIR)/DataMC.o $(ODIR)/NTupleReader.o $(ODIR)/samples.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

Data: $(ODIR)/Data.o $(ODIR)/NTupleReader.o $(ODIR)/samples.o $(ODIR)/searchBins.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

Systematic: $(ODIR)/Systematic.o $(ODIR)/NTupleReader.o $(ODIR)/samples.o $(ODIR)/searchBins.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

SignalCont: $(ODIR)/SignalCont.o $(ODIR)/NTupleReader.o $(ODIR)/samples.o $(ODIR)/searchBins.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

PrintSearchBins: $(ODIR)/PrintSearchBins.o $(ODIR)/searchBins.o
	$(LD) $^ $(LIBS) -L${CMSSW_BASE}/lib/${SCRAM_ARCH}/ -lrecipeAUXOxbridgeMT2 -o $@

clean:
	rm -f $(ODIR)/*.o $(ODIR)/*.d $(PROGRAMS) core

-include $(ODIR)/*.d