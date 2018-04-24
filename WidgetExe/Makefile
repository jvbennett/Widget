
INC_DIR :=  -I.. -I$(shell root-config --incdir) -I$(WIDGET)

LIB_DIR := -L$(WIDGET)/lib

ROOTLIBS := $(shell root-config --glibs) -lMinuit

CXX := g++
CXX_FLAGS := -O3 $(shell root-config --cflags)

LIBS := -lWidget -lElectronWidget -lHadronWidget $(ROOTLIBS) -lstdc++

SRCFILES := $(wildcard *.cc)
TARGET_EXE := $(subst .cc,_exe,$(SRCFILES))

.PHONY: default

default: $(TARGET_EXE)

# generic build rule (requires all libraries)

%_exe: %.o
	$(CXX) $(CXX_FLAGS) -o $(subst _exe,,$@) $< -O3 $(INC_DIR) $(LIB_DIR) $(LIBS)

%.o : %.cc
	$(CXX) $(CXX_FLAGS) -M -o $*.d $< $(INC_DIR); \
	cp $*.d $*.dep; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.dep; \
	rm -f $*.d $*.d.bak
	$(CXX) $(CXX_FLAGS) -c -o $@ $< $(INC_DIR)

clean:
	rm -f *.o *.dep $(subst _exe,,$(TARGET_EXE))

-include $(SRCFILES:.cc=.dep)
