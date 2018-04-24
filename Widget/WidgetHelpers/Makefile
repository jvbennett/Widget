
INC_DIR :=  -I.. -I$(shell root-config --incdir) -I../HadronWidget

SRCFILES := $(wildcard *.cc)

LIBS := $(shell root-config --glibs) -lMinuit

CXX := g++
CXX_FLAGS := -O3 $(shell root-config --cflags --libs)

.PHONY: default

default: $(LIB)

%.a: $(SRCFILES:.cc=.o)
	ar -rsv $@ $^

%.o : %.cc
	$(CXX) $(CXX_FLAGS) -M -o $*.d $< $(INC_DIR) $(LIBS); \
	cp $*.d $*.dep; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
		-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.dep; \
	rm -f $*.d $*.d.bak
	$(CXX) $(CXX_FLAGS) -c -o $@ $< $(INC_DIR) $(LIBS)

clean:
	rm *.o *.dep *.a

-include $(SRCFILES:.cc=.dep)
