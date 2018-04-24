
SRC_DIRS := HadronWidget ElectronWidget WidgetHelpers

TARGET_LIBS := $(addsuffix .a, $(addprefix lib, $(SRC_DIRS)))

.PHONY: default clean

default: lib libWidget.a

lib:
	mkdir lib

libWidget.a: $(TARGET_LIBS)
	$(foreach lib, $(TARGET_LIBS), $(shell cd lib; ar -x $(lib) ) )
	@cd lib && ar -rsv $@ *.o && ranlib $@
	@cd lib && rm -f *.o

lib%.a:
	@$(MAKE) -C $(subst lib,, $(subst .a,, $@ )) LIB=$@
	@cp $(subst lib,, $(subst .a,, $@))/$@ lib/

clean: $(addprefix clean_, $(SRC_DIRS))
	-rm -f lib/*.a

clean_%:
	@-cd $(subst clean_,, $@) && $(MAKE) clean
