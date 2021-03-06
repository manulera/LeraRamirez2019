# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

# THIS FILE SHOULD NOT BE EDITED
# ONLY EDIT FILE makefile.inc 

#include the compiler-specifications
include makefile.inc

#-----------------check the compilers and the compiler-flags--------------------

ifndef CXX
   $(error No compiler defined for $$(MACHINE)=$(MACHINE))
endif

ifndef Flags$(MODE)
   $(warning No compiler options defined for $$(MACHINE)=$(MACHINE) in mode $(MODE))
endif

ifndef LINK
  $(error No linkage-options defined for $$(MACHINE)=$(MACHINE))
endif

#command to invoke compiler:
CXXFLAGS := $(CXX) $(Flags$(MODE))

SRCDIR1 := $(addprefix src/, math base sim disp play)
SRCDIR2 := $(addprefix src/sim/, spaces hands fibers singles couples organizers)
SRCDIR  := $(SRCDIR1) $(SRCDIR2)


#command used to build the dependencies files automatically
MAKEDEP := gcc -MM $(addprefix -I, $(SRCDIR) src/SFMT)

#-----------------------SVN revision number-------------------------------------

CODE_VERSION = -D'CODE_VERSION="$(shell svnversion -n . || echo unknown)"'

COMPILER_VERSION = -D'COMPILER_VERSION="$(shell $(word 1, $(CXX)) --version | head -1)"'

#----------------------------paths----------------------------------------------

vpath %.h     $(SRCDIR)
vpath %.cc    $(SRCDIR)
vpath %.o     build/
vpath lib%.a  build/
vpath %.dep   dep/
vpath SFMT%   src/SFMT/

#----------------------------targets--------------------------------------------

# calling 'make' without arguments will make sim and play:

.PHONY: simplay
simplay: sim report play


include src/sim/makefile.inc
include src/base/makefile.inc
include src/math/makefile.inc
include src/sim/spaces/makefile.inc
include src/disp/makefile.inc
include src/play/makefile.inc

-include src/tools/makefile.inc
-include src/test/makefile.inc


# Attention: Mersenne-Twister is coded in C-language,
# and we must tell the compiler with '-x c':
SFMT.o: SFMT.c SFMT.h
	$(CXXFLAGS) -x c -DNDEBUG -DSFMT_MEXP=19937 -c $< -o build/$@


.PHONY: all bin1 bin2 bin3 alldim doc

all: sim play tools tests

alldim: bin1 bin2 bin3

bin1:
	if ! test -d bin1; then mkdir bin1; fi
	make bin1/sim bin1/play bin1/report

bin2:
	if ! test -d bin2; then mkdir bin2; fi
	make bin2/sim bin2/play bin2/report

bin3:
	if ! test -d bin3; then mkdir bin3; fi
	make bin3/sim bin3/play bin3/report

doc:
	if test -d html; then rm -rf html; fi
	mkdir html;
	doxygen src/doc/doxygen.cfg > log.txt 2> /dev/null

#------------------------------- archive ---------------------------------------
.PHONY: cytosim.tar cytosim_src.tar tar pack


cytosim.tar:
	tar cf cytosim.tar --exclude \*.cmo --exclude \*tar\* \
        --exclude bin/\* --exclude build/\*  *


cytosim.tar.bz2: cytosim.tar
	rm -f cytosim.tar.bz2;
	bzip2 cytosim.tar;


cytosim_src.tar:
	if ! test -f config.cym; then cp cym/fiber.cym config.cym; fi
	COPYFILE_DISABLE=1 tar cf cytosim_src.tar --exclude "*.o" --exclude ".*" \
	--exclude ".svn" --exclude "*~" src makefile makefile.inc *.cym cym matlab python


tar: cytosim_src.tar cytosim.tar.bz2


pack: sterile tar

#---------------------------- maintenance --------------------------------------
.PHONY: bin build clean cleaner sterile

bin:
	if ! test -d bin; then mkdir bin; fi

build:
	if ! test -d build; then mkdir build; fi


clean:
	rm -f build/*.o build/*.a


cleaner:
	rm -f *.cmo build/*.o build/*.a dep/*.dep;


sterile:
	rm -rf build/*
	rm -f  dep/*.dep
	rm -rf bin/*.dSYM
	rm -f  bin/*
	rm -f  bin1/*
	rm -f  bin2/*
	rm -f  bin3/*
	rm -f *.cmo
	rm -f log.txt;
	svn cleanup


#---------------------------- dependencies -------------------------------------
.PHONY: dep

dep:
	if ! test -d dep; then mkdir dep; else rm -f dep/*; fi
	$(foreach file, $(wildcard src/base/*.cc),  $(MAKEDEP) $(file) >> dep/part0.dep; )
	$(foreach file, $(wildcard src/math/*.cc),  $(MAKEDEP) $(file) >> dep/part1.dep; )
	$(foreach file, $(wildcard src/sim/*.cc),   $(MAKEDEP) $(file) >> dep/part2.dep; )
	$(foreach file, $(wildcard src/sim/*/*.cc), $(MAKEDEP) $(file) >> dep/part3.dep; )
	$(foreach file, $(wildcard src/disp/*.cc),  $(MAKEDEP) $(file) >> dep/part4.dep; )
	$(foreach file, $(wildcard src/play/*.cc),  $(MAKEDEP) $(file) >> dep/part5.dep; )
	$(foreach file, $(wildcard src/tools/*.cc), $(MAKEDEP) $(file) >> dep/part6.dep; )
	$(foreach file, $(wildcard src/test/*.cc),  $(MAKEDEP) $(file) >> dep/part7.dep; )


-include dep/part?.dep



