# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.




TOOLS:=frametool sieve reader report reportF analyse analyse2 analyse3


.PHONY: tools
tools: $(TOOLS)


vpath %.cc src/tools


#--------------------macros----------------------------------------------------

TOOL_DEP := libcytosim.a libcytospace.a libcytomath.a libcytobase.a

TOOL_OBJ  = $(filter %.cc, $^) $(addprefix build/, $(notdir $(filter-out %.cc, $^)))

TOOL_INC  = $(addprefix -Isrc/, SFMT math base sim sim/spaces sim/fibers sim/singles sim/organizers disp play)

TOOL_MAKE = $(CXXFLAGS) $(TOOL_INC) $(TOOL_OBJ) $(LINK) -o bin/$@


#-------------------targets----------------------------------------------------
 
 
frametool: frametool.cc
	$(TOOL_MAKE)
	$(DONE)
vpath frametool bin


sieve: sieve.cc frame_reader.o $(TOOL_DEP) 
	$(TOOL_MAKE)
	$(DONE)
vpath sieve bin


reader: reader.cc frame_reader.o $(TOOL_DEP)
	$(TOOL_MAKE)
	$(DONE)
vpath reader bin



report: report.cc frame_reader.o $(TOOL_DEP) | bin
	$(CXXFLAGS) $(TOOL_INC) $(TOOL_OBJ) $(LINK) -o bin/$@
	$(DONE)
vpath report bin


bin1/report: report.cc frame_reader.o libcytosim1D.a libcytospace1D.a libcytomath1D.a libcytobase.a
	$(CXXFLAGS) -DDIM=1 $(TOOL_INC) $(TOOL_OBJ) $(LINK) -o bin1/report
	$(DONE)
vpath bin1/report bin1


bin2/report: report.cc frame_reader.o libcytosim2D.a libcytospace2D.a libcytomath2D.a libcytobase.a
	$(CXXFLAGS) -DDIM=2 $(TOOL_INC) $(TOOL_OBJ) $(LINK) -o bin2/report
	$(DONE)
vpath bin2/report bin2


bin3/report: report.cc frame_reader.o libcytosim3D.a libcytospace3D.a libcytomath3D.a libcytobase.a
	$(CXXFLAGS) -DDIM=3 $(TOOL_INC) $(TOOL_OBJ) $(LINK) -o bin3/report
	$(DONE)
vpath bin3/report bin3




reportF: reportF.cc frame_reader.o $(TOOL_DEP)
	$(TOOL_MAKE)
	$(DONE)
vpath reportF bin

analyse: analyse.cc frame_reader.o $(TOOL_DEP)
	$(TOOL_MAKE)
	$(DONE)
vpath analyse bin

analyse2: analyse2.cc frame_reader.o $(TOOL_DEP)
	$(TOOL_MAKE)
	$(DONE)
vpath analyse2 bin

analyse3: analyse3.cc frame_reader.o $(TOOL_DEP)
	$(TOOL_MAKE)
	$(DONE)
vpath analyse3 bin


cymart: cymart.cc frame_reader.o $(TOOL_DEP)
	$(TOOL_MAKE)
	$(DONE)
vpath cymart bin


micimage: micimage.cc saveimage.cc glossary.o filepath.o tokenizer.o
micimage: stream_func.o exceptions.o random.o SFMT.o ansi_colors.o backtrace.o
	$(CXXFLAGS) $(addprefix -Isrc/, SFMT base math disp) -DNO_OPENGL -DHAS_PNG $(TOOL_OBJ) $(LIB_PNG) -o bin/$@
	$(DONE)
vpath micimage bin


#----------------------------makedep--------------------------------------------

dep/part7.dep:
	$(MAKEDEP) $(wildcard src/tools/*.cc) > $@

