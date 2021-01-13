# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.



TESTS:=test test_gillespie test_solve test_random test_math test_param test_quaternion\
       test_code test_matrix test_thread test_blas


TESTS_GL:=test_opengl test_vbo test_glut test_glapp test_platonic\
          test_rasterizer test_space test_grid test_sphere


.PHONY: tests
tests: $(TESTS) $(TESTS_GL)


vpath %.cc src/test


#--------------------macros----------------------------------------------------


TEST_OBJ=$(filter %.cc, $^) $(addprefix build/, $(notdir $(filter-out %.cc, $^)))

TEST_INC=$(addprefix -Isrc/, SFMT math base sim disp play)

TEST_MAKE=$(CXXFLAGS) $(TEST_INC) $(TEST_OBJ) $(LINK) -o bin/$@

TEST_MAKE_GL=$(CXXFLAGS) $(TEST_INC) -I src/disp $(TEST_OBJ) $(LINK) $(LINKGL) $(IMAGE_LIB) -o bin/$@


#----------------------------Targets--------------------------------------------

test: test.cc random.o SFMT.o tictoc.o backtrace.o
	$(TEST_MAKE)
	$(DONE)
vpath test bin

test_blas: test_blas.cc random.o SFMT.o backtrace.o
	$(TEST_MAKE)
	$(DONE)
vpath test_blas bin

test_code: test_code.cc random.o SFMT.o tictoc.o 
	$(CXXFLAGS) $(INC_SIM) $(TEST_OBJ) $(LINK) -o bin/$@
	$(DONE)
vpath test_code bin

test_matrix: test_matrix.cc libcytomath.a libcytobase.a 
	$(TEST_MAKE)
	$(DONE)
vpath test_matrix bin

test_param: test_param.cc glossary.o filepath.o tokenizer.o stream_func.o exceptions.o backtrace.o ansi_colors.o
	$(TEST_MAKE)
	$(DONE)
vpath test_param bin

test_quartic: test_quartic.cc quartic_solver.o random.o SFMT.o backtrace.o
	$(TEST_MAKE)
	$(DONE)
vpath test_quartic bin

test_quaternion: test_quaternion.cc libcytomath.a SFMT.o libcytobase.a 
	$(TEST_MAKE)
	$(DONE)
vpath test_quaternion bin

test_random: test_random.cc random.o SFMT.o exceptions.o filewrapper.o messages.o smath.o backtrace.o tictoc.o
	$(TEST_MAKE)
	$(DONE)
vpath test_random bin

test_simd: test_simd.cc random.o SFMT.o tictoc.o 
	$(CXXFLAGS) -Isrc/base -Isrc/math -Isrc/SFMT $^ -o bin/test_simd
	$(DONE)
vpath test_simd bin

test_sizeof: test_sizeof.cc
	$(TEST_MAKE)
	$(DONE)
vpath test_sizeof bin

test_cxx: test_cxx.cc
	$(TEST_MAKE)
	$(DONE)
vpath test_cxx bin

test_math: test_math.cc
	$(TEST_MAKE)
	$(DONE)
vpath test_math bin

test_thread: test_thread.cc
	$(TEST_MAKE)
	$(DONE)
vpath test_thread bin

test_string: test_string.cc
	$(TEST_MAKE)
	$(DONE)
vpath test_string bin

test_gillespie: test_gillespie.cc random.o SFMT.o backtrace.o
	$(TEST_MAKE)
	$(DONE)
vpath test_gillespie bin

#----------------------------Targets with Graphics------------------------------

test_opengl: test_opengl.cc
	$(CXX) $^ $(LINKGL) -o bin/$@
	$(DONE)
vpath test_opengl bin

test_vbo: test_vbo.cc
	$(CXX) $^ $(LINKGL) -o bin/$@
	$(DONE)
vpath test_vbo bin

test_glut: test_glut.cc
	$(TEST_MAKE_GL)
	$(DONE)
vpath test_glut bin

test_glapp: test_glapp.cc libcytodisp.a libcytobase.a libcytomath.a
	$(TEST_MAKE_GL)
	$(DONE)
vpath test_glapp bin

test_sphere: test_sphere.cc libcytodisp.a libcytomath.a libcytobase.a
	$(TEST_MAKE_GL)
	$(DONE)
vpath test_sphere bin

test_platonic: test_platonic.cc platonic.o libcytodisp.a libcytomath.a libcytobase.a
	$(TEST_MAKE_GL)
	$(DONE)
vpath test_platonic bin

test_rasterizer: test_rasterizer.cc rasterizerGL.o libcytomath.a libcytodisp.a libcytobase.a 
	$(TEST_MAKE_GL)
	$(DONE)
vpath test_rasterizer bin

test_space: test_space.cc libcytospace.a libcytosim.a libcytomath.a libcytodisp.a libcytobase.a 
	$(TEST_MAKE_GL)
	$(DONE)
vpath test_space bin

test_grid: test_grid.cc tictoc.o libcytomath.a libcytobase.a libcytodisp.a
	$(TEST_MAKE_GL)
	$(DONE)
vpath test_grid bin

test_solve: test_solve.cc libcytomath.a libcytobase.a
	$(TEST_MAKE_GL)
	$(DONE)
vpath test_solve bin

