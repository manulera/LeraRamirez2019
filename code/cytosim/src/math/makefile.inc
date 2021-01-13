# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.



OBJ_MATH:=smath.o vector1.o vector2.o vector3.o matrix1.o matrix2.o matrix3.o\
    	rasterizer.o project_ellipse.o platonic.o matrix.o matsparse.o matsparsesym.o\
    	matsym.o matblock.o matsparsesym1.o matsparsesym2.o\
    	matsparsesym4.o matsparsesym5.o bicgstab.o polygon.o\
    	conjgradient.o pointsonsphere.o random.o random_vector.o quartic_solver.o\


#----------------------------rules----------------------------------------------

INC_MATH:=$(addprefix -Isrc/, SFMT base math)


$(OBJ_MATH): %.o: %.cc %.h
	$(CXXFLAGS) -c $(INC_MATH) $< -o build/$@


matsparsesymblk.o: matsparsesymblk.cc matsparsesymblk.h
	$(CXXFLAGS) -c $(INC_MATH) $< -o build/$@

matsparsesymblk1D.o: matsparsesymblk.cc matsparsesymblk.h
	$(CXXFLAGS) -c $(INC_MATH) -DBSZ=1 $< -o build/$@

matsparsesymblk2D.o: matsparsesymblk.cc matsparsesymblk.h
	$(CXXFLAGS) -c $(INC_MATH) -DBSZ=2 $< -o build/$@

matsparsesymblk3D.o: matsparsesymblk.cc matsparsesymblk.h
	$(CXXFLAGS) -c $(INC_MATH) -DBSZ=3 $< -o build/$@


libcytomath.a: $(OBJ_MATH) matsparsesymblk.o
	$(MAKELIB) build/$@ $(addprefix build/, $(notdir $^))
	$(DONE)

libcytomath1D.a: $(OBJ_MATH) matsparsesymblk1D.o
	$(MAKELIB) build/$@ $(addprefix build/, $(notdir $^))
	$(DONE)

libcytomath2D.a: $(OBJ_MATH) matsparsesymblk2D.o
	$(MAKELIB) build/$@ $(addprefix build/, $(notdir $^))
	$(DONE)

libcytomath3D.a: $(OBJ_MATH) matsparsesymblk3D.o
	$(MAKELIB) build/$@ $(addprefix build/, $(notdir $^))
	$(DONE)


# build the rasterizer with graphical debugging enabled:
rasterizerGL.o: rasterizer.cc rasterizer.h
	$(CXXFLAGS) -c $(INC_MATH) -Isrc/disp -DDISPLAY $< -o build/$@

# build the grid with display functions enabled:
gridGL.o: grid.cc grid.h
	$(CXXFLAGS) -c $(INC_MATH) -Isrc/disp -DDISPLAY $< -o build/$@


