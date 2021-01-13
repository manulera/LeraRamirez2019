# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#
# File src/disp/makefile.inc

OBJ_DISPLAY=gle.o gle_color.o gle_color_list.o view.o view_prop.o \
            glapp_prop.o glapp.o


#-----------------------LIB & DEF for PNG support-------------------------------

IMAGE_LIB=
IMAGE_DEF=

ifneq ($(HAS_PNG), 0)

    IMAGE_DEF += -DHAS_PNG $(INC_PNG)
    IMAGE_LIB += $(LIB_PNG)

endif


#----------------------------targets--------------------------------------------


libcytodisp.a: $(OBJ_DISPLAY) offscreen.o saveimage.o
	$(MAKELIB) build/$@ $(addprefix build/, $(notdir $^))
	$(DONE)


saveimage.o: saveimage.cc saveimage.h
	$(CXXFLAGS) -c $(IMAGE_DEF) $< -o build/$@


offscreen.o: offscreen.cc offscreen.h offscreen_fbo.cc offscreen_glx.cc
	$(CXXFLAGS) -c $< -o build/$@


$(OBJ_DISPLAY): %.o: %.cc %.h
	$(CXXFLAGS) -c $(addprefix -Isrc/, SFMT math base) $< -o build/$@
