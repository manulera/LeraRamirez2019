# Cytosim 3.1 -  Copyright Francois Nedelec et al.  EMBL 2007-2013


#---------------- MACHINE = {mac, linux-centos, embl-cluster, linux-ubuntu, windows}
MACHINE := mac

#---------------- COMPILER = {gcc, icc, clang}
COMPILER := gcc

#---------------- MODE = {F, P, D, C}
# F = Fast
# P = Profiling
# D = Debug 
# C = Coverage
MODE := F

#---------------- INTEL Math Kernel Library
# If the Intel Math Kernel library is available on your system,
# you can set:
#     HAS_MKL = 1 for dynamic linking
#     HAS_MKL = 2 for static linking
# otherwise, use HAS_MKL = 0;

HAS_MKL := 0

#---------------- PNG image export support
# `libpng` needs to be installed to save PNG images:
#     Mac OSX:
#        www.macports.org :  port install libpng 
#        Homebrew  (http://brew.sh)
#             brew install libpng
#             brew tap homebrew/dupes
#             brew install zlib
#     CentOS:  yum install libpng-devel
#     Ubuntu:  apt-get or synaptic libpng-dev
# HAS_PNG can be:
#     0 : no PNG support
#     1 : macport installation
#     2 : homebrew installation

HAS_PNG := 2

#---------------- GIF image export support
# You need to install giflib (http://giflib.sourceforge.net) to export GIF images 
#     Mac OSX:
#         www.macports.org:  port install giflib 
#         Homebrew (http://brew.sh) :  brew install giflib
#     CentOS:  yum install giflib-devel 
#     Ubuntu:  apt-get or synaptic libgif-dev 
# HAS_GIF can be:
#     0 : no GIF support
#     1 : static linking of giflib

HAS_GIF := 0

#------------------------------------------------------------------------------------
#---------------------------  Compiler Specifications  ------------------------------
#------------------------------------------------------------------------------------

DONE     = @printf "\x1B[35m >>>>>> made %s\x1B[0m\n" $@

ALIGN   := -malign-double -falign-loops -falign-jumps -falign-functions
ALIGN   :=
WARN    := -Wall -Wno-unknown-pragmas # -Wno-unused


ifeq ($(MACHINE),mac)
    MAKELIB := libtool -static -o
else
    MAKELIB := ar rcs
endif


ifeq ($(COMPILER),gcc)

    CXX      := g++ -ffinite-math-only -ftrapping-math -Wfatal-errors
    FAST     :=-O3 -ffast-math -funroll-loops #-floop-optimize2 -ftree-vectorize -msse3
    FAST     :=-O3 -fno-tree-vectorize -ffast-math -funroll-loops -msse3 -Wno-deprecated
    DEBUG    :=-g3 -ggdb
    COVERAGE :=-fprofile-arcs -ftest-coverage
    
else ifeq ($(COMPILER),icc)

    MAKELIB  := xiar -rcs
     
    #CXX     := icpc -fp-trap=common -vec-report0
    #CXX     := icpc -bind_at_load -fvisibility=hidden -no-parallel -mkl=sequential -wd1224,161
    CXX      := icpc -vec-report0 -wd1224,161,1478
    CXX      := icpc -mkl=sequential -wd1224,161,1478
    # Intel advanced instruction sets:
    # '-xSSE4.1' for core duo (ca. 2009)
    # '-mavx' for Intel core i7 (ca. 2012)
    FAST     := -O3 -xSSE4.1
    FAST     := -O3 -Wl,-no_pie -xHost -funroll-loops -fno-math-errno
    DEBUG    := -no-pie -g -wd279,383,810,869,981,1418,1419,1572,2259
    COVERAGE := -ftest-coverage

else ifeq ($(COMPILER),clang)

    CXX      := c++ -ffinite-math-only -ftrapping-math -Wfatal-errors -mmacosx-version-min=10.6
    FAST     := -O3 -ffast-math -funroll-loops #-floop-optimize2 -msse3
    FAST     := -O4 -ffast-math -funroll-loops
    DEBUG    := -g3 -ggdb
    COVERAGE := -fprofile-arcs -ftest-coverage

endif

#------------------------------------------------------------------------------------
#------------------------------------- Linux ----------------------------------------
#------------------------------------------------------------------------------------

ifeq ($(MACHINE),embl-cluster)

    MKL_PATH := /opt/intel/mkl/lib/intel64

    FlagsD := -O0 $(DEBUG) $(WARN)
    FlagsP := $(FAST) $(ALIGN) -O0 -fno-inline -pg
    FlagsC := -O0 -fprofile-arcs -ftest-coverage
    FlagsF := $(FAST) $(ALIGN) -finline-functions
    
    ### self-compiled BLAS/LAPACK fortran library:
    # this works on the new cluster (64 bit, gcc 4.1 or higher)
    #LINK  := -static -L/g/nedelec/opt/netlib/linux_em64t -llapack -lblas -lgfortran
      
    ### 64bit Linux distribution, 2011 (distribution on EMBL cluster)
    # ln -sf /usr/lib64/liblapack.so.3 /usr/lib64/liblapack.so
    # ln -sf /usr/lib64/libblas.so.3 /usr/lib64/libblas.so
    LINK   := /usr/lib64/liblapack.so.3 /usr/lib64/libblas.so.3 -L/usr/lib64 -lpthread
    
    # linkage options for programs using graphics:
    #LINKGL: =-L/usr/lib -L/usr/lib64 -lglut -lGL -lGLU -lXt -lX11

endif


ifeq ($(MACHINE),linux-centos)

    ### 64bit Centos 5.8 March 2012
    FlagsD := -O0 $(DEBUG) $(WARN)
    FlagsP := -O0 -fno-inline -pg $(ALIGN) $(FAST)
    FlagsC := -O0 -fprofile-arcs -ftest-coverage
    FlagsF := -O3 -finline-functions $(ALIGN) $(FAST)

    ### Needs
    # blas:   yum install blas-devel
    # lapack: yum install lapack-devel
    LINK = -L/usr/lib64 -llapack -lblas -lpthread
    LINK = /usr/lib64/liblapack.so.3 /usr/lib64/libblas.so.3 -lpthread
    
    ### linkage options for programs using graphics:
    ### Needs
    ### freglut: yum install freeglut-devel
    ### libXi:   yum install libXi-devel
    ### libXmu:  yum install libXmu-devel
    ###  check for correct library yourlibname.so otherwise give full path 	
    LINKGL := -L/usr/lib64 -lglut -lGL -lGLU -lXt -lX11


    ifneq ($(HAS_PNG), 0)
      
        LIB_PNG := /usr/lib64/libpng.a /usr/lib64/libz.a
        ING_PNG :=
    
    endif

    ifneq ($(HAS_GIF), 0)

    	LIB_GIF := /usr/lib64/libgif.a
        INC_GIF :=

    endif

endif


ifeq ($(MACHINE),linux-ubuntu)

    FlagsD := -O0 $(DEBUG) $(WARN)
    FlagsP := -O0 -fno-inline -pg $(ALIGN) $(FAST)
    FlagsC := -O0 -fprofile-arcs -ftest-coverage
    FlagsF := -O3 -finline-functions $(ALIGN) $(FAST)

    ### 64bit Linux distribution Ubuntu 12.04 LTS, 2012
    ## needs : apt-get install libblas.dev
    ## needs : apt-get install liblapack.dev
    LINK:= -L/usr/lib /usr/lib/libblas.so /usr/lib/liblapack.so -lpthread

    ### needs freglut: apt-get install freeglut3-dev
    ### needs libXi:   apt-get install libXi-dev
    ### needs libXmu:  apt-get install libXmu-dev
    ###  check for correct library 
    LINKGL:= -L/usr/lib -lglut -lGL -lGLU -lXt -lX11


    ifneq ($(HAS_PNG), 0)

        LIB_PNG :=/usr/lib/libpbg.a /usr/lib/libz.a
        ING_PNG :=

    endif

    ifneq ($(HAS_GIF), 0)

        LIB_GIF :=/usr/lib/libgif.a
        ING_GIF :=

    endif
endif

#------------------------------------------------------------------------------------
#-------------------------------- Mac OSX -------------------------------------------
#------------------------------------------------------------------------------------

ifeq ($(MACHINE),mac)

    MKL_PATH := /opt/intel/mkl/lib

    FlagsC := $(COVERAGE)
    FlagsD := $(DEBUG) $(WARN)
    FlagsP := $(DEBUG) $(FAST) $(ALIGN)
    FlagsF := $(FAST) $(ALIGN)
    FlagsN := -O2 
  
    LINK   := -framework Accelerate
    LINKGL := -framework GLUT -framework OpenGL -framework AGL
    #LINKGL = -L/opt/local/lib -lglut -lgl
    
    ifeq ($(HAS_PNG), 1)
        
        # macport library:
        LIB_PNG := /opt/local/lib/libpng.a  /opt/local/lib/libz.a
        INC_PNG := -I/opt/local/include

    endif
    ifeq ($(HAS_PNG), 2)
       
        # homebrew library:
        LIB_PNG := /usr/local/lib/libpng.a /usr/local/opt/zlib/lib/libz.a
        INC_PNG := -I/usr/local/include

    endif

    ifeq ($(HAS_GIF), 1)

        # macport library:
    	LIB_GIF := /opt/local/lib/libgif.a
        INC_GIF := -I/opt/local/include

    endif
    ifeq ($(HAS_GIF), 2)
    
        # homebrew library:
        LIB_GIF := /usr/local/lib/libgif.a
        INC_GIF := -I/usr/local/include

    endif

endif

#------------------------------------------------------------------------------------
#--------------------------------- Windows ------------------------------------------
#------------------------------------------------------------------------------------

ifeq ($(MACHINE),windows)

    #--- use a home-made blas/lapack for distribution:
    # the libraries "libblas.a" and "liblapack.a" are expected to be in the
    # standard location: "/usr/lib"

    FlagsD := -g  $(WARN)
    FlagsP := -pg -O3 -march=i686 -ffast-math $(ALIGN)
    FlagsF := -O3 -finline-functions -march=i686 -ffast-math $(ALIGN)
  
    LINK   := -llapack -lblas
    LINKGL := -lglut32 -lglu32 -lopengl32
    BINEXT :=.exe

endif



#------------------------------------------------------------------------------------
#--------------------- linking with Intel Math Kernel Library -----------------------
#------------------------------------------------------------------------------------

MKL_EXT := _lp64


ifeq ($(HAS_MKL),1)

    # sequential dynamic linking:
    MKL_LIBS := -lmkl_intel$(MKL_EXT) -lmkl_sequential -lmkl_core
    
    # threaded dynamic linking:
    #MKL_LIBS := -lmkl_intel$(MKL_EXT) -lmkl_intel_thread -lmkl_core -liomp5

    # modify the linking command:
    LINK := -L$(MKL_PATH) $(MKL_LIBS) -lpthread    

endif


ifeq ($(HAS_MKL),2)

    # sequential static linking:
    MKL_LIBS := $(MKL_PATH)/libmkl_intel$(MKL_EXT).a $(MKL_PATH)/libmkl_sequential.a $(MKL_PATH)/libmkl_core.a
    
    # threaded static linking:
    #MKL_LIBS := $(MKL_PATH)/libmkl_intel$(MKL_EXT).a $(MKL_PATH)/libmkl_thread.a $(MKL_PATH)/libmkl_core.a /opt/intel/lib/libiomp5.a

    # modify the linking command:
    LINK := $(MKL_LIBS) -lpthread    

endif

