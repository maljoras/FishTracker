#!/usr/bin/make -f
# ============================================================================
#                              mexopencv Makefile
# ============================================================================
#
# The following configuration parameters are recognized:
#
# MATLABDIR             MATLAB root directory.
# MATLAB                MATLAB executable.
# MEX                   MATLAB MEX compiler frontend.
# MEXEXT                extension of MEX-files.
# DOXYGEN               Doxygen executable used to generate documentation.
# NO_CV_PKGCONFIG_HACK  Boolean. If not set, we attempt to fix the output of
#                       pkg-config for OpenCV.
# CFLAGS                Extra flags to give to the C/C++ MEX compiler.
# LDFLAGS               Extra flags to give to compiler when it invokes the
#                       linker.
#
# The above settings can be defined as shell environment variables and/or
# specified on the command line as arguments to make:
#
#    export VAR=value
#    make VAR=value
#
# The following targets are available:
#
# all      Builds mexopencv.
# contrib  Builds the extra modules. Needs opencv_contrib.
# clean    Deletes temporary and generated files.
# doc      Generates source documentation using Doxygen.
# test     Run MATLAB unit-tests.
#
# Note that the Makefile uses pkg-config to locate OpenCV, so you need to have
# the opencv.pc file accessible from the PKG_CONFIG_PATH environment variable.
#
# Required OpenCV version: 3.0
#
# ============================================================================

# programs
MATLABDIR  ?= /opt/MATLAB/R2014b
MEX        ?= $(MATLABDIR)/bin/mex
MATLAB     ?= $(MATLABDIR)/bin/matlab
DOXYGEN    ?= doxygen
MEXOPENCVDIR ?= /home/malte/work/progs/toolboxes/mexopencv


# mexopencv directories
TARGETDIR  = .
INCLUDEDIR = $(MEXOPENCVDIR)/include
LIBDIR     = $(MEXOPENCVDIR)/lib
SRCDIR     = src

# file extensions
OBJEXT     ?= o
LIBEXT     ?= a
MEXEXT     ?= $(shell $(MATLABDIR)/bin/mexext)
ifeq ($(MEXEXT),)
    $(error "MEX extension not set")
endif


# savevideo
SAVEVIDEOSRC1 = $(SRCDIR)/FrameRateCounter.cpp
SAVEVIDEOOBJ1 = $(TARGETDIR)/FrameRateCounter.$(OBJEXT)
SAVEVIDEOSRC = $(SRCDIR)/SaveVideoClass.cpp
SAVEVIDEOOBJ = $(TARGETDIR)/SaveVideoClass.$(OBJEXT)
FLYCAPINCLUDES = -I$(SRCDIR) -I/usr/include/glibmm-2.4 -I/usr/lib64/glibmm-2.4/include -I/usr/include/glib-2.0 -I/usr/lib64/glib-2.0/include/ -I/usr/include/sigc++-2.0 -I/usr/lib64/sigc++-2.0/include/ -I/usr/include/flycapture 
FLYCAPFLAGS =  -lsigc-2.0 -lglibmm-2.4 -lglib-2.0 -lstdc++ -lncurses -lflycapture

# struc
SAR = $(SRCDIR)/strucarr2strucmat.c
SARTARGET = $(TARGETDIR)/strucarr2strucmat.$(MEXEXT)

# pdist
PDIST = $(SRCDIR)/pdist2Euclidean.c
PDISTTARGET = $(TARGETDIR)/pdist2Euclidean.$(MEXEXT)

# pcenterlinedist
PCLDIST = $(SRCDIR)/pdist2CenterLine.c
PCLDISTTARGET = $(TARGETDIR)/pdist2CenterLine.$(MEXEXT)


# mexopencv files and targets
HEADERS    := $(wildcard $(INCLUDEDIR)/*.hpp) 
SRCS1      := $(SRCDIR)/FishVideoCapture_.cpp
TARGETS1   := $(TARGETDIR)/FishVideoCapture_.$(MEXEXT) 
SRCS2      := $(SRCDIR)/VideoHandler.cpp
OBJECTS    := $(TARGETDIR)/VideoHandler.$(OBJEXT) 
SRCS3      := $(SRCDIR)/FishVideoHandler_.cpp 
TARGETS2   := $(TARGETDIR)/FishVideoHandler_.$(MEXEXT)

# OpenCV flags
ifneq ($(shell pkg-config --exists --atleast-version=3 opencv; echo $$?), 0)
    $(error "OpenCV 3.0 package was not found in the pkg-config search path")
endif
CV_CFLAGS  := $(shell pkg-config --cflags opencv)
CV_LDFLAGS := $(shell pkg-config --libs opencv)
ifndef NO_CV_PKGCONFIG_HACK
LIB_SUFFIX := %.so %.dylib %.a %.la %.dll.a %.dll
CV_LDFLAGS := $(filter-out $(LIB_SUFFIX),$(CV_LDFLAGS)) \
              $(addprefix -L, \
                  $(sort $(dir $(filter $(LIB_SUFFIX),$(CV_LDFLAGS))))) \
              $(patsubst lib%, -l%, \
                  $(basename $(notdir $(filter $(LIB_SUFFIX),$(CV_LDFLAGS)))))
endif

# compiler/linker flags
override CFLAGS  +=  -I$(INCLUDEDIR)  $(CV_CFLAGS) $(FLYCAPINCLUDES) 
override LDFLAGS += -L$(LIBDIR) -lMxArray $(CV_LDFLAGS) $ $(FLYCAPFLAGS) 


# targets
all: $(SAVEVIDEOOBJ1) $(SAVEVIDEOOBJ) $(OBJECTS) $(TARGETS1) $(TARGETS2) $(SARTARGET) $(PDISTTARGET) $(PCLDISTTARGET)

$(SAVEVIDEOOBJ1): $(SAVEVIDEOSRC1)
	$(MEX) -c -cxx -largeArrayDims  $(CFLAGS)  $<

$(SAVEVIDEOOBJ): $(SAVEVIDEOSRC)
	$(MEX) -c -cxx -largeArrayDims  $(CFLAGS)  $(SAVEVIDEOOBJ1) $<

#  objects
$(OBJECTS): $(SRCS2)
	$(MEX) -c -cxx -largeArrayDims $(CFLAGS) $(SAVEVIDEOOBJ) $(SAVEVIDEOOBJ1) $<

 
# MEX-files
$(TARGETS1): $(SRCS1) 
	$(MEX) -cxx -largeArrayDims $(CFLAGS) -output $(TARGETS1) $< $(LDFLAGS)

$(TARGETS2): $(SRCS3)
	$(MEX) -cxx -largeArrayDims $(CFLAGS) $(SAVEVIDEOOBJ) $(SAVEVIDEOOBJ1) $(OBJECTS) -output $(TARGETS2) $< $(LDFLAGS)

$(SARTARGET): $(SAR)
	$(MEX) -largeArrayDims -output $(SARTARGET) $< 

$(PDISTTARGET): $(PDIST)
	$(MEX) -largeArrayDims -output $(PDISTTARGET) $< 

$(PCLDISTTARGET): $(PCLDIST)
	$(MEX) -largeArrayDims -output $(PCLDISTTARGET) $< 

clean:
	rm $(TARGETDIR)/$(TARGETS1) $(TARGETDIR)/$(TARGETS2) $(TARGETDIR)/$(OBJECTS) $(TARGETDIR)/$(SAVEVIDEOOBJ) $(TARGETDIR)/$(SAVEVIDEOOBJ1) $(TARGETDIR)/$(SARTARGET) $(TARGETDIR)/$(PDISTTARGET) $(TARGETDIR)/$(PCLDISTTARGET)

