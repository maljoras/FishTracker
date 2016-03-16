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

#flycapture paths
FLYCAPINCLUDEDIR ?= /usr/include/flycapture 

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

#aux dirs
HELPERDIR = +fish/+helper
CAPTUREDIR = +fish/+core/@FishVideoCapture/private
HANDLERDIR = +fish/+core/@FishVideoHandlerMex/private

# savevideo
SAVEVIDEOSRC1 = $(SRCDIR)/FrameRateCounter.cpp
SAVEVIDEOOBJ1 = $(TARGETDIR)/$(HANDLERDIR)/FrameRateCounter.$(OBJEXT)
SAVEVIDEOSRC = $(SRCDIR)/SaveVideoClass.cpp
SAVEVIDEOOBJ = $(TARGETDIR)/$(HANDLERDIR)/SaveVideoClass.$(OBJEXT)
FLYCAPINCLUDES = -I$(SRCDIR) $(shell pkg-config --cflags glibmm-2.4) -I$(FLYCAPINCLUDEDIR)  
FLYCAPFLAGS =  -lsigc-2.0 -lglibmm-2.4 -lglib-2.0 -lstdc++ -lncurses -lflycapture

# struc
SAR = $(SRCDIR)/strucarr2strucmat.c
SARTARGET = $(TARGETDIR)/$(HELPERDIR)/strucarr2strucmat.$(MEXEXT)

# pdist
PDIST = $(SRCDIR)/pdist2Euclidean.c
PDISTTARGET = $(TARGETDIR)/$(HELPERDIR)/pdist2Euclidean.$(MEXEXT)

# pcenterlinedist
PCLDIST = $(SRCDIR)/pdist2CenterLine.c
PCLDISTTARGET = $(TARGETDIR)/$(HELPERDIR)/pdist2CenterLine.$(MEXEXT)

# munkres
MUNKRES = $(SRCDIR)/assignDetectionsToTracks.cpp
MUNKRESTARGET = $(TARGETDIR)/$(HELPERDIR)/assignDetectionsToTracks.$(MEXEXT)


# mexopencv files and targets
HEADERS    := $(wildcard $(INCLUDEDIR)/*.hpp) 
SRCS1      := $(SRCDIR)/FishVideoCapture_.cpp
TARGETS1   := $(TARGETDIR)/$(CAPTUREDIR)/FishVideoCapture_.$(MEXEXT) 
SRCS2      := $(SRCDIR)/VideoHandler.cpp
OBJECTS    := $(TARGETDIR)/$(HANDLERDIR)/VideoHandler.$(OBJEXT) 
SRCS3      := $(SRCDIR)/FishVideoHandler_.cpp 
TARGETS2   := $(TARGETDIR)/$(HANDLERDIR)/FishVideoHandler_.$(MEXEXT)

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
all: $(SAVEVIDEOOBJ1) $(SAVEVIDEOOBJ) $(OBJECTS) $(TARGETS1) $(TARGETS2) $(SARTARGET) $(PDISTTARGET) $(PCLDISTTARGET) $(MUNKRESTARGET)

$(SAVEVIDEOOBJ1): $(SAVEVIDEOSRC1)
	$(MEX) -c -cxx -largeArrayDims  $(CFLAGS) -outdir $(TARGETDIR)/$(HANDLERDIR)/ $<

$(SAVEVIDEOOBJ): $(SAVEVIDEOSRC)
	$(MEX) -c -cxx -largeArrayDims  $(CFLAGS)  $(SAVEVIDEOOBJ1) -outdir $(TARGETDIR)/$(HANDLERDIR)/  $<

#  objects
$(OBJECTS): $(SRCS2)
	$(MEX) -c -cxx -largeArrayDims $(CFLAGS) $(SAVEVIDEOOBJ) $(SAVEVIDEOOBJ1) -outdir $(TARGETDIR)/$(HANDLERDIR)/ $<

 
# MEX-files
$(TARGETS1): $(SRCS1) 
	$(MEX) -cxx -largeArrayDims $(CFLAGS) -output $(TARGETS1) -outdir $(TARGETDIR)/$(CAPTUREDIR) $< $(LDFLAGS)

$(TARGETS2): $(SRCS3)
	$(MEX) -cxx -largeArrayDims $(CFLAGS) $(SAVEVIDEOOBJ) $(SAVEVIDEOOBJ1) $(OBJECTS) -output $(TARGETS2) -outdir $(TARGETDIR)/$(HANDLERDIR) $< $(LDFLAGS) 

$(SARTARGET): $(SAR)
	$(MEX) -largeArrayDims -outdir $(TARGETDIR)/$(HELPERDIR) -output $(SARTARGET) $< 

$(PDISTTARGET): $(PDIST)
	$(MEX) -largeArrayDims -outdir $(TARGETDIR)/$(HELPERDIR)  -output $(PDISTTARGET) $< 

$(PCLDISTTARGET): $(PCLDIST)
	$(MEX) -largeArrayDims -outdir $(TARGETDIR)/$(HELPERDIR) -output $(PCLDISTTARGET) $< 

$(MUNKRESTARGET): $(MUNKRES)
	$(MEX) -largeArrayDims -cxx CXXFLAGS='$$CXXFLAGS -std=c++11 -lstdc++' -I$(SRCDIR) -outdir $(TARGETDIR)/$(HELPERDIR) -output $(MUNKRESTARGET) $<

clean:
	rm $(TARGETS1) $(TARGETS2) $(OBJECTS) $(SAVEVIDEOOBJ) $(SAVEVIDEOOBJ1) $(SARTARGET) $(PDISTTARGET) $(PCLDISTTARGET) $(MUNKRESTARGET)

