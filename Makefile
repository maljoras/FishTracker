#!/usr/bin/make -f
# ============================================================================
#                              FishTrackerMakefile
# ============================================================================
#
# The following configuration parameters are recognized:
#
# MATLABDIR             MATLAB root directory.
# FLYCAPINCLUDEDIR      Include dir for the FlyCaptureSDK
# MEXOPENCVDIR          Path to the MEXOPENCV directory  
#
# The above settings can be defined as shell environment variables and/or
# specified on the command line as arguments to make:
#
#    export VAR=value
#    make VAR=value
#
# Note that the Makefile uses pkg-config to locate OpenCV and glibmm-2.4, so you need to have
# the opencv.pc and glibmm-2.4.pc file accessible from the PKG_CONFIG_PATH environment variable.
#
# Required OpenCV version: 3.0
#
# If no FLyCaptureSDK is available the MEX enhanced version cannot be build. Then use 
#
#   make helper
#
# to compile only the helper files
#
#
# ============================================================================

# programs
MATLABDIR  ?= /opt/MATLAB/R2014b
MATLAB     = $(MATLABDIR)/bin/matlab
MEXOPENCVDIR ?= /home/malte/work/progs/toolboxes/mexopencv
SRCDIR = src
TARGETDIR = .

#flycapture paths
FLYCAPINCLUDEDIR ?= /usr/include/flycapture 

# search whether .h can be found
ifeq ($(wildcard $(FLYCAPINCLUDEDIR)/FlyCapture2.h),)
  FLYCAPFLAG = 
  $(warning "Cannot find FlyCapture2.h. Compiling without ptGrey grabbing feature.")

  FLYCAPINCLUDES = -I$(SRCDIR) $(shell pkg-config --cflags glibmm-2.4)
  FLYCAPLIBS = $(shell pkg-config --libs glibmm-2.4) -lstdc++  
  ALLTARGET = nograb
else
  FLYCAPFLAG = -DFLYCAPTURE
  FLYCAPINCLUDES = -I$(SRCDIR) $(shell pkg-config --cflags glibmm-2.4) -I$(FLYCAPINCLUDEDIR)
  FLYCAPLIBS =  $(shell pkg-config --libs glibmm-2.4) -lstdc++ -lncurses -lflycapture
  ALLTARGET = everything
endif

# mexopencv directories
INCLUDEDIR = $(MEXOPENCVDIR)/include
LIBDIR     = $(MEXOPENCVDIR)/lib

ifeq ($(wildcard $(INCLUDEDIR)/mexopencv.hpp),)
  ALLTARGET = helper
  $(warning "Cannot find mexopencv.hpp. OpenCV/Mex-versions disabled. Compiling only helper function.")
endif

# file extensions
OBJEXT     = o
LIBEXT     = a
ifneq ($(wildcard $(MATLABDIR)/bin/mexext.bat),)
  MEXEXT     = $(shell $(MATLABDIR)/bin/mexext.bat)
  MEX        = $(MATLABDIR)/bin/mex.bat
else
  MEX        = $(MATLABDIR)/bin/mex
  MEXEXT     = $(shell $(MATLABDIR)/bin/mexext)
endif

ifeq ($(MEXEXT),)
    $(error "MEX extension not set")
endif

#aux dirs
HELPERDIR = +fish/+helper
CAPTUREDIR = +fish/+core/@FishVideoCapture/private
HANDLERDIR = +fish/+core/@FishVideoHandlerMex/private
BTRACEDIR = +fish/+core/@FishDAGraph/private
GTTDIR = +fish/@Tracker/private

# savevideo
ifneq ($(FLYCAPFLAG),)
  SAVEVIDEOSRC = $(SRCDIR)/SaveVideoClass.cpp
  SAVEVIDEOOBJ = $(TARGETDIR)/$(HANDLERDIR)/SaveVideoClass.$(OBJEXT)
endif

# struc
SAR = $(SRCDIR)/strucarr2strucmat.c
SARTARGET = $(TARGETDIR)/$(HELPERDIR)/strucarr2strucmat.$(MEXEXT)

# pdist
PDIST = $(SRCDIR)/pdist2Euclidean.c
PDISTTARGET = $(TARGETDIR)/$(HELPERDIR)/pdist2Euclidean.$(MEXEXT)

# pdist
BTRACE = $(SRCDIR)/backtrace_.c
BTRACETARGET = $(TARGETDIR)/$(BTRACEDIR)/backtrace_.$(MEXEXT)

# current tracks
GTT = $(SRCDIR)/getCurrentTracks_.c
GTTTARGET = $(TARGETDIR)/$(GTTDIR)/getCurrentTracks_.$(MEXEXT)


# pcenterlinedist
PCLDIST = $(SRCDIR)/pdist2CenterLine.c
PCLDISTTARGET = $(TARGETDIR)/$(HELPERDIR)/pdist2CenterLine.$(MEXEXT)

# munkres
MUNKRES = $(SRCDIR)/assignDetectionsToTracks.cpp
MUNKRESTARGET = $(TARGETDIR)/$(HELPERDIR)/assignDetectionsToTracks.$(MEXEXT)



ifneq ($(ALLTARGET),helper)
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
endif

# compiler/linker flags
override CFLAGS  +=  -I$(INCLUDEDIR)  $(CV_CFLAGS) $(FLYCAPINCLUDES) $(FLYCAPFLAG)
override LDFLAGS += -L$(LIBDIR) -lMxArray $(CV_LDFLAGS) $ $(FLYCAPLIBS) 


# targets
all: $(ALLTARGET) 
helper: $(SARTARGET) $(PDISTTARGET) $(PCLDISTTARGET) $(MUNKRESTARGET) $(BTRACETARGET) $(GTTTARGET)
everything: $(SAVEVIDEOOBJ) $(OBJECTS) $(TARGETS1) $(TARGETS2) helper
nograb: $(OBJECTS) $(TARGETS1) $(TARGETS2) helper


$(SAVEVIDEOOBJ): $(SAVEVIDEOSRC)
	$(MEX) -c -cxx -largeArrayDims  $(CFLAGS)  -outdir $(TARGETDIR)/$(HANDLERDIR)/  $<

#  objects
$(OBJECTS): $(SRCS2)
	$(MEX) -c -cxx -largeArrayDims $(CFLAGS) $(SAVEVIDEOOBJ) -outdir $(TARGETDIR)/$(HANDLERDIR)/ $<

 
# MEX-files
$(TARGETS1): $(SRCS1) 
	$(MEX) -cxx -largeArrayDims $(CFLAGS) -output $(TARGETS1) -outdir $(TARGETDIR)/$(CAPTUREDIR) $< $(LDFLAGS)

$(TARGETS2): $(SRCS3)
	$(MEX) -cxx -largeArrayDims $(CFLAGS) $(SAVEVIDEOOBJ) $(OBJECTS) -output $(TARGETS2) -outdir $(TARGETDIR)/$(HANDLERDIR) $< $(LDFLAGS) 

$(SARTARGET): $(SAR)
	$(MEX) -largeArrayDims -outdir $(TARGETDIR)/$(HELPERDIR) -output $(SARTARGET) $< 

$(PDISTTARGET): $(PDIST)
	$(MEX) -largeArrayDims -outdir $(TARGETDIR)/$(HELPERDIR)  -output $(PDISTTARGET) $< 

$(PCLDISTTARGET): $(PCLDIST)
	$(MEX) -largeArrayDims -outdir $(TARGETDIR)/$(HELPERDIR) -output $(PCLDISTTARGET) $< 

$(BTRACETARGET): $(BTRACE)
	$(MEX) -largeArrayDims -outdir $(TARGETDIR)/$(BTRACEDIR) -output $(BTRACETARGET) $< 

$(GTTTARGET): $(GTT)
	$(MEX) -largeArrayDims -outdir $(TARGETDIR)/$(GTTDIR) -output $(GTTTARGET) $< 

$(MUNKRESTARGET): $(MUNKRES)
	$(MEX) -largeArrayDims -cxx CXXFLAGS='$$CXXFLAGS -std=c++11 -lstdc++' -I$(SRCDIR) -outdir $(TARGETDIR)/$(HELPERDIR) -output $(MUNKRESTARGET) $<

clean:
	rm $(TARGETS1) $(TARGETS2) $(OBJECTS) $(SAVEVIDEOOBJ) $(SARTARGET) $(PDISTTARGET) $(PCLDISTTARGET) $(MUNKRESTARGET) $(BTRACETARGET) $(GTTTARGET)

