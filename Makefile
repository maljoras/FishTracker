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
SRCDIR     = .

# file extensions
OBJEXT     ?= o
LIBEXT     ?= a
MEXEXT     ?= $(shell $(MATLABDIR)/bin/mexext)
ifeq ($(MEXEXT),)
    $(error "MEX extension not set")
endif

# mexopencv files and targets
HEADERS    := $(wildcard $(INCLUDEDIR)/*.hpp)
SRCS1      := FishVideoCapture_.cpp
TARGETS1   := $(SRCS1:.cpp=.$(MEXEXT))


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
override CFLAGS  += -cxx -largeArrayDims -I$(INCLUDEDIR) $(CV_CFLAGS)
override LDFLAGS += -L$(LIBDIR) -lMxArray $(CV_LDFLAGS)

# search path for prerequisites of pattern rules
# Note that VPATH/vpath are designed to find sources, not targets!
# (http://make.mad-scientist.net/papers/how-not-to-use-vpath/)
vpath %.cpp $(SRCDIR)/$(TARGETDIR)


# targets
all: $(TARGETS1)

# MEX-files
$(TARGETDIR)/%.$(MEXEXT) \
: %.cpp $(TARGETS1)
	$(MEX) $(CFLAGS) -output ${@:.$(MEXEXT)=} $< $(LDFLAGS)

clean:
	$(TARGETDIR)/$(TARGETS1) 

doc:
	$(DOXYGEN) Doxyfile

test:
	$(MATLAB) -nodisplay -r "addpath(pwd);cd test;try,UnitTest;catch e,disp(e.getReport);end;exit;"
