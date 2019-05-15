# Makefile for building Pythia's ROOT examples.
# Axel Naumann, 2011-03-03
#
# See README file in this directory.

# A few variables used in this Makefile:
EX   := hist tree
RUNEX := $(addprefix run-,$(EX))
SRC  := $(addsuffix .cc,$(EX))
DICT := pythiaDict.cc

PYTHIA8 ?= $(PWD)/..
DICTCXXFLAGS := -I$(PYTHIA8)/include 
CXXFLAGS := $(DICTCXXFLAGS) $(shell root-config --cflags)
LDFLAGS := $(shell root-config --ldflags --glibs) -lEG -lEGPythia8 \
  -L$(PYTHIA8)/lib -lpythia8 -llhapdfdummy

# The default target: run all examples
all: $(RUNEX)

# before we run an example we need to build it
run-%: %
	LD_LIBRARY_PATH=$(PYTHIA8)/lib:$$LD_LIBRARY_PATH ./$<

# Generic rule to build all examples
$(EX): % : %.cc
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@.exe

# The tree example needs a dictionary.
tree: $(DICT)

# How to build the dictionary:
$(DICT): pythiaROOT.h pythiaLinkdef.h
	rootcint -f $(DICT) -c $(DICTCXXFLAGS) pythiaROOT.h pythiaLinkdef.h

clean:
	rm -f $(EX) pythiaDict* pytree.root
