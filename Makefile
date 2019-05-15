#
# Libraries Makefile. Some ideas from Geant4 Makefiles
#
#                  M. Kirsanov 07.04.2006

SHELL = /bin/sh

# Compilers and compiler options.
# Default.
FC   = g77
CC   = gcc
CPP  = g++
FFLAGS= -O
CCFLAGS= -O
CPPFLAGS= -g -ansi -pedantic -O -W -Wall
# Linux platform with gcc3: as default.
ifeq ($(ARCH), Linux)
 FC   = g77
 CC   = gcc
 CPP  = g++
 FFLAGS= -O
 CCFLAGS= -O
 CPPFLAGS= -g -ansi -pedantic -O -W -Wall
endif
# Linux platform with gcc4: new Fortran90 compiler.
ifeq ($(ARCH), Linux-gcc4)
 FC   = gfortran
 CC   = gcc
 CPP  = g++
 FFLAGS= -O
 CCFLAGS= -O
 CPPFLAGS= -g -ansi -pedantic -O -W -Wall
endif

# Location of directories.
TMPDIR=tmp
TOPDIR=$(shell \pwd)
INCDIR=include
SRCDIR=src
LIBDIR=lib
BINDIR=bin

# Location of libraries to be built.
targets=$(LIBDIR)/libpythia8.a
ifneq (x$(HEPMCLOCATION),x)
 targets+=$(LIBDIR)/libhepmcinterface.a
endif
ifeq (x$(PYTHIA6LOCATION),x)
 targets+=$(LIBDIR)/libpythia6.a
endif

# Main part: build Pythia8 library. 

$(TMPDIR)/%.o : $(SRCDIR)/%.cc
	@mkdir -p $(TMPDIR)
	$(CPP) $(CPPFLAGS) -c -I$(INCDIR) $< -o $@

# Creating the dependency files *.d
# The compiler with option -M is used to build the dependency strings. They
# are further edited with sed (stream editor). The first sed command adds the
# dependency for the *.d files themselves, the second one is needed because
# object files are put in the directory different from src. The last line
# removes empty *.d files produced in case of error.

$(TMPDIR)/%.d : $(SRCDIR)/%.cc
	@echo Making dependency for file $<; \
	mkdir -p $(TMPDIR); \
	$(CC) -M -I$(INCDIR) $< | \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' | \
	sed 's/$*.o/$(TMPDIR)\/$*.o/' > $@; \
	[ -s $@ ] || rm -f $@

all: $(targets)

objects := $(patsubst $(SRCDIR)/%.cc,$(TMPDIR)/%.o,$(wildcard $(SRCDIR)/*.cc))

$(LIBDIR)/libpythia8.a: $(objects)
	@mkdir -p $(LIBDIR)
	ar cru $(LIBDIR)/libpythia8.a $(objects)

deps := $(patsubst $(SRCDIR)/%.cc,$(TMPDIR)/%.d,$(wildcard $(SRCDIR)/*.cc))

# The "if" below is needed in order to avoid producing the dependency files
# when you want to just clean

ifneq ($(MAKECMDGOALS),clean)
-include $(deps)
endif

# Build Pythia6 library if a location with existing Pythia6 is not set

ifeq (x$(PYTHIA6LOCATION),x)

 $(TMPDIR)/%.o : pythia6/%.f
	@mkdir -p $(TMPDIR)
	$(FC) $(FFLAGS) -c $< -o $@

 objectsP6 := $(patsubst pythia6/%.f,$(TMPDIR)/%.o,$(wildcard pythia6/*.f))

 $(LIBDIR)/libpythia6.a : $(objectsP6)
	@mkdir -p $(LIBDIR)
	ar cru $(LIBDIR)/libpythia6.a $(objectsP6)

endif

# Build HepMC interface part if HepMC and CLHEP locations are set.

ifneq (x$(HEPMCLOCATION),x)
 ifneq (x$(CLHEPLOCATION),x)

  $(TMPDIR)/%.o : hepmcinterface/%.cc
	@mkdir -p $(TMPDIR)
	$(CPP) $(CPPFLAGS) -c -I$(INCDIR) -I$(HEPMCLOCATION)/include \
	-I$(CLHEPLOCATION)/include $< -o $@

  $(TMPDIR)/%.d : hepmcinterface/%.cc
	@echo Making dependency for file $<; \
	mkdir -p $(TMPDIR); \
	$(CC) -M -I$(INCDIR) -I$(HEPMCLOCATION)/include -I$(CLHEPLOCATION)/include $< | \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' | \
	sed 's/$*.o/$(TMPDIR)\/$*.o/' > $@; \
	[ -s $@ ] || rm -f $@

  objectsI := $(patsubst hepmcinterface/%.cc,$(TMPDIR)/%.o,$(wildcard hepmcinterface/*.cc))

  $(LIBDIR)/libhepmcinterface.a : $(objectsI)
	@mkdir -p $(LIBDIR)
	ar cru $(LIBDIR)/libhepmcinterface.a $(objectsI)

  depsI := $(patsubst hepmcinterface/%.cc,$(TMPDIR)/%.d,$(wildcard hepmcinterface/*.cc))

  ifneq ($(MAKECMDGOALS),clean)
  -include $(depsI)
  endif

 else
  $(LIBDIR)/libhepmcinterface.a : hepmcinterface/I_Pythia8.cc
	@echo ERROR, CLHEPLOCATION should be defined with HEPMCLOCATION
 endif
endif

# Clean up: remove (almost?) everything that cannot be recreated.

.PHONY: clean
clean:
	rm -f *~; rm -f \#*;
	rm -rf $(TMPDIR)
	rm -rf $(LIBDIR)
	rm -rf $(BINDIR)
	cd $(SRCDIR); rm -f *~; rm -f \#*; cd -
	cd $(INCDIR); rm -f *~; rm -f \#*; cd -
	cd doc; rm -f *~; rm -f \#*; cd -
	cd pythia6; rm -f *~; rm -f \#*; cd -
	cd examples; rm -rf *.exe; rm -f *~; rm -f \#*; rm -f core*; cd -

