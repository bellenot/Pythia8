#
# Examples Makefile.
#
#                  M. Kirsanov 07.04.2006
#                     Modified 18.11.2006

SHELL = /bin/sh

-include ../config.mk

# Location of directories.
TMPDIR=tmp
TOPDIR=$(shell \pwd)
INCDIR=include
SRCDIR=src
LIBDIR=lib
LIBDIRARCH=lib/archive
BINDIR=bin

# Pythia6 library location: in this directory or set in config file.
ifneq (x$(PYTHIA6LOCATION),x)
 ifneq (x$(PYTHIA6LIBNAME),x)
  MYPYTHIA6LOCATION=$(PYTHIA6LOCATION)
  MYPYTHIA6LIBNAME=$(PYTHIA6LIBNAME)
 else
  echo ERROR, PYTHIA6LIBNAME should be specified
  exit 1
 endif
else
 MYPYTHIA6LOCATION=../$(LIBDIRARCH)
 MYPYTHIA6LIBNAME=-lpythia6
endif

# There is no default behaviour, so remind user.
all:
	@echo "Usage: for NN = example number: make mainNN"

# Create an executable for one of the normal test programs
main01 main02 main03 main04 main05 main06 main07 main08 \
	main11 main12 main13 main14 main15 main16 main17 main18: \
	../$(LIBDIRARCH)/libpythia8.a
	@mkdir -p ../$(BINDIR)
	$(CXX) $(CXXFLAGS) -I../$(INCDIR) $@.cc -o ../$(BINDIR)/$@.exe \
	-L../$(LIBDIRARCH) -lpythia8 -L$(MYPYTHIA6LOCATION) $(MYPYTHIA6LIBNAME) \
	$(FLIBS)
	@ln -fs ../$(BINDIR)/$@.exe $@.exe

# Create an executable for the Fortran test program, as prelude to main14.cc. 
main14for: ../$(LIBDIRARCH)/libpythia8.a
	@mkdir -p ../$(BINDIR)
	$(FC) $(FFLAGS) $@.f -o ../$(BINDIR)/$@.exe \
	-L$(MYPYTHIA6LOCATION) $(MYPYTHIA6LIBNAME) $(FLIBS)
	@ln -fs ../$(BINDIR)/$@.exe $@.exe

# Create an executable linked to HepMC and CLHEP (if all goes well).
ifneq (x$(HEPMCLOCATION),x)
 ifeq ($(findstring x2,x$(HEPMCVERSION)),x2)
   main21 main22: ../$(LIBDIRARCH)/libpythia8.a ../$(LIBDIRARCH)/libhepmcinterface.a
	 @mkdir -p ../$(BINDIR)
	 $(CXX) $(CXXFLAGS) -I../$(INCDIR) -I$(HEPMCLOCATION)/include \
	 $@.cc -o ../$(BINDIR)/$@.exe \
	 -L../$(LIBDIRARCH) -lpythia8 -L$(MYPYTHIA6LOCATION) $(MYPYTHIA6LIBNAME) \
	 -lhepmcinterface \
	 -L$(HEPMCLOCATION)/lib -lHepMC \
	 $(FLIBS)
	 @ln -fs ../$(BINDIR)/$@.exe $@.exe
 else
  ifneq (x$(CLHEPLOCATION),x)
   main21 main22: ../$(LIBDIRARCH)/libpythia8.a ../$(LIBDIRARCH)/libhepmcinterface.a
	 @mkdir -p ../$(BINDIR)
	 $(CXX) $(CXXFLAGS) -I../$(INCDIR) -I$(HEPMCLOCATION)/include \
     -I$(CLHEPLOCATION)/include $@.cc -o ../$(BINDIR)/$@.exe \
	 -L../$(LIBDIRARCH) -lpythia8 -L$(MYPYTHIA6LOCATION) $(MYPYTHIA6LIBNAME) \
     -lhepmcinterface \
	 -L$(HEPMCLOCATION)/lib -lHepMC \
     -L$(CLHEPLOCATION)/lib -lCLHEP \
	 $(FLIBS)
	 @ln -fs ../$(BINDIR)/$@.exe $@.exe
  else
   main21 main22:
	@echo ERROR, apart from HepMC this target needs CLHEP, \
	variable CLHEPLOCATION
  endif
 endif
else
 main21 main22:
	@echo ERROR, this target needs HepMC, variable HEPMCLOCATION
endif

# Clean up: remove executables and outdated files.
.PHONY: clean
clean:
	rm -rf ../$(BINDIR)
	rm -rf *.exe
	rm -f *~; rm -f \#*; rm -f core*
