#
# Examples Makefile.
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
FLIBS = -lfrtbegin -lg2c
# Linux platform with gcc3: as default.
ifeq ($(ARCH), Linux)
 FC   = g77
 CC   = gcc
 CPP  = g++
 FFLAGS= -O
 CCFLAGS= -O
 CPPFLAGS= -g -ansi -pedantic -O -W -Wall
 FLIBS = -lfrtbegin -lg2c
endif
# Linux platform with gcc4: new Fortran90 compiler.
ifeq ($(ARCH), Linux-gcc4)
 FC   = gfortran
 CC   = gcc
 CPP  = g++
 FFLAGS= -O
 CCFLAGS= -O
 CPPFLAGS= -g -ansi -pedantic -O -W -Wall
 FLIBS = -lgfortran -lgfortranbegin
endif

# Location of directories.
TMPDIR=tmp
TOPDIR=$(shell \pwd)
INCDIR=include
SRCDIR=src
LIBDIR=lib
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
 MYPYTHIA6LOCATION=../$(LIBDIR)
 MYPYTHIA6LIBNAME=-lpythia6
endif

# There is no default behaviour, so remind user.
all:
	@echo "Usage: for NN = example number: make mainNN"

# Create an executable for one of the normal test programs
main01 main02 main03 main04 main05 main06 main07 main08 main09 main10 \
	main12 main13 main14:
	@mkdir -p ../$(BINDIR)
	$(CPP) $(CPPFLAGS) -I../$(INCDIR) $@.cc -o ../$(BINDIR)/$@.exe \
	-L../$(LIBDIR) -lpythia8 -L$(MYPYTHIA6LOCATION) $(MYPYTHIA6LIBNAME) \
	$(FLIBS)
	@ln -fs ../$(BINDIR)/$@.exe $@.exe

# Create an executable for the Fortran test program, as prelude to main02.cc. 
main02for:
	@mkdir -p ../$(BINDIR)
	$(FC) $(FFLAGS) $@.f -o ../$(BINDIR)/$@.exe \
	-L$(MYPYTHIA6LOCATION) $(MYPYTHIA6LIBNAME) $(FLIBS)
	@ln -fs ../$(BINDIR)/$@.exe $@.exe

# Create an executable linked to HepMC and CLHEP (if all goes well).
ifneq (x$(HEPMCLOCATION),x)
 ifneq (x$(CLHEPLOCATION),x)
  main11:
	 @mkdir -p ../$(BINDIR)
	 $(CPP) $(CPPFLAGS) -I../$(INCDIR) -I$(HEPMCLOCATION)/include \
     -I$(CLHEPLOCATION)/include $@.cc -o ../$(BINDIR)/$@.exe \
	 -L../$(LIBDIR) -lpythia8 -L$(MYPYTHIA6LOCATION) $(MYPYTHIA6LIBNAME) \
     -lhepmcinterface \
	 -L$(HEPMCLOCATION)/lib -lHepMC \
     -L$(CLHEPLOCATION)/lib -lCLHEP \
	 $(FLIBS)
	 @ln -fs ../$(BINDIR)/$@.exe $@.exe
 else
  main11:
	@echo ERROR, apart from HepMC this target needs CLHEP, \
	variable CLHEPLOCATION
 endif
else
 main11:
	@echo ERROR, this target needs HepMC, variable HEPMCLOCATION
endif

# Clean up: remove executables and outdated files.
.PHONY: clean
clean:
	rm -rf ../$(BINDIR)
	rm -rf *.exe
	rm -f *~; rm -f \#*; rm -f core*
