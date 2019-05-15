# Makefile for Pythia8 classes, and a main program.
# Still primitive; should be improved to build a proper library.
# Copyright © 2005 Torbjörn Sjöstrand

# Main program
MAIN = main01.cc

# C++ compiler
CC = g++
CFLAGS = -g -ansi -pedantic -O -W -Wall   

# f77 compiler.
FF = g77
FFLAGS = -O
FLIBS = -L/usr/lib/gcc-lib/ -lfrtbegin -lg2c

# Pythia 8 library components.
SRCS = Basics.cc Settings.cc ParticleData.cc StandardModel.cc\
	PartonDistributions.cc Event.cc Pythia.cc ProcessLevel.cc\
	PartonLevel.cc HadronLevel.cc LesHouches.cc TimeShower.cc\
	SpaceShower.cc MultipleInteractions.cc SigmaProcess.cc\
	SigmaTotal.cc Beams.cc FragmentationSystems.cc\
	FragmentationFlavZpT.cc StringFragmentation.cc\
	MiniStringFragmentation.cc ParticleDecays.cc

# Legacy Pythia 6.3 library.
FSRC = Pythia6.f

# One header file for each library file, and in addition Stdlib only header.
HDRS = $(SRCS:.cc=.h) $(SRCS:.f=.h) Stdlib.h

# One compiled file for each library file.
OBJS = $(SRCS:.cc=.o)
FOBJ = $(FSRC:.f=.o)

# Compile main program
a.out : $(MAIN) $(OBJS) $(HDRS) $(FOBJ)
	$(CC) $(CFLAGS) $(MAIN) $(OBJS) $(FOBJ) $(FLIBS)

# Compile each of the C++ classes in Pythia 8. 

Basics.o : Basics.cc Basics.h Stdlib.h
	$(CC) $(CFLAGS) -c Basics.cc

Settings.o : Settings.cc Settings.h Stdlib.h
	$(CC) $(CFLAGS) -c Settings.cc

ParticleData.o : ParticleData.cc ParticleData.h Stdlib.h Basics.h Settings.h
	$(CC) $(CFLAGS) -c ParticleData.cc

StandardModel.o : StandardModel.cc StandardModel.h Stdlib.h Basics.h\
	Settings.h ParticleData.h
	$(CC) $(CFLAGS) -c StandardModel.cc

PartonDistributions.o : PartonDistributions.cc PartonDistributions.h\
	Stdlib.h Basics.h Settings.h ParticleData.h
	$(CC) $(CFLAGS) -c PartonDistributions.cc

Beams.o : Beams.cc Beams.h Stdlib.h Basics.h Settings.h ParticleData.h\
	PartonDistributions.h Event.h FragmentationFlavZpT.h
	$(CC) $(CFLAGS) -c Beams.cc

Event.o : Event.cc Event.h Stdlib.h Basics.h Settings.h ParticleData.h
	$(CC) $(CFLAGS) -c Event.cc

Pythia.o : Pythia.cc Pythia.h Stdlib.h Basics.h Settings.h\
	ParticleData.h StandardModel.h PartonDistributions.h Beams.h Event.h\
	ProcessLevel.h PartonLevel.h HadronLevel.h LesHouches.h TimeShower.h\
	SpaceShower.h MultipleInteractions.h SigmaProcess.h SigmaTotal.h\
	FragmentationSystems.h FragmentationFlavZpT.h StringFragmentation.h\
	MiniStringFragmentation.h ParticleDecays.h
	$(CC) $(CFLAGS) -c Pythia.cc

ProcessLevel.o : ProcessLevel.cc ProcessLevel.h Stdlib.h Basics.h Settings.h\
	ParticleData.h Event.h PartonDistributions.h Beams.h LesHouches.h\
	Pythia6.h
	$(CC) $(CFLAGS) -c ProcessLevel.cc

PartonLevel.o : PartonLevel.cc PartonLevel.h Stdlib.h Basics.h Settings.h\
	ParticleData.h Event.h PartonDistributions.h Beams.h TimeShower.h\
	SpaceShower.h MultipleInteractions.h
	$(CC) $(CFLAGS) -c PartonLevel.cc

HadronLevel.o : HadronLevel.cc HadronLevel.h Stdlib.h Basics.h Settings.h\
	ParticleData.h Event.h FragmentationSystems.h FragmentationFlavZpT.h\
	StringFragmentation.h MiniStringFragmentation.h ParticleDecays.h
	$(CC) $(CFLAGS) -c HadronLevel.cc 

LesHouches.o : LesHouches.cc LesHouches.h
	$(CC) $(CFLAGS) -c LesHouches.cc

TimeShower.o : TimeShower.cc TimeShower.h Stdlib.h Basics.h\
	Settings.h ParticleData.h StandardModel.h Event.h 
	$(CC) $(CFLAGS) -c TimeShower.cc

SpaceShower.o : SpaceShower.cc SpaceShower.h Stdlib.h Basics.h\
	Settings.h ParticleData.h StandardModel.h Event.h\
	PartonDistributions.h Beams.h 
	$(CC) $(CFLAGS) -c SpaceShower.cc

MultipleInteractions.o : MultipleInteractions.cc MultipleInteractions.h\
	Stdlib.h Basics.h Settings.h ParticleData.h Event.h\
	PartonDistributions.h Beams.h SigmaTotal.h SigmaProcess.h
	$(CC) $(CFLAGS) -c MultipleInteractions.cc

SigmaProcess.o : SigmaProcess.cc SigmaProcess.h Stdlib.h Basics.h\
	Settings.h ParticleData.h StandardModel.h Event.h
	$(CC) $(CFLAGS) -c SigmaProcess.cc

SigmaTotal.o : SigmaTotal.cc SigmaTotal.h Stdlib.h Settings.h\
	ParticleData.h 
	$(CC) $(CFLAGS) -c SigmaTotal.cc

FragmentationSystems.o : FragmentationSystems.cc FragmentationSystems.h\
	Stdlib.h Basics.h Settings.h ParticleData.h Event.h
	$(CC) $(CFLAGS) -c FragmentationSystems.cc 

FragmentationFlavZpT.o : FragmentationFlavZpT.cc FragmentationFlavZpT.h\
	Stdlib.h Basics.h Settings.h ParticleData.h
	$(CC) $(CFLAGS) -c FragmentationFlavZpT.cc 

StringFragmentation.o : StringFragmentation.cc StringFragmentation.h\
	Stdlib.h Basics.h Settings.h ParticleData.h Event.h\
	 FragmentationSystems.h FragmentationFlavZpT.h
	$(CC) $(CFLAGS) -c StringFragmentation.cc 

MiniStringFragmentation.o : MiniStringFragmentation.cc\
	MiniStringFragmentation.h Stdlib.h Basics.h Settings.h\
	ParticleData.h Event.h  FragmentationSystems.h\
	FragmentationFlavZpT.h
	$(CC) $(CFLAGS) -c MiniStringFragmentation.cc 

ParticleDecays.o : ParticleDecays.cc ParticleDecays.h Stdlib.h Basics.h\
	Settings.h ParticleData.h Event.h TimeShower.h\
	FragmentationFlavZpT.h
	$(CC) $(CFLAGS) -c ParticleDecays.cc 

# Compile Pythia 6.3 legacy f77 library.
$(FOBJ) : $(FSRC)
	$(FF) $(FFLAGS) -c $(FSRC)

# Clean up current directory.
clean :
	rm -f *~
	rm -f \#*
	rm -f out*
	rm -f core*

