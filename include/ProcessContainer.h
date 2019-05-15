// This file contains the collected machinery of a process.
// ProcessContainer: contains information on a particular process.
// SetupContainers: administrates the selection/creation of processes.
// Copyright C 2007 Torbjorn Sjostrand

#ifndef Pythia8_ProcessContainer_H
#define Pythia8_ProcessContainer_H

#include "Basics.h"
#include "BeamParticle.h"
#include "Event.h"
#include "Information.h"
#include "ParticleData.h"
#include "PartonDistributions.h"
#include "PhaseSpace.h"
#include "PythiaStdlib.h"
#include "Settings.h"
#include "SigmaProcess.h"
#include "SigmaTotal.h"

namespace Pythia8 {

//**************************************************************************

// The ProcessContainer class combines pointers to matrix element and 
// phase space generator with general generation info. 

class ProcessContainer {

public:

  // Constructor. 
  ProcessContainer(SigmaProcess* sigmaProcessPtrIn = 0) 
    : sigmaProcessPtr(sigmaProcessPtrIn) {} 

  // Destructor.
  ~ProcessContainer() {delete phaseSpacePtr; delete sigmaProcessPtr;}
  
  // Store pointer to Info.
  static void setInfoPtr(Info* infoPtrIn) {infoPtr = infoPtrIn;}

  // Initialize phase space and counters.
  bool init(); 

  // Generate a trial event; accepted or not.
  bool trialProcess(); 
  
  // Give the hard subprocess.
  bool constructProcess( Event& process); 

  // Accumulate statistics after user veto.
  void accumulate() {++nAcc;}

  // Process name and code, and the number of final-state particles.
  string name()      const {return sigmaProcessPtr->name();}
  int    code()      const {return sigmaProcessPtr->code();}
  int    nFinal()    const {return sigmaProcessPtr->nFinal();}

  // Member functions for info on generation process.
  double sigmaMax()  const {return sigmaMx;}
  long   nTried()    const {return nTry;}
  long   nSelected() const {return nSel;}
  long   nAccepted() const {return nAcc;}
  double sigmaMC()  {if (nTry > nTryStat) sigmaDelta(); return sigmaFin;}
  double deltaMC()  {if (nTry > nTryStat) sigmaDelta(); return deltaFin;} 

private:

  // Static pointer to various information on the generation.
  static Info* infoPtr;

  // Constants: could only be changed in the code itself.
  static const int NSAMPLE;

  // Pointer to the subprocess matrix element.
  SigmaProcess* sigmaProcessPtr;

  // Pointer to the phase space generator.
  PhaseSpace* phaseSpacePtr;

  // Info on process.
  bool   isMinBias, isResolved, isDiffA, isDiffB, hasOctetOnium;

  // Statistics on generation process. (Long integers just in case.)
  long   nTry, nSel, nAcc, nTryStat;  
  double sigmaMx, sigmaSum, sigma2Sum, sigmaNeg, sigmaFin, deltaFin;

  // Estimate integrated cross section and its uncertainty. 
  void sigmaDelta();

};
 
//**************************************************************************

// The SetupContainers class turns the list of user-requested processes
// into a vector of ProcessContainer objects, each with a process.

class SetupContainers {

public:

  // Constructor. 
  SetupContainers() {} 
 
  // Initialization assuming all necessary data already read.
  bool init(vector<ProcessContainer*>& containerPtrs);

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_ProcessContainer_H
