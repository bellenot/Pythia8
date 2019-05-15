// This file contains the collected machinery of a process.
// ProcessContainer: contains information on a particular process.
// SetupContainers: administrates the selection/creation of processes.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_ProcessContainer_H
#define Pythia8_ProcessContainer_H

#include "Stdlib.h"
#include "Basics.h"
#include "Settings.h"
#include "ParticleData.h"
#include "Event.h"
#include "Beams.h"
#include "PartonDistributions.h"
#include "Information.h"
#include "PhaseSpace.h"
#include "SigmaHat.h"
#include "SigmaTotal.h"

namespace Pythia8 {

//**************************************************************************

// The ProcessContainer class combines pointers to matrix element and 
// phase space generator with general generation info. 

class ProcessContainer {

public:

  // Constructor. 
  ProcessContainer(SigmaHat* sigmaHatPtrIn = 0) 
    : sigmaHatPtr(sigmaHatPtrIn) {} 

  // Destructor.
  ~ProcessContainer() {delete phaseSpacePtr; delete sigmaHatPtr;}

  // Get pointer to SigmaTotal.
  void setSigmaTotalPtr(SigmaTotal* sigmaTotPtrIn) {
    sigmaTotPtr = sigmaTotPtrIn;}

  // Initialize phase space and counters.
  bool init(Info& info, PDF* pdfAPtrIn, PDF* pdfBPtrIn); 

  // Generate a trial event; accepted or not.
  bool trialProcess(); 
  
  // Give the hard subprocess.
  bool constructProcess( Info& info, Event& process); 

  // Process name and code, and the number of final-state particles.
  string name() const {return sigmaHatPtr->name();}
  int code() const {return sigmaHatPtr->code();}
  int nFinal() const {return sigmaHatPtr->nFinal();}

  // SigmaTotal is needed by elastic/diffractive processes.
  bool needsSigmaTotal() const {return !isResolved;}

  // Member functions for info on generation process.
  double sigmaMax() const {return sigmaMx;}
  int nTried() const {return nTry;}
  int nAccepted() const {return nAcc;}
  double sigmaMC() const {return (nTry == 0) ? 0. : sigmaSum / nTry;}
  double deltaMC() const {return (nTry <= 1) ? 0. :
    sqrt( max( 0., (sigma2Sum / nTry - pow2(sigmaSum / nTry)) / nTry) );} 

private:

  // Pointer to the subprocess matrix element.
  SigmaHat* sigmaHatPtr;

  // Pointer to the phase space generator.
  PhaseSpace* phaseSpacePtr;

  // Pointer to total/elastic/diffractive cross sections.
  SigmaTotal* sigmaTotPtr;

  // Info on process.
  bool isResolved, isDiffA, isDiffB;

  // Statistics on generation process.
  int nTry, nAcc;  
  double sigmaMx, sigmaSum, sigma2Sum;

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
