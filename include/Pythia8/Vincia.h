// Vincia.h is a part of the PYTHIA event generator.
// Copyright (C) 2019 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains global header information for Vincia.

#ifndef Pythia8_Vincia_H
#define Pythia8_Vincia_H

// Maths headers.
#include <limits>
#include <cmath>

// Include Pythia 8 headers.
#include "Pythia8/Event.h"
#include "Pythia8/ParticleData.h"
#include "Pythia8/PartonSystems.h"
#include "Pythia8/StandardModel.h"

// Include Vincia headers.
#include "Pythia8/VinciaAntennaFunctions.h"
#include "Pythia8/VinciaCommon.h"
#include "Pythia8/VinciaFSR.h"
#include "Pythia8/VinciaISR.h"
#include "Pythia8/VinciaQED.h"
#include "Pythia8/VinciaMG5MEs.h"

// Define namespace inside which Vincia lives.
namespace Pythia8 {

//==========================================================================

// The Vincia class. Top-level handler class for the Vincia antenna
// shower model.

class Vincia : public PhysicsBase {

public:

  // Default constructor is pure dummy.
  Vincia() {}

  // Real constructor needs a pointer to Pythia's SLHA instance (used
  // by MG5 interface).
  Vincia(SusyLesHouches* slhaPtrIn);

  // Preliminary initialization.
  void initPrel();

  // Initialize.
  bool init();

  // Method to check for and return Vincia-specific tune file. First
  // check in present working directory, then in share/Pythia8/tunes.
  string tuneFile();

  // Methods to return shower pointers.
  TimeShowerPtr  getTimesShower()    {return fsrShowerPtr;}
  TimeShowerPtr  getTimesDecShower() {return fsrShowerPtr;}
  SpaceShowerPtr getSpaceShower()    {return isrShowerPtr;}

  // Automatically set verbose level in all members.
  void setVerbose(int verboseIn);

  // Utilities for printing info and internal histograms.
  void printInfo() {fsrShowerPtr->printInfo(true);
    isrShowerPtr->printInfo(true);}
  void printHistos() {fsrShowerPtr->printHistos();}
  void writeHistos(string fileName = "vincia", string lastName = "dat") {
    fsrShowerPtr->writeHistos(fileName, lastName);}
  const Hist& getDiagnosticHistogram(string name) {
    return fsrShowerPtr->getDiagnosticHistogram(name);}

  // Public Vincia objects.
  VinciaCommon          vinCom;
  Resolution            resolution;
  shared_ptr<VinciaFSR> fsrShowerPtr;
  shared_ptr<VinciaISR> isrShowerPtr;
  QEDShower             qedShower;
  Colour                colour;
  ResScaleHook          resScaleHook;
  VinciaWeights         vinWeights;
  MECs                  mecs;

  // Auxiliary objects.
  VinciaMG5MEs          mg5mes;
  Rambo                 rambo;

  // Vectors of antenna functions.
  DGLAP         dglap;
  AntennaSetFSR antennaSetFSR;
  AntennaSetISR antennaSetISR;

  // Pointers to Pythia classes.
  SusyLesHouches* slhaPtr;

private:

  // Verbosity level.
  int verbose;
  bool isConstructed;

};

//==========================================================================

} // end Pythia8 namespace

#endif // end Pythia8_Vincia_H
