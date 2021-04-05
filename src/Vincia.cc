// VinciaCommon.cc is a part of the PYTHIA event generator.
// Copyright (C) 2019 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the headers) for the Vincia class.

#include "Pythia8/Vincia.h"

namespace Pythia8 {

//==========================================================================

// Vincia parton shower class.

//--------------------------------------------------------------------------

// Real constructor needs a pointer to Pythia's SLHA instance.

Vincia::Vincia(SusyLesHouches* slhaPtrIn) {

  // First set pointers to Pythia classes.
  isConstructed = false;
  slhaPtr       = slhaPtrIn;

  // Create VinciaFSR and VinciaISR instances.
  fsrShowerPtr = make_shared<VinciaFSR>() ;
  isrShowerPtr = make_shared<VinciaISR>() ;

}

//--------------------------------------------------------------------------

// Preliminary initialization.

void Vincia::initPrel() {

  // Verbosity level.
  setVerbose(settingsPtr->mode("Vincia:verbose"));

  // Init FSR shower pointers and default settings, beyond those set
  // by the non-virtual TimeShower::initPtr().
  fsrShowerPtr->initVinciaPtrs(&colour,isrShowerPtr,&qedShower,&mecs,
                               &resolution, &vinCom,&vinWeights);

  // Init ISR shower pointers and default settings, beyond those set
  // by the non-virtual SpaceShower::initPtr().
  isrShowerPtr->initVinciaPtrs(&colour,fsrShowerPtr,&qedShower,&mecs,
                               &resolution, &vinCom,&vinWeights);
  isConstructed = true;

}

//--------------------------------------------------------------------------

// Initialize.

bool Vincia::init() {

  // Verbosity level.
  setVerbose(settingsPtr->mode("Vincia:verbose"));
  if (verbose >= quiteloud) printOut(__METHOD_NAME__, "setting pointers...");

  // FSR and ISR antenna sets.
  antennaSetFSR.initPtr(infoPtr, &dglap);
  antennaSetISR.initPtr(infoPtr, &dglap);

  // QED Shower module.
  qedShower.initPtr(infoPtr, &vinCom);

  // Hand antenna set pointers to shower and matching objects.
  fsrShowerPtr->initAntPtr(&antennaSetFSR);
  isrShowerPtr->initAntPtr(&antennaSetISR);
  mecs.initAntPtr(&antennaSetFSR, &antennaSetISR);

  // Pass pointers on to objects that require them.
  resolution.initPtr(settingsPtr);
  rambo.initPtr(rndmPtr);
  vinCom.initPtr(infoPtr);
  mg5mes.initPtr(infoPtr, slhaPtr, &vinCom);
  mecs.initPtr(infoPtr, &mg5mes, &vinCom);
  colour.initPtr(infoPtr);
  vinWeights.initPtr(infoPtr, &vinCom);

  // Initialise classes not initialised by shower::init() methods.
  vinCom.init();
  resolution.init();
  colour.init();
  vinWeights.init();

  // MECs depend on Pythia Couplings so must be initialised after Pythia
  mecs.init();

  // Print VINCIA header and list of parameters
  bool vinciaOn = settingsPtr->mode("PartonShowers:model") == 2;
  if (verbose >= 1 && vinciaOn) fsrShowerPtr->header();

  // Verbose output
  if(verbose >= veryloud) printOut(__METHOD_NAME__, "end --------------");
  return true;

}

//--------------------------------------------------------------------------

// Method to check for and return Vincia-specific tune file.

string Vincia::tuneFile() {

  // Check if user asked for no tune.
  string file = settingsPtr->word("Vincia:tuneFile");
  if (toLower(file) == "none") return "none";

  // Now check if user-requested tune file lives in current dir (also
  // allows for absolute file name).
  if (fileExists(file)) return file;

  // Find path to PYTHIA data files (share/Pythia8/xmldoc directory).
  const char* PYTHIA8DATA = "PYTHIA8DATA";
  char* envPath = getenv(PYTHIA8DATA);
  string xmlPath = "";
  if (envPath != 0 && *envPath != '\0') {
    int i = 0;
    while (*(envPath+i) != '\0') xmlPath += *(envPath+(i++));
  } else return "none";

  // Check if we can find preset (VINCIA-specific) tune parameters.
  string presetFile = xmlPath+"/../tunes/" + file;

  // Check if this file exists, otherwise return "none"
  if (fileExists(presetFile)) return presetFile;
  else return "none";

}

//--------------------------------------------------------------------------

// Automatically set verbose level in all members.

void Vincia::setVerbose(int verboseIn) {

  verbose = verboseIn;
  vinCom.setVerbose(verboseIn);
  resolution.setVerbose(verboseIn);
  fsrShowerPtr->setVerbose(verboseIn);
  qedShower.setVerbose(verboseIn);
  isrShowerPtr->setVerbose(verboseIn);
  colour.setVerbose(verboseIn);
  mg5mes.setVerbose(verboseIn);
  mecs.setVerbose(verboseIn);

}

//==========================================================================

} // end namespace Pythia8
