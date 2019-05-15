// Header file for beam-remnants handling..
// BeamRemnants: matches the remnants between the two beams.
// Copyright C 2007 Torbjorn Sjostrand

#ifndef Pythia8_BeamRemnants_H
#define Pythia8_BeamRemnants_H

#include "Basics.h"
#include "BeamParticle.h"
#include "Event.h"
#include "FragmentationFlavZpT.h"
#include "Information.h"
#include "ParticleData.h"
#include "PartonDistributions.h"
#include "PythiaStdlib.h"
#include "Settings.h"

namespace Pythia8 {

//**************************************************************************

// This class matches the kinematics of the hard-scattering subsystems
// (with primordial kT added) to that of the two beam remnants.  

class BeamRemnants {

public:

  // Constructor.
  BeamRemnants() { }  

  // Initialize static data members.
  static void initStatic();

  // Initialization.
  bool init( Info* infoPtrIn, BeamParticle* beamAPtrIn, 
    BeamParticle* beamBPtrIn);

  // Select the flavours/kinematics/colours of the two beam remnants. 
  bool add( Event& event);

private: 

  // Static initialization data, normally only set once.
  static bool   primordialKT, doReconnect;
  static double primordialKTsoft, primordialKThard, primordialKTremnant,
                halfScaleForKT, halfMassForKT, reconnectRange, 
                pT0Ref, ecmRef, ecmPow;

  // Constants: could only be changed in the code itself.
  static const int NTRYCOLMATCH, NTRYKINMATCH;

  // Information set for events.
  int    nSys, oldSize;
  double eCM, sCM, pT0, pT20Rec;

  // Pointer to various information on the generation.
  Info* infoPtr;

  // Pointers to the two incoming beams.
  BeamParticle* beamAPtr;
  BeamParticle* beamBPtr;

  // Do the kinematics of the collision subsystems and two beam remnants. 
  bool setKinematics( Event& event);

  // Colour collapses (when one colour is mapped onto another).
  vector<int> colFrom, colTo;

  // Allow colour reconnections.
  bool reconnectColours( Event&  event);

  // Check that colours are consistent.
  bool checkColours( Event& event);

};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_BeamRemnants_H
