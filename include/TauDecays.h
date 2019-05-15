// TauDecays.h is a part of the PYTHIA event generator.
// Copyright (C) 2011 Philip Ilten, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for the TauDecays class.

#ifndef Pythia8_TauDecays_H
#define Pythia8_TauDecays_H

#include "Basics.h"
#include "Event.h"
#include "HelicityBasics.h"
#include "HelicityMatrixElements.h"
#include "PythiaComplex.h"
#include "PythiaStdlib.h"
#include "Settings.h"

namespace Pythia8 {

//==========================================================================

// TauDecays class.
// This class decays tau leptons, with helicity information.

class TauDecays {

public:

  // Constructor. 
  TauDecays() {};

  // Destructor.
  ~TauDecays() {}
  
  // Initializer.
  void init(Info* INFO, Settings* SET, ParticleData* DATA, Rndm* RNDM,
	    Couplings* SM);

  // Decay a tau or correlated tau pair.
  bool decay(int iDec, Event& event);

  // Choose a decay channel for a particle.
  vector<HelicityParticle> createChildren(HelicityParticle parent);

  // Perform an N-body isotropic decay.
  void isotropicDecay(vector<HelicityParticle>& p);

  // Write the decay to event record.
  void writeDecay(Event& event, vector<HelicityParticle>& p);

private: 

  // Flag whether a correlated tau decay should be performed.
  bool correlated;

  // Helicity matrix element pointers.
  HelicityMatrixElement* hardME;
  HelicityMatrixElement* decayME;

  // Hard process helicity matrix elements.
  HMETwoFermions2W2TwoFermions      hmeTwoFermions2W2TwoFermions;
  HMETwoFermions2Z2TwoFermions      hmeTwoFermions2Z2TwoFermions;
  HMETwoFermions2Gamma2TwoFermions  hmeTwoFermions2Gamma2TwoFermions;
  HMETwoFermions2GammaZ2TwoFermions hmeTwoFermions2GammaZ2TwoFermions;
  HMEHiggsEven2TwoFermions          hmeHiggsEven2TwoFermions;
  HMEHiggsOdd2TwoFermions           hmeHiggsOdd2TwoFermions;
  HMEHiggsCharged2TwoFermions       hmeHiggsCharged2TwoFermions;
  HMEUnpolarized                    hmeUnpolarized;

  // Tau decay helicity matrix elements.
  HMETau2Meson                    hmeTau2Meson;
  HMETau2TwoLeptons               hmeTau2TwoLeptons;
  HMETau2TwoMesonsViaVector       hmeTau2TwoMesonsViaVector;
  HMETau2TwoMesonsViaVectorScalar hmeTau2TwoMesonsViaVectorScalar;
  HMETau2ThreePions               hmeTau2ThreePions;
  HMETau2FourPions                hmeTau2FourPions;
  HMETau2PhaseSpace               hmeTau2PhaseSpace;

  // Particles of the hard process.
  HelicityParticle in1;
  HelicityParticle in2;
  HelicityParticle mediator;
  HelicityParticle out1;
  HelicityParticle out2;
  vector<HelicityParticle> particles;

  // The info pointer for the Pythia class.
  Info* info;

  // Pointer to SM coupling data.
  Couplings* smData;

  // Pointer to the particle data table.
  ParticleData* particleData;

  // Pointer to the random number generator.
  Rndm* rndm;

  // Pointer to settings database.
  Settings* settings;

  // Hardcoded constants.
  static const int NTRYCHANNEL, NTRYDECAY;

};

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_TauDecays_H
