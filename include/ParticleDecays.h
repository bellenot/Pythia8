// This file contains the classes to perform a particle decay.
// DecayHandler: base class for external handling of decays.
// ParticleDecays: decay a particle.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_ParticleDecays_H
#define Pythia8_ParticleDecays_H

#include "Basics.h"
#include "Event.h"
#include "FragmentationFlavZpT.h"
#include "Information.h"
#include "ParticleData.h"
#include "PythiaStdlib.h"
#include "Settings.h"
#include "TimeShower.h"

namespace Pythia8 {
 
//**************************************************************************

// DecayHandler is base class for the external handling of decays.
// There is only one pure virtual method, that should do the decay. 

class DecayHandler {

public:

  // A pure virtual method, wherein the derived class method does a decay.
  virtual bool decay(vector<int>& idProd, vector<double>& mProd, 
    vector<Vec4>& pProd, int iDec, const Event& event) = 0;

protected:

  // Destructor.
  virtual ~DecayHandler() {}

};
 
//**************************************************************************

// The ParticleDecays class contains the routines to decay a particle.

class ParticleDecays {

public:

  // Constructor. 
  ParticleDecays() {decayHandlePtr = 0;}

  // Possibility to pass in pointer for external handling of some decays.
  bool decayPtr( DecayHandler* decayHandlePtrIn, 
    vector<int> handledParticles);  

  // Initialize static data members.
  static void initStatic();

  // Initialize alphaStrong (needed to shower decays to partons).
  void init() {times.init();}
 
  // Perform a decay of a single particle.
  bool decay(int iDec, Event& event); 

  // Did decay result in new partons to hadronize?
  bool moreToDo() const {return moreHadronization;}

private: 

  // Static initialization data, normally only set once.
  static bool limitTau0, limitTau, limitRadius, limitCylinder, limitDecay, 
    mixB, FSRinDecays;
  static double mSafety, tau0Max, tauMax, rMax, xyMax, zMax, xBdMix, 
    xBsMix, multIncrease, multRefMass, multGoffset, colRearrange, 
    probStoU, probQandS, stopMass, sRhoDal, wRhoDal;

  // Constants: could only be changed in the code itself.
  static const int NTRYDECAY;
  static const double WTCORRECTION[11];

  // Check whether a decay is allowed, given the upcoming decay vertex.
  bool checkVertex(Particle& decayer);

  // Check for oscillations B0 <-> B0bar or B_s0 <-> B_s0bar.
  bool oscillateB(Particle& decayer);

  // Do a one-body decay.
  bool oneBody(Event& event);

  // Do a two-body decay;
  bool twoBody(Event& event);

  // Do a three-body decay;
  bool threeBody(Event& event);

  // Do a multibody decay using the M-generator algorithm.
  bool mGenerator(Event& event); 

  // Translate a partonic content into a set of actual hadrons.
  bool pickHadrons(Event& event);

  // Set colour flow and scale in a decay explicitly to partons.
  bool setColours(Event& event);

  // Pointer to a handler of external decays.
  DecayHandler* decayHandlePtr;

  // Multiplicity. Decay products positions and masses.
  int meMode, mult;
  vector<int> iProd, idProd, cols, acols, idPartons, idEnds;
  vector<double> mProd, mInv, rndmOrd;
  vector<Vec4> pInv, pProd;
  bool moreHadronization;    

  // Flavour generator; needed when required to pick hadrons.
  StringFlav flavSel;

  // Timelike parton shower, needed when explicit decay to partons.
  TimeShower times;
  
};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_ParticleDecays_H
