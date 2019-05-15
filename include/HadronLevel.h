// This file contains the main class for hadron-level generation.
// HadronLevel: handles administration of fragmentation and decay.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_HadronLevel_H
#define Pythia8_HadronLevel_H

#include "Stdlib.h"
#include "Basics.h"
#include "Settings.h"
#include "ParticleData.h"
#include "Event.h"
#include "Information.h"
#include "FragmentationSystems.h"
#include "StringFragmentation.h"
#include "MiniStringFragmentation.h"
#include "ParticleDecays.h"

namespace Pythia8 {
 
//**************************************************************************

// The HadronLevel class contains the top-level routines to generate
// the transition from the partonic to the hadronic stage of an event.

class HadronLevel {

public:

  // Constructor. 
  HadronLevel() {}

  // Possibility to pass in pointer for external handling of some decays.
  bool decayPtr( DecayHandler* decayHandlePtrIn, 
    vector<int> handledParticles) 
    { return decays.decayPtr( decayHandlePtrIn, handledParticles);}  

  // Initialize static data members.
  static void initStatic();

  // Initialize alphaStrong in ParticleDecays.
  void init() {decays.init();}
 
  // Generate the next event.
  bool next(Event& event); 

private: 

  // Static initialization data, normally only set once.
  static bool Hadronize, Decay;
  static double mStringMin, eNormJunction, mThad;

  // Constants: could only be changed in the code itself.
  static const int NTRYJNREST;
  static const double JJSTRINGM2MAX, JJSTRINGM2FRAC, CONVJNREST, MTHAD;

  // The main generator classes for hadronization and decay.
  StringFragmentation stringFrag;
  MiniStringFragmentation ministringFrag;
  ParticleDecays decays;

  // Configuration of colour-singlet systems.
  ColConfig colConfig;   
 
  // Trace colour flow in the event to form colour singlet subsystems.
  bool findSinglets(Event& event);
 
  // Trace a colour line, from a colour, from an anticolour, or in loop.
  bool traceFromCol(int indxCol, Event& event, int iJun = -1, int iCol = -1);
  bool traceFromAcol(int indxCol, Event& event, int iJun = -1, int iCol = -1); 
  bool traceInLoop(int indxCol, int indxAcol, Event& event);

  // Split junction-antijunction system into two, or simplify other way.
  bool splitJunctionPair(Event& event);

  // Colour information.
  vector<int> iColEnd, iAcolEnd, iColAndAcol, iParton, iPartonJun, 
    iPartonAntiJun, iJunLegA, iJunLegB, iJunLegC, iAntiLegA, iAntiLegB, 
    iAntiLegC, iGluLeg;
  vector<double> m2Pair; 
  
};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_HadronLevel_H
