// This file contains the classes for string fragmentation.
// StringEnd: keeps track of the fragmentation step.
// StringFragmentation: is the top-level class.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_StringFragmentation_H
#define Pythia8_StringFragmentation_H

#include "Basics.h"
#include "Event.h"
#include "Information.h"
#include "FragmentationFlavZpT.h"
#include "FragmentationSystems.h"
#include "ParticleData.h"
#include "Settings.h"
#include "Stdlib.h"

namespace Pythia8 {
 
//**************************************************************************

// The StringEnd class contains the information related to 
// one of the current endpoints of the string system.
// Only to be used inside StringFragmentation, so no private members.

class StringEnd {

public:

  // Constructor. 
  StringEnd() {}
   
  // Set up initial endpoint values from input.
  void setUp(bool fromPosIn, int iEndIn, int idOldIn, int iMaxIn,
    double pxIn, double pyIn, double GammaIn, double xPosIn, double xNegIn); 

  // Fragment off one hadron from the string system, in flavour and pT.
  void newHadron();

  // Fragment off one hadron from the string system, in momentum space,
  // by taking steps either from positive or from negative end.
  Vec4 kinematicsHadron(StringSystem& system);

  // Update string end information after a hadron has been removed.
  void update();

  // Constants: could only be changed in the code itself.
  static const double TINY;
 
  // Data members.
  bool fromPos;
  int iEnd, iMax, idOld, idNew, idHad, iPosOld, iNegOld, iPosNew, iNegNew;
  double pxOld, pyOld, pxNew, pyNew, pxHad, pyHad, mHad, mT2Had, zHad,
    GammaOld, GammaNew, xPosOld, xPosNew, xPosHad, xNegOld, xNegNew,
    xNegHad;
  Vec4 pHad, pSoFar;

  // Classes for flavour, pT and z generation.
  StringFlav flavSel;
  StringPT pTsel;
  StringZ zSel;

};
  
//**************************************************************************

// The StringFragmentation class contains the top-level routines 
// to fragment a colour singlet partonic system.

class StringFragmentation {

public:

  // Constructor. 
  StringFragmentation() {}

  // Initialize static data members.
  static void initStatic();

  // Do the fragmentation: driver routine.
  bool fragment( int iSub, ColConfig& colConfig, Event& event);

  // Find the boost matrix to the rest frame of a junction.
  static RotBstMatrix junctionRestFrame(Vec4& p0, Vec4& p1, Vec4& p2);

private: 

  // Static initialization data, normally only set once.
  static double stopMass, stopNewFlav, stopSmear, eNormJunction,
    eBothLeftJunction, eMaxLeftJunction, eMinLeftJunction, bLund;

  // Constants: could only be changed in the code itself.
  static const int NTRYFLAV, NTRYJOIN, NSTOPMASS, NTRYJNREST, 
    NTRYJNMATCH, NTRYJRFEQ;
  static const double FACSTOPMASS, CLOSEDM2MAX, CLOSEDM2FRAC, EXPMAX,
    MATCHPOSNEG, EJNWEIGHTMAX, CONVJNREST, M2MAXJRF, CONVJRFEQ;

  // Find region where to put first string break for closed gluon loop.
  vector<int> findFirstRegion(vector<int>& iPartonIn, Event& event);

  // Set flavours and momentum position for initial string endpoints. 
  void setStartEnds(int idPos, int idNeg, StringSystem systemNow);

  // Check remaining energy-momentum whether it is OK to continue.
  bool energyUsedUp(bool fromPos);

  // Produce the final two partons to complete the system.
  bool finalTwo(bool fromPos);

  // Construct a special joining region for the final two hadrons.
  StringRegion finalRegion();

  // Store the hadrons in the normal event record, ordered from one end.
  void store(Event& event);

  // Fragment off two of the string legs in to a junction. 
  bool fragmentToJunction(Event& event);

  // Temporary event record for the produced particles.
  Event hadrons;

  // Information on the system of string regions.
  StringSystem system, systemMin, systemMid;

  // Information on the two current endpoints of the fragmenting system.
  StringEnd posEnd, negEnd; 

  // Classes for flavour, pT and z generation.
  StringFlav flavSel;
  StringPT pTsel;
  StringZ zSel;

  // Data members.
  bool hasJunction, isClosed;
  int iPos, iNeg;
  double w2Rem, stopMassNow;
  Vec4 pSum, pRem, pJunctionHadrons;

  // List of partons in string system.
  vector<int> iParton;

};  
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_StringFragmentation_H
