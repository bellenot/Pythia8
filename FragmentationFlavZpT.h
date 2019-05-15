// This file contains helper classes for fragmentation.
// StringFlav is used to select quark and hadron flavours.
// StringPT is used to select transverse momenta.
// StringZ is used to sample the fragmentation function f(z).
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_FragmentationFlavZpT_H
#define Pythia8_FragmentationFlavZpT_H

#include "Stdlib.h"
#include "Basics.h"
#include "Settings.h"
#include "ParticleData.h"

namespace Pythia8 {

//**************************************************************************

// The StringFlav class is used to select quark and hadron flavours.
// Purely static, since current values are stored in the calling routines. 

class StringFlav {

public:

  // Constructor. 
  StringFlav() {}

  // Initialize static data members.
  static void initStatic();

  // Pick a new flavour (including diquarks) given an incoming one.
  static int pick(int idOld);

  // Combine two flavours (including diquarks) to produce a hadron.
  static int combine(int id1, int id2);

  // Combine two quarks to produce a diquark.
  static int makeDiquark(int id1, int id2, int idHad = 0);

private: 

  // Static initialization data, normally only set once.
  static double probQQtoQ, probStoU, probSQtoQQ, probQQ1toQQ0, 
    probQandQQ, probQandS, probQandSinQQ, probQQ1corr, probQQ1corrInv,
    probQQ1norm, mesonUspin1, mesonSspin1, mesonCspin1, mesonBspin1,
    mesonMix1[2][4], mesonMix2[2][4], suppressEta, suppressEtaPrime,
    baryonClebsch12[6], baryonClebsch32[6], baryonClebschSum[6];

};
 
//**************************************************************************

// The StringZ class is used to sample the fragmentation function f(z).
// Purely static, since current values are stored in the calling routines. 

class StringZ {

public:

  // Constructor. 
  StringZ() {}

  // Initialize static data members.
  static void initStatic();
  
  // Fragmentation function: top-level to determine parameters.
  static double zFrag( int idOld, int idNew = 0, double mT2 = 1.);

private: 

  // Static initialization data, normally only set once.
  static bool usePetersonC, usePetersonB, usePetersonH;
  static double mc2, mb2, aLund, bLund, aExtraDiquark, rFactC, rFactB, 
    rFactH, epsilonC, epsilonB, epsilonH;

  // Constants: could only be changed in the code itself.
    static const double CFROMUNITY, AFROMZERO, AFROMC, EXPMAX;

  // Fragmentation function: select z according to provided parameters.
  static double zLund( double a, double b, double c = 1.);
  static double zPeterson( double epsilon);

};
 
//**************************************************************************

// The StringPT class is used to select select transverse momenta.
// Purely static, since current values are stored in the calling routines. 

class StringPT {

public:

  // Constructor. 
  StringPT() {}

  // Initialize static data members.
  static void initStatic();

  // Return px and py separately, but really same routine.
  static double px() {return pxy();}
  static double py() {return pxy();}

private: 

  // Static initialization data, normally only set once.
  static double sigmaQ, enhancedFraction, enhancedWidth;

  // pT fragmentation spectrum.
  static double pxy();

};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_FragmentationFlavZpT_H
