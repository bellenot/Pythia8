// Header file for differential process cross sections.
// SigmaProcess: d(sigma-hat)/d(t-hat) for some 2 -> 2 processes.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_SigmaProcess_H
#define Pythia8_SigmaProcess_H

#include "Stdlib.h"
#include "Basics.h"
#include "Settings.h"
#include "ParticleData.h"
#include "StandardModel.h"
#include "Event.h"

namespace Pythia8 {
 
//**************************************************************************

// The SigmaProcess class contains simple QCD 2 -> 2 processes.

class SigmaProcess {

public:

  // Constructor.
  SigmaProcess() { }

  // Initialize static data members.
  static void initStatic();

  // Initialize alphaStrong.
  void init();

  // Choice of input parameters for calculations.
  void setupShTh(int, int, double, double);
  void setupShThAs(int, int, double, double, double);

  // Matrix elements, given inputs above.
  double gg2gg();
  double qg2qg();
  double qq2qqDiff();
  double qq2qqSame();
  double qqbar2qqbarDiff();
  double qqbar2qqbarAnti();
  double qqbar2qqbarNew();
  double qqbar2gg();
  double gg2qqbar();
  double qg2qgam();
  double qqbar2ggam();
  double qqbar2gamgam();

  // New flavour and colour flow assigned to latest process tried.
  int id3New() const {return id3;}
  int id4New() const {return id4;}
  int colFlow() const {return flow;}

  // Perform kinematics of process (in its rest frame).
  bool doKinematics(int, int, int, int, int, double, double);
  bool doKinematics();
  Particle getParton(int i) {return parton[i];}

private:

  // Static initialization data, normally only set once.
  static int alphaSorder, nQuark;
  static double alphaSvalue, alphaEM;

  // Constants: could only be changed in the code itself.
  static const double MASSMARGIN, CONVERT2MB;

  // alphaStrong calculation.
  AlphaStrong alphaScalc;

  // Data specific to current interaction.
  int id1, id2, id3, id4, idNew, flow;
  double sH, tH, uH, sH2, tH2, uH2, mNew, alphaS, alphaS2;
  Particle parton[4];

};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_SigmaProcess_H
