// This file contains the main class for multiple interactions physics.
// MultipleInteractions: generates multiple parton-parton interactions.
// Copyright © 2005 Torbjörn Sjöstrand

#ifndef Pythia8_MultipleInteractions_H
#define Pythia8_MultipleInteractions_H

#include "Stdlib.h"
#include "Basics.h"
#include "Settings.h"
#include "StandardModel.h"
#include "Event.h"
#include "Beams.h"
#include "SigmaTotal.h"
#include "SigmaProcess.h"

namespace Pythia8 {
 
//**************************************************************************

// The MultipleInteractions class contains the main methods for the 
// generation of multiple parton-parton interactions in hadronic collisions.


class MultipleInteractions {

public:

  // Constructor.
  MultipleInteractions() {sudExpPT.resize(NBINS+1);}

  // Initialize static data members.
  static void initStatic();

  // Initialize generation. Possibility to force re-initialization by hand.
  void init( BeamParticle& beamA, BeamParticle& beamB);

  // Prepare system for evolution.
  void prepare(double pTscale = 1000.);

  // Select next pT in downwards evolution.
  double pTnext( BeamParticle& beamA, BeamParticle& beamB, double pTbegAll, 
    double pTendAll);

  // Set up kinematics of acceptable interaction.
  bool scatter( BeamParticle& beamA, BeamParticle& beamB, Event& event); 

  // Statistics. (Currently dummy.)
  void statistics() {}
  
private: 

  // Static initialization data, normally only set once.
  static int alphaSorder, bProfile, nQuark, nSample;
  static double alphaSvalue, Kfactor, pT0Ref, ecmRef, ecmPow, pTmin, 
    coreRadius, coreFraction, expPow;

  // Constants: could only be changed in the code itself.
  static const int NBINS;
  static const double SIGMAFUDGE, RPT20, PT0STEP, SIGMASTEP, EXPPOWMIN,
    BSTEP, BMAX, EXPMAX, KCONVERGE, CONVERT2MB;

  // Other non-static initialization data.
  double eCM, sCM, pT0, pT20, pT2min, pTmax, pT2max, pT20R, pT20minR, 
    pT20maxR, pT20min0maxR, pT2maxmin, sigmaND, pT4dSigmaMax, 
    pT4dProbMax, dSigmaApprox, sigmaInt, zeroIntCorr, normOverlap, 
    radius2B, radius2C, fractionA, fractionB, fractionC;
  vector<double> sudExpPT;

  // Total and differential cross section calculation.
  SigmaTotal sigmaTot;
  SigmaProcess dSigmaDt1, dSigmaDt2;
  SigmaProcess* dSigmaDtSel;

  // alphaStrong calculation.
  AlphaStrong alphaS;

  // Determine constant in d(Prob)/d(pT2) < const / (pT2 + r * pT20)^2.  
  void upperEnvelope( BeamParticle& beamA, BeamParticle& beamB);

  // Integrate the parton-parton interaction cross section.
  void jetCrossSection( BeamParticle& beamA, BeamParticle& beamB);

  // Calculate the actual cross section for initialization decision.
  double sigmaPT2( BeamParticle& beamA, BeamParticle& beamB);

  // Calculate factor relating matter overlap and interaction rate.
  void overlapFactor();

  // Pick interaction rate enhancement related to impact parameter.
  void selectEnhance(double pTscale);

  // Do a quick evolution towards the next smaller pT.
  double fastPT2( double pT2beg);

  // Calculate the actual cross section to decide whether fast choice is OK.
  bool acceptPT2( BeamParticle& beamA, BeamParticle& beamB);

  // Properties specific to current system.
  int id1, id2;
  double pT2, pT2shift, x1, x2, xT, xT2, tau, y, sHat, enhanceB;

};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_MultipleInteractions_H
