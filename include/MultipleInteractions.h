// This file contains the main class for multiple interactions physics.
// MultipleInteractions: generates multiple parton-parton interactions.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_MultipleInteractions_H
#define Pythia8_MultipleInteractions_H

#include "Basics.h"
#include "Beams.h"
#include "Event.h"
#include "Information.h"
#include "PythiaStdlib.h"
#include "Settings.h"
#include "SigmaTotal.h"
#include "SigmaProcess.h"
#include "StandardModel.h"

namespace Pythia8 {
 
//**************************************************************************

// The MultipleInteractions class contains the main methods for the 
// generation of multiple parton-parton interactions in hadronic collisions.


class MultipleInteractions {

public:

  // Constructor.
  MultipleInteractions() {isInit = false; sudExpPT.resize(NBINS+1);}

  // Initialize static data members.
  static void initStatic();

  // Initialize generation. Possibility to force re-initialization by hand.
  bool init( BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn, 
    bool reInit = false);

  // Reset impact parameter choice.
  void clear() {bSetInFirst = false;}

  // Select first = hardest pT in minbias process.
  void pTfirst(); 

  // Set up kinematics for first = hardest pT in minbias process.
  void setupFirstSys( Info* infoPtr, Event& process);

  // Find whether to limit maximum scale of emissions.
  bool limitPTmax( Event& event);

  // Prepare system for evolution.
  void prepare(double pTscale = 1000.) {
    if (!bSetInFirst) overlapNext(pTscale);}

  // Select next pT in downwards evolution.
  double pTnext( double pTbegAll, double pTendAll);

  // Set up kinematics of acceptable interaction.
  void scatter( Event& event); 

  // Get some information on current interaction.
  double Q2Ren() const {return pT2Ren;}
  double alphaSH() const {return alpS;}
  double alphaEMH() const {return alpEM;}
  double x1H() const {return x1;} 
  double x2H() const {return x2;} 
  double Q2Fac() const {return pT2Fac;}
  double pdf1() const {return xPDF1now;}
  double pdf2() const {return xPDF2now;}
  double bMI() const {return bNow / bAvg;}
  double enhanceMI() const {return enhanceB / zeroIntCorr;}

  // Statistics. (Currently dummy.)
  void statistics() {}
  
private: 

  // Static initialization data, normally only set once.
  static int pTmaxMatch, alphaSorder, bProfile, nQuark, nSample;
  static double alphaSvalue, alphaEMfix, Kfactor, pT0Ref, ecmRef, ecmPow, 
    pTmin, coreRadius, coreFraction, expPow;

  // Constants: could only be changed in the code itself.
  static const bool SHIFTFACSCALE;
  static const int NBINS;
  static const double SIGMAFUDGE, RPT20, PT0STEP, SIGMASTEP, EXPPOWMIN,
    PROBATLOWB, BSTEP, BMAX, EXPMAX, KCONVERGE, CONVERT2MB;

  // Other non-static initialization data.
  double eCM, sCM, pT0, pT20, pT2min, pTmax, pT2max, pT20R, pT20minR, 
    pT20maxR, pT20min0maxR, pT2maxmin, sigmaND, pT4dSigmaMax, pT4dProbMax, 
    dSigmaApprox, sigmaInt, zeroIntCorr, normOverlap, nAvg, kNow, normPi, 
    bAvg, bDiv, probLowB, radius2B, radius2C, fracA, fracB, fracC, 
    fracAhigh, fracBhigh, fracChigh, fracABChigh, expRev, cDiv, cMax;
  vector<double> sudExpPT;
  bool isInit, lowPow;

  // Properties specific to current system.
  int id1, id2;
  double bNow, enhanceB, pT2, pT2shift, pT2Ren, pT2Fac, x1, x2, xT, xT2, 
    tau, y, sHat, tHat, uHat, alpS, alpEM, xPDF1now, xPDF2now;
  bool bSetInFirst, atLowB;

  // Pointers to the two incoming beams.
  BeamParticle* beamAPtr;
  BeamParticle* beamBPtr;

  // Total cross section parametrization.
  SigmaTotal sigmaTot;

  // Pointers to the parton-level 2 -> 2 cross sections.
  SigmaProcess *sigma2gg2ggT, *sigma2gg2ggU, *sigma2qg2qgT, *sigma2qg2qgU,
    *sigma2qq2qqSameT, *sigma2qq2qqSameU, *sigma2qqbar2qqbarSameT,
    *sigma2qqbar2qqbarSameU, *sigma2qq2qqDiffT, *sigma2qq2qqDiffU;
  SigmaProcess* dSigmaDtSel;

  // alphaStrong calculation.
  AlphaStrong alphaS;

  // Determine constant in d(Prob)/d(pT2) < const / (pT2 + r * pT20)^2.  
  void upperEnvelope();

  // Integrate the parton-parton interaction cross section.
  void jetCrossSection();

  // Evaluate "Sudakov form factor" for not having a harder interaction.
  double sudakov(double pT2sud, double enhance = 1.);

  // Do a quick evolution towards the next smaller pT.
  double fastPT2( double pT2beg);

  // Calculate the actual cross section, either for the first interaction
  // (including at initialization) or for any subsequent in the sequence. 
  double sigmaPT2(bool isFirst = false);

  // Calculate factor relating matter overlap and interaction rate.
  void overlapInit();

  // Pick impact parameter and interaction rate enhancement,
  // either before the first interaction (for minbias) or after it.
  void overlapFirst();
  void overlapNext(double pTscale);

};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_MultipleInteractions_H
