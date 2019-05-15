// Header file for the timelike final-state showers.
// TimeDipoleEnd: data on a radiating dipole end.
// TimeShower: handles the showering description.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_TimeShower_H
#define Pythia8_TimeShower_H

#include "Basics.h"
#include "Event.h"
#include "ParticleData.h"
#include "Settings.h"
#include "StandardModel.h"
#include "Stdlib.h"

namespace Pythia8 {

//**************************************************************************

// Data on radiating dipole ends; only used inside TimeShower class.

class TimeDipoleEnd {

public:

  // Constructors.
  TimeDipoleEnd() : iRadiator(-1), iRecoiler(-1), pTmax(0.), colType(0), 
    chgType(0), MEtype(0), iMEpartner(-1), MEmix(0.), MEorder(true), 
    MEsplit(true), MEgluinoDau(false) { }  
  TimeDipoleEnd(int iRadiatorIn, int iRecoilerIn, double pTmaxIn = 0., 
    int colIn = 0, int chgIn = 0, int MEtypeIn = 0, 
    int iMEpartnerIn = -1, double MEmixIn = 0., bool MEorderIn = true, 
    bool MEsplitIn = true, bool MEgluinoDauIn = false) 
    : iRadiator(iRadiatorIn), iRecoiler(iRecoilerIn), pTmax(pTmaxIn), 
    colType(colIn), chgType(chgIn), MEtype(MEtypeIn), 
    iMEpartner(iMEpartnerIn), MEmix(MEmixIn), MEorder (MEorderIn), 
    MEsplit(MEsplitIn), MEgluinoDau(MEgluinoDauIn) { }

  // Basic properties related to dipole and matrix element corrections.
  int iRadiator, iRecoiler;
  double pTmax;
  int colType, chgType, MEtype, iMEpartner;
  double MEmix;
  bool MEorder, MEsplit, MEgluinoDau;

  // Properties specific to current trial emission.
  int flavour, iAunt;
  double mRad, m2Rad, mRec, m2Rec, mDip, m2Dip, m2DipCorr, 
    pT2, m2, z, mFlavour, asymPol;   
  
} ;

//**************************************************************************

// The TimeShower class does timelike showers.

class TimeShower {

public:

  // Constructor.
  TimeShower() {dipole.reserve(20);}

  // Initialize static data members.
  static void initStatic();

  // Initialize alphaStrong and related pTmin parameters.
  void init();

  // Top-level driver routine to do a single time-like shower.
  void shower( Event& event, int iBeg, int iEnd, double pTmax);

  // Do it in several steps, for interleaved evolution.
  // Prepare system for evolution; identify ME.
  void prepare( Event& event, int iBegIn = -1, int iEndIn = -1);

  // Select next pT in downwards evolution.
  double pTnext( Event& event, double pTbegAll, double pTendAll);

  // ME corrections and kinematics that may give failure,
  bool branch( Event& event); 

  // Update dipole record if MI or ISR.
  //void update( Event& event);

  // Print dipole list; for debug mainly.
  void list(ostream& = cout);

private:

  // Static initialization data, normally only set once.
  static bool doQCDshower, doQEDshowerByQ, doQEDshowerByL, 
    doMEcorrections, doPhiPolAsym;
  static int alphaSorder, nQuark;
  static double mc, mb, mc2, mb2, alphaSvalue, alphaS2pi, pTcolCutMin,
    alphaEM,alphaEM2pi, pTchgQCut, pT2chgQCut, pTchgLCut, pT2chgLCut, 
    sin2thetaW, mZ, gammaZ;

  // Constants: could only be changed in the code itself.
  static const double SIMPLIFYROOT, XMARGIN;

  // Other non-static initialization data.
  double Lambda3flav, Lambda4flav, Lambda5flav, Lambda3flav2, Lambda4flav2, 
    Lambda5flav2, pTcolCut, pT2colCut;

  // alphaStrong calculation.
  AlphaStrong alphaS;

  // All dipole ends and a pointer to the selected hardest dipole end.
  vector<TimeDipoleEnd> dipole;
  TimeDipoleEnd* dipSel;
 
  // Evolve a QCD dipole end. 
  void pT2nextQCD( double pT2begDip, double pT2endDip, TimeDipoleEnd& dip);

  // Evolve a QED dipole end. 
  void pT2nextQED( double pT2begDip, double pT2endDip, TimeDipoleEnd& dip);

  // Find kind of QCD ME correction.
  void findMEtype( Event& event, TimeDipoleEnd& dip);

  // Find type of particle; used by findMEtype.
  int findMEparticle( int id);

  // Find mixture of V and A in gamma/Z: energy- and flavour-dependent. 
  double gammaZmix( Event& event, int iRes, int iDau1, int iDau2);

  // Set up to calculate QCD ME correction with calcMEcorr.
  double findMEcorr(TimeDipoleEnd* dip, Particle& rad, Particle& partner, 
   Particle& emt);

  // Calculate value of QCD ME correction.
  double calcMEcorr( int kind, int combiIn, double mixIn, double x1, 
    double x2, double r1, double r2);

  // Find coefficient of azimuthal asymmetry from gluon polarization.
  void findAsymPol( Event& event, TimeDipoleEnd* dip);

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_TimeShower_H
