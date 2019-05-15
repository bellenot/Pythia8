// Header file for the spacelike initial-state showers.
// SpaceDipoleEnd: radiating dipole end in ISR.
// SpaceSystem: info on one interaction (among multiple ones).
// SpaceShower: handles the showering description.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_SpaceShower_H
#define Pythia8_SpaceShower_H

#include "Basics.h"
#include "Beams.h"
#include "Event.h"
#include "Information.h"
#include "ParticleData.h"
#include "PythiaStdlib.h"
#include "Settings.h"
#include "StandardModel.h"

namespace Pythia8 {
 
//**************************************************************************

// Data on radiating dipole ends, only used inside SpaceSystem and SpaceShower.

class SpaceDipoleEnd {
  
public:

  // Constructor.
  SpaceDipoleEnd( double pTmaxIn = 0., int sideIn = 0, int colIn = 0, 
    int chgIn = 0, int MEtypeIn = 0) : pTmax(pTmaxIn), side(sideIn),
    colType(colIn), chgType(chgIn), MEtype(MEtypeIn) { }

  // Basic properties related to evolution and matrix element corrections.
  double pTmax;
  int side, colType, chgType, MEtype;
  
  // Properties specific to current trial emission.
  int idMother, idSister;  
  double pT2, z, Q2, mSister, m2Sister, pT2corr, phi;

} ;
 
//**************************************************************************

// Data on colliding system, to be used inside SpaceShower.

class SpaceSystem {

public:

  // Constructors.
  SpaceSystem() {}
  SpaceSystem( vector<int> iPartonIn) {
    for (int i = 0; i < int(iPartonIn.size()); ++i) 
    iParton.push_back( iPartonIn[i] ); }

  // List of partons belonging to dipole system.
  vector<int> iParton;
  
  // Properties of system.
  int id1, id2; 
  double m, m2, x1, x2;

  // Four dipole ends: two sides with QCD (0, 1) and QED (2, 3) separate.
  SpaceDipoleEnd dipEnd[4];

} ;

 
//**************************************************************************

// The SpaceShower class does spacelike showers.

class SpaceShower {

public:

  // Constructor.
  SpaceShower() {system.reserve(10);}

  // Initialize static data members.
  static void initStatic();

  // Initialize generation. Possibility to force re-initialization by hand.
  void init( BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn);

  // Find whether to limit maximum scale of emissions.
  bool limitPTmax( Event& event);

  // Do it in several steps, for interleaved evolution.
  // Prepare system for evolution; identify ME.
  void prepare( Event& event, bool limitPTmax, int sizeOld = 0);

  // Select next pT in downwards evolution.
  double pTnext( double pTbegAll, double pTendAll);

  // ME corrections and kinematics that may give failure,
  bool branch( Event& event); 

  // Update dipole record, if MI of FSR has occured in between.
  void update();

private: 

  // Static initialization data, normally only set once.
  static bool doQCDshower, doQEDshowerByQ, doQEDshowerByL, samePTasMI,
    doMEcorrections, doPhiPolAsym;
  static int pTmaxMatch, alphaSorder, nQuark;
  static double mc, mb, mc2, mb2,  alphaSvalue, pT0Ref, ecmRef, ecmPow, 
    pTmin, alphaEM, pTminChgQ, pTminChgL;

  // Constants: could only be changed in the code itself.
  static const double CTHRESHOLD, BTHRESHOLD, EVALPDFSTEP, TINYPDF, 
    TINYKERNELPDF, TINYPT2, HEAVYPT2EVOL, HEAVYXEVOL, EXTRASPACEQ;

  // Other non-static initialization data.
  double alphaS2pi, alphaEM2pi, Lambda3flav, Lambda4flav, Lambda5flav, 
    Lambda3flav2, Lambda4flav2, Lambda5flav2, sCM, eCM, pT0, pT20,
    pT2min, pT2minChgQ, pT2minChgL; 

  // Pointers to the two incoming beams.
  BeamParticle* beamAPtr;
  BeamParticle* beamBPtr;

  // alphaStrong calculation.
  AlphaStrong alphaS;

  // List of partons in system.
  vector<int> iParton; 

  // All (sub)systems.
  vector<SpaceSystem> system;

  // Pointers to the current and hardest (so far) system and dipole ends.
  int iSysNow;
  SpaceSystem* sysNow;
  SpaceDipoleEnd* dipEndNow; 
  int iSysSel;
  SpaceSystem* sysSel;
  SpaceDipoleEnd* dipEndSel; 
 
  // Evolve a QCD dipole end. 
  void pT2nextQCD( double pT2begDip, double pT2endDip);

  // Evolve a QCD dipole end near heavy quark threshold region. 
  void pT2nearQCDthreshold( BeamParticle& beam, double m2Massive, 
    double m2Threshold, double zMinAbs, double zMaxMassive);

  // Evolve a QED dipole end. 
  void pT2nextQED( double pT2begDip, double pT2endDip);

  // Find class of ME correction.
  void findMEtype( Event& event);

  // Provide maximum of expected ME weight; for preweighting of evolution.
  double calcMEmax( int MEtype, int idMother);

  // Provide actual ME weight for current branching.
  double calcMEcorr(int MEtype, int idMother, double M2, double z, double Q2); 

  // Find coefficient of azimuthal asymmetry from gluon polarization.
  // void findAsymPol(DipoleEnd*);

};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_SpaceShower_H
