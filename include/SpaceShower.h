// Header file for the spacelike initial-state showers.
// SpaceDipoleEnd: radiating dipole end in ISR.
// SpaceSystem: info on one interaction (among multiple ones).
// SpaceShower: handles the showering description.
// Copyright C 2007 Torbjorn Sjostrand

#ifndef Pythia8_SpaceShower_H
#define Pythia8_SpaceShower_H

#include "Basics.h"
#include "BeamParticle.h"
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
  SpaceDipoleEnd( int systemIn = 0, int sideIn = 0, int iRadiatorIn = 0, 
    int iRecoilerIn = 0, double pTmaxIn = 0., int colTypeIn = 0, 
    int chgTypeIn = 0,  int MEtypeIn = 0) : system(systemIn), side(sideIn), 
    iRadiator(iRadiatorIn), iRecoiler(iRecoilerIn), pTmax(pTmaxIn), 
    colType(colTypeIn), chgType(chgTypeIn), MEtype(MEtypeIn) { }
 
  // Store values for trial emission.
  void store( int idDaughterIn, int idMotherIn, int idSisterIn,   
    double x1In, double x2In, double m2DipIn, double pT2In, double zIn, 
    double Q2In, double mSisterIn, double m2SisterIn, double pT2corrIn, 
    double phiIn) {idDaughter = idDaughterIn; idMother = idMotherIn;
    idSister = idSisterIn; x1 = x1In; x2 = x2In; m2Dip = m2DipIn;
    pT2 = pT2In; z = zIn; Q2 = Q2In; mSister = mSisterIn; 
    m2Sister = m2SisterIn; pT2corr = pT2corrIn; phi = phiIn;}
 
  // Basic properties related to evolution and matrix element corrections.
  int    system, side, iRadiator, iRecoiler;
  double pTmax;
  int    colType, chgType, MEtype;
  
  // Properties specific to current trial emission.
  int    idDaughter, idMother, idSister;  
  double x1, x2, m2Dip, pT2, z, Q2, mSister, m2Sister, pT2corr, phi;

} ;
 
//**************************************************************************

// The SpaceShower class does spacelike showers.

class SpaceShower {

public:

  // Constructor.
  SpaceShower() {}

  // Destructor.
  virtual ~SpaceShower() {}

  // Initialize static data members.
  static void initStatic();

  // Initialize generation. Possibility to force re-initialization by hand.
  virtual void init( BeamParticle* beamAPtrIn, BeamParticle* beamBPtrIn);

  // Find whether to limit maximum scale of emissions.
  virtual bool limitPTmax( Event& event);

  // Potential enhancement factor of pTmax scale for hardest emission.
  virtual double enhancePTmax() {return pTmaxFudge;}

  // Prepare system for evolution; identify ME.
  virtual void prepare( int iSys, Event& event, bool limitPTmax = true);

  // Update dipole list after each FSR emission. Currently superfluous.
  // Usage: update( iSys, event).  
  virtual void update( int , Event& ) {}

  // Select next pT in downwards evolution.
  virtual double pTnext( Event& event, double pTbegAll, double pTendAll);

  // ME corrections and kinematics that may give failure,
  virtual bool branch( Event& event); 

  // Tell which system was the last processed one.
  int system() const {return iSysSel;} 

  // Print dipole list; for debug mainly.
  virtual void list(ostream& os = cout);

protected:

  // Static initialization data, normally only set once.
  static bool   doQCDshower, doQEDshowerByQ, doQEDshowerByL, samePTasMI,
                doMEcorrections, doPhiPolAsym;
  static int    pTmaxMatch, alphaSorder, alphaEMorder, nQuark;
  static double pTmaxFudge, mc, mb, m2c, m2b, alphaSvalue, alphaS2pi, 
                pT0Ref, ecmRef, ecmPow, pTmin, pTminChgQ, pTminChgL;

  // Constants: could only be changed in the code itself.
  static const double CTHRESHOLD, BTHRESHOLD, EVALPDFSTEP, TINYPDF, 
                      TINYKERNELPDF, TINYPT2, HEAVYPT2EVOL, HEAVYXEVOL, 
                      EXTRASPACEQ, LEPTONXMIN, LEPTONXMAX, LEPTONPT2MIN, 
                      LEPTONFUDGE;

  // Pointers to the two incoming beams.
  BeamParticle* beamAPtr;
  BeamParticle* beamBPtr;

  // Store index of last processed system.
  int iSysSel;

private: 

  // Other non-static initialization data.
  double Lambda3flav, Lambda4flav, Lambda5flav, Lambda3flav2, 
         Lambda4flav2, Lambda5flav2, sCM, eCM, pT0, pT20, pT2min, 
         pT2minChgQ, pT2minChgL; 

  // Some current values.
  int    iNow, iRec, idDaughter;
  double xDaughter, x1Now, x2Now, m2Dip;

  // alphaStrong and alphaEM calculations.
  AlphaStrong alphaS;
  AlphaEM alphaEM;

  // All dipole ends
  vector<SpaceDipoleEnd> dipEnd;

  // Pointers to the current and hardest (so far) dipole ends.
  int iDipNow, iSysNow;
  SpaceDipoleEnd* dipEndNow; 
  int iDipSel;
  SpaceDipoleEnd* dipEndSel; 
 
  // Evolve a QCD dipole end. 
  void pT2nextQCD( double pT2begDip, double pT2endDip);

  // Evolve a QCD dipole end near heavy quark threshold region. 
  void pT2nearQCDthreshold( BeamParticle& beam, double m2Massive, 
    double m2Threshold, double zMinAbs, double zMaxMassive);

  // Evolve a QED dipole end. 
  void pT2nextQED( double pT2begDip, double pT2endDip);

  // Find class of ME correction.
  int findMEtype( int iSys, Event& event);

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
