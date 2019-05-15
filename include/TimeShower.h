// Header file for the timelike final-state showers.
// TimeDipoleEnd: data on a radiating dipole end.
// TimeShower: handles the showering description.
// Copyright C 2007 Torbjorn Sjostrand

#ifndef Pythia8_TimeShower_H
#define Pythia8_TimeShower_H

#include "Basics.h"
#include "BeamParticle.h"
#include "Event.h"
#include "ParticleData.h"
#include "PythiaStdlib.h"
#include "Settings.h"
#include "StandardModel.h"

namespace Pythia8 {

//**************************************************************************

// Data on radiating dipole ends; only used inside TimeShower class.

class TimeDipoleEnd {

public:

  // Constructors.
  TimeDipoleEnd() : iRadiator(-1), iRecoiler(-1), pTmax(0.), colType(0), 
    chgType(0), gamType(0), isrType(0), system(0), MEtype(0), 
    iMEpartner(-1), MEmix(0.), MEorder(true), MEsplit(true), 
    MEgluinoDau(false) { }  
  TimeDipoleEnd(int iRadiatorIn, int iRecoilerIn, double pTmaxIn = 0., 
    int colIn = 0, int chgIn = 0, int gamIn = 0, int isrIn = 0, 
    int systemIn = 0, int MEtypeIn = 0, int iMEpartnerIn = -1, 
    double MEmixIn = 0., bool MEorderIn = true, bool MEsplitIn = true, 
    bool MEgluinoDauIn = false) : iRadiator(iRadiatorIn), 
    iRecoiler(iRecoilerIn), pTmax(pTmaxIn), colType(colIn), chgType(chgIn), 
    gamType(gamIn), isrType(isrIn), system(systemIn), MEtype(MEtypeIn), 
    iMEpartner(iMEpartnerIn), MEmix(MEmixIn), MEorder (MEorderIn), 
    MEsplit(MEsplitIn), MEgluinoDau(MEgluinoDauIn) { }

  // Basic properties related to dipole and matrix element corrections.
  int    iRadiator, iRecoiler;
  double pTmax;
  int    colType, chgType, gamType, isrType, system, MEtype, iMEpartner;
  double MEmix;
  bool   MEorder, MEsplit, MEgluinoDau;

  // Properties specific to current trial emission.
  int    flavour, iAunt;
  double mRad, m2Rad, mRec, m2Rec, mDip, m2Dip, m2DipCorr, 
         pT2, m2, z, mFlavour, asymPol;   
  
} ;

//**************************************************************************

// The TimeShower class does timelike showers.

class TimeShower {

public:

  // Constructor.
  TimeShower() {}

  // Destructor.
  virtual ~TimeShower() {}

  // Initialize static data members.
  static void initStatic();

  // Initialize alphaStrong and related pTmin parameters.
  virtual void init( BeamParticle* beamAPtrIn = 0, 
    BeamParticle* beamBPtrIn = 0);

  // Potential enhancement factor of pTmax scale for hardest emission.
  virtual double enhancePTmax() {return pTmaxFudge;}

  // Top-level routine to do a full time-like shower in resonance decay.
  virtual int shower( int iBeg, int iEnd, Event& event, double pTmax);

  // Prepare system for evolution after each new interaction; identify ME.
  virtual void prepare( int iSys, Event& event);

  // Update dipole list after each ISR emission.  
  virtual void update( int iSys, Event& event);

  // Select next pT in downwards evolution.
  virtual double pTnext( Event& event, double pTbegAll, double pTendAll);

  // ME corrections and kinematics that may give failure,
  virtual bool branch( Event& event); 

  // Tell which system was the last processed one.
  int system() const {return iSysSel;}; 

  // Print dipole list; for debug mainly.
  virtual void list( ostream& os = cout);

protected:

  // Static initialization data, normally only set once.
  static bool   doQCDshower, doQEDshowerByQ, doQEDshowerByL, 
                doQEDshowerByGamma, doMEcorrections, doPhiPolAsym,
                allowBeamRecoil;
  static int    alphaSorder, alphaEMorder, nGluonToQuark, nGammaToQuark, 
                nGammaToLepton;
  static double pTmaxFudge, mc, mb, m2c, m2b, alphaSvalue, alphaS2pi, 
                pTcolCutMin, pTchgQCut, pT2chgQCut, pTchgLCut, pT2chgLCut, 
                mMaxGamma, m2MaxGamma, octetOniumFraction, mZ, gammaZ, 
                thetaWRat;

  // Constants: could only be changed in the code itself.
  static const double SIMPLIFYROOT, XMARGIN, TINYPDF, LARGEM2;

  // Pointers to the two incoming beams.
  BeamParticle* beamAPtr;
  BeamParticle* beamBPtr;

  // Store index of last processed system.
  int iSysSel;

private:

  // Other non-static initialization data.
  double Lambda3flav, Lambda4flav, Lambda5flav, Lambda3flav2, Lambda4flav2, 
         Lambda5flav2, pTcolCut, pT2colCut;

  // alphaStrong and alphaEM calculations.
  AlphaStrong alphaS;
  AlphaEM     alphaEM;

  // All dipole ends and a pointer to the selected hardest dipole end.
  vector<TimeDipoleEnd> dipEnd;
  TimeDipoleEnd* dipSel;

  // Setup a dipole end, either QCD or QED/photon one.
  void setupQCDdip( int iSys, int i, int colTag, int colSign, Event& event);
  void setupQEDdip( int iSys, int i, int chgType, int gamType, Event& event); 

  // Evolve a QCD dipole end. 
  void pT2nextQCD( double pT2begDip, double pT2sel, TimeDipoleEnd& dip,
    Event& event);

  // Evolve a QED dipole end (except photon). 
  void pT2nextQED( double pT2begDip, double pT2sel, TimeDipoleEnd& dip,
    Event& event);

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
