// ResonanceProperties.h is a part of the PYTHIA event generator.
// Copyright (C) 2007 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for resonance properties: dynamical widths etc. 
// ResonanceProperties: base class for all resonances.
// ResonanceGmZ, ...: derived classes for individual resonances.

#ifndef Pythia8_ResonanceProperties_H
#define Pythia8_ResonanceProperties_H

#include "Basics.h"
#include "Event.h"
#include "ParticleData.h"
#include "PythiaStdlib.h"
#include "Settings.h"
#include "StandardModel.h"

namespace Pythia8 {
  
//**************************************************************************

// The ResonanceProperties is the base class.

class ResonanceProperties {

public: 

  // Destructor.
  virtual ~ResonanceProperties() {}

  // Initialize static data members.
  static void initStatic();
 
  // Calculate and store partial and total widths at the nominal mass. 
  static void widthInit() {} 

  // For Higgs only: return pretabulated width for particular channel.
  // Usage: widthChan( mHat, idAbs1, idAbs2).
  static double widthChan(double, int = 0, int = 0) {return 1.;}  

  // Return fraction of width open for particle and antiparticle.
  static double openFrac( int , int ) {return 1.;}
 
  // Calculate the total/open width for given mass, charge and instate.
  virtual double width(double mHat, int idResSgn = 1, int idIn = 0, 
    bool openOnly = false, bool setBR = false) = 0; 

  // Select decay channel for given mass, charge and instate.
  virtual DecayChannel& dynamicDecay( double mHat, int idResSgn, int idIn) = 0;
  
  // Evaluate weight for angles in sequential, process-independent decays.
  // Usage: weightDecayAngles( process, iResBeg, iResEnd), where 
  // iResBeg <= i < iResEnd is range of sister partons to test decays of.
  static double weightDecayAngles( Event&, int, int);

protected:

  // Static initialization data, normally only set once.
  static int    alphaSorder, alphaEMorder;
  static double alphaSvalue;

  // Static alphaStrong and alphaElectromagnetic calculation.
  static AlphaStrong alphaS;
  static AlphaEM     alphaEM;

  // Constants: could only be changed in the code itself.
  static const int    NPOINT;
  static const double MASSMARGIN;

  // Simple routines for matrix-element integration over Breit-Wigners.
  static double numInt1BW(double mHat, double m1, double Gamma1, 
    double mMin1, double m2, int meMode = 0);
  static double numInt2BW(double mHat, double m1, double Gamma1, 
    double mMin1, double m2, double Gamma2, double mMin2, int meMode = 0);

};
  
//**************************************************************************

// The ResonanceGmZ class handles the gamma*/Z0 resonance.

class ResonanceGmZ : public ResonanceProperties {

public:

  // Constructor. 
  ResonanceGmZ() {} 
 
  // Initialize static data members.
  static void initStatic();
 
  // Calculate and store partial and total widths at the nominal mass. 
  static void widthInit(); 

  // Return fraction of width open for particle and antiparticle.
  static double openFrac( int idResSgn1, int idResSgn2 = 0);
 
  // Calculate the total/open width for given mass, charge and instate.
  virtual double width(double mHat, int = 23, int idIn = 0, 
    bool openOnly = false, bool setBR = false); 

  // Select decay channel for given mass, charge and instate.
  virtual DecayChannel& dynamicDecay( double mHat, int idResSgn, int idIn) {
    width( mHat, idResSgn, idIn, true, true);
    return particlePtr->decay.dynamicPick(idResSgn);}

private: 

  // Identify particle species.
  static int    idRes;

  // Static properties and couplings.
  static int    gmZmode;
  static double mRes, GammaRes, m2Res, GamMRat, thetaWRat, openPos, openNeg;

  // Pointer to properties of the particle species.
  static ParticleDataEntry* particlePtr;

};
  
//**************************************************************************

// The ResonanceW class handles the W+- resonance.

class ResonanceW : public ResonanceProperties {

public:

  // Constructor. 
  ResonanceW() {} 
 
  // Initialize static data members.
  static void initStatic();
 
  // Calculate and store partial and total widths at the nominal mass. 
  static void widthInit(); 

  // Return fraction of width open for particle and antiparticle.
  static double openFrac( int idResSgn1, int idResSgn2 = 0);
 
  // Calculate the total/open width for given mass, charge and instate.
  virtual double width(double mHat, int idResSgn = 24, int = 0, 
    bool openOnly = false, bool setBR = false); 

  // Select decay channel for given mass, charge and instate.
  virtual DecayChannel& dynamicDecay( double mHat, int idResSgn, int idIn) {
    width( mHat, idResSgn, idIn, true, true);
    return particlePtr->decay.dynamicPick(idResSgn);}

private: 

  // Identify particle species.
  static int    idRes;

  // Static properties and couplings.
  static double mRes, GammaRes, m2Res, GamMRat, thetaWRat, openPos, openNeg;

  // Pointer to properties of the particle species.
  static ParticleDataEntry* particlePtr;

};
  
//**************************************************************************

// The ResonanceTop class handles the top/antitop resonance.

class ResonanceTop : public ResonanceProperties {

public:

  // Constructor. 
  ResonanceTop() {} 
 
  // Initialize static data members.
  static void initStatic();
 
  // Calculate and store partial and total widths at the nominal mass. 
  static void widthInit(); 

  // Return fraction of width open for particle and antiparticle.
  static double openFrac( int idResSgn1, int idResSgn2 = 0);
 
  // Calculate the total/open width for given mass, charge and instate.
  virtual double width(double mHat, int idResSgn = 24, int = 0, 
    bool openOnly = false, bool setBR = false); 

  // Select decay channel for given mass, charge and instate.
  virtual DecayChannel& dynamicDecay( double mHat, int idResSgn, int idIn) {
    width( mHat, idResSgn, idIn, true, true);
    return particlePtr->decay.dynamicPick(idResSgn);}

  // Evaluate weight for W decay distribution in t -> W b -> f fbar b.
  static double weightDecayAngles( Event& process, int iResBeg, 
    int iResEnd);

private: 

  // Identify particle species.
  static int    idRes;

  // Static properties and couplings.
  static double mRes, GammaRes, m2Res, GamMRat, thetaWRat, m2W, 
                openPos, openNeg;

  // Pointer to properties of the particle species.
  static ParticleDataEntry* particlePtr;

};
  
//**************************************************************************

// The ResonanceSMH class handles the SM Higgs resonance.

class ResonanceSMH : public ResonanceProperties {

public:

  // Constructor. 
  ResonanceSMH() {} 
 
  // Initialize static data members.
  static void initStatic();
 
  // Calculate and store partial and total widths at the nominal mass. 
  static void widthInit(); 

  // Return pretabulated width for particular channel.
  // Usage: widthChan( mHat, idAbs1, idAbs2).
  static double widthChan(double mHat, int id1Abs = 0, int = 0);  

  // Return fraction of width open for particle and antiparticle.
  static double openFrac( int idResSgn1, int idResSgn2 = 0);
 
  // Calculate the total/open width for given mass, charge and instate.
  virtual double width(double mHat, int = 25, int = 0, 
    bool openOnly = false, bool setBR = false); 

  // Select decay channel for given mass, charge and instate.
  virtual DecayChannel& dynamicDecay( double mHat, int idResSgn, int idIn) {
    width( mHat, idResSgn, idIn, true, true);
    return particlePtr->decay.dynamicPick(idResSgn);}

  // Evaluate weight for Z0/W+- decay distributions in H -> Z0/W+ Z0/W- -> 4f.
  static double weightDecayAngles( Event& process, int iResBeg, 
    int iResEnd);

private: 

  // Identify particle species.
  static int    idRes;

  // Static properties and couplings.
  static bool   linearWidthWWZZ; 
  static int    SMHiggsParity;
  static double mRes, GammaRes, m2Res, GamMRat, sin2tW, cos2tW, mZ, mW,
                GammaZ, GammaW, GammaT, widTable[25], SMHiggsEta, 
                openPos, openNeg;

  // Pointer to properties of the particle species.
  static ParticleDataEntry* particlePtr;

  // Constants: could only be changed in the code itself.
  static const double MASSMIN, GAMMAMARGIN;

  // Sum up quark loop contributions in H0 -> g + g.
  static double eta2gg(double mHat);

  // Sum up quark loop contributions in H0 -> gamma + gamma.
  static double eta2gaga(double mHat);

  // Sum up quark loop contributions in H0 -> gamma + Z0.
  static double eta2gaZ(double mHat);

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_ResonanceProperties_H
