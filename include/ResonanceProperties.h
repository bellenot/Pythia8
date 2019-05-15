// Header file for resonance properties: dynamical widths etc. 
// ResonanceProperties: base class for all resonances.
// ResonanceGmZ, ...: derived classes for individual resonances.
// Copyright C 2007 Torbjorn Sjostrand

#ifndef Pythia8_ResonanceProperties_H
#define Pythia8_ResonanceProperties_H

#include "Basics.h"
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
  static void widthInit(); 

  // Calculate the in-state widths: dummy implementation.
  // Usage: widthIn( mHat, idAbsIn).
  virtual double widthIn(double, int = 0) {return 1.;}  
 
  // Calculate the total width at a given energy and/or (relative) outwidths.
  virtual double width(double mHat, int idIn = 0, bool openOnly = false, 
    bool setBR = false) = 0; 

  // Select decay channel for given mass.
  virtual DecayChannel& dynamicDecay( double mHat, int idIn = 0) = 0;

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
 
  // Calculate the total width at a given energy.
  virtual double width(double mHat, int idIn = 0, bool openOnly = false, 
    bool setBR = false); 

  // Select decay channel for given mass and incoming flavour.
  virtual DecayChannel& dynamicDecay( double mHat, int idIn = 0) {
    width( mHat, idIn, true, true);
    return particlePtr->decay.dynamicPick(idRes);}

private: 

  // Identify particle species.
  static int idRes;

  // Static properties and couplings.
  static int    gmZmode;
  static double mRes, GammaRes, m2Res, GamMRat, thetaWRat;

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

  // Calculate the in-state widths for use in cross sections.
  virtual double widthIn(double, int = 0);  
 
  // Calculate the total width at a given energy.
  virtual double width(double mHat, int , bool openOnly = false, 
    bool setBR = false); 

  // Select decay channel for given mass.
  virtual DecayChannel& dynamicDecay( double mHat, int idIn = 0) {
    width( mHat, idIn, true, true);
    return particlePtr->decay.dynamicPick(idRes);}

private: 

  // Identify particle species.
  static int idRes;

  // Static properties and couplings.
  static double mRes, GammaRes, m2Res, GamMRat, thetaWRat;

  // Pointer to properties of the particle species.
  static ParticleDataEntry* particlePtr;

};
  
//**************************************************************************

// The ResonanceH class handles the SM Higgs resonance.

class ResonanceH : public ResonanceProperties {

public:

  // Constructor. 
  ResonanceH() {} 
 
  // Initialize static data members.
  static void initStatic();
 
  // Calculate and store partial and total widths at the nominal mass. 
  static void widthInit(); 

  // Calculate the in-state widths for use in cross sections.
  virtual double widthIn(double mHat, int idAbs = 0);  
 
  // Calculate the total or open width at a given energy.
  virtual double width(double mHat, int , bool openOnly = false, 
    bool setBR = false); 

  // Select decay channel for given mass.
  virtual DecayChannel& dynamicDecay( double mHat, int idIn = 0) {
    width( mHat, idIn, true, true);
    return particlePtr->decay.dynamicPick(idRes);}

private: 

  // Identify particle species.
  static int idRes;

  // Static properties and couplings.
  static bool   linearWidthWWZZ; 
  static double mRes, GammaRes, m2Res, GamMRat, sin2tW, cos2tW, mZ, mW,
                GammaZ, GammaW, GammaT, widIn[25];

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
