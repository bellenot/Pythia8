// Header file for resonance properties: dynamical widths etc. 
// ResonanceProperties: base class for all resonances.
// ResonanceGmZ, ...: derived classes for individual resonances.
// Copyright C 2006 Torbjorn Sjostrand

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

protected:

  // Static initialization data, normally only set once.
  static int alphaSorder;
  static double alphaSvalue, alphaEMfix;

  // Static alphaStrong and alphaElectromagnetic calculation.
  static AlphaStrong alphaScalc;
  static AlphaEM alphaEMcalc;

  // Constants: could only be changed in the code itself.
  static const double MASSMARGIN;

};
  
//**************************************************************************

// The ResonanceGmZ class handles the gamma*/Z0 resonance.

class ResonanceGmZ : public ResonanceProperties {

public:

  // Constructor. 
  ResonanceGmZ() {} 
 
  // Initialize static data members.
  static void initStatic();
 
  // Calculate the total width at a given energy.
  double width(double mH); 

  // Cross section for given incoming flavour to any (allowed) final flavour. 
  // Presupposes that width(mH) has already been called.
  double sigma(int id);

  // Select decay channel for given mass and incoming flavour.
  DecayChannel& dynamicDecay( double mH, int idIn);

  // Select cosine of decay angle for given incoming and outgoing flavour.
  double cosTheta( double mH, int idIn, int idOut);

private: 

  // Identify particle species.
  static int idRes;

  // Pointer to properties of the particle species.
  static ParticleDataEntry* particlePtr;

  // Static properties and couplings.
  static int gmZmode;
  static double mRes, GammaRes, m2Res, GamMRat, thetaWRat;

  // Current values.
  double widSum[4], widNorm,sigNorm, gamNorm, intNorm, resNorm;

};
  
//**************************************************************************

// The ResonanceW class handles the W+- resonance.

class ResonanceW : public ResonanceProperties {

public:

  // Constructor. 
  ResonanceW() {} 
 
  // Initialize static data members.
  static void initStatic();
 
  // Calculate the total width at a given energy.
  double width(double mH); 

  // Cross section for an incoming flavour to any (allowed) final flavour. 
  // Presupposes that width(mH) has already been called.
  double sigma(int id);

  // Select decay channel for given mass (incoming flavour dummy).
  DecayChannel& dynamicDecay( double mH, int idIn = 0);

  // Select cosine of decay angle for given incoming and outgoing flavour.
  double cosTheta( double mH, int idIn, int idOut1, int idOut2);

private: 

  // Identify particle species.
  static int idRes;

  // Pointer to properties of the particle species.
  static ParticleDataEntry* particlePtr;

  // Static properties and couplings.
  static double mRes, GammaRes, m2Res, GamMRat, thetaWRat;

  // Current values.
  double widSum, sigOpen;

  // Dummy to avoid compiler warnings.
  int idSave;

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_ResonanceProperties_H
