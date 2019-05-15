// ResonanceWidths.h is a part of the PYTHIA event generator.
// Copyright (C) 2007 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for resonance properties: dynamical widths etc. 
// ResonanceWidths: base class for all resonances.
// ResonanceGmZ, ...: derived classes for individual resonances.

#ifndef Pythia8_ResonanceWidths_H
#define Pythia8_ResonanceWidths_H

#include "Basics.h"
#include "ParticleData.h"
#include "PythiaStdlib.h"
#include "Settings.h"
#include "StandardModel.h"

namespace Pythia8 {

//**************************************************************************

// Forward references to ParticleData and StandardModel classes.
class DecayChannel;
class DecayTable;
class ParticleDataEntry;
class AlphaStrong;
class AlphaEM;
  
//**************************************************************************

// The ResonanceWidths is the base class.

class ResonanceWidths {

public: 

  // Destructor.
  virtual ~ResonanceWidths() {}

  // Initialize static data members in base class.
  static void initStatic();
 
  // Calculate and store partial and total widths at the nominal mass. 
  virtual void init() {} 
 
  // Calculate the total/open width for given mass, charge and instate.
  virtual double width(int , double mHat, int idIn = 0, 
    bool openOnly = false, bool setBR = false) = 0; 

  // Return fraction of width open for particle and antiparticle.
  double openFrac(int idSgn) {return (idSgn > 0) ? openPos : openNeg;}

  // For Higgs only: return pretabulated width for particular channel.
  // Usage: widthChan( mHat, idAbs1, idAbs2).
  virtual double widthChan(double, int = 0, int = 0) {return 1.;} 

  // Return forced rescaling factor of resonance width.
  double widthRescaleFactor() {return forceFactor;} 

protected:

  // Constructor.
  ResonanceWidths() {}

  // Static initialization data, normally only set once.
  static int    alphaSorder, alphaEMorder;
  static double alphaSvalue;

  // Static alphaStrong and alphaElectromagnetic calculation.
  static AlphaStrong alphaS;
  static AlphaEM     alphaEM;

  // Constants: could only be changed in the code itself.
  static const int    NPOINT;
  static const double MASSMARGIN, MINWIDTH, BETAMIN;

  // Pointer to properties of the particle species.
  ParticleDataEntry* particlePtr;

  // Particle properties always locally present.
  int    idRes;
  bool   doForceWidth;
  double mRes, GammaRes, m2Res, GamMRat, openPos, openNeg, forceFactor;

  // Set up standard properties.
  void initBasic(ParticleDataEntry* particlePtrIn = 0);

  // Simple routines for matrix-element integration over Breit-Wigners.
  double numInt1BW(double mHat, double m1, double Gamma1, double mMin1, 
    double m2, int meMode = 1);
  double numInt2BW(double mHat, double m1, double Gamma1, double mMin1, 
    double m2, double Gamma2, double mMin2, int meMode = 1);

};
  
//**************************************************************************

// The ResonanceNeutral class handles a generic resonance that is  
// its own antiparticle, and therefore is neutral.

class ResonanceNeutral : public ResonanceWidths {

public:

  // Constructor. 
  ResonanceNeutral(ParticleDataEntry* particlePtrIn) 
    {initBasic(particlePtrIn);} 
 
  // Calculate and store partial and total widths at the nominal mass. 
  virtual void init(); 
 
  // Calculate the total/open width for given mass, charge and instate.
  virtual double width(int , double mHat, int = 0, 
    bool openOnly = false, bool setBR = false); 

private: 

};
  
//**************************************************************************

// The ResonanceCharged class handles a generic resonance that has a
// nonidentical antiparticle, usually but not necessarily charged.

class ResonanceCharged : public ResonanceWidths {

public:

  // Constructor. 
  ResonanceCharged(ParticleDataEntry* particlePtrIn) 
    {initBasic(particlePtrIn);} 
 
  // Calculate and store partial and total widths at the nominal mass. 
  virtual void init(); 
 
  // Calculate the total/open width for given mass, charge and instate.
  virtual double width(int , double mHat, int = 0, 
    bool openOnly = false, bool setBR = false); 

private: 

};
  
//**************************************************************************

// The ResonanceGmZ class handles the gamma*/Z0 resonance.

class ResonanceGmZ : public ResonanceWidths {

public:

  // Constructor. 
  ResonanceGmZ(ParticleDataEntry* particlePtrIn) 
    {initBasic(particlePtrIn);} 
 
  // Calculate and store partial and total widths at the nominal mass. 
  virtual void init(); 
 
  // Calculate the total/open width for given mass, charge and instate.
  virtual double width(int , double mHat, int idIn = 0, 
    bool openOnly = false, bool setBR = false); 

private: 

  // Locally stored properties and couplings.
  int    gmZmode;
  double thetaWRat;

};
  
//**************************************************************************

// The ResonanceW class handles the W+- resonance.

class ResonanceW : public ResonanceWidths {

public:

  // Constructor. 
  ResonanceW(ParticleDataEntry* particlePtrIn) 
    {initBasic(particlePtrIn);} 
 
  // Calculate and store partial and total widths at the nominal mass. 
  virtual void init(); 
 
  // Calculate the total/open width for given mass, charge and instate.
  virtual double width(int idSgn, double mHat, int = 0, 
    bool openOnly = false, bool setBR = false); 

private: 

  // Locally stored properties and couplings.
  double thetaWRat;

};
  
//**************************************************************************

// The ResonanceTop class handles the top/antitop resonance.

class ResonanceTop : public ResonanceWidths {

public:

  // Constructor. 
  ResonanceTop(ParticleDataEntry* particlePtrIn) 
    {initBasic(particlePtrIn);} 
 
  // Calculate and store partial and total widths at the nominal mass. 
  virtual void init(); 
 
  // Calculate the total/open width for given mass, charge and instate.
  virtual double width(int idSgn, double mHat, int = 0, 
    bool openOnly = false, bool setBR = false); 

private: 

  // Locally stored properties and couplings.
  double thetaWRat, m2W;

};
  
//**************************************************************************

// The ResonanceFour class handles fourth-generation resonances.

class ResonanceFour : public ResonanceWidths {

public:

  // Constructor. 
  ResonanceFour(ParticleDataEntry* particlePtrIn) 
    {initBasic(particlePtrIn);} 
 
  // Calculate and store partial and total widths at the nominal mass. 
  virtual void init(); 
 
  // Calculate the total/open width for given mass, charge and instate.
  virtual double width(int idSgn, double mHat, int = 0, 
    bool openOnly = false, bool setBR = false); 

private: 

  // Locally stored properties and couplings.
  double thetaWRat, m2W;

};
  
//**************************************************************************

// The ResonanceH class handles the SM and BSM Higgs resonance.
// higgsType = 0 : SM H; = 1: h^0/H_1; = 2 : H^0/H_2; = 3 : A^0/A_3.

class ResonanceH : public ResonanceWidths {

public:

  // Constructor. 
  ResonanceH(int higgsTypeIn ,ParticleDataEntry* particlePtrIn) 
    : higgsType(higgsTypeIn) {initBasic(particlePtrIn);} 
 
  // Calculate and store partial and total widths at the nominal mass. 
  virtual void init(); 
 
  // Calculate the total/open width for given mass, charge and instate.
  virtual double width(int , double mHat, int = 0, 
    bool openOnly = false, bool setBR = false); 

  // Return pretabulated width for particular channel.
  // Usage: widthChan( mHat, idAbs1, idAbs2).
  virtual double widthChan(double mHat, int id1Abs = 0, int = 0);  

private: 

  // Constants: could only be changed in the code itself.
  static const double MASSMIN, GAMMAMARGIN;

  // Higgs type in current instance.
  int    higgsType;

  // Locally stored properties and couplings.
  bool   useCubicWidth; 
  double sin2tW, cos2tW, mZ, mW, mHchg, GammaZ, GammaW, GammaT, 
         coup2d, coup2u, coup2l, coup2Z, coup2W, coup2Hchg, coup2H1H1, 
         coup2A3A3, coup2H1Z, coup2A3Z, coup2A3H1, coup2HchgW,
         widTable[25];

  // Sum up loop contributions in Higgs -> g + g.
  double eta2gg(double mHat);

  // Sum up loop contributions in Higgs -> gamma + gamma.
  double eta2gaga(double mHat);

  // Sum up loop contributions in Higgs -> gamma + Z0.
  double eta2gaZ(double mHat);

};
  
//**************************************************************************

// The ResonanceHchg class handles the H+- resonance.

class ResonanceHchg : public ResonanceWidths {

public:

  // Constructor. 
  ResonanceHchg(ParticleDataEntry* particlePtrIn) 
    {initBasic(particlePtrIn);} 
 
  // Calculate and store partial and total widths at the nominal mass. 
  virtual void init(); 
 
  // Calculate the total/open width for given mass, charge and instate.
  virtual double width(int idSgn, double mHat, int = 0, 
    bool openOnly = false, bool setBR = false); 

private: 

  // Locally stored properties and couplings.
  bool   useCubicWidth; 
  double thetaWRat, mW, tanBeta, tan2Beta, coup2H1W;

};
  
//**************************************************************************

// The ResonanceGraviton class handles the excited Graviton resonance.

class ResonanceGraviton : public ResonanceWidths {

public:

  // Constructor. 
  ResonanceGraviton(ParticleDataEntry* particlePtrIn) 
    {initBasic(particlePtrIn);} 
 
  // Calculate and store partial and total widths at the nominal mass. 
  virtual void init(); 
 
  // Calculate the total/open width for given mass, charge and instate.
  virtual double width(int , double mHat, int = 0, 
    bool openOnly = false, bool setBR = false); 

private: 

  // Locally stored properties and couplings.
  double kappaMG;

};

//**************************************************************************

// The ResonanceLeptoquark class handles the LQ/LQbar resonance.

class ResonanceLeptoquark : public ResonanceWidths {

public:

  // Constructor. 
  ResonanceLeptoquark(ParticleDataEntry* particlePtrIn) 
    {initBasic(particlePtrIn);} 
 
  // Calculate and store partial and total widths at the nominal mass. 
  virtual void init(); 
 
  // Calculate the total/open width for given mass, charge and instate.
  virtual double width(int idSgn, double mHat, int = 0, 
    bool openOnly = false, bool setBR = false); 

private: 

  // Locally stored properties and couplings.
  double kCoup;

};
  
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_ResonanceWidths_H
