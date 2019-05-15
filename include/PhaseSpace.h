// PhaseSpace.h is a part of the PYTHIA event generator.
// Copyright (C) 2007 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for phase space generators in kinematics selection.
// PhaseSpace: base class for phase space generators.
// Base class for derived classes> 2 ->1 , 2 -> 2 (+el/diff/minbias), 2 -> 3. 

#ifndef Pythia8_PhaseSpace_H
#define Pythia8_PhaseSpace_H

#include "Basics.h"
#include "BeamParticle.h"
#include "Information.h"
#include "MultipleInteractions.h"
#include "ParticleData.h"
#include "PartonDistributions.h"
#include "PythiaStdlib.h"
#include "SigmaProcess.h"
#include "SigmaTotal.h"
#include "Settings.h"
#include "UserHooks.h"

namespace Pythia8 {

//**************************************************************************

// Forward reference to the UserHooks class.
class UserHooks;
 
//**************************************************************************

// PhaseSpace is a base class for  phase space generators 
// used in the selection of hard-process kinematics.

class PhaseSpace {

public:

  // Destructor.
  virtual ~PhaseSpace() {}

  // Initialize static data members.
  static void initStatic();

  // Store static pointers to beams, SigmaTotal and user hooks.
  static void setStaticPtrs( BeamParticle* beamAPtrIn, 
    BeamParticle* beamBPtrIn, SigmaTotal* sigmaTotPtrIn,
    UserHooks* userHooksPtrIn = 0);

  // Give in pointer to cross section and cm energy.
  void initInfo(SigmaProcess* sigmaProcessPtrIn, double eCMIn);

  // A pure virtual method, wherein an optimization procedure
  // is used to determine how phase space should be sampled.
  virtual bool setupSampling() = 0; 

  // A pure virtual method, wherein a trial event kinematics 
  // is to be selected in the derived class
  virtual bool trialKin(bool inEvent = true) = 0; 

  // A pure virtual method, wherein the accepted event kinematics 
  // is to be constructed in the derived class
  virtual bool finalKin() = 0; 

  // Allow for nonisotropic decays when ME's available.
  void   decayKinematics( Event& process);

  // Give back current or maximum cross section, or set latter.
  double sigmaNow() const {return sigmaNw;}
  double sigmaMax() const {return sigmaMx;}
  bool   newSigmaMax() const {return newSigmaMx;}
  void   setSigmaMax(double sigmaMaxIn) {sigmaMx = sigmaMaxIn;}

  // Give back constructed four-vectors and known masses.
  Vec4   p(int i)   const {return pH[i];} 
  double m(int i)   const {return mH[i];}

  // Give back other event properties.
  double ecm()      const {return eCM;}
  double x1()       const {return x1H;}
  double x2()       const {return x2H;}
  double sHat()     const {return sH;}
  double tHat()     const {return tH;}
  double uHat()     const {return uH;}
  double pTHat()    const {return pTH;}
  double thetaHat() const {return theta;}
  double phiHat()   const {return phi;}
  double runBW3()   const {return runBW3H;}
  double runBW4()   const {return runBW4H;}
  double runBW5()   const {return runBW5H;}

  // Inform whether beam particles are resolved in partons or scatter directly.
  virtual bool isResolved() const {return true;}

protected:

  // Constructor.
  PhaseSpace() {}

  // Static initialization data, normally only set once.
  static bool   useBreitWigners, showSearch, showViolation;
  static int    gmZmodeGlobal;
  static double mHatGlobalMin, mHatGlobalMax, pTHatGlobalMin, pTHatGlobalMax, 
                pTHatMinDiverge, minWidthBreitWigners;

  // Constants: could only be changed in the code itself.
  static const int    NMAXTRY, NTRY3BODY;
  static const double SAFETYMARGIN, TINY, EVENFRAC, SAMESIGMA, WIDTHMARGIN, 
                      SAMEMASS, MASSMARGIN, EXTRABWWTMAX, THRESHOLDSIZE, 
                      THRESHOLDSTEP, YRANGEMARGIN, LEPTONXMIN, LEPTONXMAX, 
                      LEPTONXLOGMIN, LEPTONXLOGMAX, LEPTONTAUMIN;

  // Static information on incoming beams.
  static BeamParticle* beamAPtr;
  static BeamParticle* beamBPtr;
  static int    idA, idB;
  static double mA, mB; 
  static bool   hasLeptonBeams, hasPointLeptons;
  
  // Static pointer to the total/elastic/diffractive cross section object.
  static SigmaTotal* sigmaTotPtr;

  // Static pointer to userHooks object for user interaction with program.
  static UserHooks* userHooksPtr;
  static bool canModifySigma;

  // Center-of-mass energy.
  double eCM, s; 

  // Cross section information.
  bool   newSigmaMx;
  int    gmZmode;
  double wtBW, sigmaNw, sigmaMx, sigmaNeg;

  // Pointer to cross section. 
  SigmaProcess* sigmaProcessPtr; 

  // Process-specific kinematics properties, almost always available.
  double mHatMin, mHatMax, sHatMin, sHatMax, pTHatMin, pTHatMax, 
         pT2HatMin, pT2HatMax; 

  // Event-specific kinematics properties, almost always available.
  double x1H, x2H, m3, m4, m5, s3, s4, s5, mHat, sH, tH, uH, pAbs, p2Abs, 
         pTH, theta, phi, betaZ;
  Vec4   pH[6];
  double mH[6];

  // Much common code for normal 2 -> 1, 2 -> 2 and 2 -> 3 cases:

  // Determine how phase space should be sampled.
  void setup3Body();
  bool setupSampling123(bool is2, bool is3, ostream& os = cout); 

  // Select a trial kinematics phase space point.
  bool trialKin123(bool is2, bool is3, bool inEvent = true, ostream& os = cout); 

  // Presence and properties of any s-channel resonances.
  int    idResA, idResB;
  double mResA, mResB, GammaResA, GammaResB, tauResA, tauResB, widResA, 
         widResB;
  bool   sameResMass;

  // Kinematics properties specific to 2 -> 1/2/3.
  bool   useMirrorWeight; 
  double tau, y, z, tauMin, tauMax, yMax, zMin, zMax, ratio34, unity34, 
         zNeg, zPos, wtTau, wtY, wtZ, wt3Body, runBW3H, runBW4H, runBW5H, 
         intTau0, intTau1, intTau2, intTau3, intTau4, intTau5, intTau6, 
         intY01, intY2, intY34, mTchan1, sTchan1, mTchan2, sTchan2, 
         frac3Flat, frac3Pow1, frac3Pow2; 
  Vec4   p3cm, p4cm, p5cm;

  // Coefficients for optimized selection in 2 -> 1/2/3.
  int    nTau, nY, nZ;
  double tauCoef[8], yCoef[8], zCoef[8], tauCoefSum[8], yCoefSum[8], 
         zCoefSum[8];

  // Calculate kinematical limits for 2 -> 1/2/3.
  bool limitTau(bool is2, bool is3);
  bool limitY();
  bool limitZ();

  // Select kinematical variable between defined limits for 2 -> 1/2/3.
  void selectTau(int iTau, double tauVal, bool is2);
  void selectY(int iY, double yVal);
  void selectZ(int iZ, double zVal);
  bool select3Body();

  // Solve equation system for better phase space coefficients in 2 -> 1/2/3.
  void solveSys( int n, int bin[8], double vec[8], double mat[8][8],
    double coef[8], ostream& os = cout); 

  // Properties specific to resonance mass selection in 2 -> 2 and 2 -> 3.
  bool   useBW[6]; 
  int    idMass[6];
  double mPeak[6], sPeak[6], mWidth[6], mMin[6], mMax[6], mw[6], wmRat[6], 
         mLower[6], mUpper[6], sLower[6], sUpper[6], fracFlat[6], fracInv[6], 
         fracInv2[6], atanLower[6], atanUpper[6], intBW[6], intFlat[6], 
         intInv[6], intInv2[6]; 

  // Setup mass selection for one resonance at a time. Split in two parts.
  void   setupMass1(int iM);
  void   setupMass2(int iM, double distToThresh);

  // Do mass selection and find the associated weight.
  void   trialMass(int iM);
  double weightMass(int iM);

};
 
//**************************************************************************

// A derived class with 2 -> 1 kinematics set up in tau, y.

class PhaseSpace2to1tauy : public PhaseSpace {

public:

  // Constructor.
  PhaseSpace2to1tauy() {}

  // Optimize subsequent kinematics selection.
  virtual bool setupSampling() {if (!setupMass()) return false;
    return setupSampling123(false, false);} 

  // Construct the trial kinematics.
  virtual bool trialKin(bool inEvent = true) {wtBW = 1.; 
    return trialKin123(false, false, inEvent);}

  // Construct the final event kinematics.
  virtual bool finalKin();

private:

  // Set up allowed mass range.
  bool setupMass();

};
 
//**************************************************************************

// A derived class with 2 -> 2 kinematics set up in tau, y, z = cos(theta).

class PhaseSpace2to2tauyz : public PhaseSpace {

public:

  // Constructor.
  PhaseSpace2to2tauyz() {}

  // Optimize subsequent kinematics selection.
  virtual bool setupSampling() {if (!setupMasses()) return false; 
    return setupSampling123(true, false);} 

  // Construct the trial kinematics.
  virtual bool trialKin(bool inEvent = true) {
    if (!trialMasses()) return false; 
    return trialKin123(true, false, inEvent);}

  // Construct the final event kinematics.
  virtual bool finalKin();

private:

  // Set up for fixed or Breit-Wigner mass selection.
  bool setupMasses();

  // Select fixed or Breit-Wigner-distributed masses.
  bool trialMasses();

  // Pick off-shell initialization masses when on-shell not allowed.
  bool constrainedM3M4();
  bool constrainedM3();
  bool constrainedM4();

};
 
//**************************************************************************

// A derived class with 2 -> 2 kinematics set up for elastic or 
// diffractive scattering.

class PhaseSpace2to2eldiff : public PhaseSpace {

public:

  // Constructor.
  PhaseSpace2to2eldiff(bool diffAin = false, bool diffBin = false)
    : diffA(diffAin), diffB(diffBin) {}

  // Construct the trial or final event kinematics.
  virtual bool setupSampling(); 
  virtual bool trialKin(bool inEvent = true); 
  virtual bool finalKin(); 

  // Are beam particles resolved in partons or scatter directly?
  virtual bool isResolved() const {return false;}

private:

  // Constants: could only be changed in the code itself.
  static const int    NTRY;
  static const double EXPMAX, DIFFMASSMAX;

  // Kinematics properties specific to 2 -> 2.
  bool   diffA, diffB;
  double m3ElDiff, m4ElDiff, cRes, sResXB, sResAX, sProton,
         s1, s2, bMin, lambda12, lambda34, tLow, tUpp, tAux;

};
 
//**************************************************************************

// A derived class for minumum bias event. Hardly does anything, since
// the real action is taken care of by the MultipleInteractions class.

class PhaseSpace2to2minbias : public PhaseSpace {

public:

  // Constructor.
  PhaseSpace2to2minbias() {}

  // Construct the trial or final event kinematics.
  virtual bool setupSampling() {sigmaNw = sigmaProcessPtr->sigmaHat();
    sigmaMx = sigmaNw; return true;}
  virtual bool trialKin( bool ) {return true;}  
  virtual bool finalKin() {return true;}

private:

};
 
//**************************************************************************

// A derived class with 2 -> 3 kinematics 1 + 2 -> 3 + 4 + 5 set up in 
// tau, y, pT2_4, pT2_5, phi_4, phi_5 and y_3 (partial cylindrical symmetry).

class PhaseSpace2to3tauycyl : public PhaseSpace {

public:

  // Constructor.
  PhaseSpace2to3tauycyl() {}

  // Optimize subsequent kinematics selection.
  virtual bool setupSampling() {if (!setupMasses()) return false; 
    setup3Body(); return setupSampling123(false, true);} 

  // Construct the trial kinematics.
  virtual bool trialKin(bool inEvent = true) {
    if (!trialMasses()) return false; 
    return trialKin123(false, true, inEvent);}

  // Construct the final event kinematics.
  virtual bool finalKin();

private:

  // Constants: could only be changed in the code itself.
  static const int    NITERNR;

  // Set up for fixed or Breit-Wigner mass selection.
  bool setupMasses();

  // Select fixed or Breit-Wigner-distributed masses.
  bool trialMasses();

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_PhaseSpace_H
 
