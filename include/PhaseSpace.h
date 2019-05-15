// Header file for phase space generators in kinematics selection.
// PhaseSpace: base class for phase space generators.
// PhaseSpace2to2tauyz: 2 -> 2 processes in tau, y, z = cos(theta).
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_PhaseSpace_H
#define Pythia8_PhaseSpace_H

#include "Basics.h"
#include "Beams.h"
#include "Information.h"
#include "MultipleInteractions.h"
#include "ParticleData.h"
#include "PartonDistributions.h"
#include "PythiaStdlib.h"
#include "SigmaProcess.h"
#include "SigmaTotal.h"
#include "Settings.h"

namespace Pythia8 {
 
//**************************************************************************

// PhaseSpace is a base class for  phase space generators 
// used in the selection of hard-process kinematics.

class PhaseSpace {

public:

  // Destructor.
  virtual ~PhaseSpace() {}

  // Initialize static data members.
  static void initStatic();

  // Store pointer to SigmaTotal and info on beams.
  void setSigmaTotalPtr(SigmaTotal* sigmaTotPtrIn, int idAIn, int idBIn,
    double mAIn, double mBIn) { sigmaTotPtr = sigmaTotPtrIn; idA = idAIn;
    idB = idBIn; mA = mAIn; mB = mBIn;}

  // Give in pointer to cross section and cm energy.
  void initInfo(SigmaProcess* sigmaProcessPtrIn, double eCMIn);

  // A pure virtual method, wherein an optimization procedure
  // is used to determine how phase space should be sampled.
  virtual bool setupSampling() = 0; 

  // A pure virtual method, wherein a trial event kinematics 
  // is to be selected in the derived class
  virtual bool trialKin() = 0; 

  // A pure virtual method, wherein the accepted event kinematics 
  // is to be constructed in the derived class
  virtual bool finalKin() = 0; 

  // Give back current or maximum cross section.
  double sigmaNow() const {return sigmaNw;}
  double sigmaMax() const {return sigmaMx;}

  // Give back constructed four-vectors and known masses.
  Vec4 p(int i) const {return pH[i];} 
  double m(int i) const {return mH[i];}

  // Give back other event properties.
  double x1() const {return x1H;}
  double x2() const {return x2H;}
  double sHat() const {return sH;}
  double tHat() const {return tH;}
  double uHat() const {return uH;}
  double pTHat() const {return pTH;}
  double thetaHat() const {return theta;}
  double phiHat() const {return phi;}

  // Are beam particles resolved in partons or scatter directly?
  virtual bool isResolved() const {return true;}

protected:

  // Constructor.
  PhaseSpace() {}

  // Static initialization data, normally only set once.
  static double mHatMin, mHatMax, pTHatMin, pTHatMax, m3Min, m3Max,
    m4Min, m4Max, sHatMin, sHatMax, pT2HatMin, pT2HatMax, m3SMin, 
    m3SMax, m4SMin, m4SMax;
  static bool showSearch;

  // Constants: could only be changed in the code itself.
  static const double SAFETYMARGIN, TINY, EVENFRAC, SAMESIGMA, 
    WIDTHMARGIN, SAMEMASS;

  // Center-of-mass energy.
  double eCM, s; 

  // Cross section information.
  double sigmaNw, sigmaMx;

  // Pointer to cross section. 
  SigmaProcess* sigmaProcessPtr; 

  // Event-specific kinematics properties, almost always available.
  double x1H, x2H, m3, m4, m3S, m4S, mHat, sH, tH, uH, 
    pAbs, p2Abs, pTH, theta, phi, betaZ;
  Vec4 pH[6];
  double mH[6];

  // Much common code for normal 2 -> 1 and 2 -> 2 (2 -> 1/2) machineries:

  // Determine how phase space should be sampled.
  bool setupSampling1or2(bool is2); 

  // Select a trial kinematics phase space point.
  bool trialKin1or2(bool is2); 

  // Presence and properties of any s-channel resonances.
  int idResA, idResB;
  double mResA, mResB, GammaResA, GammaResB, tauResA, tauResB, widResA, 
    widResB, tRatA, tRatB, aUppA, aUppB, aLowA, aLowB;
  bool sameResMass;

  // Kinematics properties specific to 2 -> 1/2.
  double tau, y, z, tauMin, tauMax, yMax, zMin, zMax, ratio34, unity34,
   zNeg, zPos, wtTau, wtY, wtZ, intTau0, intTau1, intTau2, intTau3, 
   intTau4, intTau5, intY01, intY2; 

  // Coefficients for optimized selection in 2 -> 1/2.
  double tauCoef[8], yCoef[8], zCoef[8];

  // Calculate kinematical limits for 2 -> 1/2.
  bool limitTau(bool is2);
  bool limitY();
  bool limitZ();

  // Select kinematical variable between defined limits for 2 -> 1/2.
  void selectTau(int iTau, double tauVal, bool is2);
  void selectY(int iY, double yVal);
  void selectZ(int iZ, double zVal);

  // Solve equation system for better phase space coefficients in 2 -> 1/2.
  void solveSys( int n, int bin[8], double vec[8], double mat[8][8],
    double coef[8]); 

  // For elastic/diffractive/minbias need information on incoming beams.
  int idA, idB;
  double mA, mB; 
  
  // And also pointer to the total/elastic/diffractive cross section.
  SigmaTotal* sigmaTotPtr;

};
 
//**************************************************************************

// A derived class with 2 -> 1 kinematics set up in tau, y.

class PhaseSpace2to1tauy : public PhaseSpace {

public:

  // Constructor.
  PhaseSpace2to1tauy() {}

  // Construct the trial kinematics, using methods in base class.
  virtual bool setupSampling() {return setupSampling1or2(false);} 
  virtual bool trialKin() {return trialKin1or2(false);}

  // Construct the final event kinematics.
  virtual bool finalKin();

private:

};
 
//**************************************************************************

// A derived class with 2 -> 2 kinematics set up in tau, y, z = cos(theta).

class PhaseSpace2to2tauyz : public PhaseSpace {

public:

  // Constructor.
  PhaseSpace2to2tauyz() {}

  // Construct the trial kinematics, using methods in base class.
  virtual bool setupSampling() {return setupSampling1or2(true);} 
  virtual bool trialKin() {return trialKin1or2(true);}

  // Construct the final event kinematics.
  virtual bool finalKin();

private:

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
  virtual bool trialKin(); 
  virtual bool finalKin(); 

  // Are beam particles resolved in partons or scatter directly?
  virtual bool isResolved() const {return false;}

private:

  // Constants: could only be changed in the code itself.
  static const int NTRY;
  static const double EXPMAX, DIFFMASSMAX;

  // Kinematics properties specific to 2 -> 2.
  bool diffA, diffB;
  double m3ElDiff, m4ElDiff, cRes, sResXB, sResAX, sProton,
    s1, s2, s3, s4, bMin, lambda12, lambda34, tLow, tUpp, tAux;

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
  virtual bool trialKin() {return true;}  
  virtual bool finalKin() {return true;}

private:

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_PhaseSpace_H
 
