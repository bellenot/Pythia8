// Header file for hard-process differential cross sections.
// SigmaProcess: base class for cross sections.
// Sigma0Process: base class for unresolved processes, 
// derived from SigmaProcess.
// Sigma2Process: base class for 2 -> 2 processes, 
// derived from SigmaProcess.
// Actual physics processes are found in separate files:
// SigmaQCD for QCD processes;
// SigmaEW for electroweak processes (including photon production);
// SigmaSUSY for supersymmetric production.
// Copyright C 2007 Torbjorn Sjostrand

#ifndef Pythia8_SigmaProcess_H
#define Pythia8_SigmaProcess_H

#include "Basics.h"
#include "Event.h"
#include "InFlux.h"
#include "Information.h"
#include "ParticleData.h"
#include "PartonDistributions.h"
#include "PythiaStdlib.h"
#include "ResonanceProperties.h"
#include "Settings.h"
#include "SigmaTotal.h"
#include "StandardModel.h"
#include "SusyLesHouches.h"

namespace Pythia8 {
 
//**************************************************************************

// SigmaProcess is the base class for cross section calculations.

class SigmaProcess {

public:

  // Destructor.
  virtual ~SigmaProcess() {}

  // Initialize (some) static data members.
  static void initStatic();

  // Store pointer to SigmaTotal and info on beams.
  static void setSigmaTotalPtr(SigmaTotal* sigmaTotPtrIn, int idAIn, 
    int idBIn, double mAIn, double mBIn) { sigmaTotPtr = sigmaTotPtrIn; 
    idA = idAIn; idB = idBIn; mA = mAIn; mB = mBIn; 
    hasLeptonBeams = ( ParticleDataTable::isLepton(idA) 
    || ParticleDataTable::isLepton(idB) ); }

  // Store pointer to SUSY Les Houches Accord.
  static void setSlhaPtr(SusyLesHouches* slhaIn) {slha = slhaIn;}

  // Initialize process. Only used for some processes.
  virtual void initProc() {}

  // Initialize incoming flux. Not used for multiple interactions. 
  virtual void initFlux() {}

  // Evaluate sigma for unresolved, sigmaHat(sHat) for 2 -> 1 processes, 
  // and d(sigmaHat)/d(tHat) for (resolved) 2 -> 2 processes. 
  virtual double sigmaHat() {return 0.;}

  // Convolute above with parton flux. 
  virtual double sigmaPDF() {return 0.;}

  // Input and complement kinematics for resolved 2 -> 1 process.
  virtual bool set1Kin( double x1in, double x2in, double sHin)
    {return (x1in + x2in + sHin > 0.);} 

  // Input and complement kinematics for resolved 2 -> 2 process.
  virtual bool set2Kin( double x1in, double x2in, double sHin, double tHin,
    double m3in, double m4in)
    {return (x1in + x2in + sHin + tHin + m3in + m4in > 0.);} 

  // Ditto, but for Multiple Interactions applications, so different input.
  virtual bool set2KinMI( int id1in, int id2in, double x1in, double x2in,
    double sHin, double tHin, double uHin, double alpSin, double alpEMin,
    bool needMasses, double m3in, double m4in) {
    return ( (id1in + id2in + x1in + x2in + sHin + tHin + uHin + alpSin 
    + alpEMin +m3in + m4in > 0.) && needMasses);} 

  // Select incoming parton channel and extract parton densities (resolved).
  virtual void pickInState() {}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol() = 0;

  // Perform kinematics for a Multiple Interaction, in its rest frame.
  virtual bool final2KinMI() {return true;}

  // Process name and code, and the number of final-state particles.
  virtual string name() const {return "unnamed process";}
  virtual int code() const {return 0;}
  virtual int nFinal() const {return 2;};

  // Special treatment needed for elastic and diffractive processes.
  virtual bool isMinBias() const {return false;};
  virtual bool isResolved() const {return true;};
  virtual bool isDiffA() const {return false;};
  virtual bool isDiffB() const {return false;};

  // Flavours in 2 -> 2 processes where masses needed from beginning. 
  // (For a light quark masses will be used in the final kinematics,
  // but not at the matrix-element level. For a gluon no masses at all.) 
  virtual int id3Mass() const {return 0;}
  virtual int id4Mass() const {return 0;}

  // Special treatment needed if process contains an s-channel resonance.
  virtual int resonanceA() const {return 0;}
  virtual int resonanceB() const {return 0;}

  // Tell whether tHat and uHat are swapped (= same as swap 3 and 4).
  bool swappedTU() const {return swapTU;}
  
  // Give back particle properties: flavours, colours, masses, or all.
  int id(int i) const {return idH[i];}
  int col(int i) const {return colH[i];} 
  int acol(int i) const {return acolH[i];}
  double m(int i) const {return mH[i];}
  Particle getParton(int i) {return parton[i];}

  // Give back couplings and parton densities. Not all known for minbias.
  double Q2Ren() const {return Q2RenH;}
  double alphaEMH() const {return alpEM;}
  double alphaSH() const {return alpS;}
  double Q2Fac() const {return Q2FacH;}
  double pdf1() const {return pdf1H;}
  double pdf2() const {return pdf2H;}

  // Give back angles; relevant only for multipe-interactions processes.
  double thetaMI() const {return atan2( sinTheta, cosTheta);}
  double phiMI() const {return phi;}
  double sHBetaMI() const {return sHBeta;}
  double pT2MI() const {return pT2Mass;}

  // Remove empty inFlux channels and optionally list remaining ones.
  void checkChannels() {
    if (inFluxPtr != 0) inFluxPtr->checkChannels(name());}

protected:

  // Constructor.
  SigmaProcess() {inFluxPtr = 0; for (int i = 0; i < 6; ++i) mH[i] = 0.;}

  // Static initialization data, normally only set once.
  static int alphaSorder,alphaEMorder, nQuark;
  static double alphaSvalue;

  // Static alphaStrong and alphaElectromagnetic calculation.
  static AlphaStrong alphaS;
  static AlphaEM alphaEM;

  // Constants: could only be changed in the code itself.
  static const double CONVERT2MB, MASSMARGIN;

  // Static information on incoming beams.
  static int idA, idB;
  static double mA, mB; 
  static bool hasLeptonBeams;
  
  // Static pointer to the total/elastic/diffractive cross section object.
  static SigmaTotal* sigmaTotPtr;

  // Static pointer to the SLHA object.
  static SusyLesHouches* slha;

  // Store Q2 renormalization and factorization scales, and related values.
  double Q2RenH, alpEM, alpS, Q2FacH, pdf1H, pdf2H;

  // Store flavour, colour, anticolour, mass, angles and the whole particle.
  int idH[6], colH[6], acolH[6];
  double mH[6], cosTheta, sinTheta, phi, sHMass, sHBeta, pT2Mass;
  Particle parton[6];

  // Pointer to the parton-flux object.
  InFlux* inFluxPtr;

  // Store whether tHat and uHat are swapped (= same as swap 3 and 4).
  bool swapTU;

  // Set flavour, colour and anticolour.
  void setId( int id1in = 0, int id2in = 0, int id3in = 0, int id4in = 0,
    int id5in = 0) {idH[1] = id1in; idH[2] = id2in; idH[3] = id3in; 
    idH[4] = id4in; idH[5] = id5in;}
  void swapIdM34() { swap(idH[3], idH[4]); swap(mH[3], mH[4]);}
  void setColAcol( int col1 = 0, int acol1 = 0, 
    int col2 = 0, int acol2 = 0, int col3 = 0, int acol3 = 0, 
    int col4 = 0, int acol4 = 0, int col5 = 0, int acol5 = 0) {
    colH[1] = col1; acolH[1] = acol1; colH[2] = col2; acolH[2] = acol2;
    colH[3] = col3; acolH[3] = acol3; colH[4] = col4; acolH[4] = acol4;
    colH[5] = col5; acolH[5] = acol5; }
  void swapColAcol() { swap(colH[1], acolH[1]); 
    swap(colH[2], acolH[2]); swap(colH[3], acolH[3]); 
    swap(colH[4], acolH[4]); swap(colH[5], acolH[5]);}
  void swapCol1234() { swap(colH[1], colH[2]); swap(colH[3], colH[4]); 
    swap(acolH[1], acolH[2]); swap(acolH[3], acolH[4]);}
  void swapCol12() { swap(colH[1], colH[2]); swap(acolH[1], acolH[2]);}
  void swapCol34() { swap(colH[3], colH[4]); swap(acolH[3], acolH[4]);}

};
 
//**************************************************************************

// Sigma0Process is the base class for unresolved and minimum-bias processes. 
// It is derived from SigmaProcess.

class Sigma0Process : public SigmaProcess {

public:

  // Destructor.
  virtual ~Sigma0Process() {}

  // Number of final-state particles.
  virtual int nFinal() const {return 2;};

  // Evaluate sigma for unresolved processes. 
  virtual double sigmaHat() {return 0.;}

  // Since no PDF's there is no difference from above. 
  double sigmaPDF() {return sigmaHat();}

protected:

  // Constructor.
  Sigma0Process() {}

};
 
//**************************************************************************

// Sigma1Process is the base class for 2 -> 1 processes.
// It is derived from SigmaProcess.

class Sigma1Process : public SigmaProcess {

public:

  // Destructor.
  virtual ~Sigma1Process() {if (inFluxPtr !=0) delete inFluxPtr;}

  // Number of final-state particles.
  virtual int nFinal() const {return 1;};

  // Evaluate sigmaHat(sHat) for resolved 2 -> 1 processes. 
  virtual double sigmaHat() {return 0.;}

  // Convolute sigmaHat(sHat) with parton flux. 
  virtual double sigmaPDF() {
    if (inFluxPtr == 0) return 0.;
    // Note that order is important: sigmaHat may set info for inFlux.
    return sigmaHat() * inFluxPtr->flux( x1, x2, Q2FacH);}

  // Input and complement kinematics for resolved 2 -> 2 process.
  virtual bool set1Kin( double x1in, double x2in, double sHin); 

  // Select incoming parton channel and extract parton densities.
  void pickInState() {inFluxPtr->pick(); id1 = inFluxPtr->id1();
    id2 = inFluxPtr->id2(); pdf1H = inFluxPtr->pdf1();
    pdf2H = inFluxPtr->pdf2();} 

protected:

  // Constructor.
  Sigma1Process() {}

  // Store subprocess kinematics quantities.
  int id1, id2, id3;
  double x1, x2, sH, sH2;

};
 
//**************************************************************************

// Sigma2Process is the base class for 2 -> 2 processes.
// It is derived from SigmaProcess.

class Sigma2Process : public SigmaProcess {

public:

  // Destructor.
  virtual ~Sigma2Process() {if (inFluxPtr !=0) delete inFluxPtr;}

  // Number of final-state particles.
  virtual int nFinal() const {return 2;};

  // Evaluate d(sigmaHat)/d(tHat) for resolved 2 -> 2 processes. 
  virtual double sigmaHat() {return 0.;}

  // Convolute d(sigmaHat)/d(tHat) with parton flux. 
  virtual double sigmaPDF() {
    if (inFluxPtr == 0) return 0.;
    // Note that order is important: sigmaHat may set info for inFlux.
    return sigmaHat() * inFluxPtr->flux( x1, x2, Q2FacH);}

  // Input and complement kinematics for resolved 2 -> 2 process.
  virtual bool set2Kin( double x1in, double x2in, double sHin, double tHin,
    double m3in, double m4in); 

  // Ditto, but for Multiple Interactions applications, so different input.
  virtual bool set2KinMI( int id1in, int id2in, double x1in, double x2in,
    double sHin, double tHin, double uHin, double alpSin, double alpEMin,
    bool needMasses, double m3in, double m4in);

  // Select incoming parton channel and extract parton densities.
  void pickInState() {inFluxPtr->pick(); id1 = inFluxPtr->id1();
    id2 = inFluxPtr->id2(); pdf1H = inFluxPtr->pdf1();
    pdf2H = inFluxPtr->pdf2();} 

  // Perform kinematics for a Multiple Interaction, in its rest frame.
  virtual bool final2KinMI();

protected:

  // Constructor.
  Sigma2Process() {}

  // Store subprocess kinematics quantities.
  bool id12IsSet;
  int id1, id2, id3, id4;
  double x1, x2, sH, tH, uH, sH2, tH2, uH2, m3, s3, m4, s4, pT2;

};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_SigmaProcess_H
 
