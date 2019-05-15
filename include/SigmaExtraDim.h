// SigmaExtraDim.h is a part of the PYTHIA event generator.
// Copyright (C) 2009 Torbjorn Sjostrand.
// Copyright (C) 2009 Stefan Ask for the *LED* routines.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for extra-dimensional-process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma(1/2)Process.

#ifndef Pythia8_SigmaExtraDim_H
#define Pythia8_SigmaExtraDim_H

#include "SigmaProcess.h"

namespace Pythia8 {
 
//**************************************************************************

// A derived class for g g -> G^* (excited graviton state).

class Sigma1gg2GravitonStar : public Sigma1Process {

public:

  // Constructor.
  Sigma1gg2GravitonStar() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for G* decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()       const {return "g g -> G*";}
  virtual int    code()       const {return 5001;}
  virtual string inFlux()     const {return "gg";}
  virtual int    resonanceA() const {return idGstar;}

private:

  // Parameters set at initialization or for current kinematics. 
  int    idGstar;
  double mRes, GammaRes, m2Res, GamMRat, kappaMG, sigma;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* gStarPtr;

};
 
//**************************************************************************

// A derived class for f fbar -> G^* (excited graviton state).

class Sigma1ffbar2GravitonStar : public Sigma1Process {

public:

  // Constructor.
  Sigma1ffbar2GravitonStar() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat() {return (abs(id1) < 9) ? sigma0 / 3. : sigma0;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for G* decay angle.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()       const {return "f fbar -> G*";}
  virtual int    code()       const {return 5002;}
  virtual string inFlux()     const {return "ffbarSame";}
  virtual int    resonanceA() const {return idGstar;}

private:

  // Parameters set at initialization or for current kinematics. 
  int    idGstar;
  double mRes, GammaRes, m2Res, GamMRat, kappaMG, sigma0;

  // Pointer to properties of the particle species, to access decay channels.
  ParticleDataEntry* gStarPtr;

};
 
//**************************************************************************

// A derived class for g g -> G^* g (excited graviton state).

class Sigma2gg2GravitonStarg : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2GravitonStarg() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight: currently isotropic (except secondary top decay)..
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()    const {return "g g -> G* g";}
  virtual int    code()    const {return 5003;}
  virtual string inFlux()  const {return "gg";}
  virtual int    id3Mass() const {return idGstar;}

private:

  // Parameters set at initialization or for current kinematics. 
  int    idGstar;
  double mRes, GammaRes, m2Res, GamMRat, kappaMG, openFrac, sigma;

};
 
//**************************************************************************

// A derived class for q g -> G^* q (excited graviton state).

class Sigma2qg2GravitonStarq : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2GravitonStarq() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight: currently isotropic (except secondary top decay)..
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()    const {return "q g -> G* q";}
  virtual int    code()    const {return 5004;}
  virtual string inFlux()  const {return "qg";}
  virtual int    id3Mass() const {return idGstar;}

private:

  // Parameters set at initialization or for current kinematics. 
  int    idGstar;
  double mRes, GammaRes, m2Res, GamMRat, kappaMG, openFrac, sigma;

};
 
//**************************************************************************

// A derived class for q qbar -> G^* g (excited graviton state).

class Sigma2qqbar2GravitonStarg : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2GravitonStarg() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight: currently isotropic (except secondary top decay)..
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()    const {return "q qbar -> G* g";}
  virtual int    code()    const {return 5005;}
  virtual string inFlux()  const {return "qqbarSame";}
  virtual int    id3Mass() const {return idGstar;}

private:

  // Parameters set at initialization or for current kinematics. 
  int    idGstar;
  double mRes, GammaRes, m2Res, GamMRat, kappaMG, openFrac, sigma;

};
 
//**************************************************************************

// A derived class for g g -> G g (real graviton emission in 
// large extra dimensions). 

class Sigma2gg2LEDGravitong : public Sigma2Process {

public:

  //+++ Constructor.
  Sigma2gg2LEDGravitong() {}

  //+++ Initialize process. 
  virtual void initProc(); 

  //+++ Calculate flavour-independent parts of cross section; 
  //+++ first step when inflavours unknown. 
  virtual void sigmaKin();

  //+++ Evaluate sigmaHat(sHat); second step for given inflavours. 
  virtual double sigmaHat();

  //+++ Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  //+++ Info on the subprocess.
  virtual string name()       const {return "g g -> G g";}
  virtual int    code()       const {return 5021;}
  virtual string inFlux()     const {return "gg";}
  virtual int    id3Mass()    const {return 5000039;}
  virtual int    id4Mass()    const {return 21;}

private:

  /** 
   * Function GAMMA(n/2), n = integer.
   * n = 2k   ==> Gamma = (k-1)!
   * n = 2k+1 ==> Gamma = sqrt(PI) * ProdSum _(l=0)^(k-1) [l + 1/2]
   */
  double funcGammaIntHalf( int );

private:

  bool m_trunc;
  int  m_nGrav, m_idG;
  double mG, mGS, m_sigma0, m_MD, m_constantTerm;

};

//**************************************************************************

// A derived class for q g -> G q (real graviton emission in 
// large extra dimensions). 

class Sigma2qg2LEDGravitonq : public Sigma2Process {

public:

  //+++ Constructor.
  Sigma2qg2LEDGravitonq() {} 

  //+++ Initialize process. 
  virtual void initProc(); 

  //+++ Calculate flavour-independent parts of cross section; 
  // first step when inflavours unknown. 
  virtual void sigmaKin();

  //+++ Evaluate sigmaHat(sHat); second step for given inflavours. 
  virtual double sigmaHat();

  //+++ Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  //+++ Info on the subprocess.
  virtual string name()       const {return "q g -> G q";}
  virtual int    code()       const {return 5022;}
  virtual string inFlux()     const {return "qg";}
  virtual int    id3Mass()    const {return 5000039;}
  //  virtual int    id4Mass()    const {return 21;}

private:

  /** 
   * Function GAMMA(n/2), n = integer.
   * n = 2k   ==> Gamma = (k-1)!
   * n = 2k+1 ==> Gamma = sqrt(PI) * ProdSum _(l=0)^(k-1) [l + 1/2]
   */
  double funcGammaIntHalf( int );

private: 

  bool m_trunc;
  int  m_nGrav, m_idG;
  double mG, mGS, m_sigma0, m_MD, m_constantTerm;
  
};

//**************************************************************************

// A derived class for q qbar -> G g (real graviton emission in 
// large extra dimensions). 

class Sigma2qqbar2LEDGravitong : public Sigma2Process {

public:

  //+++ Constructor.
  Sigma2qqbar2LEDGravitong( ) {}

  //+++ Initialize process. 
  virtual void initProc(); 

  //+++ Calculate flavour-independent parts of cross section; 
  //+++ first step when inflavours unknown. 
  virtual void sigmaKin();

  //+++ Evaluate sigmaHat(sHat); second step for given inflavours. 
  virtual double sigmaHat();

  //+++ Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  //+++ Info on the subprocess.
  virtual string name()       const {return "q qbar -> G g";}
  virtual int    code()       const {return 5023;}
  virtual string inFlux()     const {return "qqbarSame";}
  virtual int    id3Mass()    const {return 5000039;}
  virtual int    id4Mass()    const {return 21;}

private:

  /** 
   * Function GAMMA(n/2), n = integer.
   * n = 2k   ==> Gamma = (k-1)!
   * n = 2k+1 ==> Gamma = sqrt(PI) * ProdSum _(l=0)^(k-1) [l + 1/2]
   */
  double funcGammaIntHalf( int );

private: 

  bool m_trunc;
  int  m_nGrav, m_idG;
  double mG, mGS, m_sigma0, m_MD, m_constantTerm;
 
};

//**************************************************************************

// A derived class for f fbar -> U/G g (real graviton emission in 
// large extra dimensions or unparticle emission). 

class Sigma2ffbar2LEDUnparticleZ : public Sigma2Process {

public:

  /**
   * Constructor:
   * bool Graviton  = true, to use graviton specific settings
   */
  Sigma2ffbar2LEDUnparticleZ( bool ); 

  //+++ Initialize process. 
  virtual void initProc(); 

  //+++ Calculate flavour-independent parts of cross section; 
  //+++ first step when inflavours unknown. 
  virtual void sigmaKin();

  //+++ Evaluate sigmaHat(sHat); second step for given inflavours.
  virtual double sigmaHat();

  //+++ Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  //+++ Info on the subprocess.
  virtual string name()       const {return 
    (m_graviton ? "f fbar -> G Z" : "f fbar -> U Z") ;}
  virtual int    code()       const {return (m_graviton ? 5024 : 5041);}
  virtual string inFlux()     const {return "ffbarSame";}
  virtual int    id3Mass()    const {return 5000039;} 
  virtual int    id4Mass()    const {return 23;}
  virtual int    resonanceA() const {return 23;}
  virtual int    gmZmode()    const {return 2;}

private:

  /** 
   * Function GAMMA(n/2), n = integer. 
   * n = 2k   ==> Gamma = (k-1)!
   * n = 2k+1 ==> Gamma = sqrt(PI) * ProdSum _(l=0)^(k-1) [l + 1/2]
   */
  double funcGammaIntHalf( int );

  /** 
   * Function GAMMA(x, n), 
   * x = real value, 
   * n = terms in product sum. 
   * For x != negative integer or 0 (singularities).
   * http://en.wikipedia.org/wiki/Gamma_function
   */
  double funcGammaReal( double, int );

private:

  // Constants: could only be changed in the code itself.
  static const double FIXRATIO; 

  int    m_spin, m_idG;
  bool   m_graviton, m_trunc;
  double m_dU, m_LambdaU, m_lambda, m_ratio, m_lambdaPrime, 
         m_constantTerm, m_nGrav;
  double sHS, tHS, uHS, tHC, uHC, tHQ, uHQ, tHuH, mU, mUS, mZ, widZ, 
         mZS, mwZS, thetaWRat, m_sigma0;

};

//**************************************************************************

// A derived class for f fbar -> U/G gamma (real graviton emission in 
// large extra dimensions or unparticle emission). 

class Sigma2ffbar2LEDUnparticlegamma : public Sigma2Process {

public:

  /**
   * Constructor:
   * bool Graviton  = true, to use graviton specific settings
   */
  Sigma2ffbar2LEDUnparticlegamma( bool ); 

  //+++ Initialize process. 
  virtual void initProc(); 

  //+++ Calculate flavour-independent parts of cross section; 
  //+++ first step when inflavours unknown. 
  virtual void sigmaKin();

  //+++ Evaluate sigmaHat(sHat); second step for given inflavours. 
  virtual double sigmaHat();

  //+++ Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  //+++ Info on the subprocess.
  virtual string name()       const {return 
    (m_graviton ? "f fbar -> G gamma" : "f fbar -> U gamma") ;}
  virtual int    code()       const {return (m_graviton ? 5025 : 5042);}
  virtual string inFlux()     const {return "ffbarSame";}
  virtual int    id3Mass()    const {return 5000039;}
  virtual int    id4Mass()    const {return 22;}

private:

  /** 
   * Function GAMMA(n/2), n = integer.
   * n = 2k   ==> Gamma = (k-1)!
   * n = 2k+1 ==> Gamma = sqrt(PI) * ProdSum _(l=0)^(k-1) [l + 1/2]
   */
  double funcGammaIntHalf( int );

  /** 
   * Function GAMMA(x, n), 
   * x = real value, 
   * n = terms in product sum. 
   * For x != negative integer or 0 (singularities).
   * http://en.wikipedia.org/wiki/Gamma_function
   */
  double funcGammaReal( double, int );

private:

  // Constants: could only be changed in the code itself.
  static const double FIXRATIO; 

  int    m_spin, m_idG;
  bool   m_graviton;
  bool   m_trunc;
  double m_dU, m_LambdaU, m_lambda, m_ratio, m_lambdaPrime, 
         m_constantTerm, m_nGrav;  
  double sHS, tHS, uHS, tHC, uHC, tHQ, uHQ, tHuH, mU, mUS, mZ, widZ, 
         mZS, mwZS, thetaWRat, m_sigma0;

};

//**************************************************************************

// A derived class for f fbar -> (LED G*/U*) -> gamma gamma 
// (virtual graviton/unparticle exchange). 

class Sigma2ffbar2LEDgammagamma : public Sigma2Process {

public:

  /**
   * Constructor:
   * bool Graviton  = true, to use LED graviton settings
   */
  Sigma2ffbar2LEDgammagamma( bool );

  //+++ Initialize process. 
  virtual void initProc(); 

  //+++ Calculate flavour-independent parts of cross section; 
  //+++ first step when inflavours unknown. 
  virtual void sigmaKin();

  //+++ Evaluate sigmaHat(sHat); second step for given inflavours. 
  virtual double sigmaHat();

  //+++ Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  //+++ Info on the subprocess.
  virtual string name()       const {return 
    (m_graviton ? "f fbar -> (LED G*) -> gamma gamma" : "f fbar -> (U*) -> gamma gamma") ;}
  virtual int    code()       const {return (m_graviton ? 5026 : 5043);}
  virtual string inFlux()     const {return "ffbarSame";}

private:

  /** 
   * Function GAMMA(x, n), 
   * x = real value, 
   * n = terms in product sum. 
   * For x != negative integer or 0 (singularities).
   * http://en.wikipedia.org/wiki/Gamma_function
   */
  double funcGammaReal( double, int );

private:

  int    m_spin;
  bool   m_graviton;
  double m_dU, m_LambdaU, m_lambda, m_lambda2chi, 
         m_term1, m_term2, m_term3;
  
};

//**************************************************************************

// A derived class for g g -> (LED G*/U*) -> gamma gamma 
// (virtual graviton/unparticle exchange). 

class Sigma2gg2LEDgammagamma : public Sigma2Process {

public:

  /**
   * Constructor:
   * bool Graviton  = true, to use LED graviton settings
   */
  Sigma2gg2LEDgammagamma( bool );

  //+++ Initialize process. 
  virtual void initProc(); 

  //+++ Calculate flavour-independent parts of cross section; 
  //+++ first step when inflavours unknown. 
  virtual void sigmaKin();

  //+++ Evaluate sigmaHat(sHat); second step for given inflavours. 
  virtual double sigmaHat();

  //+++ Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  //+++ Info on the subprocess.
  virtual string name()       const {return 
    (m_graviton ? "g g -> (LED G*) -> gamma gamma" : "g g -> (U*) -> gamma gamma") ;}
  virtual int    code()       const {return (m_graviton ? 5027 : 5044);}
  virtual string inFlux()     const {return "gg";}

private:

  /** 
   * Function GAMMA(x, n), 
   * x = real value, 
   * n = terms in product sum. 
   * For x != negative integer or 0 (singularities).
   * http://en.wikipedia.org/wiki/Gamma_function
   */
  double funcGammaReal( double, int );

private:

  int    m_spin;
  bool   m_graviton;
  double m_dU, m_LambdaU, m_lambda, m_lambda2chi, m_sigma0;
  
};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_SigmaExtraDim_H
