// Header file for Higgs process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma2Process.
// Copyright C 2007 Torbjorn Sjostrand

#ifndef Pythia8_SigmaHiggs_H
#define Pythia8_SigmaHiggs_H

#include "SigmaProcess.h"

namespace Pythia8 {
 
//**************************************************************************

// A derived class for f fbar -> H0 (Standard Model Higgs).

class Sigma1ffbar2H : public Sigma1Process {

public:

  // Constructor.
  Sigma1ffbar2H() {}

  // Initialize process. 
  virtual void initProc(); 

  // Initialize parton-flux object. 
  virtual void initFlux(); 

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()    const {return "f fbar -> H0 (SM)";}
  virtual int code()       const {return 801;}
  virtual int resonanceA() const {return 25;}

private:

  // A H0 resonance object provides coupling and propagator expressions.
  ResonanceH HRes;
  double mRes, GammaRes, m2Res, GamMRat;

};
 
//**************************************************************************

// A derived class for g g -> H0 (Standard Model Higgs).

class Sigma1gg2H : public Sigma1Process {

public:

  // Constructor.
  Sigma1gg2H() {}

  // Initialize process. 
  virtual void initProc(); 

  // Initialize parton-flux object. 
  virtual void initFlux(); 

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()    const {return "g g -> H0 (SM)";}
  virtual int code()       const {return 802;}
  virtual int resonanceA() const {return 25;}

private:

  // A H0 resonance object provides coupling and propagator expressions.
  ResonanceH HRes;
  double mRes, GammaRes, m2Res, GamMRat;

};
 
//**************************************************************************

// A derived class for gamma gamma -> H0 (Standard Model Higgs).

class Sigma1gmgm2H : public Sigma1Process {

public:

  // Constructor.
  Sigma1gmgm2H() {}

  // Initialize process. 
  virtual void initProc(); 

  // Initialize parton-flux object. 
  virtual void initFlux(); 

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()    const {return "gamma gamma -> H0 (SM)";}
  virtual int code()       const {return 803;}
  virtual int resonanceA() const {return 25;}

private:

  // A H0 resonance object provides coupling and propagator expressions.
  ResonanceH HRes;
  double mRes, GammaRes, m2Res, GamMRat;

};
 
//**************************************************************************

// A derived class for f fbar -> H0 Z0 (Standard Model Higgs).

class Sigma2ffbar2HZ : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2HZ() {}

  // Initialize process. 
  virtual void initProc(); 

  // Initialize parton-flux object. 
  virtual void initFlux(); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()    const {return "f fbar -> H0 Z0 (SM)";}
  virtual int code()       const {return 804;}
  virtual int id3Mass()    const {return 25;}
  virtual int id4Mass()    const {return 23;}
  virtual int resonanceA() const {return 23;}
  virtual int gmZmode()    const {return 2;}


private:

  // Store Z0 mass and width.
  double mZ, widZ, mZS, mwZS, thetaWRat;

};
 
//**************************************************************************

// A derived class for f fbar -> H0 W+- (Standard Model Higgs).

class Sigma2ffbar2HW : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2HW() {}

  // Initialize process. 
  virtual void initProc(); 

  // Initialize parton-flux object. 
  virtual void initFlux(); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()    const {return "f fbar -> H0 W+- (SM)";}
  virtual int code()       const {return 805;}
  virtual int id3Mass()    const {return 25;}
  virtual int id4Mass()    const {return 24;}
  virtual int resonanceA() const {return 24;}

private:

  // Store W+- mass and width, and couplings.
  double mW, widW, mWS, mwWS, thetaWRat;

};
 
//**************************************************************************

// A derived class for q g -> H0 q (Standard Model Higgs).

class Sigma2qg2Hq : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2Hq() {}

  // Initialize process. 
  virtual void initProc(); 

  // Initialize parton-flux object. 
  virtual void initFlux(); 

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Evaluate weight for decay angles.
  virtual double weightDecay( Event& process, int iResBeg, int iResEnd); 

  // Info on the subprocess.
  virtual string name()    const {return "q g -> H0 q (SM)";}
  virtual int code()       const {return 806;}
  virtual int id3Mass()    const {return 25;}

private:

  // Store standard prefactor.
  double m2W, thetaWRat;

};
  
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_SigmaHiggs_H
