// PartonDistributions.h is a part of the PYTHIA event generator.
// Copyright (C) 2009 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for parton densities.
// PDF: base class.
// GRV94L: derived class for the GRV 94L parton densities.
// CTEQ5L: derived class for the CTEQ 5L parton densities.
// LHAPDF: derived class for interface to the LHAPDF library.
// GRVpiL: derived class for the GRV LO pion parton densities.
// PomPDF: derived class for Q2-independent Pomeron parton densities.
// Lepton: derived class for parton densities inside a lepton.
// LeptonPoint: derived class for unresolved lepton (mainly dummy).

#ifndef Pythia8_PartonDistributions_H
#define Pythia8_PartonDistributions_H

#include "Basics.h"
#include "Info.h"
#include "ParticleData.h"
#include "PythiaStdlib.h"

namespace Pythia8 {

//**************************************************************************

// Base class for parton distribution functions.

class PDF {

public:

  // Constructor.
  PDF(int idBeamIn = 2212) {idBeam = idBeamIn; idBeamAbs = abs(idBeam);
    setValenceContent(); idSav = 9; xSav = -1.; Q2Sav = -1.; 
    isSet = true; isInit = false;}

  // Destructor.
  virtual ~PDF() {}

  // Confirm that PDF has been set up (important for LHAPDF).
  bool isSetup() {return isSet;}

  // Dynamic choice of meson valence flavours for pi0, K0S, K0L, Pomeron.
  void newValenceContent(int idVal1In, int idVal2In) {
    idVal1 = idVal1In; idVal2 = idVal2In;}

  // Allow extrapolation beyond boundaries. This is optional.
  virtual void setExtrapolate(bool) {}

  // Read out parton density
  double xf(int id, double x, double Q2);

  // Read out valence and sea part of parton densities.
  double xfVal(int id, double x, double Q2);
  double xfSea(int id, double x, double Q2);
  
protected:

  // Store relevant quantities.
  int    idBeam, idBeamAbs, idSav, idVal1, idVal2;
  double xSav, Q2Sav;
  double xu, xd, xubar, xdbar, xs, xc, xb, xg, xlepton, xgamma,
         xuVal, xuSea, xdVal, xdSea;
  bool   isSet, isInit, hasLimits; 

  // Resolve valence content for assumed meson. Possibly modified later.
  void setValenceContent();

  // Update parton densities.
  virtual void xfUpdate(int id, double x, double Q2) = 0; 

};
 
//**************************************************************************

// Gives the GRV 94 L (leading order) parton distribution function set
// in parametrized form. Authors: M. Glueck, E. Reya and A. Vogt.

class GRV94L : public PDF {

public:

  // Constructor.
  GRV94L(int idBeamIn = 2212) : PDF(idBeamIn) {}

private:

  // Update PDF values.
  void xfUpdate(int id, double x, double Q2);

  // Auxiliary routines used during the updating.
  double grvv (double x, double n, double ak, double bk, double a, 
    double b, double c, double d);
  double grvw (double x, double s, double al, double be, double ak, 
    double bk, double a, double b, double c, double d, double e, double es);
  double grvs (double x, double s, double sth, double al, double be, 
    double ak, double ag, double b, double d, double e, double es);

};
 
//**************************************************************************

// Gives the CTEQ 5 L (leading order) parton distribution function set
// in parametrized form. Parametrization by J. Pumplin. Authors: CTEQ.

class CTEQ5L : public PDF {

public:

  // Constructor.
  CTEQ5L(int idBeamIn = 2212) : PDF(idBeamIn) {}

private:

  // Update PDF values.
  void xfUpdate(int id, double x, double Q2);

};
 
//**************************************************************************

// Provide interface to the LHAPDF library of parton densities.

class LHAPDF : public PDF {

public:

  // Constructor.
  LHAPDF(int idBeamIn, string setName, int member,  int nSetIn = 1, 
    Info* infoPtr = 0) : PDF(idBeamIn), nSet(nSetIn) 
    {init( setName, member, infoPtr);} 

  // Allow extrapolation beyond boundaries. This is optional.
  void setExtrapolate(bool extrapol); 

private:

  // Initialization of PDF set.
  void init(string setName, int member, Info* infoPtr);

  // Update all PDF values.
  void xfUpdate(int , double x, double Q2);

  // Current set and pdf values.
  int    nSet;
  double xfArray[13];
  double xPhoton;

  // Keep track of latest initialized PDF, so does not have to repeat.
  static string latestSetName;
  static int    latestMember, latestNSet;   

};
 
//**************************************************************************

// Gives the GRV 1992 pi+ (leading order) parton distribution function set
// in parametrized form. Authors: Glueck, Reya and Vogt.

class GRVpiL : public PDF {

public:

  // Constructor.
  GRVpiL(int idBeamIn = 221) : PDF(idBeamIn) {}

private:

  // Update PDF values.
  void xfUpdate(int id, double x, double Q2);

};

//**************************************************************************

// Gives generic Q2-independent Pomeron PDF.

class PomPDF : public PDF {

public:

  // Constructor.
  PomPDF(int idBeamIn = 990, double PomGluonAIn = 0., 
  double PomGluonBIn = 0., double PomQuarkAIn = 0., double PomQuarkBIn = 0., 
  double PomQuarkFracIn = 0., double PomStrangeSuppIn = 0.) : PDF(idBeamIn), 
  PomGluonA(PomGluonAIn), PomGluonB(PomGluonBIn), PomQuarkA(PomQuarkAIn), 
  PomQuarkB(PomQuarkBIn), PomQuarkFrac(PomQuarkFracIn), 
  PomStrangeSupp(PomStrangeSuppIn) {init();}

private:

  // Stored value for PDF choice.
  double PomGluonA, PomGluonB, PomQuarkA, PomQuarkB, PomQuarkFrac, 
         PomStrangeSupp, normGluon, normQuark;

  // Initialization of some constants.
  void init();

  // Update PDF values.
  void xfUpdate(int id, double x, double);

};
 
//**************************************************************************

// Gives electron (or muon, or tau) parton distribution.
 
class Lepton : public PDF {

public:

  // Constructor.
  Lepton(int idBeamIn = 11) : PDF(idBeamIn) {}

private:

  // Constants: could only be changed in the code itself.
  static const double ALPHAEM;

  // Update PDF values.
  void xfUpdate(int id, double x, double Q2);

  // The lepton mass, set at initialization.
  double m2Lep;

};
 
//**************************************************************************

// Gives electron (or other lepton) parton distribution when unresolved.
 
class LeptonPoint : public PDF {

public:

  // Constructor.
  LeptonPoint(int idBeamIn = 11) : PDF(idBeamIn) {}

private:

  // Update PDF values in trivial way. 
  void xfUpdate(int , double , double ) {xlepton = 1; xgamma = 0.;}

};

} // end namespace Pythia8
 
//**************************************************************************

#endif // Pythia8_PartonDistributions_H
