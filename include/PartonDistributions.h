// Header file for parton densities.
// PDF: base class.
// GRV94L: derived class for the GRV 94L parton densities.
// CTEQ5L: derived class for the CTEQ 5L parton densities.
// Lepton: derived class for parton densities inside a lepton.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_PartonDistributions_H
#define Pythia8_PartonDistributions_H

#include "Basics.h"
#include "ParticleData.h"
#include "Stdlib.h"

namespace Pythia8 {

//**************************************************************************

// Base class for parton distribution functions.

class PDF {

public:

  // Constructor.
  PDF(int idBeamIn = 2212) {idBeam = idBeamIn; idSav = 9; xSav = -1.; 
    Q2Sav = -1.; isInit = false;}

  // Destructor.
  virtual ~PDF() {}

  // Read out parton density
  double xf(int, double, double);

  // Read out valence and sea part of parton densities.
  double xfVal(int, double, double);
  double xfSea(int, double, double);
  
protected:

  // Store relevant quantities.
  int idBeam, idSav;
  double xSav, Q2Sav;
  double xu, xd, xubar, xdbar, xs, xc, xb, xg, xlepton, xgamma,
    xuVal, xuSea, xdVal, xdSea;
  bool isInit; 

  // Update parton densities.
  virtual void xfUpdate(int id, double x, double Q2) = 0; 

};

//*********

// Gives the GRV 94 L (leading order) parton distribution function set
// in parametrized form. Authors: M. Glueck, E. Reya and A. Vogt.

class GRV94L : public PDF {

public:

  GRV94L(int idBeamIn = 2212) : PDF(idBeamIn) {;}

private:

  void xfUpdate(int id, double x, double Q2);
  double grvv (double x, double n, double ak, double bk, double a, 
    double b, double c, double d);
  double grvw (double x, double s, double al, double be, double ak, 
    double bk, double a, double b, double c, double d, double e, double es);
  double grvs (double x, double s, double sth, double al, double be, 
    double ak, double ag, double b, double d, double e, double es);

};

//*********

// Gives the GRV 94 L (leading order) parton distribution function set
// in parametrized form. Authors: M. Glueck, E. Reya and A. Vogt.

class CTEQ5L : public PDF {

public:

  CTEQ5L(int idBeamIn = 2212) : PDF(idBeamIn) {;}

private:

  void xfUpdate(int id, double x, double Q2);

};

//*********

// Gives electron (or muon, or tau) parton distribution.
 
class Lepton : public PDF {

public:

  Lepton(int idBeamIn = 11) : PDF(idBeamIn) {;}

private:

  void xfUpdate(int id, double x, double Q2);

};

} // end namespace Pythia8
 
//**************************************************************************

#endif // Pythia8_PartonDistributions_H
