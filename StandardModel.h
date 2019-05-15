// This file gives access to some Standard Model parameters.
// AlphaStrong: fix or first- or second-order running alpha_strong.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_StandardModel_H
#define Pythia8_StandardModel_H

#include "Stdlib.h"
#include "ParticleData.h"

namespace Pythia8 {

//**************************************************************************

// The AlphaStrong class calculates the alpha_strong value at an arbitrary 
// scale, given the value at m_Z, to zeroth, first or second order.

class AlphaStrong {

public:

  // Constructors.
  AlphaStrong() : isInit(false) {}
  AlphaStrong(double valueIn, int orderIn = 1) { 
    init( valueIn, orderIn) ;}

  // Initialization for given value at M_Z and given order.
  void init(double valueIn = 0.12, int orderIn = 1);

  // alpha_S value and Lambda values.
  double alphaS(double scale2);
  double alphaS1Ord(double scale2);
  double alphaS2OrdCorr(double scale2);
  double Lambda3() const { return Lambda3Save; }
  double Lambda4() const { return Lambda4Save; }
  double Lambda5() const { return Lambda5Save; }

private:

  // Constants: could only be changed in the code itself.
  static const int NITER;

  // Data members.
  bool isInit, lastCallToFull;
  int order;
  double valueRef, valueNow, scale2Now, Lambda3Save, Lambda4Save, 
    Lambda5Save, Lambda3Save2, Lambda4Save2, Lambda5Save2, mc, mb, mZ, 
    mc2, mb2;

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_StandardModel_H
