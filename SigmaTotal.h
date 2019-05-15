// This file contains the class for cross section parametrizations.
// SigmaTotal: total and partial cross section in hadron-hadron collisions.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_SigmaTotal_H
#define Pythia8_SigmaTotal_H

#include "Stdlib.h"
#include "Settings.h"
#include "ParticleData.h"

namespace Pythia8 {
 
//**************************************************************************

// The SigmaTotal class contains parametrizations of total, elastic and 
// diffractive cross sections, and of the respective slope parameter.

class SigmaTotal {

public:

  // Constructor, also with incoming beams and CM energy.
  SigmaTotal() {};
  SigmaTotal(int idA, int idB, double eCM) { init(idA, idB, eCM) ;}

  // Initialize static data members.
  static void initStatic();

  // Calculate, or recalculate for new beams or new energy.
  bool init(int, int, double); 

  // Read out total and partial cross sections.
  double sigmaTot() const {return sigTot;}
  double sigmaEl() const {return sigEl;}
  double sigmaXB() const {return sigXB;}
  double sigmaAX() const {return sigAX;}
  double sigmaXX() const {return sigXX;}
  double sigmaND() const {return sigND;}

  // Read out slope b in exp(b*t) dependence.
  double bSlopeEl() const {return bEl;}
  double bSlopeXB(double sX) const { return 2.*bB + alP2 * log(s/sX) ;}
  double bSlopeAX(double sX) const { return 2.*bA + alP2 * log(s/sX) ;} 
  double bSlopeXX(double sX1, double sX2) const { 
    return alP2 * log( exp(4.) + s * s0 / (sX1 * sX2) ) ;}   

private:

  // Static initialization data, normally only set once.
  static bool setOwn;
  static double sigTotOwn, sigElOwn, sigXBOwn, sigAXOwn, sigXXOwn;

  // Constants: could only be changed in the code itself.
  static const int IHADATABLE[], IHADBTABLE[], ISDTABLE[], IDDTABLE[];
  static const double MMIN, EPSILON, ETA, X[], Y[], BETA0[], BHAD[],
    ALPHAPRIME, CONVERTEL, CONVERTSD, CONVERTDD, MMIN0, CRES, MRES0, 
    CSD[10][8], CDD[10][9], SPROTON;

  // Store values found by init.
  double sigTot, sigEl, sigXB, sigAX, sigXX, sigND, bEl, s, bA, bB,
    alP2, s0, exp4;

};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_SigmaTotal_H
