// SigmaSUSY.h is a part of the PYTHIA event generator.
// Copyright (C) 2008 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for Supersymmetric process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma2Process.

#ifndef Pythia8_SigmaSUSY_H
#define Pythia8_SigmaSUSY_H

#include "PythiaComplex.h"
#include "SigmaProcess.h"

namespace Pythia8 {
 
//**************************************************************************

// A derived class for q qbar -> neutralino_i neutralino_j.

class Sigma2qqbar2chi0chi0 : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2chi0chi0(int id3chiIn, int id4chiIn, int codeIn) 
   : id3chi(id3chiIn), id4chi(id4chiIn), codeSave(codeIn) {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "qqbarSame";}
  virtual int    id3Mass() const {return id3;}
  virtual int    id4Mass() const {return id4;}

private:

  // Values stored for later use.  
  int     id3chi, id4chi, codeSave;
  string  nameSave;
  double  sigma0, ui, uj, ti, tj, sz, d;
  complex propZ;

  // Couplings.
  // Shorthand for sin2thetaW, mZ, and GammaZ.
  double  sin2W, mZ, wZ;      
  // qqZ couplings.
  double  LqqZ[10], RqqZ[10]; 
  // qsqchi_i couplings.
  complex LsqXi[10][10], RsqXi[10][10];
  // qsqchi_j couplings.
  complex LsqXj[10][10], RsqXj[10][10];
  complex OL, OR;

};
  
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_SigmaSUSY_H
