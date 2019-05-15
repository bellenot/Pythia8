// SigmaSUSY.h is a part of the PYTHIA event generator.
// Copyright (C) 2009 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for Supersymmetric process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma2Process.

#ifndef Pythia8_SigmaSUSY_H
#define Pythia8_SigmaSUSY_H

#include "PhaseSpace.h"
#include "PythiaComplex.h"
#include "SigmaProcess.h"

namespace Pythia8 {
 
//**************************************************************************

// CoupSUSY
// Auxiliary class to compute and store various SM and SUSY couplings.

class CoupSUSY {

public:

  // Constructor
  CoupSUSY() {isInit=false;};

  // Initialize
  static void initStatic(SusyLesHouches*);

  // Status flag
  static bool isInit;

  // Z and W pole masses and widths
  static double mWpole, wWpole, mZpole, wZpole;

  // Running masses and weak mixing angle 
  // (default to pole values if no running available)
  static double mW, mZ, sin2W, sinW, cosW;

  // Tanbeta
  static double tanb, cosb, sinb;

  // ~qq~g couplings
  static complex Lsddg[7][4], Rsddg[7][4];
  static complex Lsuug[7][4], Rsuug[7][4];

  // ~chi0~chi0Z couplings
  static complex OLpp[6][6], ORpp[6][6];

  // ~chi+~chi-Z couplings
  static complex OLp[3][3], ORp[3][3];

  // ~chi0~chi+W couplings
  static complex OL[6][3], OR[6][3];

  // qqZ couplings 
  static complex LqqZ[7], RqqZ[7]; 

  // ~q~qZ couplings 
  static complex LsdsdZ[7][7], RsdsdZ[7][7]; 
  static complex LsusuZ[7][7], RsusuZ[7][7]; 

  // udW couplings
  static complex LudW[4][4], RudW[4][4];

  // ~u~dW couplings
  static complex LsusdW[7][7], RsusdW[7][7];

  // ~qq~chi0 couplings
  static complex LsddX[7][4][6], RsddX[7][4][6];
  static complex LsuuX[7][4][6], RsuuX[7][4][6];

  // ~du~chi+ couplings
  static complex LsduX[7][4][3], RsduX[7][4][3];

  // ~ud~chi+ couplings
  static complex LsudX[7][4][3], RsudX[7][4][3];

  // Function to return neutralino and chargino flavour codes
  static int idNeut(int idChi) {
    int id = 0;
    if      (idChi == 1) id = 1000022; 
    else if (idChi == 2) id = 1000023; 
    else if (idChi == 3) id = 1000025; 
    else if (idChi == 4) id = 1000035; 
    else if (idChi == 5) id = 1000045; 
    return id;
  }

  // Function to return neutralino and chargino flavour codes
  static int idChar(int idChi) {
    int id = 0;
    if      (idChi ==  1) id =  1000024; 
    else if (idChi == -1) id = -1000024;     
    else if (idChi ==  2) id =  1000037; 
    else if (idChi == -2) id = -1000037; 
    return id;
  }

private:

  // Pointer to SLHA instance
  static SusyLesHouches* slhaPtr;

  // Debug flag
  static const bool   DEBUG;

};

//**************************************************************************

// A derived class for q qbar -> neutralino_i neutralino_j.

class Sigma2qqbar2chi0chi0 : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2chi0chi0() { };

  // Constructor.
  Sigma2qqbar2chi0chi0(int id3chiIn, int id4chiIn, int codeIn) { 

    // Save ordering indices and process code
    id3chi   = id3chiIn; 
    id4chi   = id4chiIn; 
    codeSave = codeIn; 

    // Construct id codes from ordering indices.
    id3                  = 1000022; 
    if (id3chi == 2) id3 = 1000023; 
    if (id3chi == 3) id3 = 1000025; 
    if (id3chi == 4) id3 = 1000035; 
    if (id3chi == 5) id3 = 1000045; 
    id4                  = 1000022; 
    if (id4chi == 2) id4 = 1000023; 
    if (id4chi == 3) id4 = 1000025; 
    if (id4chi == 4) id4 = 1000035; 
    if (id4chi == 5) id4 = 1000045; 
    
  }

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
  virtual string inFlux()  const {return "qq";}
  virtual int    id3Mass() const {return abs(id3);}
  virtual int    id4Mass() const {return abs(id4);}
  virtual int    resonanceA() const {return 23;}

 protected:  

  // Values stored for later use.  
  int     id3chi, id4chi, codeSave;
  string  nameSave;
  double  sigma0, ui, uj, ti, tj;
  complex propZ;

};

//**************************************************************************

// A derived class for q qbar -> neutralino_i chargino_j.

class Sigma2qqbar2charchi0 : public Sigma2qqbar2chi0chi0 {

public:

  // Constructor.
  Sigma2qqbar2charchi0(int id3chiIn, int id4chiIn, int codeIn) {
    
    // Save ordering indices and process code
    id3chi   = id3chiIn; 
    id4chi   = id4chiIn; 
    codeSave = codeIn; 

    // Construct id codes from ordering indices.
    id3 = (abs(id3chi) == 2) ? 1000037 : 1000024; 
    if (id3chi < 0)  id3 = -id3;

    id4                  = 1000022; 
    if (id4chi == 2) id4 = 1000023; 
    if (id4chi == 3) id4 = 1000025; 
    if (id4chi == 4) id4 = 1000035; 
    if (id4chi == 5) id4 = 1000045; 
    
  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();
  
  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  virtual int    resonanceA() const {return 24;}

protected :

  complex propW;

};

//**************************************************************************

// A derived class for q qbar -> chargino+_i chargino-_j.

class Sigma2qqbar2charchar : public Sigma2qqbar2chi0chi0 {

public:

  // Constructor.
  Sigma2qqbar2charchar(int id3chiIn, int id4chiIn, int codeIn) {

    // Save ordering indices and process code
    id3chi   = id3chiIn; 
    id4chi   = id4chiIn; 
    codeSave = codeIn; 

    // Construct id codes from ordering indices.
    id3 = (abs(id3chi) == 2) ?  1000037 :  1000024; 
    id4 = (abs(id4chi) == 2) ? -1000037 : -1000024; 

  }

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();
  
  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

};

//**************************************************************************

// A derived class for q g -> neutralino_i squark_j (and cc)

class Sigma2qg2chi0squark : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2chi0squark() { };

  // Constructor.
  Sigma2qg2chi0squark(int id3chiIn, int id4sqIn, bool isUp, int codeIn) { 

    // Save ordering indices and process code
    id3chi   = id3chiIn; 
    id4sq    = id4sqIn; 
    codeSave = codeIn; 

    // Construct id codes from ordering indices.
    id3                  = 1000022; 
    if (id3chi == 2) id3 = 1000023; 
    if (id3chi == 3) id3 = 1000025; 
    if (id3chi == 4) id3 = 1000035; 
    if (id3chi == 5) id3 = 1000045; 
    id4                  = 1000001 + (isUp ? 1 : 0); 
    if (id4sq  == 2) id4 = 1000003 + (isUp ? 1 : 0); 
    if (id4sq  == 3) id4 = 1000005 + (isUp ? 1 : 0);
    if (id4sq  == 4) id4 = 2000001 + (isUp ? 1 : 0); 
    if (id4sq  == 5) id4 = 2000003 + (isUp ? 1 : 0); 
    if (id4sq  == 6) id4 = 2000005 + (isUp ? 1 : 0); 
    
  }

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
  virtual string inFlux()  const {return "qg";}
  virtual int    id3Mass() const {return abs(id3);}
  virtual int    id4Mass() const {return abs(id4);}

 protected:  

  // Values stored for later use.  
  int     id3chi, id4sq, codeSave;
  string  nameSave;
  double  sigma0, ui, uj, ti, tj;

};

//**************************************************************************

// A derived class for q g -> chargino_i squark_j (incl cc)

class Sigma2qg2charsquark : public Sigma2qg2chi0squark {

public:

  // Constructor.
  Sigma2qg2charsquark(int id3chiIn, int id4sqIn, bool isUp, int codeIn) {

    // Save ordering indices and process code
    id3chi   = id3chiIn; 
    id4sq    = id4sqIn; 
    codeSave = codeIn; 

    // Construct id codes from ordering indices.
    id3Sav                       = 1000024; 
    if (abs(id3chi) == 2) id3Sav = 1000037; 
    if (id3chi < 0) id3Sav       = -id3Sav;
    id4Sav                       = 1000001 + (isUp ? 1 : 0); 
    if (id4sq  == 2) id4Sav      = 1000003 + (isUp ? 1 : 0); 
    if (id4sq  == 3) id4Sav      = 1000005 + (isUp ? 1 : 0);
    if (id4sq  == 4) id4Sav      = 2000001 + (isUp ? 1 : 0); 
    if (id4sq  == 5) id4Sav      = 2000003 + (isUp ? 1 : 0); 
    if (id4sq  == 6) id4Sav      = 2000005 + (isUp ? 1 : 0); 

    // Initial values, can be swapped to charge conjugates event by event.
    id3 = id3Sav;
    id4 = id4Sav;

  }

  // Initialize process. 
  virtual void initProc(); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

 private:
  int id3Sav, id4Sav;

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_SigmaSUSY_H

