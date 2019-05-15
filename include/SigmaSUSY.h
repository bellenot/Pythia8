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
  static complex LsddG[7][4], RsddG[7][4];
  static complex LsuuG[7][4], RsuuG[7][4];
  static complex getLsqqG(int iGenSq, int idQ) {
    // Assume generation index for Squark. Translate if PDG code instead.
    if (abs(iGenSq) > 1000000) 
      iGenSq =  3*(abs(iGenSq)/2000000) + (abs(iGenSq)%10+1)/2;
    return (abs(idQ)%2 == 0) ? LsuuG[iGenSq][abs(idQ)/2]
      : LsddG[iGenSq][(abs(idQ)+1)/2] ;
  }
  static complex getRsqqG(int iGenSq, int idQ) {
    // Assume generation index for Squark. Translate if PDG code instead.
    if (abs(iGenSq) > 1000000) 
      iGenSq =  3*(abs(iGenSq)/2000000) + (abs(iGenSq)%10+1)/2;
    return (abs(idQ)%2 == 0) ? RsuuG[iGenSq][abs(idQ)/2]
      : RsddG[iGenSq][(abs(idQ)+1)/2] ;
  }

  // ~chi0~chi0Z couplings
  static complex OLpp[6][6], ORpp[6][6];

  // ~chi+~chi-Z couplings
  static complex OLp[3][3], ORp[3][3];

  // ~chi0~chi+W couplings
  static complex OL[6][3], OR[6][3];

  // qqZ couplings 
  static double LqqZ[7], RqqZ[7]; 

  // ~q~qZ couplings 
  static complex LsdsdZ[7][7], RsdsdZ[7][7]; 
  static complex LsusuZ[7][7], RsusuZ[7][7]; 
  static complex getLsqsqZ(int idSq1, int idSq2) {    
    if (abs(idSq1)%2 != abs(idSq2)%2) return complex(0.0,0.0);
    int iGen1 = 3*(abs(idSq1)/2000000) + (abs(idSq1)%10+1)/2;
    int iGen2 = 3*(abs(idSq2)/2000000) + (abs(idSq2)%10+1)/2;
    return (abs(idSq1)%2 == 0) ? LsusuZ[iGen1][iGen2]
      : LsdsdZ[iGen1][iGen2];
  }
  static complex getRsqsqZ(int idSq1, int idSq2) {    
    if (abs(idSq1)%2 != abs(idSq2)%2) return complex(0.0,0.0);
    int iGen1 = 3*(abs(idSq1)/2000000) + (abs(idSq1)%10+1)/2;
    int iGen2 = 3*(abs(idSq2)/2000000) + (abs(idSq2)%10+1)/2;
    return (abs(idSq1)%2 == 0) ? RsusuZ[iGen1][iGen2]
      : RsdsdZ[iGen1][iGen2];
  }

  // udW couplings
  static complex LudW[4][4], RudW[4][4];

  // ~u~dW couplings
  static complex LsusdW[7][7], RsusdW[7][7];

  // ~qq~chi0 couplings
  static complex LsddX[7][4][6], RsddX[7][4][6];
  static complex LsuuX[7][4][6], RsuuX[7][4][6];
  static complex getLsqqX(int iSq, int idQ, int iNeut) {
    return (abs(idQ)%2 == 0) ? LsuuX[iSq][abs(idQ)/2][iNeut] 
      : LsddX[iSq][(abs(idQ)+1)/2][iNeut] ;
  }
  static complex getRsqqX(int iSq, int idQ, int iNeut) {
    return (abs(idQ)%2 == 0) ? RsuuX[iSq][abs(idQ)/2][iNeut] 
      : RsddX[iSq][(abs(idQ)+1)/2][iNeut] ;
  }
 
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

  // Function to return sup flavour codes
  static int idSup(int iSup) {
    int id = 0;
    int sgn = ( iSup > 0 ) ? 1 : -1;
    iSup = abs(iSup);
    if      (iSup ==  1) id =  1000002; 
    else if (iSup ==  2) id =  1000004;     
    else if (iSup ==  3) id =  1000006; 
    else if (iSup ==  4) id =  2000002; 
    else if (iSup ==  5) id =  2000004;     
    else if (iSup ==  6) id =  2000006; 
    return sgn*id;
  }

  // Function to return sdown flavour codes
  static int idSdown(int iSdown) {
    int id = 0;
    int sgn = ( iSdown > 0 ) ? 1 : -1;
    iSdown = abs(iSdown);
    if      (iSdown ==  1) id =  1000001; 
    else if (iSdown ==  2) id =  1000003;     
    else if (iSdown ==  3) id =  1000005; 
    else if (iSdown ==  4) id =  2000001; 
    else if (iSdown ==  5) id =  2000003;     
    else if (iSdown ==  6) id =  2000005; 
    return sgn*id;
  }

  // Function to return sdown flavour codes
  static int idSlep(int iSlep) {
    int id = 0;
    int sgn = ( iSlep > 0 ) ? 1 : -1;
    iSlep = abs(iSlep);
    if      (iSlep ==  1) id =  1000011; 
    else if (iSlep ==  2) id =  1000013;     
    else if (iSlep ==  3) id =  1000015; 
    else if (iSlep ==  4) id =  2000011; 
    else if (iSlep ==  5) id =  2000013;     
    else if (iSlep ==  6) id =  2000015; 
    return sgn*id;
  }

  // Function to return a particle name, given pdg code
  static string getName(int pdgCode) {    

    // Absolute value and corresponding SM code
    int codeA = abs(pdgCode);
    int idSM  = codeA % 1000000;

    // Name
    string name;

    // Flag to indicate whether SLHA1 or SLHA2
    bool slha1 = false;

    // SM particles
    if (codeA == idSM) {
      // Neutrinos
      if (!slhaPtr->upmns.exists()) slha1=true;
      if (codeA == 12) name = (slha1) ? "nu_e" : "nu_1";
      if (codeA == 14) name = (slha1) ? "nu_mu" : "nu_2";
      if (codeA == 16) name = (slha1) ? "nu_tau" : "nu_3";
    }

    // Squarks
    else if ( idSM <= 6) {
      // up squarks
      if (idSM % 2 == 0) {
	// If SLHA1, return old PDG names
	if (slhaPtr->stopmix.exists()) slha1 = true;
	if (codeA == 1000002) name = (slha1) ? "~u_L" : "~u_1";
	if (codeA == 1000004) name = (slha1) ? "~c_L" : "~u_2";
	if (codeA == 1000006) name = (slha1) ? "~t_1" : "~u_3";
	if (codeA == 2000002) name = (slha1) ? "~u_R" : "~u_4";
	if (codeA == 2000004) name = (slha1) ? "~c_R" : "~u_5";
	if (codeA == 2000006) name = (slha1) ? "~t_2" : "~u_6";
      } 
      // down squarks
      else {
	// If SLHA1, return old PDG names
	if (slhaPtr->sbotmix.exists()) slha1 = true;
	if (codeA == 1000001) name = (slha1) ? "~d_L" : "~d_1";
	if (codeA == 1000003) name = (slha1) ? "~s_L" : "~d_2";
	if (codeA == 1000005) name = (slha1) ? "~b_1" : "~d_3";
	if (codeA == 2000001) name = (slha1) ? "~d_R" : "~d_4";
	if (codeA == 2000003) name = (slha1) ? "~s_R" : "~d_5";
	if (codeA == 2000005) name = (slha1) ? "~b_2" : "~d_6";
      }
      if (pdgCode < 0) name += "bar";
    } 

    // Sleptons
    else if ( idSM <= 19 ) {

      // Sneutrinos
      if (idSM % 2 == 0) {
	// If SLHA1, return old PDG names
	if (slhaPtr->staumix.exists()) slha1 = true;
	if (codeA == 1000012) name = (slha1) ? "~nu_eL" : "~nu_1";
	if (codeA == 1000014) name = (slha1) ? "~nu_muL" : "~nu_2";
	if (codeA == 1000016) name = (slha1) ? "~nu_tauL" : "~nu_3";
	if (codeA == 2000012) name = (slha1) ? "~nu_eR" : "~nu_4";
	if (codeA == 2000014) name = (slha1) ? "~nu_muR" : "~nu_5";
	if (codeA == 2000016) name = (slha1) ? "~nu_tauR" : "~nu_6";
	if (pdgCode < 0) name += "bar";
      }
      // charged sleptons
      else {
	// If SLHA1, return old PDG names
	if (slhaPtr->staumix.exists()) slha1 = true;
	if (codeA == 1000011) name = (slha1) ? "~e_L" : "~e_1";
	if (codeA == 1000013) name = (slha1) ? "~mu_L" : "~e_2";
	if (codeA == 1000015) name = (slha1) ? "~tau_1" : "~e_3";
	if (codeA == 2000011) name = (slha1) ? "~e_R" : "~e_4";
	if (codeA == 2000013) name = (slha1) ? "~mu_R" : "~e_5";
	if (codeA == 2000015) name = (slha1) ? "~tau_2" : "~e_6";
	if (pdgCode < 0) name += "-";
	else name += "+";
      }
    }

    else if (codeA == 1000021) name = "~g";
    else if (codeA == 1000022) name = "~chi_10";
    else if (codeA == 1000023) name = "~chi_20";
    else if (codeA == 1000024) name = (pdgCode > 0) ? "~chi_1+" : "~chi_1-";
    else if (codeA == 1000025) name = "~chi_30";
    else if (codeA == 1000035) name = "~chi_40";
    else if (codeA == 1000037) name = (pdgCode > 0) ? "~chi_2+" : "~chi_2-";

    return name;

  }

  // is NMSSM 
  static bool isNMSSM;

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

  // Basic process information
  int     id3chi, id4chi, codeSave;
  string  nameSave;

  // Values stored for later use
  double  sigma0, ui, uj, ti, tj, openFracPair;
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

  // Basic process information
  int     id3chi, id4sq, codeSave;
  string  nameSave;

  // Values stored for later use
  double  sigma0, ui, uj, ti, tj, openFracPair;

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
    if (isUp)             id3Sav = -id3Sav;
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

  // Basic process information
  int id3Sav, id4Sav;

};

//**************************************************************************

// A derived class for q q' -> ~q_i ~q_j 

class Sigma2qq2squarksquark : public Sigma2Process {

public:

  // Constructor.
  Sigma2qq2squarksquark() {}

  // Constructor.
  Sigma2qq2squarksquark(int id3In, int id4In, int codeIn) { 

    // Save ordering indices and process code
    id3Sav = id3In;
    id4Sav = id4In;
    codeSave = codeIn; 
    // Initial values (flipped for c.c.)
    id3    = id3Sav;
    id4    = id4Sav;

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
  virtual int    id3Mass() const {return abs(id3Sav);}
  virtual int    id4Mass() const {return abs(id4Sav);}

private:

  // Basic process information
  int     id3Sav, id4Sav, codeSave, iGen3, iGen4, nNeut;
  string  nameSave;
  bool    isUD;  

  // Storage of mass squares
  double m2Glu;
  vector<double> m2Neut, m2Char;

  // Flavor-independent prefactors.
  double sigmaChar, sigmaNeut, sigmaGlu;
  double sigmaCharNeut, sigmaCharGlu, sigmaNeutGlu;
  double openFracPair;

  // Point-by-point info
  double tGlu, uGlu;
  vector<double> tNeut, uNeut, tChar, uChar;
  double sumCt, sumCu, sumNt, sumNu, sumGt, sumGu, sumInterference;

};

//**************************************************************************

// A derived class for q qbar' -> ~q_i ~q*_j 

class Sigma2qqbar2squarkantisquark : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2squarkantisquark() {}

  // Constructor.
  Sigma2qqbar2squarkantisquark(int id3In, int id4In, int codeIn) { 

    // Save ordering indices and process code
    // (always store squark first, antisquark second)
    id3Sav = abs(id3In);
    id4Sav = -abs(id4In);
    codeSave = codeIn; 
    // Initial values 
    id3    = id3Sav;
    id4    = id4Sav;

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
  virtual int    id3Mass() const {return abs(id3Sav);}
  virtual int    id4Mass() const {return abs(id4Sav);}

private:

  // Basic process information
  int     id3Sav, id4Sav, codeSave, iGen3, iGen4, nNeut;
  string  nameSave;
  bool    isUD, isCC;

  // Storage of mass squares
  double m2Glu;
  vector<double> m2Neut;

  // Flavor-independent prefactors: EW, strong, and interference
  double xW;
  double openFracPair;
  double sigmaEW, sigmaGlu, sigmaEWG;

  // Point-by-point info
  double tGlu, uGlu;
  vector<double> tNeut, uNeut;
  complex propZW; 
  double sumColS, sumColT, sumColSCC, sumColTCC, sumInterference;

};

//**************************************************************************

// A derived class for g g -> ~q ~q*

class Sigma2gg2squarkantisquark : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2squarkantisquark() {}

  // Constructor.
  Sigma2gg2squarkantisquark(int id34In, int codeIn) { 

    // Save ordering indices and process code
    // (always store squark first, antisquark second)
    id3Sav = abs(id34In);
    id4Sav = -abs(id34In);
    codeSave = codeIn; 
    // Initial values 
    id3    = id3Sav;
    id4    = id4Sav;

  }

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).  
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return nameSave;}
  virtual int    code()    const {return codeSave;}
  virtual string inFlux()  const {return "gg";}
  virtual int    id3Mass() const {return abs(id3Sav);}
  virtual int    id4Mass() const {return abs(id4Sav);}

private:

  // Basic process information
  int     id3Sav, id4Sav, codeSave;
  string  nameSave;
  double sigma, m2Sq, openFracPair;

  // Color flow info
  double sumColT, sumColU, sumInterference;

};

//**************************************************************************

// A derived class for q g -> ~q ~g

class Sigma2qg2squarkgluino : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2squarkgluino() {}

  // Constructor.
  Sigma2qg2squarkgluino(int id3In, int codeIn) { 

    // Save ordering indices and process code
    id3Sav = abs(id3In);
    codeSave = codeIn; 
    // Initial values 
    id3    = id3Sav;
    id4    = 1000021;

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
  virtual int    id3Mass() const {return abs(id3Sav);}
  virtual int    id4Mass() const {return 1000021;}

private:

  // Basic process information
  int     id3Sav, codeSave;
  string  nameSave;
  double sigmaA, sigmaB, comFac, m2Glu, m2Sq, openFracPair;

};

//**************************************************************************

// A derived class for g g -> gluino gluino.

class Sigma2gg2gluinogluino : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2gluinogluino() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).  
  virtual double sigmaHat() {return sigma;}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return "g g -> gluino gluino";}
  virtual int    code()    const {return 1201;}
  virtual string inFlux()  const {return "gg";}
  virtual int    id3Mass() const {return 1000021;}
  virtual int    id4Mass() const {return 1000021;}

private:

  // Values stored for process type and colour flow selection.
  double sigTS, sigUS, sigTU, sigSum, sigma, openFracPair;

};

//**************************************************************************

// A derived class for q qbar -> gluino gluino.

class Sigma2qqbar2gluinogluino : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2gluinogluino() {}

  // Initialize process. 
  virtual void initProc(); 

  // Calculate flavour-independent parts of cross section.
  virtual void sigmaKin();

  // Evaluate d(sigmaHat)/d(tHat).  
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name()    const {return "q qbar -> gluino gluino";}
  virtual int    code()    const {return 1202;}
  virtual string inFlux()  const {return "qqbarSame";}
  virtual int    id3Mass() const {return 1000021;}
  virtual int    id4Mass() const {return 1000021;}

private:

  // Values stored for process type and colour flow selection.
  double openFracPair, s34Avg, sigS, tHG, uHG, tHG2, uHG2;

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_SigmaSUSY_H

