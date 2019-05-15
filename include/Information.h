// This file contains classes that keep track of generic event info.
// Info: contains information on the generation process.
// ErrorMessages: table with all warnings and errors encountered.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_Information_H
#define Pythia8_Information_H

#include "Stdlib.h"
#include "Settings.h"

namespace Pythia8 {
 
//**************************************************************************

// The Info class contains a mixed bag of information on the event
// generation activity, especially on the current subprocess properties. 
// This is used by the generation machinery, but can also be read
// by the user.

class Info {

public:

  // Constructor. 
  Info() {} 

  // Listing of most available information.
  void list(ostream& os = cout);
  
  // Beam particles (in rest frame). CM energy of event.
  int idA() const {return idAM;}
  int idB() const {return idBM;}
  double pzA() const {return pzAM;}
  double pzB() const {return pzBM;}
  double eA() const {return eAM;}
  double eB() const {return eBM;}
  double mA() const {return mAM;}
  double mB() const {return mBM;}
  double eCM() const {return eCMM;}
  double s() const {return sM;}

  // Process name and code, and the number of final-state particles.
  string name() const {return nameSave;}
  int code() const {return codeSave;}    
  int nFinal() const {return nFinalSave;}

  // Are beam particles resolved, with pdf's? Are they diffractive?
  bool isResolved() const {return isRes;}
  bool isDiffractiveA() const {return isDiffA;} 
  bool isDiffractiveB() const {return isDiffB;} 

  // Incoming parton flavours, x values and densities at given Q2 scale.
  int id1() const {return id1H;}
  int id2() const {return id2H;}
  double x1() const {return x1H;}
  double x2() const {return x2H;}
  double pdf1() const {return pdf1H;}
  double pdf2() const {return pdf2H;}
  double Q2pdf() const {return Q2pdfH;}

  // Mandelstam variables (notation as if subcollision).
  double mHat() const {return sqrt(sH);}   
  double sHat() const {return sH;}   
  double tHat() const {return tH;}   
  double uHat() const {return uH;}   
  double pTHat() const {return pTH;} 
  double pT2Hat() const {return pTH*pTH;} 
  double m3Hat() const {return m3H;}   
  double m4Hat() const {return m4H;} 
  double thetaHat() const {return thetaH;}   

private:

  // Set info on the two incoming beams: only from Pythia class!
  friend class Pythia;
  void setBeamA( int idAin, double pzAin, double eAin, double mAin) {
    idAM = idAin; pzAM = pzAin; eAM = eAin; mAM = mAin;}
  void setBeamB( int idBin, double pzBin, double eBin, double mBin) {
    idBM = idBin; pzBM = pzBin; eBM = eBin; mBM = mBin;}
  void setECM( double eCMin) {eCMM = eCMin; sM = eCMM * eCMM;}

  // Set info on the (sub)process: only from ProcessContainer class.
  friend class ProcessContainer;
  void setType( string nameIn, int codeIn, int nFinalIn, 
    bool isResolvedIn = true, bool isDiffractiveAin = false, 
    bool isDiffractiveBin = false) {nameSave = nameIn; 
    codeSave = codeIn; nFinalSave = nFinalIn; isRes = isResolvedIn; 
    isDiffA = isDiffractiveAin; isDiffB = isDiffractiveBin;
    nTotal = 2 + nFinalSave;}
  void setPDF( int id1In, int id2In,  double pdf1In, double pdf2In, 
    double Q2pdfIn) {id1H = id1In; id2H = id2In; pdf1H = pdf1In; 
    pdf2H = pdf2In; Q2pdfH = Q2pdfIn;}
  void setKin( double x1In, double x2In, double sHatIn, double tHatIn, 
    double uHatIn, double pTHatIn, double m3HatIn, double m4HatIn, 
    double thetaHatIn, double phiHatIn) {x1H = x1In; x2H = x2In; 
    sH = sHatIn; tH = tHatIn; uH = uHatIn; pTH = pTHatIn; m3H = m3HatIn; 
    m4H = m4HatIn; thetaH = thetaHatIn; phiH = phiHatIn;}

  // Store global quantities. 
  int idAM, idBM;
  double pzAM, eAM,mAM, pzBM, eBM, mBM, eCMM, sM;

  // Store current-event quantities.
  string nameSave;
  int codeSave, nFinalSave, nTotal, id1H, id2H;
  bool isRes, isDiffA, isDiffB;  
  double x1H, x2H, pdf1H, pdf2H, Q2pdfH, sH, tH, uH, pTH, m3H, m4H, 
    thetaH, phiH;

};

//**************************************************************************

// This class holds info on all messages received and how many times.

class ErrorMessages {

public:

  // Constructor.
  ErrorMessages() {}

  // Normal init not needed for now. reInit resets to empty map.
  static void init() {}
  static void reInit() {messages.clear();}

  // Initialize static data members.
  static void initStatic();
  
  // Print a message the first few times. Insert in database.
  static void message(string messageIn, string extraIn = " ");

  // Print statistics on errors/warnings.
  static void statistics(ostream& os = cout);

private:

  // Static initialization data, normally only set once.
  static int timesToPrint;

  // Map for all error messages.
  static map<string, int> messages;
 
};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_Information_H
