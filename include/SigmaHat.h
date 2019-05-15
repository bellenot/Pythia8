// Header file for hard-process differential cross sections.
// InFlux: base class for combinations of incoming partons.
// InFluxgg, InFluxqqAnti, InFluxqg: derived classes.
// SigmaHat: base class for cross sections.
// SigmaHgg2gg, ... : derived classes for specific cross sections.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_SigmaHat_H
#define Pythia8_SigmaHat_H

#include "Stdlib.h"
#include "Basics.h"
#include "Settings.h"
#include "ParticleData.h"
#include "StandardModel.h"
#include "PartonDistributions.h"
#include "Information.h"
#include "SigmaTotal.h"

namespace Pythia8 {

//**************************************************************************

// InFlux is the base class for the combined incoming parton flux.

class InFlux {

public:

  // Destructor.
  virtual ~InFlux() {}

  // Initialize static data members.
  static void initStatic();

  // Initialization: store parton-density pointers. Construct channels.
  void init(PDF* pdfAPtrIn, PDF* pdfBPtrIn) {pdfAPtr = pdfAPtrIn; 
    pdfBPtr = pdfBPtrIn; initChannels();} 

  // Initialization of process-specific allowed combinations. 
  virtual void initChannels() = 0; 

  // Add even powers of e to weight.
  void weightCharge(int ePow = 0);

  // Calculate products of parton densities for allowed combinations.
  double flux(double x1, double x2, double Q2);

  // Information on combination of required partons from the two beams.
  int nAB() const {return weightAB.size();}
  int idA(int i) const {return idPartonPairA[i];} 
  int idB(int i) const {return idPartonPairB[i];} 
  double wtAB(int i) const {return weightAB[i];} 
  double fluxwtAB(int i) const {return fluxweightAB[i];} 

  // Pick one of the possible channels according to their weight.
  void pick();
  int id1() const {return idNow1;}
  int id2() const {return idNow2;}
  double pdf1() const {return pdfNow1;}
  double pdf2() const {return pdfNow2;}

protected:

  // Static initialization data, normally only set once.
  static int nQuark;

  // Pointers to parton densities. 
  PDF* pdfAPtr; 
  PDF* pdfBPtr;  

  // Partons in beams and their allowed combinations with weights.
  vector<int> idPartonA, idPartonB, idPartonPairA, idPartonPairB;
  vector<double> pdfA, pdfB, pdfPairA, pdfPairB, weightAB, fluxweightAB;
  int idNow1, idNow2;
  double pdfNow1, pdfNow2, fluxwtSum;

};
 
//**************************************************************************

// A derived class for dummy cases, e.g. elastic/diffractive scattering.

class InFluxNone : public InFlux {

public:

  // Constructor.
  InFluxNone() {}

private:

  // Initialize values.
  virtual void initChannels();  

};
 
//**************************************************************************

// A derived class for g g incoming state.

class InFluxgg : public InFlux {

public:

  // Constructor.
  InFluxgg() {}

private:

  // Initialize values.
  virtual void initChannels();  

};
 
//**************************************************************************

// A derived class for q g incoming state.

class InFluxqg : public InFlux {

public:

  // Constructor.
  InFluxqg() {}

private:

  // Initialize values.
  virtual void initChannels();  

};
 
//**************************************************************************

// A derived class for q qbbar' or q q' (qbar qbar') incoming states, 
// with q' != q.

class InFluxqqbarqqDiff : public InFlux {

public:

  // Constructor.
  InFluxqqbarqqDiff() {}

private:

  // Initialize values.
  virtual void initChannels();  

};
 
//**************************************************************************

// A derived class for q q' (qbar qbar') incoming states, with q' != q.

class InFluxqqDiff : public InFlux {

public:

  // Constructor.
  InFluxqqDiff() {}

private:

  // Initialize values.
  virtual void initChannels();  

};
 
//**************************************************************************

// A derived class for q q (qbar qbar) incoming states.

class InFluxqqSame : public InFlux {

public:

  // Constructor.
  InFluxqqSame() {}

private:

  // Initialize values.
  virtual void initChannels();  

};
 
//**************************************************************************

// A derived class for q qbar'  incoming state.

class InFluxqqbarDiff : public InFlux {

public:

  // Constructor.
  InFluxqqbarDiff() {}

private:

  // Initialize values.
  virtual void initChannels();  

};
 
//**************************************************************************

// A derived class for q qbar antiparticle incoming state.

class InFluxqqbarSame : public InFlux {

public:

  // Constructor.
  InFluxqqbarSame() {}

private:

  // Initialize values.
  virtual void initChannels();  

};
 
//**************************************************************************

// SigmaHat is the base class for cross section calculations.

class SigmaHat {

public:

  // Destructor.
  virtual ~SigmaHat() {if (inFluxPtr !=0) delete inFluxPtr;}

  // Initialize static data members.
  static void initStatic();

  // Get pointer to SigmaTotal.
  void setSigmaTotalPtr(SigmaTotal* sigmaTotPtrIn) {
    sigmaTotPtr = sigmaTotPtrIn;}

  // Store info of use during the generation. 
  void initInfo(Info& info);

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr) = 0; 

  // Evaluate integrated sigma for elastic/diffractive processes. 
  virtual double sigma0() {return 0.;} 

  // Evaluate sigmaHat(sHat) for resolved 2 -> 1 processes. 
  virtual double sigma1( double x1, double x2, double sH) 
    {return x1 + x2 + sH;} 

  // Evaluate d(sigmaHat)/d(tHat) for resolved 2 -> 1 processes. 
  virtual double sigma2( double x1, double x2, double sH, double tH) 
    {return x1 + x2 + sH + tH;} 

  // Select incoming parton channel and extract parton densities.
  void pickInState() {inFluxPtr->pick();} 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol() = 0;

  // Process name and code, and the number of final-state particles.
  virtual string name() const {return "unnamed process";}
  virtual int code() const {return 0;}
  virtual int nFinal() const {return 2;};

  // Special treatment needed for elastic and diffractive processes
  virtual bool isResolved() const {return true;};
  virtual bool isDiffA() const {return false;};
  virtual bool isDiffB() const {return false;};

  // Give back flavours, colours and PDF values.
  int id(int i) const {return idH[i];}
  int col(int i) const {return colH[i];} 
  int acol(int i) const {return acolH[i];} 
  double pdf1() const {return inFluxPtr->pdf1();}
  double pdf2() const {return inFluxPtr->pdf2();}
  double Q2pdf() const {return q2pdf;}

protected:

  // Constructor.
  SigmaHat() {}

  // Static initialization data, normally only set once.
  static int alphaSorder, nQuark;
  static double alphaSvalue, alphaEM;

  // Static alphaStrong calculation.
  static AlphaStrong alphaScalc;

  // Constants: could only be changed in the code itself.
  static const double CONVERT2MB;

  // Information on incoming beams.
  int idA, idB;
  double mA, mB; 

  // Current data; Q2 scale for PDF's.
  double q2pdf;
  
  // Pointer to the parton density flux object.
  InFlux* inFluxPtr;
  
  // Pointer to the total/elastic/diffractive cross section  object.
  SigmaTotal* sigmaTotPtr;

  // Set flavour, colour and anticolour.
  void setId( int id1 = 0, int id2 = 0, int id3 = 0, int id4 = 0,
    int id5 = 0) {idH[1] = id1; idH[2] = id2; idH[3] = id3; 
    idH[4] = id4; idH[5] = id5;}
  void setColAcol( int col1 = 0, int acol1 = 0, 
    int col2 = 0, int acol2 = 0, int col3 = 0, int acol3 = 0, 
    int col4 = 0, int acol4 = 0, int col5 = 0, int acol5 = 0) {
    colH[1] = col1; acolH[1] = acol1; colH[2] = col2; acolH[2] = acol2;
    colH[3] = col3; acolH[3] = acol3; colH[4] = col4; acolH[4] = acol4;
    colH[5] = col5; acolH[5] = acol5; }
  void swapColAcol() { swap(colH[1], acolH[1]); 
    swap(colH[2], acolH[2]); swap(colH[3], acolH[3]); 
    swap(colH[4], acolH[4]); swap(colH[5], acolH[5]);}
  void swap1234() { swap(colH[1], colH[2]); swap(colH[3], colH[4]); 
    swap(acolH[1], acolH[2]); swap(acolH[3], acolH[4]);}

  // Store flavour, colour and anticolour.
  int idH[6], colH[6], acolH[6];

};
 
//**************************************************************************

// A derived class for elastic scattering A B -> A B.

class SigmaHAB2AB : public SigmaHat {

public:

  // Constructor.
  SigmaHAB2AB() {}

  // Initialize dummy parton flux object.
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

  // Evaluate sigma. 
  virtual double sigma0() { return sigmaTotPtr->sigmaEl();} 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "A B -> A B elastic";}
  virtual int code() const {return 102;}
  virtual int nFinal() const {return 2;}
  virtual bool isResolved() const {return false;}

private:

};
 
//**************************************************************************

// A derived class for single diffractive scattering A B -> X B.

class SigmaHAB2XB : public SigmaHat {

public:

  // Constructor.
  SigmaHAB2XB() {}

  // Initialize dummy parton flux object.
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

  // Evaluate sigma. 
  virtual double sigma0() { return sigmaTotPtr->sigmaXB();} 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "A B -> X B single diffractive";}
  virtual int code() const {return 103;}
  virtual int nFinal() const {return 2;}
  virtual bool isResolved() const {return false;}
  virtual bool isDiffA() const {return true;};

private:

};
 
//**************************************************************************

// A derived class for single diffractive scattering A B -> A X.

class SigmaHAB2AX : public SigmaHat {

public:

  // Constructor.
  SigmaHAB2AX() {}

  // Initialize dummy parton flux object.
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

  // Evaluate sigma. 
  virtual double sigma0() { return sigmaTotPtr->sigmaAX();} 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "A B -> A X single diffractive";}
  virtual int code() const {return 104;}
  virtual int nFinal() const {return 2;}
  virtual bool isResolved() const {return false;}
  virtual bool isDiffB() const {return true;};

private:

};
 
//**************************************************************************

// A derived class for double diffractive scattering A B -> X X.

class SigmaHAB2XX : public SigmaHat {

public:

  // Constructor.
  SigmaHAB2XX() {}

  // Initialize dummy parton flux object.
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

  // Evaluate sigma. 
  virtual double sigma0() { return sigmaTotPtr->sigmaXX();} 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "A B -> X X double diffractive";}
  virtual int code() const {return 105;}
  virtual int nFinal() const {return 2;}
  virtual bool isResolved() const {return false;}
  virtual bool isDiffA() const {return true;};
  virtual bool isDiffB() const {return true;};

private:

};
 
//**************************************************************************

// A derived class for g g -> g g.

class SigmaHgg2gg : public SigmaHat {

public:

  // Constructor.
  SigmaHgg2gg() {}

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigma2( double x1, double x2, double sH, double tH); 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "g g -> g g";}
  virtual int code() const {return 111;}
  virtual int nFinal() const {return 2;}

private:

  // Values stored for colour flow selection.
  double sigTS, sigUS, sigTU, sigSum;

};
 
//**************************************************************************

// A derived class for g g -> q qbar (q = u, d, s, i.e. almost massless).

class SigmaHgg2qqbar : public SigmaHat {

public:

  // Constructor.
  SigmaHgg2qqbar() {}

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigma2( double x1, double x2, double sH, double tH); 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "g g -> q qbar (uds)";}
  virtual int code() const {return 112;}
  virtual int nFinal() const {return 2;}

private:

  // Values stored for colour flow selection.
  int idNew;
  double mNew, m2New, sigTS, sigUS, sigSum;

};
 
//**************************************************************************

// A derived class for q g -> q g (q = u, d, s, c, b).
// Use massless approximation also for Q since no alternative.

class SigmaHqg2qg : public SigmaHat {

public:

  // Constructor.
  SigmaHqg2qg() {}

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigma2( double x1, double x2, double sH, double tH); 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "q g -> q g (udscb)";}
  virtual int code() const {return 113;}
  virtual int nFinal() const {return 2;}

private:

  // Values stored for colour flow selection.
  double mNew, m2New, sigTS, sigTU, sigSum;

};
 
//**************************************************************************

// A derived class for q qbar' -> q qbar' or q q' -> q q' 
// (qbar qbar' -> qbar qbar'), q' != q.

class SigmaHqq2qqDiff : public SigmaHat {

public:

  // Constructor.
  SigmaHqq2qqDiff() {}

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigma2( double x1, double x2, double sH, double tH); 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "q q(bar)' -> q q(bar)'";}
  virtual int code() const {return 114;}
  virtual int nFinal() const {return 2;}

 private:

  // Values stored for colour flow selection.
  double sigT;

};
 
//**************************************************************************

// A derived class for q q -> q q (qbar qbar -> qbar qbar).

class SigmaHqq2qqSame : public SigmaHat {

public:

  // Constructor.
  SigmaHqq2qqSame() {}

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigma2( double x1, double x2, double sH, double tH); 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "q q -> q q";}
  virtual int code() const {return 115;}
  virtual int nFinal() const {return 2;}

 private:

  // Values stored for colour flow selection.
  double sigT, sigU, sigTU, sigSum;

};
 
//**************************************************************************

// A derived class for q qbar -> q qbar.

class SigmaHqqbar2qqbarSame : public SigmaHat {

public:

  // Constructor.
  SigmaHqqbar2qqbarSame() {}

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigma2( double x1, double x2, double sH, double tH); 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "q qbar -> q qbar";}
  virtual int code() const {return 116;}
  virtual int nFinal() const {return 2;}

 private:

  // Values stored for colour flow selection.
  double sigT;

};
 
//**************************************************************************

// A derived class for q qbar -> q' qbar'.

class SigmaHqqbar2qqbarNew : public SigmaHat {

public:

  // Constructor.
  SigmaHqqbar2qqbarNew() {}

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigma2( double x1, double x2, double sH, double tH); 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "q qbar -> q' qbar' (uds)";}
  virtual int code() const {return 117;}
  virtual int nFinal() const {return 2;}

 private:

  // Values stored for colour flow selection.
  int idNew;
  double mNew, m2New, sigS;

};
 
//**************************************************************************

// A derived class for q qbar -> g g.

class SigmaHqqbar2gg : public SigmaHat {

public:

  // Constructor.
  SigmaHqqbar2gg() {}

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigma2( double x1, double x2, double sH, double tH); 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "q qbar -> g g";}
  virtual int code() const {return 118;}
  virtual int nFinal() const {return 2;}

 private:

  // Values stored for colour flow selection.
  double sigTS, sigUS, sigSum;

};
 
//**************************************************************************

// A derived class for q qbar -> Q Qbar (Q = c or b).

class SigmaHqqbar2QQbar : public SigmaHat {

public:

  // Constructor.
  SigmaHqqbar2QQbar(int idIn, int codeIn, string nameIn) 
    : idNew(idIn), codeSave(codeIn), nameSave(nameIn) {}

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigma2( double x1, double x2, double sH, double tH); 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return nameSave;}
  virtual int code() const {return codeSave;}
  virtual int nFinal() const {return 2;}

 private:

  // Values stored for colour flow selection.
  int idNew, codeSave;
  string nameSave;
  double mNew, m2New, sigS;

};
 
//**************************************************************************

// A derived class for q g -> q gamma (q = u, d, s, c, b).
// Use massless approximation also for Q since no alternative.

class SigmaHqg2qgamma : public SigmaHat {

public:

  // Constructor.
  SigmaHqg2qgamma() {}

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigma2( double x1, double x2, double sH, double tH); 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "q g -> q gamma (udscb)";}
  virtual int code() const {return 131;}
  virtual int nFinal() const {return 2;}

private:

  // Values stored for colour flow selection.
  double mNew, m2New, sigUS;

};
 
//**************************************************************************

// A derived class for q qbar -> g gamma.

class SigmaHqqbar2ggamma : public SigmaHat {

public:

  // Constructor.
  SigmaHqqbar2ggamma() {}

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigma2( double x1, double x2, double sH, double tH); 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "q qbar -> g gamma";}
  virtual int code() const {return 132;}
  virtual int nFinal() const {return 2;}

private:

};
 
//**************************************************************************

// A derived class for g g -> g gamma.

class SigmaHgg2ggamma : public SigmaHat {

public:

  // Constructor.
  SigmaHgg2ggamma() {}

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigma2( double x1, double x2, double sH, double tH); 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "g g -> g gamma";}
  virtual int code() const {return 133;}
  virtual int nFinal() const {return 2;}

private:
  
  double chargeSum;

};
 
//**************************************************************************

// A derived class for q qbar -> gamma gamma.

class SigmaHqqbar2gammagamma : public SigmaHat {

public:

  // Constructor.
  SigmaHqqbar2gammagamma() {}

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigma2( double x1, double x2, double sH, double tH); 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "q qbar -> gamma gamma";}
  virtual int code() const {return 134;}
  virtual int nFinal() const {return 2;}

private:

  // Values stored for colour flow selection.
  double sigTU;

};
 
//**************************************************************************

// A derived class for g g -> gamma gamma.

class SigmaHgg2gammagamma : public SigmaHat {

public:

  // Constructor.
  SigmaHgg2gammagamma() {}

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigma2( double x1, double x2, double sH, double tH); 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "g g -> gamma gamma";}
  virtual int code() const {return 135;}
  virtual int nFinal() const {return 2;}

private:
  
  double charge2Sum;

};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_SigmaHat_H
 
