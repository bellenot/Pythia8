// Header file for hard-process differential cross sections.
// InFlux: base class for combinations of incoming partons.
// InFluxgg, InFluxqqAnti, InFluxqg, ...: derived classes.
// SigmaProcess: base class for cross sections.
// Sigma0Process: base class for unresolved processes, 
// derived from SigmaProcess.
// Sigma2Process: base class for 2 -> 2 processes, 
// derived from SigmaProcess.
// Sigma0MinBias, Sigma2gg2gg, ... : derived classes 
// for specific cross sections.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_SigmaProcess_H
#define Pythia8_SigmaProcess_H

#include "Basics.h"
#include "Event.h"
#include "Information.h"
#include "ParticleData.h"
#include "PartonDistributions.h"
#include "Settings.h"
#include "SigmaTotal.h"
#include "StandardModel.h"
#include "Stdlib.h"

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

// A derived class for g g incoming state.

class InFluxgg : public InFlux {

public:

  // Constructor.
  InFluxgg() {}

  // Destructor.
  ~InFluxgg() {}

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

// SigmaProcess is the base class for cross section calculations.

class SigmaProcess {

public:

  // Destructor.
  virtual ~SigmaProcess() {}

  // Initialize static data members.
  static void initStatic();

  // Store pointer to SigmaTotal and info on beams (unresolved only).
  void setSigmaTotalPtr(SigmaTotal* sigmaTotPtrIn, int idAIn, int idBIn,
    double mAIn, double mBIn) { sigmaTotPtr = sigmaTotPtrIn; idA = idAIn;
    idB = idBIn; mA = mAIn; mB = mBIn;}

  // Evaluate sigma for unresolved, sigmaHat(sHat) for 2 -> 1 processes, 
  // and d(sigmaHat)/d(tHat) for (resolved) 2 -> 2 processes. 
  virtual double sigmaHat() {return 0.;}

  // Convolute above with parton flux. 
  virtual double sigmaPDF() {return 0.;}

  // Initialize parton flux object for resolved processes (dummy here). 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {
    pdfAPtrDummy = pdfAPtr; pdfBPtrDummy = pdfBPtr;}

  // Input and complement kinematics for resolved 2 -> 2 process.
  virtual bool set2Kin( double x1in, double x2in, double sHin, double tHin)
    {return (x1in + x2in + sHin + tHin > 0.) ? true : false;} 

  // Ditto, but for Multiple Interactions applications, so different input.
  virtual bool set2KinMI( int id1in, int id2in, double x1in, double x2in,
    double sHin, double tHin, double uHin, double alpSin) {
    return (id1in + id2in + x1in + x2in + sHin + tHin + uHin + alpSin > 0.) 
    ? true : false;} 

  // Select incoming parton channel and extract parton densities (resolved).
  virtual void pickInState() {}

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol() = 0;

  // Perform kinematics for a Multiple Interaction, in its rest frame.
  virtual bool final2KinMI() {return true;}

  // Process name and code, and the number of final-state particles.
  virtual string name() const {return "unnamed process";}
  virtual int code() const {return 0;}
  virtual int nFinal() const {return 2;};

  // Special treatment needed for elastic and diffractive processes
  virtual bool hasSigmaTot() const {return false;}
  virtual bool isMinBias() const {return false;};
  virtual bool isResolved() const {return true;};
  virtual bool isDiffA() const {return false;};
  virtual bool isDiffB() const {return false;};

  // Give back particle properties: flavours, colours, masses, or all.
  int id(int i) const {return idH[i];}
  int col(int i) const {return colH[i];} 
  int acol(int i) const {return acolH[i];}
  double m(int i) const {return mH[i];}
  Particle getParton(int i) {return parton[i];}

  // Give back couplings and parton densities. Not all known for minbias.
  double Q2Ren() const {return Q2RenH;}
  double alphaEM() const {return alpEM;}
  double alphaS() const {return alpS;}
  double Q2Fac() const {return Q2FacH;}
  double pdf1() const {return pdf1H;}
  double pdf2() const {return pdf2H;}

  // Give back angles; relevant only for minimum-bias process.
  double theta() const { return thetaH;}
  double phi() const { return phiH;}

protected:

  // Constructor.
  SigmaProcess() {setDefaults();}

  // Static initialization data, normally only set once.
  static int alphaSorder, nQuark;
  static double alphaSvalue, alphaEMfix;

  // Static alphaStrong calculation.
  static AlphaStrong alphaScalc;

  // Constants: could only be changed in the code itself.
  static const double CONVERT2MB, MASSMARGIN;

  // Store Q2 renormalization and factorization scales, and related values.
  double Q2RenH, alpEM, alpS, Q2FacH, pdf1H, pdf2H;

  // Store flavour, colour, anticolour, mass, angles and the whole particle.
  int idH[6], colH[6], acolH[6];
  double mH[6], thetaH, phiH;
  Particle parton[6];

  // Set some default values at instantiation.
  void setDefaults() {alpEM = alphaEMfix; alpS = alphaSvalue;
    for (int i = 0; i < 6; ++i) mH[i] = 0.;}

  // Set flavour, colour and anticolour.
  void setId( int id1in = 0, int id2in = 0, int id3in = 0, int id4in = 0,
    int id5in = 0) {idH[1] = id1in; idH[2] = id2in; idH[3] = id3in; 
    idH[4] = id4in; idH[5] = id5in;}
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

  // The following are only used by unresolved processes.

  // Information on incoming beams.
  int idA, idB;
  double mA, mB; 
  
  // Pointer to the total/elastic/diffractive cross section object.
  SigmaTotal* sigmaTotPtr;

  // The following are only used by resolved processes.
  
  // Pointers to the incoming parton densities.
  PDF* pdfAPtrDummy;
  PDF* pdfBPtrDummy;

};
 
//**************************************************************************

// Sigma0Process is the base class for unresolved and minimum-bias processes. 
// It is derived from SigmaProcess.

class Sigma0Process : public SigmaProcess {

public:

  // Destructor.
  virtual ~Sigma0Process() {}

  // Evaluate sigma for unresolved processes. 
  virtual double sigmaHat() {return 0.;}

  // Since no PDF's there is no difference from above. 
  double sigmaPDF() {return sigmaHat();}

protected:

  // Constructor.
  Sigma0Process() {}

};
 
//**************************************************************************

// Sigma2Process is the base class for 2 -> 2 processes.
// It is derived from SigmaProcess.

class Sigma2Process : public SigmaProcess {

public:

  // Destructor.
  virtual ~Sigma2Process() {if (inFluxPtr !=0) delete inFluxPtr;}

  // Evaluate d(sigmaHat)/d(tHat) for resolved 2 -> 2 processes. 
  virtual double sigmaHat() {return 0.;}

  // Convolute d(sigmaHat)/d(tHat) with parton flux. 
  virtual double sigmaPDF() {
    return inFluxPtr->flux( x1, x2, Q2FacH) * sigmaHat();}

  // Initialize parton-flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr) = 0;

  // Input and complement kinematics for resolved 2 -> 2 process.
  bool set2Kin( double x1in, double x2in, double sHin, double tHin); 

  // Ditto, but for Multiple Interactions applications, so different input.
  virtual bool set2KinMI( int id1in, int id2in, double x1in, double x2in,
    double sHin, double tHin, double uHin, double alpSin);

  // Select incoming parton channel and extract parton densities.
  void pickInState() {inFluxPtr->pick(); id1 = inFluxPtr->id1();
    id2 = inFluxPtr->id2(); pdf1H = inFluxPtr->pdf1();
    pdf2H = inFluxPtr->pdf2();} 

  // Perform kinematics for a Multiple Interaction, in its rest frame.
  virtual bool final2KinMI();

protected:

  // Constructor.
  Sigma2Process() {masslessKin = true;}

  // Store subprocess kinematics quantities.
  int id1, id2, id3, id4;
  double x1, x2, sH, tH, uH, sH2, tH2, uH2, m3, m3S, m4, m4S, pT2;
  bool masslessKin;

  // Pointer to the parton-flux object.
  InFlux* inFluxPtr;

};
 
//**************************************************************************

// A derived class for minimum-bias (inelastic, nondiffractive) events.

class Sigma0minBias : public Sigma0Process {

public:

  // Constructor.
  Sigma0minBias() {}

  // Evaluate sigma. 
  virtual double sigmaHat() {return sigmaTotPtr->sigmaND();} 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol() {}

  // Info on the subprocess.
  virtual string name() const {return "minimum bias";}
  virtual int code() const {return 101;}
  virtual int nFinal() const {return 2;}
  virtual bool hasSigmaTot() const {return true;}
  virtual bool isMinBias() const {return true;}

private:

};
 
//**************************************************************************

// A derived class for elastic scattering A B -> A B.

class Sigma0AB2AB : public Sigma0Process {

public:

  // Constructor.
  Sigma0AB2AB() {}

  // Evaluate sigma. 
  virtual double sigmaHat() { return sigmaTotPtr->sigmaEl();} 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "A B -> A B elastic";}
  virtual int code() const {return 102;}
  virtual int nFinal() const {return 2;}
  virtual bool hasSigmaTot() const {return true;}
  virtual bool isResolved() const {return false;}

private:

};
 
//**************************************************************************

// A derived class for single diffractive scattering A B -> X B.

class Sigma0AB2XB : public Sigma0Process {

public:

  // Constructor.
  Sigma0AB2XB() {}

  // Evaluate sigma. 
  virtual double sigmaHat() { return sigmaTotPtr->sigmaXB();} 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "A B -> X B single diffractive";}
  virtual int code() const {return 103;}
  virtual int nFinal() const {return 2;}
  virtual bool hasSigmaTot() const {return true;}
  virtual bool isResolved() const {return false;}
  virtual bool isDiffA() const {return true;};

private:

};
 
//**************************************************************************

// A derived class for single diffractive scattering A B -> A X.

class Sigma0AB2AX : public Sigma0Process {

public:

  // Constructor.
  Sigma0AB2AX() {}

  // Evaluate sigma. 
  virtual double sigmaHat() { return sigmaTotPtr->sigmaAX();} 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "A B -> A X single diffractive";}
  virtual int code() const {return 104;}
  virtual int nFinal() const {return 2;}
  virtual bool hasSigmaTot() const {return true;}
  virtual bool isResolved() const {return false;}
  virtual bool isDiffB() const {return true;};

private:

};
 
//**************************************************************************

// A derived class for double diffractive scattering A B -> X X.

class Sigma0AB2XX : public Sigma0Process {

public:

  // Constructor.
  Sigma0AB2XX() {}

  // Evaluate sigma. 
  virtual double sigmaHat() { return sigmaTotPtr->sigmaXX();} 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "A B -> X X double diffractive";}
  virtual int code() const {return 105;}
  virtual int nFinal() const {return 2;}
  virtual bool hasSigmaTot() const {return true;}
  virtual bool isResolved() const {return false;}
  virtual bool isDiffA() const {return true;};
  virtual bool isDiffB() const {return true;};

private:

};
 
//**************************************************************************

// A derived class for g g -> g g.

class Sigma2gg2gg : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2gg() {}

  // Destructor.
  ~Sigma2gg2gg() {}

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

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

class Sigma2gg2qqbar : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2qqbar() {}

  // Destructor.
  ~Sigma2gg2qqbar() {}

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

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

class Sigma2qg2qg : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2qg() {}

  // Destructor.
  ~Sigma2qg2qg() {}

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

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

class Sigma2qq2qqDiff : public Sigma2Process {

public:

  // Constructor.
  Sigma2qq2qqDiff() {}

  // Destructor.
  ~Sigma2qq2qqDiff() {}

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

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

class Sigma2qq2qqSame : public Sigma2Process {

public:

  // Constructor.
  Sigma2qq2qqSame() {}

  // Destructor.
  ~Sigma2qq2qqSame() {}

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

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

class Sigma2qqbar2qqbarSame : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2qqbarSame() {}

  // Destructor.
  ~Sigma2qqbar2qqbarSame() {}

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

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

class Sigma2qqbar2qqbarNew : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2qqbarNew() {}

  // Destructor.
  ~Sigma2qqbar2qqbarNew() {}

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

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

class Sigma2qqbar2gg : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2gg() {}

  // Destructor.
  ~Sigma2qqbar2gg() {}

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

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

// A derived class for g g -> Q Qbar (Q = c or b).

class Sigma2gg2QQbar : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2QQbar(int idIn, int codeIn, string nameIn) 
    : idNew(idIn), codeSave(codeIn), nameSave(nameIn) {masslessKin = false;}

  // Destructor.
  ~Sigma2gg2QQbar() {}

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

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
  double mNew, m2New, sigTS, sigUS, sigSum;

};
 
//**************************************************************************

// A derived class for q qbar -> Q Qbar (Q = c or b).

class Sigma2qqbar2QQbar : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2QQbar(int idIn, int codeIn, string nameIn) 
    : idNew(idIn), codeSave(codeIn), nameSave(nameIn) {masslessKin = false;}

  // Destructor.
  ~Sigma2qqbar2QQbar() {}

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

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

class Sigma2qg2qgamma : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2qgamma() {}

  // Destructor.
  ~Sigma2qg2qgamma() {}

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

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

class Sigma2qqbar2ggamma : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2ggamma() {}

  // Destructor.
  ~Sigma2qqbar2ggamma() {}

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

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

class Sigma2gg2ggamma : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2ggamma() {}

  // Destructor.
  ~Sigma2gg2ggamma() {}

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

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

class Sigma2qqbar2gammagamma : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2gammagamma() {}

  // Destructor.
  ~Sigma2qqbar2gammagamma() {}

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

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

class Sigma2gg2gammagamma : public Sigma2Process {

public:

  // Constructor.
  Sigma2gg2gammagamma() {}

  // Destructor.
  ~Sigma2gg2gammagamma() {}

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Initialize parton flux object. 
  virtual void initFlux( PDF* pdfAPtr, PDF* pdfBPtr); 

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

#endif // Pythia8_SigmaProcess_H
 
