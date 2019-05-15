// Header file for electroweak process differential cross sections.
// Contains classes derived from SigmaProcess via Sigma2Process.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_SigmaEW_H
#define Pythia8_SigmaEW_H

#include "SigmaProcess.h"

namespace Pythia8 {

 
//**************************************************************************

// A derived class for q g -> q gamma (q = u, d, s, c, b).
// Use massless approximation also for Q since no alternative.

class Sigma2qg2qgamma : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2qgamma() {}

  // Destructor.
  ~Sigma2qg2qgamma() {}

  // Initialize process, especially parton-flux object. 
  virtual void initProc(); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "q g -> q gamma (udscb)";}
  virtual int code() const {return 131;}

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

  // Initialize process, especially parton-flux object. 
  virtual void initProc(); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "q qbar -> g gamma";}
  virtual int code() const {return 132;}

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

  // Initialize process, especially parton-flux object. 
  virtual void initProc(); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "g g -> g gamma";}
  virtual int code() const {return 133;}

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

  // Initialize process, especially parton-flux object. 
  virtual void initProc(); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "q qbar -> gamma gamma";}
  virtual int code() const {return 134;}

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

  // Initialize process, especially parton-flux object. 
  virtual void initProc(); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "g g -> gamma gamma";}
  virtual int code() const {return 135;}

private:
  
  double charge2Sum;

};
 
//**************************************************************************

// A derived class for f f' -> f f' via t-channel gamma*/Z0 exchange.

class Sigma2ff2ff9gmZ : public Sigma2Process {

public:

  // Constructor.
  Sigma2ff2ff9gmZ() {}

  // Destructor.
  ~Sigma2ff2ff9gmZ() {}

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat();

  // Initialize process, especially parton-flux object. 
  virtual void initProc(); 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "f f' -> f f' (t-channel gamma*/Z0)";}
  virtual int code() const {return 141;}

private:

  //  Z parameters for propagator.
  int gmZmode;
  double mZ, mZS, thetaWRat;

};
 
//**************************************************************************

// A derived class for f_1 f_2 -> f_3 f_4 via t-channel W+- exchange.

class Sigma2ff2ff9W : public Sigma2Process {

public:

  // Constructor.
  Sigma2ff2ff9W() {}

  // Destructor.
  ~Sigma2ff2ff9W() {}

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat();

  // Initialize process, especially parton-flux object. 
  virtual void initProc(); 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "f_1 f_2 -> f_3 f_4 (t-channel W+-)";}
  virtual int code() const {return 142;}

private:

  //  W parameters for propagator.
  double mW, mWS, thetaWRat;

};
 
//**************************************************************************

// A derived class for q q' -> Q q" via t-channel W+- exchange.
// Related to Sigma2ff2ffViaW class, but with massive matrix elements.

class Sigma2qq2Qq9W : public Sigma2Process {

public:

  // Constructor.
  Sigma2qq2Qq9W(int idIn, int codeIn, string nameIn) 
    : idNew(idIn), codeSave(codeIn), nameSave(nameIn) {}

  // Destructor.
  ~Sigma2qq2Qq9W() {}

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat();

  // Initialize process, especially parton-flux object. 
  virtual void initProc(); 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return nameSave;}
  virtual int code() const {return codeSave;}
  virtual int id3Mass() const {return idNew;}

private:

  // Values stored for process type. W parameters for propagator.
  int idNew, codeSave;
  string nameSave;
  double mW, mWS, thetaWRat;

};
 
//**************************************************************************

// A derived class for f fbar -> gamma*/Z0.

class Sigma1ffbar2gmZ : public Sigma1Process {

public:

  // Constructor.
  Sigma1ffbar2gmZ() {}

  // Destructor.
  ~Sigma1ffbar2gmZ() {}

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat();

  // Initialize process, especially parton-flux object. 
  virtual void initProc(); 

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "f fbar' -> gamma*/Z0";}
  virtual int code() const {return 151;}
  virtual int resonanceA() const {return 23;}

private:

  // A Z0 resonance object provides coupling and propagator expressions.
  ResonanceGmZ GmZRes;

};
 
//**************************************************************************

// A derived class for f fbar' -> W+-.

class Sigma1ffbar2W : public Sigma1Process {

public:

  // Constructor.
  Sigma1ffbar2W() {}

  // Destructor.
  ~Sigma1ffbar2W() {}

  // Initialize process, especially parton-flux object. 
  virtual void initProc(); 

  // Evaluate sigmaHat(sHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "f fbar' -> W+-";}
  virtual int code() const {return 152;}
  virtual int resonanceA() const {return 24;}

private:

  // A W+- resonance object provides coupling and propagator expressions.
  ResonanceW WRes;

};
 
//**************************************************************************

// A derived class for f fbar' -> Z0 W+-. (Here pure Z0, unfortunately.)

class Sigma2ffbar2ZW : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2ZW() {}

  // Destructor.
  ~Sigma2ffbar2ZW() {}

  // Initialize process, especially parton-flux object. 
  virtual void initProc(); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "f fbar' -> Z0 W+- (no gamma*!)";}
  virtual int code() const {return 162;}
  virtual int id3Mass() const {return 23;}
  virtual int id4Mass() const {return 24;}
  virtual int resonanceA() const {return 24;}

private:

  // Store W+- mass and width, and couplings.
  double mW, widW, mWS, mwWS, thetaWRat, thetaWpt, thetaWmm, lu, ld;

};
 
//**************************************************************************

// A derived class for f fbar -> W+ W-.

class Sigma2ffbar2WW : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2WW() {}

  // Destructor.
  ~Sigma2ffbar2WW() {}

  // Initialize process, especially parton-flux object. 
  virtual void initProc(); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "f fbar -> W+ W-";}
  virtual int code() const {return 163;}
  virtual int id3Mass() const {return 24;}
  virtual int id4Mass() const {return -24;}
  virtual int resonanceA() const {return 23;}

private:

  // Store Z0 mass and width.
  double mZ, widZ, mZS, mwZS, thetaWRat;

};
 
//**************************************************************************

// A derived class for q qbar' -> W+- g.

class Sigma2qqbar2Wg : public Sigma2Process {

public:

  // Constructor.
  Sigma2qqbar2Wg() {}

  // Destructor.
  ~Sigma2qqbar2Wg() {}

  // Initialize process, especially parton-flux object. 
  virtual void initProc(); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "q qbar' -> W+- g";}
  virtual int code() const {return 176;}
  virtual int id3Mass() const {return 24;}

private:

};
 
//**************************************************************************

// A derived class for q g -> W+- q'.

class Sigma2qg2Wq : public Sigma2Process {

public:

  // Constructor.
  Sigma2qg2Wq() {}

  // Destructor.
  ~Sigma2qg2Wq() {}

  // Initialize process, especially parton-flux object. 
  virtual void initProc(); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "q g-> W+- q'";}
  virtual int code() const {return 177;}
  virtual int id3Mass() const {return 24;}

private:

};
 
//**************************************************************************

// A derived class for f fbar' -> W+- gamma.

class Sigma2ffbar2Wgm : public Sigma2Process {

public:

  // Constructor.
  Sigma2ffbar2Wgm() {}

  // Destructor.
  ~Sigma2ffbar2Wgm() {}

  // Initialize process, especially parton-flux object. 
  virtual void initProc(); 

  // Evaluate d(sigmaHat)/d(tHat). 
  virtual double sigmaHat();

  // Select flavour, colour and anticolour.
  virtual void setIdColAcol();

  // Info on the subprocess.
  virtual string name() const {return "f fbar' -> W+- gamma";}
  virtual int code() const {return 178;}
  virtual int id3Mass() const {return 24;}

private:

};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_SigmaEW_H
