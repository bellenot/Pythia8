// Header file for Les Houches Accord user process information.
// LHAinit: base class for initialization information.
// LHEevnt: Base class for event information. 
// LHAinitLHEF: derived class for initilization from Les Houches Event File.
// LHAevntLHEF: derived class for events from Les Houches Evewnt File.
// Code for interfacing with Fortran commonblocks is found in LHAFortran.h.
// Copyright C 2007 Torbjorn Sjostrand

#ifndef Pythia8_LesHouches_H
#define Pythia8_LesHouches_H

#include "PythiaStdlib.h"

namespace Pythia8 {

//**************************************************************************

// LHAinit is base class for initialization information 
// from an external parton-level generator.

class LHAinit {

public:

  // A pure virtual method set, wherein all initialization information 
  // is supposed to be set in the derived class. Can do this by reading a 
  // file or some other way, as desired. Returns false if it did not work. 
  virtual bool set() = 0; 

  // Give back info on beams.
  int    idBeamA() const {return idBeamAsave;}
  int    idBeamB() const {return idBeamBsave;}
  double eBeamA() const {return eBeamAsave;}
  double eBeamB() const {return eBeamBsave;}
  int    pdfGroupBeamA() const {return pdfGroupBeamAsave;}
  int    pdfGroupBeamB() const {return pdfGroupBeamBsave;}
  int    pdfSetBeamA() const {return pdfSetBeamAsave;}
  int    pdfSetBeamB() const {return pdfSetBeamBsave;}
    
  // Give back weight strategy.
  int    strategy() const {return strategySave;}

  // Give back info on processes.
  int    size() const {return processes.size();} 
  int    idProcess(int proc) const {return processes[proc].idPr;} 
  double xSec(int proc) const {return processes[proc].xSecPr;}    
  double xErr(int proc) const {return processes[proc].xErrPr;}    
  double xMax(int proc) const {return processes[proc].xMaxPr;} 
   
  // Print the info; useful to check that setting it worked.
  void   list(ostream& os = cout);  

protected:

  // Constructor. Sets default to be that events come with unit weight.
  LHAinit(int strategyIn = 3) : strategySave(strategyIn) 
    { processes.reserve(10);} 

  // Destructor.
  virtual ~LHAinit() {}

  // Input beam info.
  void beamA(int idIn, double eIn, int pdfGroupIn = 0, int pdfSetIn = 0) 
    { idBeamAsave = idIn; eBeamAsave = eIn; pdfGroupBeamAsave = pdfGroupIn;  
    pdfSetBeamAsave = pdfSetIn;} 
  void beamB(int idIn, double eIn, int pdfGroupIn = 0, int pdfSetIn = 0) 
    { idBeamBsave = idIn; eBeamBsave = eIn; pdfGroupBeamBsave = pdfGroupIn;  
    pdfSetBeamBsave = pdfSetIn;} 

  // Input process weight strategy.
  void strategy(int strategyIn) {strategySave = strategyIn;} 

  // Input process info.
  void process(int idIn, double xSecIn = 1., double xErrIn = 0., 
    double xMaxIn = 1.) 
    { processes.push_back( Process( idIn, xSecIn, xErrIn, xMaxIn)); }

private:

  // Event weighting and mixing strategy.
  int strategySave;

  // Beam particle properties.
  int idBeamAsave, idBeamBsave;
  double eBeamAsave, eBeamBsave;
  int pdfGroupBeamAsave, pdfGroupBeamBsave, pdfSetBeamAsave, pdfSetBeamBsave;

  // A nested class for processes...
  class Process {
  public:
    Process() : idPr(0), xSecPr(0.), xErrPr(0.), xMaxPr(0.) { }
    Process(int idIn, double xSecIn, double xErrIn, double xMaxIn) :
      idPr(idIn), xSecPr(xSecIn), xErrPr(xErrIn), xMaxPr(xMaxIn) { }
    int idPr;
    double xSecPr, xErrPr, xMaxPr;
  } ;

  // ...so that the process list can be kept as a vector.
  vector<Process> processes;
};

//**************************************************************************

// LHAevnt is base class for event information 
// from an external parton-level generator.

class LHAevnt {

public:

  // A pure virtual method set, wherein all information on the next event
  // is supposed to be set in the derived class. Can do this by reading a 
  // file or some other way, as desired. Returns false if it did not work. 
  virtual bool set() = 0; 

  // Give back process number, weight, scale, alpha_em, alpha_s.
  int    idProc() const {return idPr;} 
  double weight() const {return weightPr;} 
  double scale() const {return scalePr;} 
  double alphaQED() const {return alphaQEDPr;} 
  double alphaQCD() const {return alphaQCDPr;} 

  // Give back info on separate particle.
  int    size() const {return particles.size();}
  int    id(int part) const {return particles[part].idPa;}
  int    status(int part) const {return particles[part].statusPa;}
  int    mother1(int part) const {return particles[part].mother1Pa;}
  int    mother2(int part) const {return particles[part].mother2Pa;}
  int    col1(int part) const {return particles[part].col1Pa;}
  int    col2(int part) const {return particles[part].col2Pa;}
  double px(int part) const {return particles[part].pxPa;}
  double py(int part) const {return particles[part].pyPa;}
  double pz(int part) const {return particles[part].pzPa;}
  double e(int part) const {return particles[part].ePa;}
  double m(int part) const {return particles[part].mPa;}
  double tau(int part) const {return particles[part].tauPa;}
  double spin(int part) const {return particles[part].spinPa;}

  // Optional: give back info on parton density values of event.
  bool   pdfIsSet() const {return pdfIsSetSv;}
  int    id1() const {return id1Sv;}
  int    id2() const {return id2Sv;}
  double x1() const {return x1Sv;}
  double x2() const {return x2Sv;}
  double scalePDF() const {return scalePDFSv;}
  double xpdf1() const {return xpdf1Sv;}
  double xpdf2() const {return xpdf2Sv;}

  // Print the info; useful to check that reading an event worked.
  void   list(ostream& os = cout);  

protected:

  // Constructor.
  LHAevnt() { particles.reserve(20); }

  // Destructor.
  virtual ~LHAevnt() {}
 
  // Input info on the selected process.
  void process(int idProcIn = 0, double weightIn = 1., double scaleIn = 0.,
    double alphaQEDIn = 0.0073, double alphaQCDIn = 0.12) 
    { idPr = idProcIn; weightPr = weightIn; scalePr = scaleIn; 
    alphaQEDPr = alphaQEDIn; alphaQCDPr = alphaQCDIn; 
    // Clear particle list. Add empty zeroth particle for correct indices.
    particles.clear(); particle(0); pdfIsSetSv = false;}

  // Input particle info, one particle at the time.
  void particle(int idIn, int statusIn = 0, int mother1In = 0, 
    int mother2In = 0, int col1In = 0, int col2In = 0, double pxIn = 0., 
    double pyIn = 0., double pzIn = 0., double eIn = 0., double mIn = 0., 
    double tauIn = 0., double spinIn = 9.) { 
    particles.push_back( Particle( idIn, statusIn, mother1In, mother2In, 
    col1In, col2In, pxIn, pyIn, pzIn, eIn, mIn, tauIn, spinIn)); }

  // Optionally input info on parton density values of event.
  void pdf(int id1In, int id2In, double x1In, double x2In, 
    double scalePDFIn, double xpdf1In, double xpdf2In) 
    { id1Sv = id1In; id2Sv = id2In; x1Sv = x1In; x2Sv = x2In;
    scalePDFSv = scalePDFIn; xpdf1Sv = xpdf1In; xpdf2Sv = xpdf2In;
    pdfIsSetSv = true;}

private:

  // Store info on the selected process. 
  int idPr;
  double weightPr, scalePr, alphaQEDPr, alphaQCDPr;

  // A nested class for particles...
  class Particle {
  public:
    Particle() : idPa(0), statusPa(0), mother1Pa(0), mother2Pa(0),
      col1Pa(0), col2Pa(0), pxPa(0.), pyPa(0.), pzPa(0.), ePa(0.),
      mPa(0.), tauPa(0.), spinPa(9.) { }
    Particle(int idIn, int statusIn, int mother1In, int mother2In,
      int col1In, int col2In, double pxIn, double pyIn, double pzIn, 
      double eIn, double mIn, double tauIn, double spinIn) :
      idPa(idIn), statusPa(statusIn), mother1Pa(mother1In), 
      mother2Pa(mother2In), col1Pa(col1In), col2Pa(col2In), pxPa(pxIn), 
      pyPa(pyIn), pzPa(pzIn), ePa(eIn), mPa(mIn), tauPa(tauIn), 
      spinPa(spinIn) { }
    int idPa, statusPa, mother1Pa, mother2Pa, col1Pa, col2Pa ;
    double pxPa, pyPa, pzPa, ePa, mPa, tauPa, spinPa ;
  } ;

  // ...so that the particle list can be kept as a vector.
  vector<Particle> particles;

  // Optional info on parton density values of event.
  bool   pdfIsSetSv;
  int    id1Sv, id2Sv;
  double x1Sv, x2Sv, scalePDFSv, xpdf1Sv, xpdf2Sv;

};

//**************************************************************************

// A derived class with initialization information from Les Houches Event File.

class LHAinitLHEF : public LHAinit {

public:

  // Constructor.
  LHAinitLHEF(const char* fileIn) : is(fileIn) {}

  // Routine for doing the job of reading and setting initialization info.  
  bool set(); 

private:
 
  // File from which to read.
  ifstream is;

};

//**************************************************************************

// A derived class with event information from Les Houches Event File.

class LHAevntLHEF : public LHAevnt {

public:

  // Constructor.
  LHAevntLHEF(const char* fileIn) : is(fileIn) {}

  // Routine for doing the job of reading and setting info on next event.  
  bool set(); 

private:
 
  // File from which to read.
  ifstream is;

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_LesHouches_H
