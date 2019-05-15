// Header file for Les Houches Accord user process information.
// LHAinit: base class for initialization information.
// LHEevnt: Base class for event information. 
// LHAinitFortran: derived class with the HEPRUP Fortran initialization info.
// LHAevntFortran: derived class with the HEPEUP Fortran event info.
// LHAinitPythia6: derived class reading initilization file from Pythia 6.
// LHAevntPythia6: derived class reading event file from Pythia 6.
// Copyright © 2005 Torbjörn Sjöstrand

#ifndef LesHouches_H
#define LesHouches_H

// Use Stdlib facilities for input and output.
#include <iostream>
#include <fstream>
#include <iomanip>
using std::istream; 
using std::ifstream; 
using std::ostream; 
using std::ios; 
using std::fixed; 
using std::scientific; 
using std::setprecision; 
using std::setw; 

// Use vector template class from Stdlib.
#include <vector>
using std::vector; 

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
  int idBeamA() const {return idBeamAsave;}
  int idBeamB() const {return idBeamBsave;}
  double eBeamA() const {return eBeamAsave;}
  double eBeamB() const {return eBeamBsave;}
  int pdfGroupBeamA() const {return pdfGroupBeamAsave;}
  int pdfGroupBeamB() const {return pdfGroupBeamBsave;}
  int pdfSetBeamA() const {return pdfSetBeamAsave;}
  int pdfSetBeamB() const {return pdfSetBeamBsave;}
    
  // Give back weight strategy.
  int strategy() const {return strategySave;}

  // Give back info on processes.
  int size() const {return processes.size();} 
  int idProcess(int proc) const {return processes[proc].idPr;} 
  double xSec(int proc) const {return processes[proc].xSecPr;}    
  double xErr(int proc) const {return processes[proc].xErrPr;}    
  double xMax(int proc) const {return processes[proc].xMaxPr;} 
   
  // Print the info; useful to check that setting it worked.
  friend ostream& operator<<(ostream&, const LHAinit&) ;

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
  int idProc() const {return idPr;} 
  double weight() const {return weightPr;} 
  double scale() const {return scalePr;} 
  double alphaQED() const {return alphaQEDPr;} 
  double alphaQCD() const {return alphaQCDPr;} 

  // Give back info on separate particle.
  int size() const {return particles.size();}
  int id(int part) const {return particles[part].idPa;}
  int status(int part) const {return particles[part].statusPa;}
  int mother1(int part) const {return particles[part].mother1Pa;}
  int mother2(int part) const {return particles[part].mother2Pa;}
  int col1(int part) const {return particles[part].col1Pa;}
  int col2(int part) const {return particles[part].col2Pa;}
  double px(int part) const {return particles[part].pxPa;}
  double py(int part) const {return particles[part].pyPa;}
  double pz(int part) const {return particles[part].pzPa;}
  double e(int part) const {return particles[part].ePa;}
  double m(int part) const {return particles[part].mPa;}
  double tau(int part) const {return particles[part].tauPa;}
  double spin(int part) const {return particles[part].spinPa;}

  // Print the info; useful to check that reading an event worked.
  friend ostream& operator<<(ostream&, const LHAevnt&) ;

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
    particles.clear(); particle(0); }

  // Input particle info, one particle at the time.
  void particle(int idIn, int statusIn = 0, int mother1In = 0, 
    int mother2In = 0, int col1In = 0, int col2In = 0, double pxIn = 0., 
    double pyIn = 0., double pzIn = 0., double eIn = 0., double mIn = 0., 
    double tauIn = 0., double spinIn = 9.) { 
    particles.push_back( Particle( idIn, statusIn, mother1In, mother2In, 
    col1In, col2In, pxIn, pyIn, pzIn, eIn, mIn, tauIn, spinIn)); }
  
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
};

//**************************************************************************

// A derived class with initialization information from 
// the HEPRUP Fortran commonblock.

class LHAinitFortran : public LHAinit {

public:

  // Constructor.
  LHAinitFortran() {}

  // Routine for doing the job of setting initialization info.  
  bool set(); 

};

//**************************************************************************

// A derived class with event information from 
// the HEPEUP Fortran commonblock.

class LHAevntFortran: public LHAevnt {

public:

  // Constructor.
  LHAevntFortran() {}

  // Routine for doing the job of setting info on next event.  
  bool set(); 

};

//**************************************************************************

// A derived class with initialization information from Pythia 6.3.

class LHAinitPythia6 : public LHAinit {

public:

  // Constructor.
  LHAinitPythia6(const char* fileIn) : is(fileIn) {}

  // Routine for doing the job of reading and setting initialization info.  
  bool set(); 

private:
 
  // File from which to read.
  ifstream is;

};

//**************************************************************************

// A derived class with event information from Pythia 6.3.

class LHAevntPythia6 : public LHAevnt {

public:

  // Constructor.
  LHAevntPythia6(const char* fileIn) : is(fileIn) {}

  // Routine for doing the job of reading and setting info on next event.  
  bool set(); 

private:
 
  // File from which to read.
  ifstream is;

};

//**************************************************************************

#endif // LesHouches_H
