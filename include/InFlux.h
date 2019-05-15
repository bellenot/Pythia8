// Header file for incoming parton flux. Used by SigmaProcess.
// InBeam, InPair: simple helper classes.
// InFlux: base class for combinations of incoming partons.
// InFluxgg, InFluxqqAnti, InFluxqg, ...: derived classes.
// Copyright C 2007 Torbjorn Sjostrand

#ifndef Pythia8_InFlux_H
#define Pythia8_InFlux_H

#include "Basics.h"
#include "Beams.h"
#include "PartonDistributions.h"
#include "PythiaStdlib.h"
#include "Settings.h"
#include "StandardModel.h"

namespace Pythia8 {

//**************************************************************************

// InBeam is a simple helper class for partons and their flux in a beam.

class InBeam {

public:

  // Constructor.
  InBeam( int idIn = 0) : id(idIn), pdf(0.) {}

  // Values.
  int id;
  double pdf;

};

//**************************************************************************

// InPair is a simple helper class for colliding parton pairs and their flux.

class InPair {

public:

  // Constructor.
  InPair( int idAIn = 0, int idBIn = 0) : idA(idAIn), idB(idBIn), 
    fixWeight(1.), varWeight(1.), pdfA(0.), pdfB(0.), fluxWeight(0.) {}

  // Values.
  int idA, idB;
  double fixWeight, varWeight, pdfA, pdfB, fluxWeight;

};

//**************************************************************************

// InFlux is the base class for the combined incoming parton flux.

class InFlux {

public:

  // Destructor.
  virtual ~InFlux() {}

  // Initialize static data members.
  static void initStatic();

  // Initialization: store parton-density pointers. 
  static void setBeamPtr( BeamParticle* beamAPtrIn, 
    BeamParticle* beamBPtrIn) {beamAPtr = beamAPtrIn; 
    beamBPtr = beamBPtrIn;} 

  // Initialization of process-specific allowed combinations. 
  virtual void initChannels() = 0; 

  // Multiply fixed weight by squared charge.
  void weightCharge2();

  // Multiply fixed weight by squared CKM matrix elements.
  void weightCKM2();

  // Multiply fixed weight by sum of squared CKM matrix elements.
  void weightCKM2sum(int mode = 1, int idQ = 0);

  // Multiply by colour factor 1/3 or 1/8 for colour annihilation graphs.
  void weightInvCol();

  // Multiply fixed weight by spin factor 2 for neutrinos.
  void weightNeutrinoSpin();

  // Multiply fixed weight by flavour-specific factor (catch-all).
  void weightFixed(double nowWeight = 1.) {
    for (int i = 0; i < sizePair(); ++i) inPair[i].fixWeight *= nowWeight;} 
  void weightFixed(int id1, int id2, double nowWeight, 
    bool flipSide = true, bool conjugate = true, bool allGen = true);

  // Remove empty channels and optionally list remaining ones.
  void checkChannels(string processName, ostream& os = cout);

  // Introduce flavour-specific event-by-event weight factor.
  void weightInState(double nowWeight = 1.) {
    for (int i = 0; i < sizePair(); ++i) inPair[i].varWeight = nowWeight;} 
  void weightInState(int id1, int id2, double nowWeight, 
    bool flipSide = true, bool conjugate = true, bool allGen = true);

  // Calculate products of parton densities for allowed combinations.
  double flux(double x1, double x2, double Q2);

  // Information on combination of required partons from the two beams.
  int nAB() const {return sizePair();}
  int idA(int i) const {return inPair[i].idA;} 
  int idB(int i) const {return inPair[i].idB;} 
  double fixWeightAB(int i) const {return inPair[i].fixWeight;} 
  double varWeightAB(int i) const {return inPair[i].varWeight;} 
  double fluxWeightAB(int i) const {return inPair[i].fluxWeight;} 

  // Pick one of the possible channels according to their weight.
  void pick();
  int id1() const {return idNow1;}
  int id2() const {return idNow2;}
  double pdf1() const {return pdfNow1;}
  double pdf2() const {return pdfNow2;}

protected:

  // Constructor.
  InFlux() {}

  // Static initialization data, normally only set once.
  static int nQuark;
  static bool showChannels;

  // Static pointers to beams.
  static BeamParticle* beamAPtr;
  static BeamParticle* beamBPtr;

  // Constants: could only be changed in the code itself.
  static const int MAPTOFIRST[40] ;

  // Partons in beams, with pdf's.
  vector<InBeam> inBeamA;
  vector<InBeam> inBeamB;
  void addBeamA(int id) {inBeamA.push_back(InBeam(id));}
  void addBeamB(int id) {inBeamB.push_back(InBeam(id));}
  int sizeBeamA() const {return inBeamA.size();}
  int sizeBeamB() const {return inBeamB.size();}

  // Allowed colliding parton pairs, with pdf's.
  vector<InPair> inPair;
  void addPair(int idA, int idB) {inPair.push_back(InPair(idA, idB));}
  int sizePair() const {return inPair.size();}

  // Currently picked incoming channel, and flux sum of all channels.
  int idNow1, idNow2;
  double pdfNow1, pdfNow2, fluxwtSum;

};
 
//**************************************************************************

// A derived class for g g incoming state.

class InFluxgg : public InFlux {

public:

  // Constructor.
  InFluxgg() {initChannels();}

private:

  // Initialize values.
  virtual void initChannels();  

};
 
//**************************************************************************

// A derived class for q g incoming state.

class InFluxqg : public InFlux {

public:

  // Constructor.
  InFluxqg() {initChannels();}

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
  InFluxqqbarqqDiff() {initChannels();}

private:

  // Initialize values.
  virtual void initChannels();  

};
 
//**************************************************************************

// A derived class for q q' (qbar qbar') incoming states, with q' != q.

class InFluxqqDiff : public InFlux {

public:

  // Constructor.
  InFluxqqDiff() {initChannels();}

private:

  // Initialize values.
  virtual void initChannels();  

};
 
//**************************************************************************

// A derived class for q q (qbar qbar) incoming states.

class InFluxqqSame : public InFlux {

public:

  // Constructor.
  InFluxqqSame() {initChannels();}

private:

  // Initialize values.
  virtual void initChannels();  

};
 
//**************************************************************************

// A derived class for q qbar'  incoming state.

class InFluxqqbarDiff : public InFlux {

public:

  // Constructor.
  InFluxqqbarDiff() {initChannels();}

private:

  // Initialize values.
  virtual void initChannels();  

};
 
//**************************************************************************

// A derived class for q qbar antiparticle incoming state.

class InFluxqqbarSame : public InFlux {

public:

  // Constructor.
  InFluxqqbarSame() {initChannels();}

private:

  // Initialize values.
  virtual void initChannels();  

};
 
//**************************************************************************

// A derived class for f f', f fbar', fbar fbar' incoming state,
// where f an f' may be same or different.

class InFluxff : public InFlux {

public:

  // Constructor.
  InFluxff() {initChannels();}

private:

  // Initialize values.
  virtual void initChannels();  

};
 
//**************************************************************************

// A derived class for f fbar antiparticle incoming state.

class InFluxffbarSame : public InFlux {

public:

  // Constructor.
  InFluxffbarSame() {initChannels();}

private:

  // Initialize values.
  virtual void initChannels();  

};
 
//**************************************************************************

// A derived class for f fbar' incoming states with net charge +-1.

class InFluxffbarChg : public InFlux {

public:

  // Constructor.
  InFluxffbarChg() {initChannels();}

private:

  // Initialize values.
  virtual void initChannels();  

};
  
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_InFlux_H
 
