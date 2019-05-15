// Header file for incoming parton flux. Used by SigmaProcess.
// InFlux: base class for combinations of incoming partons.
// InFluxgg, InFluxqqAnti, InFluxqg, ...: derived classes.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_InFlux_H
#define Pythia8_InFlux_H

#include "Basics.h"
#include "PartonDistributions.h"
#include "PythiaStdlib.h"
#include "Settings.h"
#include "StandardModel.h"

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

  // Add even powers of CKM matrix element to weight.
  void weightCKM(int ckmPow = 0);

  // Multiply by colour factor 1/3 or 1/8 for colour annihilation graphs.
  void weightInvCol();

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

// A derived class for f fbar' incoming states with net charge +-1.

class InFluxffbarChg : public InFlux {

public:

  // Constructor.
  InFluxffbarChg() {}

private:

  // Initialize values.
  virtual void initChannels();  

};
  
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_InFlux_H
 
