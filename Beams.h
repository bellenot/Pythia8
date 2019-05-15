// Header file for information on incoming beams.
// ResolvedParton: an initiator or remnant in beam.
// BeamParticle: contains partons, parton densities, etc.
// BeamRemnants: matches the remnants between the two beams.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_Beams_H
#define Pythia8_Beams_H

#include "Stdlib.h"
#include "Basics.h"
#include "Settings.h"
#include "ParticleData.h"
#include "PartonDistributions.h"
#include "Event.h"
#include "FragmentationFlavZpT.h"

namespace Pythia8 {

//**************************************************************************

// This class holds info on a parton resolved inside the incoming beam,
// i.e. either an initiator (part a hard scattering or multiple interaction)
// or a remnant (part of the beam remnant treatment).

// The companion code is -1 from onset and for g, is -2 for an unmatched 
// sea quark, is >= 0 for a matched sea quark, with the number giving the 
// companion position, and is -3 for a valence quark.

class ResolvedParton {

public:

  // Constructor.
  ResolvedParton( int lineIn = 0, int idIn = 0, double xIn = 0., 
    int companionIn = -1) : lineRes(lineIn), idRes(idIn), xRes(xIn), 
    companionRes(companionIn), xqCompRes(0.), mRes(0.), colRes(0),
    acolRes(0) { } 

  // Set info on initiator or remnant parton.
  void line( int lineIn) {lineRes = lineIn;} 
  void id( int idIn) {idRes = idIn;} 
  void x( double xIn) {xRes = xIn;} 
  void lineidx( int lineIn, int idIn, double xIn) {lineRes = lineIn;
    idRes = idIn; xRes = xIn;} 
  void companion( int companionIn) {companionRes = companionIn;} 
  void xqCompanion( double xqCompIn) {xqCompRes = xqCompIn;} 
  void p(Vec4 pIn) {pRes = pIn;}
  void px(double pxIn) {pRes.px(pxIn);}
  void py(double pyIn) {pRes.py(pyIn);}
  void pz(double pzIn) {pRes.pz(pzIn);}
  void e(double eIn) {pRes.e(eIn);}
  void m(double mIn) {mRes = mIn;}
  void col(int colIn) {colRes = colIn;}
  void acol(int acolIn) {acolRes = acolIn;}
  void cols(int colIn = 0,int acolIn = 0) 
    {colRes = colIn; acolRes = acolIn;}  

  // Get info on initiator or remnant parton.
  int line() const {return lineRes;} 
  int id() const {return idRes;} 
  double x() const {return xRes;} 
  int companion() const {return companionRes;} 
  bool isValence() const {return (companionRes == -3) ? true : false;}
  bool isUnmatched() const {return (companionRes == -2) ? true : false;}
  bool isCompanion() const {return (companionRes >= 0) ? true : false;}
  double xqCompanion() const {return xqCompRes;} 
  Vec4 p() const {return pRes;}
  double px() const {return pRes.px();}
  double py() const {return pRes.py();}
  double pz() const {return pRes.pz();}
  double e() const {return pRes.e();}
  double m() const {return mRes;}
  double pT() const {return pRes.pT();}
  double mT2() const {return mRes*mRes + pRes.pT2();}
  int col() const {return colRes;}
  int acol() const {return acolRes;}
 
private:

  // Properties of a resolved parton. 
  int lineRes, idRes;
  double xRes;
  // Companion code and distribution value, if any.
  int companionRes; 
  double xqCompRes;
  // Four-momentum and mass; for remnant kinematics construction.
  Vec4 pRes;
  double mRes;
  // Colour codes.
  int colRes, acolRes;

};

//**************************************************************************

// This class holds info on a beam particle in the evolution of 
// initial-state radiation and multiple interactions.

class BeamParticle {

public:

  // Constructor.
  BeamParticle() { }  

  // Initialize static data members.
  static void initStatic();

  // Initialize. Possibility to force re-initialization by hand.
  void init( int idIn, double pzIn, double eIn, double mIn, PDF* pdfInPtr);

  // Member functions for detailed input.
  void id(int idIn) {idBeam = idIn; initBeamKind();}
  void p(Vec4 pIn) {pBeam = pIn;}
  void px(double pxIn) {pBeam.px(pxIn);}
  void py(double pyIn) {pBeam.py(pyIn);}
  void pz(double pzIn) {pBeam.pz(pzIn);}
  void e(double eIn) {pBeam.e(eIn);}
  void m(double mIn) {mBeam = mIn;}
  void pdf(PDF* pdfInPtr) {pdfBeamPtr = pdfInPtr;}

  // Member functions for output.
  int id() const {return idBeam;}
  Vec4 p() const {return pBeam;}
  double px() const {return pBeam.px();}
  double py() const {return pBeam.py();}
  double pz() const {return pBeam.pz();}
  double e() const {return pBeam.e();}
  double m() const {return mBeam;}
  PDF* pdf() const {return pdfBeamPtr;}
  bool isLepton() const {return isLeptonBeam;}
  // As hadrons here we only count those we know how to handle remnants for.
  bool isHadron() const {return isHadronBeam;}
  bool isMeson() const {return isMesonBeam;}
  bool isBaryon() const {return isBaryonBeam;}

  // Maximum x remaining after previous MI and ISR, plus safety margin.
  double xMax(int iSkip = -1);
 
  // Standard parton distributions, used for hard interactions.
  double xf(int id, double x, double Q2) 
    {return pdfBeamPtr->xf(id, x, Q2);}

  // Ditto, split into valence and sea parts (where gluon counts as sea).
  double xfVal(int id, double x, double Q2) 
    {return pdfBeamPtr->xfVal(id, x, Q2);}
  double xfSea(int id, double x, double Q2) 
    {return pdfBeamPtr->xfSea(id, x, Q2);}

  // Rescaled parton distributions, as needed for MI and ISR.
  // For ISR also allow split valence/sea, and only return relevant part.
  double xfMI(int id, double x, double Q2) 
    {return xfModified(-1, id, x, Q2);}
  double xfISR(int indexMI, int id, double x, double Q2) 
    {return xfModified( indexMI, id, x, Q2);}

  // Decide whether chosen quark is valence, sea or companion.
  void pickValSeaComp();

  // Initialize kind of incoming beam particle.
  void initBeamKind();

  // Overload index operator to access a resolved parton from the list.
  ResolvedParton& operator[](int i) {return resolved[i];}

  // Total number of partons extracted from beam, and initiators only.
  int size() const {return resolved.size();}
  int sizeInit() const {return nInit;}

  // Clear list of resolved partons. 
  void clear() {resolved.resize(0);}

  // Add a resolved parton to list. 
  int append( int line, int id, double x, int companion = -1)
    {resolved.push_back( ResolvedParton( line, id, x, companion) );
    return resolved.size() - 1;}

  // Print extracted parton list; for debug mainly.
  void list(ostream& os = cout); 

  // How many different flavours, and how many quarks of given flavour.
  int nValenceKinds() const {return nValKinds;}
  int nValence(int id) const {for (int i = 0; i < nValKinds; ++i) 
    if (id == idVal[i]) return nVal[i]; return 0;}

  // Add extra remnant flavours to make valence and sea come out right. 
  bool remnantFlavours(); 

  // Correlate all initiators and remnants to make a colour singlet. 
  bool remnantColours(Event& event, vector<int>& colFrom,
    vector<int>& colTo); 

  // Tell whether a junction has been resolved, and its junction colours.
  bool hasJunction() const {return hasJunctionBeam;}  
  int junctionCol(int i) const {return junCol[i];}
  void junctionCol(int i, int col) {junCol[i] = col;}

private: 

  // Static initialization data, normally only set once.
  static int maxValQuark, companionPower;
  static bool allowJunction;

  // Basic properties of a beam particle.
  int idBeam, idBeamAbs;  
  Vec4 pBeam;
  double mBeam;
  PDF* pdfBeamPtr;

  // Beam kind. Valence flavour content for hadrons.
  bool isLeptonBeam, isHadronBeam, isMesonBeam, isBaryonBeam;
  int nValKinds, idVal[3], nVal[3];

  // Current parton density, by valence, sea and companion.
  int idSave, iSkipSave, nValLeft[3]; 
  double xqgTot, xqVal, xqgSea, xqCompSum;

  // The list of resolved partons.
  vector<ResolvedParton> resolved; 

  // Put vectors here rather than in members to avoid segmentation fault.
  vector<int> iVal, iGlu, iGluRndm, colList, acolList;        

  // Status after all initiators have been accounted for. Junction content.
  int nInit;
  bool hasJunctionBeam;
  int junCol[3];

  // Routine to calculate pdf's given previous interactions.
  double xfModified( int iSkip, int id, double x, double Q2); 

  // Fraction of hadron momentum sitting in a valence quark distribution.
  double xValFrac(int j, double Q2);
  double Q2ValFracSav, uValInt, dValInt;

  // Fraction of hadron momentum sitting in a companion quark distribution.
  double xCompFrac(double xs);

  // Value of companion quark PDF, also given the sea quark x.
  double xCompDist(double xc, double xs);

};

//**************************************************************************

// This class matches the kinematics of the hard-scattering subsystems
// (with primordial kT added) to that of the two beam remnants.  

class BeamRemnants {

public:

  // Constructor.
  BeamRemnants() { }  

  // Initialize static data members.
  static void initStatic();

  // Do the matching/kinematics of the two beam remnants. 
  bool add( BeamParticle& beamA, BeamParticle& beamB, Event& event);

private: 

  // Static initialization data, normally only set once.
  static double primordialKTwidth, valencePowerMeson, valencePowerUinP,
    valencePowerDinP, valenceDiqEnhance;
  static int companionPower;

  // Constants: could only be changed in the code itself.
  static const int NTRYCOLMATCH, NTRYKINMATCH;

  // Put vectors here rather than in functions to avoid segmentation fault.
  vector<int> colSave, acolSave, colFrom, colTo, colList, acolList;

  // Pick unrescaled x of remnant parton (valence or sea).
  double xRemnant(int i, BeamParticle& beam);

  // Check that colours are consistent.
  bool checkColours( BeamParticle& beamA, BeamParticle& beamB, Event& event);

};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_Beams_H
