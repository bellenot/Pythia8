// Header file for the Sphericity and CellJet classes.
// Sphericity: sphericity analysis of the event.
// CellJet: calorimetric cone jet finder. 
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_Analysis_H
#define Pythia8_Analysis

#include "Basics.h"
#include "Event.h"
#include "Information.h"
#include "PythiaStdlib.h"
#include "Settings.h"

namespace Pythia8 {

//**************************************************************************

// Sphericity class.
// This class performs (optionally modified) sphericity analysis on an event.

class Sphericity {

public: 

  // Constructor.
  Sphericity(double powerIn = 2., int selectIn = 2) 
    : power(powerIn), select(selectIn) {powerInt = 0;
    if (abs(power - 1.) < 0.01) powerInt = 1;
    if (abs(power - 2.) < 0.01) powerInt = 2;
    powerMod = 0.5 * power - 1.;}
  
  // Analyze event.
  bool analyze(Event& event);

  // Return info on results of analysis.
  double sph() const {return 1.5 * (eVal2 + eVal3);}
  double apl() const {return 1.5 * eVal3;}
  double eigenValue(int i) const {return (i < 2) ? eVal1 :
    ( (i < 3) ? eVal2 : eVal3 ) ;}
  Vec4 eigenVector(int i) const {return (i < 2) ? eVec1 :
    ( (i < 3) ? eVec2 : eVec3 ) ;}

  // Provide a listing of the info.
  void list(ostream& = cout);

private: 

  // Constants: could only be changed in the code itself.
  static const int NSTUDYMIN;
  static const double P2MIN, EIGENVALUEMIN;

  // Properties of analysis.
  double power;
  int select, powerInt; 
  double powerMod;

  // Outcome of analysis.
  double eVal1, eVal2, eVal3; 
  Vec4 eVec1, eVec2, eVec3; 

};  

//**************************************************************************

// SingleCell class.
// Simple helper class to CellJet for a cell and its contents. 

class SingleCell {

public:

  // Constructor.
  SingleCell(int iCellIn = 0, double etaCellIn = 0., double phiCellIn = 0., 
    double eTcellIn = 0., int multiplicityIn = 0) : iCell(iCellIn), 
    etaCell(etaCellIn), phiCell(phiCellIn), eTcell(eTcellIn), 
    multiplicity(multiplicityIn), canBeSeed(true), isUsed(false),
    isAssigned(false) {}

  // Properties of cell.
  int iCell;
  double etaCell, phiCell, eTcell;
  int multiplicity;
  bool canBeSeed, isUsed, isAssigned;

} ;

//**************************************************************************

// SingleCellJet class.
// Simple helper class to CellJet for a jet and its contents. 

class SingleCellJet {

public:

  // Constructor.
  SingleCellJet(double eTjetIn = 0., double etaCenterIn = 0., 
    double phiCenterIn = 0., double etaWeightedIn = 0.,
    double phiWeightedIn = 0., int multiplicityIn = 0,
    Vec4 pMassiveIn = 0.) : eTjet(eTjetIn), etaCenter(etaCenterIn), 
    phiCenter(phiCenterIn), etaWeighted(etaWeightedIn), 
    phiWeighted(phiWeightedIn), multiplicity(multiplicityIn),
    pMassive(pMassiveIn) {}

  // Properties of jet.
  double eTjet, etaCenter, phiCenter, etaWeighted, phiWeighted;
  int multiplicity;
  Vec4 pMassive;  

} ;

//**************************************************************************

// CellJet class.
// This class performs a cone jet search in (eta, phi, E_T) space.

class CellJet {

public: 

  // Constructor.
  CellJet(double eTjetMinIn = 20., double coneRadiusIn = 0.7, 
    int selectIn = 2, double etaMaxIn = 5., int nEtaIn = 50,
    int nPhiIn = 32, double eTseedIn = 1.5, int smearIn = 0,
    double resolutionIn = 0.5, double upperCutIn = 2.,
    double thresholdIn = 0.) : eTjetMin(eTjetMinIn), 
    coneRadius(coneRadiusIn), select(selectIn), etaMax(etaMaxIn), 
    nEta(nEtaIn), nPhi(nPhiIn), eTseed(eTseedIn), smear(smearIn),
    resolution(resolutionIn), upperCut(upperCutIn), 
    threshold(thresholdIn) { }
  
  // Analyze event.
  bool analyze(Event& event);

  // Return info on results of analysis.
  int size() const {return jets.size();}
  double eT(int i) const {return jets[i].eTjet;}
  double etaCenter(int i) const {return jets[i].etaCenter;}
  double phiCenter(int i) const {return jets[i].phiCenter;}
  double etaWeighted(int i) const {return jets[i].etaWeighted;}
  double phiWeighted(int i) const {return jets[i].phiWeighted;}
  double multiplicity(int i) const {return jets[i].multiplicity;}
  Vec4 pMassless(int i) const {return jets[i].eTjet * Vec4(
    cos(jets[i].phiWeighted), sin(jets[i].phiWeighted),
    sinh(jets[i].etaWeighted), cosh(jets[i].etaWeighted) );}
  Vec4 pMassive(int i) const {return jets[i].pMassive;}
  double m(int i) const {return jets[i].pMassive.mCalc();}

  // Provide a listing of the info.
  void list(ostream& = cout);

private: 

  // Properties of analysis.
  double eTjetMin, coneRadius; 
  int select; 
  double etaMax; 
  int nEta, nPhi; 
  double eTseed;
  int smear;
  double resolution, upperCut, threshold;

  // Outcome of analysis: ET-ordered list of jets. 
  vector<SingleCellJet> jets;

};  

//**************************************************************************

} // end namespace Pythia8

#endif // end Pythia8_Analysis_H

