// This file contains auxiliary classes in the fragmentation process.
// ColSinglet contains info on an individual singlet.
// ColConfig describes the colour configuration of the whole event.
// StringRegion keeps track on string momenta and directions.
// StringSystem contains all the StringRegions of the colour singlet.
// Copyright © 2005 Torbjörn Sjöstrand

#ifndef Pythia8_FragmentationSystems_H
#define Pythia8_FragmentationSystems_H

#include "Stdlib.h"
#include "Basics.h"
#include "Settings.h"
#include "ParticleData.h"
#include "Event.h"
#include "FragmentationFlavZpT.h"

namespace Pythia8 {
 
//**************************************************************************

// The ColSinglet class contains info on an individual singlet.
// Only to be used inside ColConfig, so no private members. 

class ColSinglet {
  
public:

  // Constructors.
  ColSinglet() : iBegin(0), iEnd(0), pSum(0., 0., 0.,0.), mass(0.),
    massExcess(0.), hasJunction(false), isClosed(false), 
    isCollected(false) {}
  ColSinglet(int iBeginIn, int iEndIn, Vec4 pSumIn, double massIn, 
    double massExcessIn, bool hasJunctionIn = false,
    bool isClosedIn = false, bool isCollectedIn = false) 
    : iBegin(iBeginIn), iEnd(iEndIn), pSum(pSumIn), mass(massIn), 
    massExcess(massExcessIn), hasJunction(hasJunctionIn),
    isClosed(isClosedIn), isCollected(isCollectedIn) {}

  // Methods for information extraction.
  int size() const {return iEnd - iBegin;}

  // Stored quantities. iBegin and iEnd refer to iParton(Save) list.
  // (Workaround for segmentation faults for vector inside vector.)
  int iBegin, iEnd;
  Vec4 pSum;
  double mass, massExcess;
  bool hasJunction, isClosed, isCollected;
  
};
 
//**************************************************************************

// The ColConfig class describes the colour configuration of the whole event. 
// (Keep track of which belongs together, e.g. from same Z0??)

class ColConfig {

public:

  // Constructor.
  ColConfig() {}

  // Initialize static data members.
  static void initStatic();

  // Number of colour singlets.
  int size() const {return singlets.size();}

  // Overload index operator to access separate colour singlets.
  ColSinglet& operator[](int iSub) {return singlets[iSub];}

  // Partons inside a singlet.
  int iPartonSize(int iSub) const {return singlets[iSub].size();}
  int& iParton(int iSub, int i) 
    {return iPartonSave[ singlets[iSub].iBegin + i];}  
  int front(int iSub) const {return iPartonSave[ singlets[iSub].iBegin];}  
  int back(int iSub) const {return iPartonSave[ singlets[iSub].iEnd - 1];}  

  // Clear contents.
  void clear() {singlets.resize(0); iPartonSave.resize(0);} 

  // Insert a new colour singlet system in ascending mass order. 
  // Calculate its properties. Join nearby partons.
  void insert( vector<int>& iPartonIn, Event& event); 

  // Collect all partons of singlet to be consecutively ordered.
  void collect(int iSub, Event& event); 

  // List all currently identified singlets.
  void list(ostream& = cout);

private:

  // Static initialization data, normally only set once.
  static double mJoin, mJoinJunction, mStringMin;
 
  // List of all separate colour singlets.
  vector<ColSinglet> singlets;

  // Colour-ordered list of all parton positions.
  vector<int> iPartonSave;

  // Join two legs of junction to a diquark for small invariant masses.
  bool joinJunction( vector<int>& iPartonIn, Event& event, 
    double massExcessIn); 

};
 
//**************************************************************************

// The StringRegion class contains the information related to 
// one string section in the evolution of a multiparton system. 
// Only to be used inside StringFragmentation and MiniStringFragmentation,
// so no private members.

// Currently a number of simplifications, in particular ??
// 1) No popcorn baryon production.
// 2) Simplified treatment of pT in stepping and joining.

class StringRegion {

public:

  // Constructor. 
  StringRegion() : isSetUp(false), isEmpty(true) {}

  // Initialize static data members.
  static void initStatic();

  // Initialization data, normally only set once.
  static double mJoin, m2Join;

  // Constants: could only be changed in the code itself.
  static const double TINY;

  // Data members.
  bool isSetUp, isEmpty;
  Vec4 pPos, pNeg, eX, eY;
  double w2;
  double xPosProj, xNegProj, pxProj, pyProj;

  // Set up four-vectors for longitudinal and transverse directions.
  void setUp(Vec4 p1, Vec4 p2, bool isMassless = false);

  // Construct a four-momentum from (x+, x-, px, py).
  Vec4 pHad( double xPos, double xNeg, double px, double py) 
    { return xPos * pPos + xNeg * pNeg + px * eX + py * eY; }

  // Project a four-momentum onto (x+, x-, px, py). Read out projection.
  void project(Vec4 pIn);
  void project( double px, double py, double pz, double e) 
    { project( Vec4( px, py, pz, e) ); }
  double xPos() const {return xPosProj;} 
  double xNeg() const {return xNegProj;} 
  double px() const {return pxProj;} 
  double py() const {return pyProj;} 

};
 
//**************************************************************************

// The StringSystem class contains the complete set of all string regions.
// Only to be used inside StringFragmentation, so no private members.

class StringSystem {

public:

  // Constructor. 
  StringSystem() {}

  // Set up system from parton list.
  void setUp(vector<int>& iSys, Event& event);

  // Calculate string region from (iPos, iNeg) pair.
  int iReg( int iPos, int iNeg) const 
    {return (iPos * (indxReg - iPos)) / 2 + iNeg;}

  // Reference to string region specified by (iPos, iNeg) pair.
  StringRegion& region(int iPos, int iNeg) {return system[iReg(iPos, iNeg)];} 

  // Reference to low string region specified either by iPos or iNeg.
  StringRegion& regionLowPos(int iPos) {return system[iReg(iPos, iMax - iPos)];} 
  StringRegion& regionLowNeg(int iNeg) {return system[iReg(iMax - iNeg, iNeg)];} 

  // Main content: a vector with all the string regions of the system. 
  vector<StringRegion> system;

  // Other data members.
  int sizePartons, sizeStrings, sizeRegions, indxReg, iMax; 

};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_FragmentationSystems_H
