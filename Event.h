// Header file for the Particle and Event classes.
// Particle: information on an instance of a particle.
// Junction: information on a junction between three colours.
// Event: list of particles in the current event.
// Copyright © 2005 Torbjörn Sjöstrand

#ifndef Pythia8_Event_H
#define Pythia8_Event_H

#include "Stdlib.h"
#include "Basics.h"
#include "Settings.h"
#include "ParticleData.h"

namespace Pythia8 {

//**************************************************************************

// Particle class.
// This class holds info on a particle in general.

class Particle {

public:

  // Constructors.
  Particle() : idSave(0), statusSave(0), mother1Save(0), mother2Save(0), 
    daughter1Save(0), daughter2Save(0), colSave(0), acolSave(0), 
    pSave(Vec4(0.,0.,0.,0.)), mSave(0.), scaleSave(0.), 
    hasVertexSave(false), vProdSave(Vec4(0.,0.,0.,0.)), tauSave(0.), 
    particlePtr(0) { }  
  Particle(int idIn, int statusIn = 0, int mother1In = 0, int mother2In = 0, 
    int daughter1In = 0, int daughter2In = 0, int colIn = 0, int acolIn = 0, 
    Vec4 pIn = Vec4(0.,0.,0.,0.), double mIn = 0., double scaleIn = 0.) 
    : idSave(idIn), statusSave(statusIn), mother1Save(mother1In), 
    mother2Save(mother2In), daughter1Save(daughter1In), 
    daughter2Save(daughter2In), colSave(colIn), acolSave(acolIn), 
    pSave(pIn), mSave(mIn), scaleSave(scaleIn), hasVertexSave(false), 
    vProdSave(Vec4(0.,0.,0.,0.)), tauSave(0.),
    particlePtr(ParticleDataTable::particleDataPtr(idIn)) { }  
  Particle(int idIn, int statusIn = 0, int mother1In = 0, 
    int mother2In = 0, int daughter1In = 0, int daughter2In = 0,
    int colIn = 0, int acolIn = 0, double pxIn = 0., double pyIn = 0., 
    double pzIn = 0., double eIn = 0., double mIn = 0., double scaleIn = 0.) 
    : idSave(idIn), statusSave(statusIn), mother1Save(mother1In), 
    mother2Save(mother2In), daughter1Save(daughter1In), 
    daughter2Save(daughter2In), colSave(colIn), acolSave(acolIn), 
    pSave(Vec4(pxIn, pyIn, pzIn, eIn)), mSave(mIn), scaleSave(scaleIn), 
    hasVertexSave(false), vProdSave(Vec4(0.,0.,0.,0.)), tauSave(0.), 
    particlePtr(ParticleDataTable::particleDataPtr(idIn)) { }  
  Particle(const Particle& pt) : idSave(pt.idSave), statusSave(pt.statusSave), 
    mother1Save(pt.mother1Save), mother2Save(pt.mother2Save), 
    daughter1Save(pt.daughter1Save), daughter2Save(pt.daughter2Save), 
    colSave(pt.colSave), acolSave(pt.acolSave), pSave(pt.pSave), 
    mSave(pt.mSave), scaleSave(pt.scaleSave), hasVertexSave(pt.hasVertexSave), 
    vProdSave(pt.vProdSave), tauSave(pt.tauSave), particlePtr(pt.particlePtr) 
    { } 
  Particle& operator=(const Particle& pt) {if (this != &pt) {
    idSave = pt.idSave; statusSave = pt.statusSave; 
    mother1Save = pt.mother1Save; mother2Save = pt.mother2Save; 
    daughter1Save = pt.daughter1Save; daughter2Save = pt.daughter2Save; 
    colSave = pt.colSave; acolSave = pt.acolSave; pSave = pt.pSave; 
    mSave = pt.mSave; scaleSave = pt.scaleSave; 
    hasVertexSave = pt.hasVertexSave; vProdSave = pt.vProdSave; 
    tauSave = pt.tauSave; particlePtr = pt.particlePtr; } return *this; } 
      
  // Member functions for input.
  void id(int idIn) {idSave = idIn; 
    particlePtr = ParticleDataTable::particleDataPtr(idIn);}
  void status(int statusIn) {statusSave = statusIn;}
  void statusPos() {statusSave = abs(statusSave);}
  void statusNeg() {statusSave = -abs(statusSave);}
  void statusCode(int statusIn) {statusSave = 
    (statusSave > 0) ? abs(statusIn) : -abs(statusIn);}
  void mother1(int mother1In) {mother1Save = mother1In;}
  void mother2(int mother2In) {mother2Save = mother2In;}
  void mothers(int mother1In = 0, int mother2In = 0) 
    {mother1Save = mother1In; mother2Save = mother2In;}
  void daughter1(int daughter1In) {daughter1Save = daughter1In;}
  void daughter2(int daughter2In) {daughter2Save = daughter2In;}
  void daughters(int daughter1In = 0, int daughter2In = 0) 
    {daughter1Save = daughter1In; daughter2Save = daughter2In;}  
  void col(int colIn) {colSave = colIn;}
  void acol(int acolIn) {acolSave = acolIn;}
  void cols(int colIn = 0,int acolIn = 0) {colSave = colIn; 
    acolSave = acolIn;}  
  void p(Vec4 pIn) {pSave = pIn;}
  void p(double pxIn, double pyIn, double pzIn, double eIn) 
    {pSave.p(pxIn, pyIn, pzIn, eIn);}
  void px(double pxIn) {pSave.px(pxIn);}
  void py(double pyIn) {pSave.py(pyIn);}
  void pz(double pzIn) {pSave.pz(pzIn);}
  void e(double eIn) {pSave.e(eIn);}
  void m(double mIn) {mSave = mIn;}
  void scale(double scaleIn) {scaleSave = scaleIn;}
  void vProd(Vec4 vProdIn) {vProdSave = vProdIn; hasVertexSave = true;}
  void vProd(double xProdIn, double yProdIn, double zProdIn, double tProdIn)
    {vProdSave.p(xProdIn, yProdIn, zProdIn, tProdIn); hasVertexSave = true;}
  void xProd(double xProdIn) {vProdSave.px(xProdIn); hasVertexSave = true;} 
  void yProd(double yProdIn) {vProdSave.py(yProdIn); hasVertexSave = true;} 
  void zProd(double zProdIn) {vProdSave.pz(zProdIn); hasVertexSave = true;} 
  void tProd(double tProdIn) {vProdSave.e(tProdIn); hasVertexSave = true;} 
  void tau(double tauIn) {tauSave = tauIn;} 

  // Member functions for output.
  int id() const {return idSave;}
  int status() const {return statusSave;}
  int mother1() const {return mother1Save;}
  int mother2() const {return mother2Save;}
  int daughter1() const {return daughter1Save;}
  int daughter2() const {return daughter2Save;}
  int col() const {return colSave;}
  int acol() const {return acolSave;}
  Vec4 p() const {return pSave;}
  double px() const {return pSave.px();}
  double py() const {return pSave.py();}
  double pz() const {return pSave.pz();}
  double e() const {return pSave.e();}
  double m() const {return mSave;}
  double scale() const {return scaleSave;}
  bool hasVertex() const {return hasVertexSave;}
  Vec4 vProd() const {return vProdSave;}
  double xProd() const {return vProdSave.px();}
  double yProd() const {return vProdSave.py();}
  double zProd() const {return vProdSave.pz();}
  double tProd() const {return vProdSave.e();}
  double tau() const {return tauSave;}

  // Member functions for output; derived int and bool quantities.
  int statusAbs() const {return abs(statusSave);}
  bool remains() const {return (statusSave > 0) ? true : false;}
  bool isQ() const {return (abs(idSave) < 9 && idSave != 0) ? true : false;}
  bool isNotQ() const {return (abs(idSave) > 9 || idSave == 0) ? true : false;}
  bool isG() const {return (idSave == 21) ? true : false;}
  bool isNotG() const {return (idSave == 21) ? false : true;}
  bool isQorG() const {return ( isQ() || isG() ) ? true : false;}
  bool isQQ() const {int idAbs = abs(idSave); return ( idAbs > 1000 
    && idAbs < 10000 && (idAbs/10)%10 == 0) ? true : false;}
  bool isQorQQ() const {return ( isQ() || isQQ() ) ? true : false;} 
  bool isL() const {return (abs(idSave) > 10 && abs(idSave) < 19) ? true : false;}
  bool hasCol() const {return (colSave > 0 || acolSave > 0) ? true : false;}

  // Member functions for output; derived double quantities.
  double m2() const {return mSave*mSave;}
  double mCalc() const {return pSave.mCalc();}
  double m2Calc() const {return pSave.m2Calc();}
  double eCalc() const {return sqrt(mSave*mSave + pSave.pAbs2());}
  double pT() const {return pSave.pT();}
  double pT2() const {return pSave.pT2();}
  double mT() const {return sqrt(mSave*mSave + pSave.pT2());}
  double mT2() const {return mSave*mSave + pSave.pT2();}
  double pAbs() const {return pSave.pAbs();}
  double pAbs2() const {return pSave.pAbs2();}
  double theta() const {return pSave.theta();}
  double phi() const {return pSave.phi();}
  double thetaXZ() const {return pSave.thetaXZ();}
  double pPlus() const {return pSave.pPlus();}
  double pMinus() const {return pSave.pMinus();}
  double y() const;
  double eta() const; 
  Vec4 vDec() const {return vProdSave + tauSave * pSave / mSave;}
  double xDec() const {return vProdSave.px() + tauSave * pSave.px() / mSave;}
  double yDec() const {return vProdSave.py() + tauSave * pSave.py() / mSave;}
  double zDec() const {return vProdSave.pz() + tauSave * pSave.pz() / mSave;}
  double tDec() const {return vProdSave.e() + tauSave * pSave.e() / mSave;}

  // Further output, based on a pointer to a ParticleDataEntry object.
  string name() const {return particlePtr->name(idSave);}
  string nameWithStatus() const {
    if (statusSave > 0) return particlePtr->name(idSave);
    return "(" + particlePtr->name(idSave) + ")"; }
  double m0() const {return particlePtr->m0();}
  double mass() const {return particlePtr->mass();}
  double constituentMass() const {return particlePtr->constituentMass();}
  double tau0() const {return particlePtr->tau0();}
  int colType() const {return particlePtr->colType(idSave);}
  double charge() const {return  particlePtr->charge(idSave);}
  int icharge() const {return particlePtr->charge3(idSave);}
  bool isCharged() const {return (particlePtr->charge3(idSave) == 0) 
    ? false : true;}
  bool isNeutral() const {return (particlePtr->charge3(idSave) == 0) 
    ? true : false;}
  int spinType() const {return particlePtr->spinType();}
  bool canDecay() const {return particlePtr->canDecay();}
  bool mayDecay() const {return particlePtr->mayDecay();}
  ParticleDataEntry& particleData() const {return *particlePtr;}

  // Member functions that perform operations.
  void rescale3(double fac) {pSave.rescale3(fac);}
  void rescale4(double fac) {pSave.rescale4(fac);}
  void rescale5(double fac) {pSave.rescale4(fac); mSave *= fac;}
  void rot(double theta, double phi) {pSave.rot(theta, phi);
    if (hasVertexSave) vProdSave.rot(theta, phi);} 
  void bst(double betaX, double betaY, double betaZ) {
    pSave.bst(betaX, betaY, betaZ);
    if (hasVertexSave) vProdSave.bst(betaX, betaY, betaZ);}
  void bst(double betaX, double betaY, double betaZ, double gamma) {
    pSave.bst(betaX, betaY, betaZ, gamma);
    if (hasVertexSave) vProdSave.bst(betaX, betaY, betaZ, gamma);}
  void bst(const Vec4& vec) {pSave.bst(vec);
    if (hasVertexSave) vProdSave.bst(vec);}
  void rotbst(const RotBstMatrix& M) {pSave.rotbst(M);
    if (hasVertexSave) vProdSave.rotbst(M);} 

  // Print a particle
  friend ostream& operator<<(ostream&, const Particle&) ;

private:

  // Constants: could only be changed in the code itself.
  static const double TINY;

  // Properties of the current particle.
  int idSave, statusSave, mother1Save, mother2Save, daughter1Save, 
    daughter2Save, colSave, acolSave;
  Vec4 pSave;
  double mSave, scaleSave;
  bool hasVertexSave;
  Vec4 vProdSave;
  double tauSave;

  // Pointer to properties of the particle species.
  ParticleDataEntry* particlePtr;

};

// Invariant mass of a pair and its square.
// (Not part of class proper, but tightly linked.)
double m(const Particle&, const Particle&); 
double m2(const Particle&, const Particle&); 

//**************************************************************************

// The juction class stores what kind of junction it is, the colour indices 
// of the legs at the junction and as far out as legs have been traced,
// and the status codes assigned for fragmentation of each leg.

class Junction {

public:

  // Constructors.
  Junction() : remainsSave(true), kindSave(0) { 
    for (int j = 0; j < 3; ++j) {
    colSave[j] = 0; endColSave[j] = 0; statusSave[j] = 0; } }
  Junction( int kindIn, int col0In, int col1In, int col2In) 
    : remainsSave(true), kindSave(kindIn) {colSave[0] = col0In; 
    colSave[1] = col1In; colSave[2] = col2In; 
    for (int j = 0; j < 3; ++j) {
    endColSave[j] = colSave[j]; statusSave[j] = 0; } }
  Junction(const Junction& ju) : remainsSave(ju.remainsSave), 
    kindSave(ju.kindSave) { for (int j = 0; j < 3; ++j) {
    colSave[j] = ju.colSave[j]; endColSave[j] = ju.endColSave[j]; 
    statusSave[j] = ju.statusSave[j]; } }
  Junction& operator=(const Junction& ju) {if (this != &ju) { 
    remainsSave = ju.remainsSave; kindSave =  ju.kindSave; 
    for (int j = 0; j < 3; ++j) { colSave[j] = ju.colSave[j]; 
    endColSave[j] = ju.endColSave[j]; statusSave[j] = ju.statusSave[j]; } } 
    return *this; }  

  // Set values.
  void remains(bool remainsIn) {remainsSave = remainsIn;}
  void col(int j, int colIn) {colSave[j] = colIn; endColSave[j] = colIn;}
  void endCol(int j, int endColIn) {endColSave[j] = endColIn;}
  void status(int j, int statusIn) {statusSave[j] = statusIn;}

  // Read out value.
  bool remains() const {return remainsSave;}
  int kind() const {return kindSave;}
  int col(int j) const {return colSave[j];}
  int endCol(int j) const {return endColSave[j];}
  int status(int j) const {return statusSave[j];}
 
private:

  // Kind, positions of the three ends and their status codes.
  bool remainsSave;
  int kindSave, colSave[3], endColSave[3], statusSave[3];

};

//**************************************************************************

// The Event class holds all info on the generated event.

class Event {
    
public:

  // Constructor.
  Event(int capacity = 100) {entry.reserve(capacity);
    headerList = "----------------------------------------";}

  // Initialize static data members.
  static void initStatic();

  // Overload index operator to access element of event record.
  Particle& operator[](int i) {return entry[i];}
  const Particle& operator[](int i) const {return entry[i];}

  // Event record size.
  int size() const {return entry.size();}

  // Put a new particle at the end of the event record; return index.
  int append(Particle entryIn) {    
    entry.push_back(entryIn); 
    if (entryIn.col() > maxColTag) maxColTag = entryIn.col();   
    if (entryIn.acol() > maxColTag) maxColTag = entryIn.acol();
    return entry.size() - 1;
  }
  int append(int id, int status = 0, int mother1 = 0, int mother2 = 0, 
    int daughter1 = 0, int daughter2 = 0, int col = 0, int acol = 0, 
    Vec4 p = Vec4(0.,0.,0.,0.), double m = 0., double scale = 0.) {    
    entry.push_back( Particle(id, status, mother1, mother2, daughter1,
    daughter2, col, acol, p, m, scale) ); 
    if (col > maxColTag) maxColTag = col;   
    if (acol > maxColTag) maxColTag = acol;
    return entry.size() - 1;
  }

  // Add a copy of an existing particle at the end of the event record.
  int copy(int iCopy, int newStatus = 0);

  // Implement reference "back" to access last element.
  Particle& back() {return entry.back();}

  // List the particles in an event.
  void list(ostream& = cout);  

  // Set header specification for event listing.
  void header( string headerIn) {
    headerList.replace(0, headerIn.length() + 2, headerIn + "  ");}

  // Clear event record, or remove last n entries.
  void clear() {entry.resize(0); maxColTag = startColTag; junction.resize(0);}
  void popBack(int nRemove = 1) { if (nRemove ==1) entry.pop_back();
    else {int newSize = max( 0, size() - nRemove); entry.resize(newSize);} } 

  // Save or restore the size of the event record (throwing at the end).
  void saveSize() {savedSize = entry.size();}
  void restoreSize() {entry.resize(savedSize);}   

  // Initialize and access colour tag information.
  void initColTag(int colTag = 0) {maxColTag = max( colTag,startColTag);}
  int lastColTag() const {return maxColTag;}
  int nextColTag() {return ++maxColTag;}

  // Access scale for which event as a whole is defined.
  void scale( double scaleIn) {scaleSave = scaleIn;}
  double scale() const {return scaleSave;}

  // Find complete list of daughters and mothers.
  vector<int> motherList(int i);
  vector<int> daughterList(int i);
 
  // Trace the first and last copy of one and the same particle.
  int iTopCopy(int i);
  int iBotCopy(int i);

  // Trace the first and last copy of a particle, using flavour match.
  int iTopCopyId(int i);
  int iBotCopyId(int i);

  // Find list of sisters, also tracking up and down identical copies.
  vector<int> sisterList(int i);
  vector<int> sisterListTopBot(int i);

  // Check whether two particles have a direct mother-daughter relation.
  bool isAncestor(int i, int iAncestor);

  // Add a junction to the list, study it or extra input.
  void appendJunction( int kind, int col0, int col1, int col2)  
    { junction.push_back( Junction( kind, col0, col1, col2) );} 
  void appendJunction(Junction junctionIn) {junction.push_back(junctionIn);} 
  int sizeJunction() const {return junction.size();}
  bool remainsJunction(int i) const {return junction[i].remains();}
  void remainsJunction(int i, bool remainsIn) {junction[i].remains(remainsIn);}
  int kindJunction(int i) const {return junction[i].kind();}
  int colJunction( int i, int j) const {return junction[i].col(j);}
  void colJunction( int i, int j, int colIn) {junction[i].col(j, colIn);}
  int endColJunction( int i, int j) const {return junction[i].endCol(j);}
  void endColJunction( int i, int j, int endColIn) 
    {junction[i].endCol(j, endColIn);}
  int statusJunction( int i, int j) const {return junction[i].status(j);}
  void statusJunction( int i, int j, int statusIn) 
    {junction[i].status(j, statusIn);}
  Junction& getJunction(int i) {return junction[i];}
  void eraseJunction(int i);

  // Save or restore the size of the junction list (throwing at the end).
  void saveJunctionSize() {savedJunctionSize = junction.size();}
  void restoreJunctionSize() {junction.resize(savedJunctionSize);}   

  // Member functions for rotations and boosts of an event.
  void rot(double theta, double phi) 
    {for (int i = 0; i < size(); ++i) entry[i].rot(theta, phi);} 
  void bst(double betaX, double betaY, double betaZ) 
    {for (int i = 0; i < size(); ++i) entry[i].bst(betaX, betaY, betaZ);}
  void bst(double betaX, double betaY, double betaZ, double gamma) 
    {for (int i = 0; i < size(); ++i) entry[i].bst(betaX, betaY, betaZ, 
    gamma);}
  void bst(const Vec4& vec) 
    {for (int i = 0; i < size(); ++i) entry[i].bst(vec);}
  void rotbst(const RotBstMatrix& M) 
    {for (int i = 0; i < size(); ++i) entry[i].rotbst(M);}
 
private: 

  // Static initialization data, normally only set once.
  static int startColTag;
  static bool listFinalOnly, listScaleAndVertex, listMothersAndDaughters, 
    extraBlankLine, listJunctions; 

  // Constants: could only be changed in the code itself.
  static const int IPERLINE;

  // The event: a vector containing all particles (entries).
  vector<Particle> entry;

  // The list of junctions.
  vector<Junction> junction;

  // The maximum colour tag of the event so far.
  int maxColTag;

  // Saved entry and junction list sizes, for simple restoration.
  int savedSize, savedJunctionSize;

  // The scale of the event; linear quantity in GeV.
  double scaleSave;

  // Header specification in event listing (at most 40 characters wide).
  string headerList;
  
};

//**************************************************************************

} // end namespace Pythia8

#endif // end Pythia8_Event_H
