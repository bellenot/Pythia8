// Header file for the classes containing particle data.
// DecayChannel contains info on a single decay channel.
// DecayTable contains all decay channels of a particle.
// ParticleDataEntry contains info on a single particle species.
// ParticleDataTable  collects info on all particles as a map.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_ParticleData_H
#define Pythia8_ParticleData_H

#include "Basics.h"
#include "PythiaStdlib.h"
#include "Settings.h"

namespace Pythia8 {

//**************************************************************************

// This class holds info on a single decay channel.

class DecayChannel {

public:
  // Constructors.
  DecayChannel(double bratIn = 0., int modeIn = 0, int prod0 = 0, 
    int prod1 = 0, int prod2 = 0, int prod3 = 0, int prod4 = 0, 
    int prod5 = 0, int prod6 = 0, int prod7 = 0) : brat(bratIn), 
    mode(modeIn), nProd(0), hasChangedSave(true) {prod[0] = prod0; 
    prod[1] = prod1; prod[2] = prod2; prod[3] = prod3; prod[4] = prod4; 
    prod[5] = prod5; prod[6] = prod6; prod[7] = prod7; 
    for (int i = 0; i < 8; ++i) if (prod[i] != 0 && i == nProd) ++nProd; }  
  DecayChannel(const DecayChannel& chan) : brat(chan.brat), mode(chan.mode),
    nProd(chan.nProd), hasChangedSave(true) 
    { for (int i = 0; i < 8; ++i) prod[i] = chan.prod[i]; } 
  DecayChannel& operator=(const DecayChannel& chan) { if (this != &chan)
    { brat = chan.brat; mode = chan.mode; nProd = chan.nProd; 
    for (int i = 0; i < 8; ++i) prod[i] = chan.prod[i]; 
    hasChangedSave =true;} return *this; }  

  // Member functions for input.
  void branchingRatio(double bratIn) {brat = bratIn; hasChangedSave = true;}
  void rescaleBR(double fac) {brat *= fac; hasChangedSave = true;} 
  void modeME(int modeIn) {mode = modeIn; hasChangedSave = true;} 
  void multiplicity(int multIn)  {nProd = multIn; hasChangedSave = true;} 
  void product(int i, int prodIn) {prod[i] = prodIn; hasChangedSave = true;}
  void setHasChanged(bool hasChangedIn) {hasChangedSave = hasChangedIn;}

  // Member functions for output.
  double branchingRatio() const {return brat;}
  int modeME() const {return mode;}
  int multiplicity() const {return nProd;} 
  int product(int i) const {return (i >= 0 && i < nProd) ? prod[i] : 0;} 
  bool hasChanged() const { return hasChangedSave;}

private:

  // Decay channel info.
  double brat;
  int mode, nProd, prod[8];
  bool hasChangedSave;

};

//**************************************************************************

// This class holds info on all decay channels of a particle.

class DecayTable {

public:

  // Constructor.
  DecayTable() {}

  // Overload index operator to access a channel in the decay table.
  DecayChannel& operator[](int i){return channel[i];}
  const DecayChannel& operator[](int i) const {return channel[i];}

  // Add a decay channel to the decay table.
  void addChannel(double brat = 0., int mode = 0, int prod0 = 0, 
    int prod1 = 0, int prod2 = 0, int prod3 = 0, int prod4 = 0, 
    int prod5 = 0, int prod6 = 0, int prod7 = 0) { 
    channel.push_back( DecayChannel( brat, mode, prod0, prod1, prod2, 
    prod3, prod4, prod5, prod6, prod7) ); }

  // Decay table size.
  int size() const {return channel.size();}

  // Random choice of decay mode according to branching ratios.
  DecayChannel& pick();

  // Rescale sum of branching ratios to unity.
  void rescaleBR(double newSumBR = 1.);

private:

  // A vector containing all the decay channels of the particle.
  vector<DecayChannel> channel;

};

//**************************************************************************

// This class holds info on a single particle species.

class ParticleDataEntry {

public:

  // Constructors: for antiparticle exists or not.
  ParticleDataEntry(int idIn = 0, string nameIn = " ", int charge3In = 0, 
    int colTypeIn = 0, double m0In = 0., double widthIn = 0., 
    double rangeIn = 0., double tau0In = 0., bool mayDecayIn = false,
    bool isResonanceIn = false) 
    : idAbs(abs(idIn)), nameSave(nameIn), antiNameSave("void"),  
    hasAntiSave(false), mayDecaySave(mayDecayIn), 
    isResonanceSave(isResonanceIn), externalDecaySave(false),
    charge3Save(charge3In), colTypeSave(colTypeIn), m0Save(m0In), 
    widthSave (widthIn), rangeSave(rangeIn), tau0Save(tau0In), 
    hasChangedSave(true) { constituentMassCalc(); }   
  ParticleDataEntry(int idIn, string nameIn, string antiNameIn, 
    int charge3In = 0, int colTypeIn = 0, double m0In = 0.,  
    double widthIn = 0., double rangeIn = 0., double tau0In = 0., 
    bool mayDecayIn = false, bool isResonanceIn = false) 
    : idAbs(abs(idIn)), nameSave(nameIn), antiNameSave(antiNameIn), 
    hasAntiSave(true), mayDecaySave(mayDecayIn), 
    isResonanceSave(isResonanceIn), externalDecaySave(false), 
    charge3Save(charge3In), colTypeSave(colTypeIn), m0Save(m0In), 
    widthSave (widthIn), rangeSave(rangeIn), tau0Save(tau0In), 
    hasChangedSave(true) { constituentMassCalc(); 
    if (tolower(antiNameIn) == "void") { hasAntiSave = false; } }

  // Initialize static data members.
  static void initStatic();

  // Change current values (or set if not set before).
  void setName(string nameIn) {nameSave = nameIn; hasChangedSave = true;}
  void setAntiName(string antiNameIn) {antiNameSave = antiNameIn; 
    hasChangedSave = true;}
  void setNames(string nameIn, string antiNameIn) {nameSave = nameIn; 
    antiNameSave = antiNameIn; hasAntiSave = true; if (tolower(antiNameIn) 
    == "void") hasAntiSave = false; hasChangedSave = true;}
  void setCharge3(int charge3In) {charge3Save = charge3In; 
    hasChangedSave = true;}
  void setColType(int colTypeIn) {colTypeSave = colTypeIn; 
    hasChangedSave = true;}
  void setM0(double m0In) {m0Save = m0In; constituentMassCalc(); 
    hasChangedSave = true;}
  void setWidth(double widthIn) {widthSave = widthIn; hasChangedSave = true;}
  void setRange(double rangeIn) {rangeSave = rangeIn; hasChangedSave = true;}
  void setTau0(double tau0In) {tau0Save = tau0In; hasChangedSave = true;}
  void setMayDecay(bool mayDecayIn) {mayDecaySave = mayDecayIn; 
    hasChangedSave = true;}
  void setIsResonance(bool isResonanceIn) {isResonanceSave = isResonanceIn; 
    hasChangedSave = true;}
  void setAll(string nameIn, string antiNameIn, int charge3In = 0, 
    int colTypeIn = 0, double m0In = 0., double widthIn = 0., 
    double rangeIn = 0., double tau0In = 0., bool mayDecayIn = false,
    bool isResonanceIn = false) 
    {nameSave = nameIn; antiNameSave = antiNameIn; hasAntiSave = true; 
    if (tolower(antiNameIn) == "void") hasAntiSave = false;
    charge3Save = charge3In; colTypeSave = colTypeIn; m0Save = m0In; 
    constituentMassCalc(); widthSave = widthIn; rangeSave = rangeIn;
    tau0Save = tau0In; mayDecaySave = mayDecayIn; 
    isResonanceSave = isResonanceIn; hasChangedSave = true;}
  void setExternalDecay(bool externalDecayIn) 
    {externalDecaySave = externalDecayIn;}
  void setHasChanged(bool hasChangedIn) {hasChangedSave = hasChangedIn;}
  void setHasChangedAll(bool hasChangedIn) {hasChangedSave = hasChangedIn;
    for (int i = 0; i < decay.size(); ++i) 
    decay[i].setHasChanged(hasChangedIn);}
  void rescaleBR(double newSumBR = 1.) {decay.rescaleBR(newSumBR);}

  // Give back current values. 
  int id() const { return idAbs; }
  bool hasAnti() const { return hasAntiSave; } 
  string name(int idIn = 1) const { 
    return (idIn > 0) ? nameSave : antiNameSave; } 
  int charge3(int idIn = 1) const { 
    return (idIn > 0) ? charge3Save : -charge3Save; } 
  double charge(int idIn = 1) const { 
    return (idIn > 0) ? charge3Save / 3. : -charge3Save / 3.; } 
  int colType(int idIn = 1) const { if (colTypeSave == 2) return colTypeSave;
    return (idIn > 0) ? colTypeSave : -colTypeSave; } 
  double m0() const { return m0Save; } 
  double constituentMass() const { return constituentMassSave; } 
  double width() const { return widthSave; } 
  double range() const { return rangeSave; } 
  double tau0() const { return tau0Save; } 
  bool isVisible() const;
  bool isInvisible() const { return !isVisible();}
  int spinType() const;
  bool mayDecay() const { return mayDecaySave; } 
  bool isResonance() const { return isResonanceSave; } 
  bool externalDecay() const { return externalDecaySave; } 
  bool hasChanged() const { return hasChangedSave;}
  bool hasChangedAny() const { if (hasChangedSave) return true;
    for (int i = 0; i < decay.size(); ++i) 
    if (decay[i].hasChanged()) return true; return false;}

  // Derived quantity.
  bool canDecay() const { return (decay.size() > 0) ? true : false;} 

  // Calculation of mass, picked according to Breit-Wigner.
  double mass(); 

  // Calculation of constituent masses, hardcoded in a special routine.
  void constituentMassCalc();

  // The decay table.
  DecayTable decay;

private:

  // Static initialization data, normally only set once.
  static int modeBreitWigner;

  // Constants: could only be changed in the code itself.
  static const double NARROWMASS;

  // Particle data.
  int idAbs;
  string nameSave, antiNameSave;
  bool hasAntiSave, mayDecaySave, isResonanceSave, externalDecaySave;
  int charge3Save, colTypeSave;
  double m0Save, constituentMassSave, widthSave, rangeSave, tau0Save;
  bool hasChangedSave;

};

//**************************************************************************

// This class holds a map of all ParticleDataEntries.

class ParticleDataTable {

public:

  // Constructor.
  ParticleDataTable() {}
 
  // Read in database from specific file.
  static bool init(string startFile = "../doc/ParticleData.xml") ;

  // Overwrite existing database by reading from specific file.
  static bool reInit(string startFile) ;

  // Read in one update from a single line.
  static bool readString(string lineIn, bool warn = true) ; 
 
  // Read in updates from user-defined file.
  static bool readFile(string updateFile, bool warn = true) ;

  // Print out table of whole database, or of only part of it.
  static void listAll(ostream& os = cout) {list(false, os);} 
  static void listChanged(ostream& os = cout) {list(true, os);} 
  static void list(bool changedOnly = false, ostream& os = cout) ; 
  static void list(int idList, ostream& os = cout) {
    vector<int> idListTemp; idListTemp.push_back(idList); 
    list( idListTemp, os);} 
  static void list(vector<int> idList, ostream& os = cout) ; 

  // Query existence of an entry.
  static bool isParticle(int idIn) {
    if (pdt.find(abs(idIn)) == pdt.end()) return false;
    if (idIn > 0 || pdt[abs(idIn)].hasAnti()) return true;
    return false; }
 
  // Add new entry.
  static void addParticle(int idIn, string nameIn = " ", int charge3In = 0, 
    int colTypeIn = 0, double m0In = 0., double widthIn = 0., 
    double rangeIn = 0., double tau0In = 0., bool mayDecayIn = false,
    bool isResonanceIn = false) 
    { pdt[abs(idIn)] = ParticleDataEntry(idIn, nameIn, charge3In, 
    colTypeIn, m0In, widthIn, rangeIn, tau0In, mayDecayIn, isResonanceIn); }  
  static void addParticle(int idIn, string nameIn, string antiNameIn, 
    int charge3In = 0, int colTypeIn = 0, double m0In = 0.,
    double widthIn = 0., double rangeIn = 0., double tau0In = 0., 
    bool mayDecayIn = false, bool isResonanceIn = false) 
    { pdt[abs(idIn)] = ParticleDataEntry(idIn, nameIn, antiNameIn, charge3In, 
    colTypeIn, m0In, widthIn, rangeIn, tau0In, mayDecayIn, isResonanceIn); }  
  
  // Return pointer to entry.
  static ParticleDataEntry* particleDataPtr(int idIn) {
    return (isParticle(idIn)) ? &pdt[abs(idIn)] : &pdt[0]; }

  // Change current values (or set if not set before).
  static void name(int idIn, string nameIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setName(nameIn); }
  static void antiName(int idIn, string antiNameIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setAntiName(antiNameIn); }
  static void names(int idIn, string nameIn, string antiNameIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setNames(nameIn, antiNameIn); }
  static void charge3(int idIn, int charge3In) {
    if (isParticle(idIn)) pdt[abs(idIn)].setCharge3(charge3In); }
  static void colType(int idIn, int colTypeIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setColType(colTypeIn); }
  static void m0(int idIn, int m0In) {
    if (isParticle(idIn)) pdt[abs(idIn)].setM0(m0In); }
  static void width(int idIn, int widthIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setWidth(widthIn); }
  static void range(int idIn, int rangeIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setRange(rangeIn); }
  static void tau0(int idIn, int tau0In) {
    if (isParticle(idIn)) pdt[abs(idIn)].setTau0(tau0In); }
  static void mayDecay(int idIn, bool mayDecayIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setMayDecay(mayDecayIn); }
  static void isResonance(int idIn, bool isResonanceIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setIsResonance(isResonanceIn); }
  static void allParticle(int idIn, string nameIn, string antiNameIn, 
    int charge3In = 0, int colTypeIn = 0, double m0In = 0.,
    double widthIn = 0., double rangeIn = 0., double tau0In = 0., 
    bool mayDecayIn = false, bool isResonanceIn = false) 
    { if (isParticle(idIn)) pdt[abs(idIn)].setAll( nameIn, antiNameIn, 
    charge3In, colTypeIn, m0In, widthIn, rangeIn, tau0In, mayDecayIn,
    isResonanceIn); }  
  static void externalDecay(int idIn, bool externalDecayIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setExternalDecay(externalDecayIn); }
  static void hasChanged(int idIn, bool hasChangedIn) {
    if (isParticle(idIn)) pdt[abs(idIn)].setHasChanged(hasChangedIn); }
  static void rescaleBR(int idIn, double newSumBR = 1.) {
    if (isParticle(idIn)) pdt[abs(idIn)].rescaleBR(newSumBR); }
 
  // Give back current values. 
  static bool hasAnti(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].hasAnti() : false ; } 
  static string name(int idIn) {
    return (isParticle(abs(idIn))) ? pdt[abs(idIn)].name(idIn) : " "; }
  static int charge3(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].charge3(idIn) : 0 ; } 
  static double charge(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].charge(idIn) : 0 ; } 
  static int colType(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].colType(idIn) : 0 ; } 
  static double m0(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].m0() : 0. ; } 
  static double constituentMass(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].constituentMass() : 0. ; } 
  static double width(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].width() : 0. ; } 
  static double range(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].range() : 0. ; } 
  static double tau0(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].tau0() : 0. ; } 
  static bool mayDecay(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].mayDecay() : false ; } 
  static bool isResonance(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].isResonance() : false ; } 
  static bool externalDecay(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].externalDecay() : false ; } 
  static bool hasChanged(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].hasChanged() : false ; } 

  // Calculate a mass, picked according to Breit-Wigner.
  static double mass(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].mass() : 0. ; } 
  static int spinType(int idIn) {
    return isParticle(idIn) ? pdt[abs(idIn)].spinType() : 0 ; } 

private:

  // All particle data stored in a map.
  static map<int, ParticleDataEntry> pdt;

  // Flag that initialization has been performed.
  static bool isInit;

  // Methods to read a <particle> and <channel> line. 
  static bool readParticle(string line);
  static bool readChannel(string line);
  static ParticleDataEntry* particlePtr;

};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_ParticleData_H
