// HISubCollisionModel.h is a part of the PYTHIA event generator.
// Copyright (C) 2023 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the definition of the ImpactParmeterGenerator,
// SubCollision, and SubCollisionModel classes, as well as a set of
// subclasses of SubCollisionModel.
//
// ImpactParameterGenerator: distributes nuclei in impact parameter space.
// SubCollision: a collision between a projectile and a target Nucleon.
// SubCollisionModel: Models the collision probabilities of nucleons.
// BlackSubCollisionModel: A very simple SubCollisionModel.
// NaiveSubCollisionModel: A very simple SubCollisionModel.
// DoubleStrikmanSubCollisionModel: A more advanced SubCollisionModel.

#ifndef Pythia8_HISubCollisionModel_H
#define Pythia8_HISubCollisionModel_H

#include "Pythia8/Pythia.h"
#include "Pythia8/HINucleusModel.h"

namespace Pythia8 {

//==========================================================================

// SubCollision represents a possible collision between a projectile
// and a target nucleon.

class SubCollision {

public:

  // This defines the type of a binary nucleon collison.
  enum CollisionType {
    NONE,       // This is not a collision.
    ELASTIC,    // This is an elastic scattering
    SDEP,       // The projectile is diffractively excited.
    SDET,       // The target is diffractively excited.
    DDE,        // Both projectile and target are diffractively excited.
    CDE,        // Both excited but with central diffraction.
    ABS         // This is an absorptive (non-diffractive) collision.
  };

  // Constructor with configuration.
  SubCollision(Nucleon & projIn, Nucleon & targIn,
               double bIn, double bpIn, CollisionType typeIn)
    : proj(&projIn), targ(&targIn), b(bIn), bp(bpIn), type(typeIn) {}

  // Default constructor.
  SubCollision()
    : proj(0), targ(0), b(0.0), bp(0.0), type(NONE) {}

  // Used to order sub-collisions in a set.
  bool operator< (const SubCollision& s) const { return b < s.b; }

  // Return 0 if neither proj or target are neutrons, 1 if target is
  // neutron, 2 if projectile is neutron, and 3 if both are neutrons.
  int nucleons() const {return ( abs(targ->id()) == 2112? 1: 0 ) +
      ( abs(proj->id()) == 2112? 2: 0 );}

  // The projectile nucleon.
  Nucleon* proj;

  // The target nucleon.
  Nucleon* targ;

  // The impact parameter distance between the nucleons in femtometer.
  double b;

  // The impact parameter distance between the nucleons scaled like
  // in Pythia to have unit average for non-diffractive collisions.
  double bp;

  // The type of collision.
  CollisionType type;

};

//==========================================================================

// The SubCollisionSet gives a set of subcollisions between two nuclei.

class SubCollisionSet {

public:

  // Default constructor.
  SubCollisionSet() = default;

  // Constructor with subcollisions.
  SubCollisionSet(multiset<SubCollision> subCollisionsIn, double TIn)
    : subCollisionsSave(subCollisionsIn), TSave(TIn) {}

  // Reset the subcollisions.
  bool empty() const { return subCollisionsSave.empty(); }


  double T() const { return TSave; }

  // Iterators over the subcollisions.
  multiset<SubCollision>::const_iterator begin() const {
    return subCollisionsSave.begin(); }
  multiset<SubCollision>::const_iterator end() const {
    return subCollisionsSave.end(); }

private:

  // Saved subcollisions.
  multiset<SubCollision> subCollisionsSave;
  double TSave;

};

//==========================================================================

// The SubCollisionModel is is able to model the collision between two
// nucleons to tell which type of collision has occurred. The model
// may manipulate the corresponding state of the nucleons.

class SubCollisionModel {

public:

  // Internal class to report cross section estimates.
  struct SigEst {
    // The cross sections (tot, nd, dd, sdp, sdt, cd, el, bslope).
    vector<double> sig;

    // The estimated error (squared)
    vector<double> dsig2;

    // Which cross sections were actually fitted
    vector<bool> fsig;

    // The estimate of the average (and squared error) impact
    // parameter for inelastic non-diffractive collisions.
    double avNDb, davNDb2;

    // Constructor for zeros.
    SigEst(): sig(8, 0.0), dsig2(8, 0.0), fsig(8, false),
              avNDb(0.0), davNDb2(0.0) {}

  };

  // The default constructor is empty.
  SubCollisionModel(int nParm): sigTarg(8, 0.0), sigErr(8, 0.05),
    parmSave(nParm),
    NInt(100000), NPop(20), sigFuzz(0.2), impactFudge(1),
    fitPrint(true), avNDb(1.0*femtometer),
    projPtr(), targPtr(), sigTotPtr(), settingsPtr(), infoPtr(), rndmPtr() {}

  // Virtual destructor.
  virtual ~SubCollisionModel() {}

  // Create a new SubCollisionModel of the given model.
  static shared_ptr<SubCollisionModel> create(int model);

  // Virtual init method.
  virtual bool init(double eCMIn);

  // Initialize the pointers.
  void initPtr(NucleusModel & projIn, NucleusModel & targIn,
               SigmaTotal & sigTotIn, Settings & settingsIn,
               Info & infoIn, Rndm & rndmIn) {
    projPtr = &projIn;
    targPtr = &targIn;
    sigTotPtr = &sigTotIn;
    settingsPtr = &settingsIn;
    infoPtr = &infoIn;
    rndmPtr = &rndmIn;
    loggerPtr = infoIn.loggerPtr;
  }

  // Take two vectors of nucleons and an impact parameter vector and
  // produce the corrsponding sub-collisions. Note that states of the
  // nucleons may be changed. The function in this abstract base
  // class will reset the nucleon states for convenience. The
  // sub-collisions are ordered in the impact parameter distance
  // between the nucleons. The T-variable will be set to the summed
  // elastic amplitude.
  virtual SubCollisionSet getCollisions(Nucleus& proj, Nucleus& targ) = 0;

  // Access the nucleon-nucleon cross sections assumed
  // for this model.

  // The total cross section.
  double sigTot() const { return sigTarg[0]; }

  // The elastic cross section.
  double sigEl() const { return sigTarg[6]; }

  // The central diffractive excitation cross section.
  double sigCDE() const { return sigTarg[5]; }

  // The single diffractive excitation cross section (both sides summed).
  double sigSDE() const { return sigTarg[3] + sigTarg[4]; }

  // The single diffractive excitation cross section (excited projectile).
  double sigSDEP() const { return sigTarg[3]; }

  // The single diffractive excitation cross section (excited target).
  double sigSDET() const { return sigTarg[4]; }

  // The double diffractive excitation cross section.
  double sigDDE() const { return sigTarg[2]; }

  // The non-diffractive (absorptive) cross section.
  double sigND() const { return sigTarg[1]; }

  // The elastic b-slope parameter.
  double bSlope() const { return sigTarg[7]; }

  // Update internally stored cross sections.
  void updateSig();

  // Calculate the cross sections for the given set of parameters.
  virtual SigEst getSig() const { return SigEst(); }

  // Return the average non-diffractive impact parameter.
  double avNDB() const { return avNDb; }

  // Calculate the Chi2 for the given cross section estimates.
  double Chi2(const SigEst & sigs, int npar) const;

  // Use a simplified genetic algorithm to fit the parameters.
  virtual bool evolve(int nGenerations, double eCM);

  // Get the number of free parameters for the model.
  int nParms() const { return parmSave.size(); }

  // Set the parameters of this model.
  void setParm(const vector<double>& parmIn) {
    for (size_t i = 0; i < parmSave.size(); ++i)
      parmSave[i] = parmIn[i];
  }

  // Get the current parameters of this model.
  vector<double> getParm() const {
    return parmSave;
  }

  // Get the minimum allowed parameter values for this model.
  virtual vector<double> minParm() const = 0;

  // Get the maximum allowed parameter values for this model.
  virtual vector<double> maxParm() const = 0;

  // Set beam kinematics.
  void setKinematics(double eCMIn);

private:

  // The nucleon-nucleon cross sections targets for this model
  // (tot, nd, dd, sdp, sdt, cd, el, bslope) and the required precision.
  vector<double> sigTarg, sigErr;

  // Generate parameters based on run settings and the evolutionary algorithm.
  bool genParms();

  // Save/load parameter configuration from disk.
  bool saveParms(string fileName) const;
  bool loadParms(string fileName);

protected:

  // Saved parameters.
  vector<double> parmSave;

  // The parameters stearing the fitting of internal parameters to
  // the different nucleon-nucleon cross sections.
  int NInt, NPop;
  double sigFuzz;
  double impactFudge;
  bool fitPrint;

  // The estimated average impact parameter distance (in femtometer)
  // for absorptive collisions.
  double avNDb;

  // Info from the controlling HeavyIons object.
  NucleusModel* projPtr;
  NucleusModel* targPtr;
  SigmaTotal* sigTotPtr;
  Settings* settingsPtr;
  Info* infoPtr;
  Rndm* rndmPtr;
  Logger* loggerPtr;

  // For variable energies.
  bool doVarECM;
  double eMin{}, eMax{};
  int eCMPts;
  vector<LogInterpolator> subCollParms;

};

//==========================================================================

// The most naive sub-collision model, assuming static nucleons and
// an absorptive cross section equal to the total inelastic. No
// fluctuations, meaning no diffraction.

class BlackSubCollisionModel : public SubCollisionModel {

public:

  // The default constructor simply lists the nucleon-nucleon cross sections.
  BlackSubCollisionModel() : SubCollisionModel(0) {}

  // Virtual destructor.
  virtual ~BlackSubCollisionModel() {}

  // Get the minimum and maximum allowed parameter values for this model.
  vector<double> minParm() const override { return vector<double>(); }
  vector<double> maxParm() const override { return vector<double>(); }

  // Take two vectors of Nucleons and an impact parameter vector and
  // produce the corrsponding sub-collisions. Note that states of the
  // nucleons may be changed.
  virtual SubCollisionSet getCollisions(Nucleus& proj, Nucleus& targ) override;

};

//==========================================================================

// A very simple sub-collision model, assuming static nucleons and
// just assuring that the individual nucleon-nucleon cross sections
// are preserved.

class NaiveSubCollisionModel : public SubCollisionModel {

public:

  // The default constructor simply lists the nucleon-nucleon cross sections.
  NaiveSubCollisionModel() : SubCollisionModel(0) {}

  // Virtual destructor.
  virtual ~NaiveSubCollisionModel() {}

  // Get the minimum and maximum allowed parameter values for this model.
  vector<double> minParm() const override { return vector<double>(); }
  vector<double> maxParm() const override { return vector<double>(); }

  // Take two vectors of Nucleons and an impact parameter vector and
  // produce the corrsponding sub-collisions. Note that states of the
  // nucleons may be changed.
  virtual SubCollisionSet getCollisions(Nucleus& proj, Nucleus& targ) override;

};

//==========================================================================

// A sub-collision model where each nucleon has a fluctuating
// "radius" according to a Strikman-inspired distribution.

class DoubleStrikmanSubCollisionModel : public SubCollisionModel {

public:

  // The default constructor simply lists the nucleon-nucleon cross sections.
  DoubleStrikmanSubCollisionModel(int modein = 0) : SubCollisionModel(3),
    sigd(parmSave[0]), k0(parmSave[1]), alpha(parmSave[2]),
    opacityMode(modein) {}

  // Virtual destructor.
  virtual ~DoubleStrikmanSubCollisionModel() {}

  // Take two vectors of Nucleons and an impact parameter vector and
  // produce the corrsponding sub-collisions. Note that states of the
  // nucleons may be changed.
  virtual SubCollisionSet getCollisions(Nucleus& proj, Nucleus& targ) override;

  // Calculate the cross sections for the given set of parameters.
  SigEst getSig() const override;

  // Get the minimum and maximum allowed parameter values for this model.
  vector<double> minParm() const override { return {  1.0,  0.01, 0.0 }; }
  vector<double> maxParm() const override { return { 60.0, 60.00, 20.0 }; }

private:

  // Saturation scale of the nucleus.
  double& sigd;

  // The power in the Gamma distribution.
  double& k0;

  // Power of the saturation scale
  double& alpha;

  // Optional mode for opacity.
  int opacityMode;

  // Return the average radius deduced from other parameters and
  // the toal cross section.
  double r0() const {
    return sqrt(sigTot() / (M_PI * (2.0 * k0 + 4.0 * k0 * k0)));
  }

  // The opacity of the collision at a given sigma.
  double opacity(double sig) const {
    sig /= sigd;
    if ( opacityMode == 1 ) sig = 1.0/sig;
    return sig > numeric_limits<double>::epsilon() ?
      pow(-expm1(-1.0/sig), alpha) : 1.0;
  }

  // Return the elastic amplitude for a projectile and target state
  // and the impact parameter between the corresponding nucleons.
  double Tpt(const Nucleon::State & p,
             const Nucleon::State & t, double b) const {
    double sig = M_PI*pow2(p[0] + t[0]);
    double grey = opacity(sig);
    return sig/grey > b*b*2.0*M_PI? grey: 0.0;
  }

  // Helper functions.
  static void shuffle(double PND1, double PND2,
                      double & PW1, double & PW2);
  static void shuffle(double & PEL11, double P11,
                      double P12, double P21, double P22);
  static double pnw(double PWp, double PWt, double PND) {
    return ( 1.0 - PWp <= 0.0 || 1.0 - PWt <= 0.0 )?
      0.0: (1.0 - PWp)*(1.0 - PWt)/(1.0 - PND);
  }

};

//==========================================================================

// ImpactParameterGenerator is able to generate a specific impact
// parameter together with a weight such that aweighted average over
// any quantity X(b) corresponds to the infinite integral over d^2b
// X(b). This base class gives a Gaussian profile, d^2b exp(-b^2/2w^2).

class ImpactParameterGenerator {

public:

  // The default constructor takes a general width (in femtometers) as
  // argument.
  ImpactParameterGenerator()
    : widthSave(0.0), collPtr(0), projPtr(0), targPtr(0),
      settingsPtr(0), rndmPtr(0) {}

  // Virtual destructor.
  virtual ~ImpactParameterGenerator() {}

  // Virtual init method.
  virtual bool init();
  void initPtr(Info & infoIn, SubCollisionModel & collIn,
    NucleusModel & projIn, NucleusModel & targIn);

  // Return a new impact parameter and set the corresponding weight provided.
  virtual Vec4 generate(double & weight) const;

  // Set the width (in femtometers).
  void width(double widthIn) { widthSave = widthIn; }

  // Get the width.
  double width() const { return widthSave; }

  // Update width based on the associated subcollision and nucleus models.
  void updateWidth();

private:

  // The width of a distribution.
  double widthSave;

protected:

  // Pointers from the controlling HeavyIons object.
  Info* infoPtr;
  SubCollisionModel* collPtr;
  NucleusModel* projPtr;
  NucleusModel* targPtr;
  Settings* settingsPtr;
  Rndm* rndmPtr;
  Logger* loggerPtr;

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_HISubCollisionModel_H
