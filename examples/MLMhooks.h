// MLMhooks.h is a part of the PYTHIA event generator.
// Copyright (C) 2012 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Richard Corke (richard.corke@thep.lu.se)
// This file provides the MLMhooks class to perform MLM merging.
// Example usage is shown in main32.cc, and further details
// can be found in the 'Alpgen and MLM Merging' manual page.

#ifndef MLMHOOKS_H
#define MLMHOOKS_H

// Includes
#include "Pythia.h"
using namespace Pythia8;

//==========================================================================

// Preprocessor settings

// Debug flag to give verbose information at all stages of the merging.
//#define MLMDEBUG
// Debug flag to enable some extra checks in the code
//#define MLMCHECK

//==========================================================================

// Declaration of main MLMhooks class to perform MLM matching.
// Note that it is defined with virtual inheritance, so that it can
// be combined with other UserHooks classes, see e.g. main32.cc.

class MLMhooks : virtual public UserHooks {

public:

  // Constructor and destructor
  MLMhooks() : cellJet(NULL), slowJet(NULL) {}
  ~MLMhooks() {
    if (cellJet) delete cellJet;
    if (slowJet) delete slowJet; 
  }

  // Initialisation
  virtual bool initAfterBeams();

  // Process level vetos
  virtual bool canVetoProcessLevel();
  virtual bool doVetoProcessLevel(Event &);

  // Parton level vetos (before beam remnants and resonance decays)
  virtual bool canVetoPartonLevelEarly();
  virtual bool doVetoPartonLevelEarly(const Event &);

private:

  // Different steps of the MLM matching algorithm
  void sortIncomingProcess(const Event &);
  void jetAlgorithmInput(const Event &, int);
  void runJetAlgorithm();
  bool matchPartonsToJets(int);
  int  matchPartonsToJetsLight();
  int  matchPartonsToJetsHeavy();

  // DeltaR between two 4-vectors (eta and y variants)
  inline double Vec4eta(const Vec4 &pIn) {
    return -log(tan(pIn.theta() / 2.));
  }
  inline double Vec4y(const Vec4 &pIn) {
    return 0.5 * log((pIn.e() + pIn.pz()) / (pIn.e() - pIn.pz()));
  }
  inline double deltaReta(const Vec4 &p1, const Vec4 &p2) {
    double dEta = abs(Vec4eta(p1) - Vec4eta(p2));
    double dPhi = abs(p1.phi() - p2.phi());
    if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
    return sqrt(dEta*dEta + dPhi*dPhi);
  }
  inline double deltaRy(const Vec4 &p1, const Vec4 &p2) {
    double dy   = abs(Vec4y(p1) - Vec4y(p2));
    double dPhi = abs(p1.phi() - p2.phi());
    if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
    return sqrt(dy*dy + dPhi*dPhi);
  }

  // Function to sort typeIdx vectors into descending eT/pT order.
  // Uses a selection sort, as number of partons generally small
  // and so efficiency not a worry.
  void sortTypeIdx(vector < int > &vecIn) {
    for (size_t i = 0; i < vecIn.size(); i++) {
      size_t jMax = i;
      double vMax = (jetAlgorithm == 1) ?
                    eventProcess[vecIn[i]].eT() :
                    eventProcess[vecIn[i]].pT();
      for (size_t j = i + 1; j < vecIn.size(); j++) {
        double vNow = (jetAlgorithm == 1) ?
                      eventProcess[vecIn[j]].eT() :
                      eventProcess[vecIn[j]].pT();
        if (vNow > vMax) {
          vMax = vNow;
          jMax = j;
        }
      }
      if (jMax != i) swap(vecIn[i], vecIn[jMax]);
    }
  }

  // Master switch for merging
  bool   doMerge;

  // Maximum and current number of jets
  int    nJetMax, nJet;

  // Jet algorithm parameters
  int    jetAlgorithm;
  double eTjetMin, coneRadius, etaJetMax, etaJetMaxAlgo;

  // CellJet specific
  int    nEta, nPhi;
  double eTseed, eTthreshold;

  // SlowJet specific
  int    slowJetPower;

  // Merging procedure parameters
  int    jetAllow, jetMatch, exclusiveMode;
  double coneMatchLight, coneRadiusHeavy, coneMatchHeavy;
  bool   exclusive;

  // Event records to store original incoming process, final-state of the
  // incoming process and what will be passed to the jet algorithm.
  // Not completely necessary to store all steps, but makes tracking the
  // steps of the algorithm a lot easier.
  Event eventProcessOrig, eventProcess, workEventJet;

  // Internal jet algorithms
  CellJet *cellJet;
  SlowJet *slowJet;

  // Sort final-state of incoming process into light/heavy jets and 'other'
  vector < int > typeIdx[3];
  set    < int > typeSet[3];

  // Momenta output of jet algorithm (to provide same output regardless of
  // the selected jet algorithm)
  vector < Vec4 > jetMomenta;

  // Store the minimum eT/pT of matched light jets
  double eTpTlightMin;

  // Constants
  static const double GHOSTENERGY, ZEROTHRESHOLD;
};

//--------------------------------------------------------------------------

// Main implementation of MLMhooks class. This may be split out to a
// separate C++ file if desired, but currently included here for ease
// of use.

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// The energy of ghost particles. For technical reasons, this cannot be
// set arbitrarily low, see 'Particle::TINY' in 'Event.cc' for details.
const double MLMhooks::GHOSTENERGY   = 1e-15;

// A zero threshold value for double comparisons.
const double MLMhooks::ZEROTHRESHOLD = 1e-10;

//--------------------------------------------------------------------------

// Initialisation routine automatically called from Pythia::init().
// Setup all parts needed for the merging.

bool MLMhooks::initAfterBeams() {

  // Read in parameters
  doMerge         = settingsPtr->flag("MLM:merge");
  nJet            = settingsPtr->mode("MLM:nJet");
  nJetMax         = settingsPtr->mode("MLM:nJetMax");
  jetAlgorithm    = settingsPtr->mode("MLM:jetAlgorithm");
  eTjetMin        = settingsPtr->parm("MLM:eTjetMin");
  coneRadius      = settingsPtr->parm("MLM:coneRadius");
  etaJetMax       = settingsPtr->parm("MLM:etaJetMax");
  // Use etaJetMax + coneRadius in input to jet algorithms
  etaJetMaxAlgo   = etaJetMax + coneRadius;
  // CellJet specific
  nEta            = settingsPtr->mode("MLM:nEta");
  nPhi            = settingsPtr->mode("MLM:nPhi");
  eTseed          = settingsPtr->parm("MLM:eTseed");
  eTthreshold     = settingsPtr->parm("MLM:eTthreshold");
  // SlowJet specific
  slowJetPower    = settingsPtr->mode("MLM:slowJetPower");
  // Matching procedure
  jetAllow        = settingsPtr->mode("MLM:jetAllow");
  jetMatch        = settingsPtr->mode("MLM:jetMatch");
  coneMatchLight  = settingsPtr->parm("MLM:coneMatchLight");
  coneRadiusHeavy = settingsPtr->parm("MLM:coneRadiusHeavy");
  if (coneRadiusHeavy < 0.) coneRadiusHeavy = coneRadius;
  coneMatchHeavy  = settingsPtr->parm("MLM:coneMatchHeavy");
  exclusiveMode   = settingsPtr->mode("MLM:exclusive");

  // If not merging, then done
  if (!doMerge) return true;

  // Exclusive mode; if set to 2, then set based on nJet/nJetMax
  if (exclusiveMode == 2) {

    // No nJet or nJetMax, so default to exclusive mode
    if (nJet < 0 || nJetMax < 0) {
      infoPtr->errorMsg("Warning in MLMhooks:init: "
          "missing jet multiplicity information; running in exclusive mode");
      exclusive = true;

    // Inclusive if nJet == nJetMax, exclusive otherwise
    } else {
      exclusive = (nJet == nJetMax) ? false : true;
    }

  // Otherwise, just set as given
  } else {
    exclusive = (exclusiveMode == 0) ? false : true;
  }

  // Initialise chosen jet algorithm. CellJet.
  if (jetAlgorithm == 1) {

    // Extra options for CellJet. nSel = 1 means that all final-state
    // particles are taken and we retain control of what to select.
    // smear/resolution/upperCut are not used and are set to default values.
    int    nSel = 2, smear = 0;
    double resolution = 0.5, upperCut = 2.;
    cellJet = new CellJet(etaJetMaxAlgo, nEta, nPhi, nSel,
                          smear, resolution, upperCut, eTthreshold);

  // SlowJet
  } else if (jetAlgorithm == 2) {
    slowJet = new SlowJet(slowJetPower, coneRadius, eTjetMin, etaJetMaxAlgo);
  }

  // Check the jetMatch parameter; option 2 only works with SlowJet
  if (jetAlgorithm == 1 && jetMatch == 2) {
    infoPtr->errorMsg("Warning in MLMhooks:init: "
        "jetMatch = 2 only valid with SlowJet algorithm. "
        "Reverting to jetMatch = 1.");
    jetMatch = 1;
  }

  // Setup local event records
  eventProcessOrig.init("(eventProcessOrig)", particleDataPtr);
  eventProcess.init("(eventProcess)", particleDataPtr);
  workEventJet.init("(workEventJet)", particleDataPtr);

  // Print information
  string jetStr  = (jetAlgorithm ==  1) ? "CellJet" :
                   (slowJetPower == -1) ? "anti-kT" :
                   (slowJetPower ==  0) ? "C/A"     :
                   (slowJetPower ==  1) ? "kT"      : "unknown";
  string modeStr = (exclusive)         ? "exclusive" : "inclusive";
  stringstream nJetStr, nJetMaxStr;
  if (nJet >= 0)    nJetStr    << nJet;    else nJetStr    << "unknown";
  if (nJetMax >= 0) nJetMaxStr << nJetMax; else nJetMaxStr << "unknown";
  cout << endl
       << " *-------  MLM matching parameters  -------*" << endl
       << " |  nJet                |  " << setw(14)
       << nJetStr.str() << "  |" << endl
       << " |  nJetMax             |  " << setw(14)
       << nJetMaxStr.str() << "  |" << endl
       << " |  Jet algorithm       |  " << setw(14)
       << jetStr << "  |" << endl
       << " |  eTjetMin            |  " << setw(14)
       << eTjetMin << "  |" << endl
       << " |  coneRadius          |  " << setw(14)
       << coneRadius << "  |" << endl
       << " |  etaJetMax           |  " << setw(14)
       << etaJetMax << "  |" << endl
       << " |  jetAllow            |  " << setw(14)
       << jetAllow << "  |" << endl
       << " |  jetMatch            |  " << setw(14)
       << jetMatch << "  |" << endl
       << " |  coneMatchLight      |  " << setw(14)
       << coneMatchLight << "  |" << endl
       << " |  coneRadiusHeavy     |  " << setw(14)
       << coneRadiusHeavy << "  |" << endl
       << " |  coneMatchHeavy      |  " << setw(14)
       << coneMatchHeavy << "  |" << endl
       << " |  Mode                |  " << setw(14)
       << modeStr << "  |" << endl
       << " *-----------------------------------------*" << endl;

  return true;
}

//--------------------------------------------------------------------------

// Process level veto. Stores incoming event for later.

bool MLMhooks::canVetoProcessLevel() { return doMerge; }

bool MLMhooks::doVetoProcessLevel(Event& process) { 

  // Copy incoming process
  eventProcessOrig = process;
  return false;
}

//--------------------------------------------------------------------------

// Early parton level veto (before beam remnants and resonance showers)

bool MLMhooks::canVetoPartonLevelEarly() { return doMerge; }

bool MLMhooks::doVetoPartonLevelEarly(const Event& event) {

  // 1) Sort the original incoming process. After this step is performed,
  //    the following assignments have been made:
  //      eventProcessOrig - the original incoming process
  //      eventProcess     - the final-state of the incoming process with
  //                         resonance decays removed (and resonances
  //                         themselves now with positive status code)
  //      typeIdx[0/1/2]   - indices into 'eventProcess' of
  //                         light jets/heavy jets/other
  //      typeSet[0/1/2]   - indices into 'event' of
  //                         light jets/heavy jets/other
  //      workEvent        - partons from the hardest subsystem
  //                         + ISR + FSR only
  sortIncomingProcess(event);

  // Debug
#ifdef MLMDEBUG
  // Begin
  cout << endl << "---------- Begin MLM Debug ----------" << endl;

  // Original incoming process
  cout << endl << "Original incoming process:";
  eventProcessOrig.list();

  // Final-state of original incoming process
  cout << endl << "Final-state incoming process:";
  eventProcess.list();

  // List categories of sorted particles
  for (size_t i = 0; i < typeIdx[0].size(); i++) 
    cout << ((i == 0) ? "Light jets: " : ", ")
         << setw(3) << typeIdx[0][i];
  for (size_t i = 0; i < typeIdx[1].size(); i++) 
    cout << ((i == 0) ? "\nHeavy jets: " : ", ")
         << setw(3) << typeIdx[1][i];
  for (size_t i = 0; i < typeIdx[2].size(); i++) 
    cout << ((i == 0) ? "\nOther:      " : ", ")
         << setw(3) << typeIdx[2][i];

  // Full event at this stage
  cout << endl << endl << "Event:";
  event.list();

  // Work event (partons from hardest subsystem + ISR + FSR)
  cout << endl << "Work event:";
  workEvent.list();
#endif

  // 2) Light/heavy jets: iType = 0 (light jets), 1 (heavy jets)
  int iTypeEnd = (typeIdx[1].empty()) ? 1 : 2;
  for (int iType = 0; iType < iTypeEnd; iType++) {

    // 2a) Find particles which will be passed from the jet algorithm.
    //     Input from 'workEvent' and output in 'workEventJet'.
    jetAlgorithmInput(event, iType);

    // Debug
#ifdef MLMDEBUG
    // Jet algorithm event
    cout << endl << "Jet algorithm event (iType = " << iType << "):";
    workEventJet.list();
#endif

    // 2b) Run jet algorithm on 'workEventJet'.
    //     Output is stored in jetMomenta.
    runJetAlgorithm();

    // 2c) Match partons to jets and decide if veto is necessary
    if (matchPartonsToJets(iType) == true) {
      // Debug
#ifdef MLMDEBUG
      cout << endl << "Event vetoed" << endl
           << "----------  End MLM Debug  ----------" << endl;
#endif
      return true;
    }
  }

  // Debug
#ifdef MLMDEBUG
  cout << endl << "Event accepted" << endl
       << "----------  End MLM Debug  ----------" << endl;
#endif

  // If we reached here, then no veto
  return false;
}

//--------------------------------------------------------------------------

// Step (1): sort the incoming particles

void MLMhooks::sortIncomingProcess(const Event &event) {

  // Remove resonance decays from original process and keep only final
  // state. Resonances will have positive status code after this step.
  omitResonanceDecays(eventProcessOrig, true);
  eventProcess = workEvent;

  // Sort original process final state into light/heavy jets and 'other'.
  // Criteria:
  //   1 <= ID <= 5 and massless, or ID == 21 --> light jet (typeIdx[0])
  //   4 <= ID <= 6 and massive               --> heavy jet (typeIdx[1])
  //   All else                               --> other     (typeIdx[2])
  // Note that 'typeIdx' stores indices into 'eventProcess' (after resonance
  // decays are omitted), while 'typeSet' stores indices into the original
  // process record, 'eventProcessOrig', but these indices are also valid
  // in 'event'.
  for (int i = 0; i < 3; i++) {
    typeIdx[i].clear();
    typeSet[i].clear();
  }
  for (int i = 0; i < eventProcess.size(); i++) {
    // Ignore nonfinal and default to 'other'
    if (!eventProcess[i].isFinal()) continue;
    int idx = 2;

    // Light jets
    if (eventProcess[i].id() == 21 || (eventProcess[i].idAbs() <= 5 &&
        abs(eventProcess[i].m()) < ZEROTHRESHOLD))
      idx = 0;

    // Heavy jets
    else if (eventProcess[i].idAbs() >= 4 && eventProcess[i].idAbs() <= 6)
      idx = 1;

    // Store
    typeIdx[idx].push_back(i);
    typeSet[idx].insert(eventProcess[i].daughter1());
  }

  // Extract partons from hardest subsystem + ISR + FSR only into
  // workEvent. Note no resonance showers or MPIs.
  subEvent(event);
}

//--------------------------------------------------------------------------

// Step (2a): pick which particles to pass to the jet algorithm

void MLMhooks::jetAlgorithmInput(const Event &event, int iType) {

  // Take input from 'workEvent' and put output in 'workEventJet'
  workEventJet = workEvent;

  // Loop over particles and decide what to pass to the jet algorithm
  for (int i = 0; i < workEventJet.size(); ++i) {
    if (!workEventJet[i].isFinal()) continue;

    // jetAllow option to disallow certain particle types
    if (jetAllow == 1) {

      // Original AG+Py6 algorithm explicitly excludes tops,
      // leptons and photons.
      int id = workEventJet[i].idAbs();
      if ((id >= 11 && id <= 16) || id == 6 || id == 22) {
        workEventJet[i].statusNeg();
        continue;
      }
    }

    // Get the index of this particle in original event
    int idx = workEventJet[i].daughter1();

    // Start with particle idx, and afterwards track mothers
    while (true) {

      // Light jets
      if (iType == 0) {

        // Do not include if originates from heavy jet or 'other'
        if (typeSet[1].find(idx) != typeSet[1].end() ||
            typeSet[2].find(idx) != typeSet[2].end()) {
          workEventJet[i].statusNeg();
          break;
        }

        // Made it to start of event record so done
        if (idx == 0) break;
        // Otherwise next mother and continue
        idx = event[idx].mother1();

      // Heavy jets
      } else if (iType == 1) {

        // Only include if originates from heavy jet
        if (typeSet[1].find(idx) != typeSet[1].end()) break;

        // Made it to start of event record with no heavy jet mother,
        // so DO NOT include particle
        if (idx == 0) {
          workEventJet[i].statusNeg();
          break;
        }

        // Otherwise next mother and continue
        idx = event[idx].mother1();

      } // if (iType)
    } // while (true)
  } // for (i)

  // For jetMatch = 2, insert ghost particles corresponding to
  // each hard parton in the original process
  if (jetMatch == 2) {
    for (int i = 0; i < int(typeIdx[iType].size()); i++) {
      // Get y/phi of the parton
      Vec4   pIn = eventProcess[typeIdx[iType][i]].p();
      double y   = Vec4y(pIn);
      double phi = pIn.phi();

      // Create a ghost particle and add to the workEventJet
      double e   = GHOSTENERGY;
      double e2y = exp(2. * y);
      double pz  = e * (e2y - 1.) / (e2y + 1.);
      double pt  = sqrt(e*e - pz*pz);
      double px  = pt * cos(phi);
      double py  = pt * sin(phi);
      workEventJet.append(Particle(21, 99, 0, 0, 0, 0, 0, 0,
                                px, py, pz, e, 0., 0, 9.));

      // Extra check on reconstructed y/phi values. If many warnings
      // of this type, GHOSTENERGY may be set too low.
#ifdef MLMCHECK
      int lastIdx = workEventJet.size() - 1;
      if (abs(y   - workEventJet[lastIdx].y())   > ZEROTHRESHOLD ||
          abs(phi - workEventJet[lastIdx].phi()) > ZEROTHRESHOLD)
        infoPtr->errorMsg("Warning in MLMhooks:jetAlgorithmInput: "
            "ghost particle y/phi mismatch");
#endif

    } // for (i)
  } // if (jetMatch == 2)
}

//--------------------------------------------------------------------------

// Step (2b): run jet algorithm and provide common output

void MLMhooks::runJetAlgorithm() {

  // Run the jet clustering algorithm
  if (jetAlgorithm == 1)
    cellJet->analyze(workEventJet, eTjetMin, coneRadius, eTseed);
  else
    slowJet->analyze(workEventJet);

  // Extract four-momenta of jets with |eta| < etaJetMax and
  // put into jetMomenta. Note that this is done backwards as
  // jets are removed with SlowJet.
  jetMomenta.clear();
  int iJet = (jetAlgorithm == 1) ? cellJet->size() - 1:
                                   slowJet->sizeJet() - 1;
  for (int i = iJet; i > -1; i--) {
    Vec4 jetMom = (jetAlgorithm == 1) ? cellJet->pMassive(i) :
                                        slowJet->p(i);
    double eta = Vec4eta(jetMom);

    if (abs(eta) > etaJetMax) {
      if (jetAlgorithm == 2) slowJet->removeJet(i);
      continue;
    }
    jetMomenta.push_back(jetMom);
  }

  // Reverse jetMomenta to restore eT/pT ordering
  reverse(jetMomenta.begin(), jetMomenta.end());
}

//--------------------------------------------------------------------------

// Step (2c): veto decision (returning true vetoes the event)

bool MLMhooks::matchPartonsToJets(int iType) {

  // Use two different routines for light/heavy jets as
  // different veto conditions and for clarity
  if (iType == 0) return (matchPartonsToJetsLight() > 0);
  else            return (matchPartonsToJetsHeavy() > 0);
}

//--------------------------------------------------------------------------

// Step(2c): light jets
// Return codes are given indicating the reason for a veto.
// Although not currently used, they are a useful debugging tool:
//   0 = no veto
//   1 = veto as number of jets less than number of partons
//   2 = veto as exclusive mode and number of jets greater than
//       number of partons
//   3 = veto as inclusive mode and there would be an extra jet
//       that is harder than any matched soft jet
//   4 = veto as there is a parton which does not match a jet

int MLMhooks::matchPartonsToJetsLight() {

  // Always veto if number of jets is less than original number of jets
  if (jetMomenta.size() < typeIdx[0].size()) return 1;
  // Veto if in exclusive mode and number of jets bigger than original
  if (exclusive && jetMomenta.size() > typeIdx[0].size()) return 2;

  // Sort partons by eT/pT
  sortTypeIdx(typeIdx[0]);

  // Number of hard partons
  int nParton = typeIdx[0].size();

  // Keep track of which jets have been assigned a hard parton
  vector < bool > jetAssigned;
  jetAssigned.assign(jetMomenta.size(), false);

  // Jet matching procedure: (1) deltaR between partons and jets
  if (jetMatch == 1) {

    // Loop over light hard partons and get 4-momentum
    for (int i = 0; i < nParton; i++) {
      Vec4 p1 = eventProcess[typeIdx[0][i]].p();

      // Track which jet has the minimal dR measure with this parton
      int    jMin  = -1;
      double dRmin = 0.;

      // Loop over all jets (skipping those already assigned).
      for (int j = 0; j < int(jetMomenta.size()); j++) {
        if (jetAssigned[j]) continue;

        // DeltaR between parton/jet and store if minimum
        double dR = (jetAlgorithm == 1) ?
            deltaReta(p1, jetMomenta[j]) : deltaRy(p1, jetMomenta[j]);
        if (jMin < 0 || dR < dRmin) {
          dRmin = dR;
          jMin  = j;
        }
      } // for (j)

      // Check for jet-parton match
      if (jMin >= 0 && dRmin < coneRadius * coneMatchLight) {

        // If the matched jet is not one of the nParton hardest jets,
        // the extra left over jet would be harder than some of the
        // matched jets. This is disallowed, so veto.
        if (jMin >= nParton) return 3;

        // Mark jet as assigned.
        jetAssigned[jMin] = true;

      // If no match, then event will be vetoed in all cases
      } else return 4;

    } // for (i)

  // Jet matching procedure: (2) ghost particles in SlowJet
  } else {

    // Loop over added 'ghost' particles and find if assigned to a jet
    for (int i = workEventJet.size() - nParton;
        i < workEventJet.size(); i++) {
      int jMin = slowJet->jetAssignment(i);

      // Veto if:
      //  1) not one of nParton hardest jets
      //  2) not assigned to a jet
      //  3) jet has already been assigned
      if (jMin >= nParton)               return 3;
      if (jMin < 0 || jetAssigned[jMin]) return 4;

      // Mark jet as assigned
      jetAssigned[jMin] = true;

    } // for (i)
  } // if (jetMatch)

  // Minimal eT/pT (CellJet/SlowJet) of matched light jets. Needed
  // later for heavy jet vetos in inclusive mode.
  if (nParton > 0)
    eTpTlightMin = (jetAlgorithm == 1) ? jetMomenta[nParton - 1].eT()
                                       : jetMomenta[nParton - 1].pT();
  else
    eTpTlightMin = -1.;

  // No veto
  return 0;
}

//--------------------------------------------------------------------------

// Step(2c): heavy jets
// Return codes are given indicating the reason for a veto.
// Although not currently used, they are a useful debugging tool:
//   0 = no veto as there are no extra jets present
//   1 = veto as in exclusive mode and extra jets present
//   2 = veto as in inclusive mode and extra jets were harder
//       than any matched light jet

int MLMhooks::matchPartonsToJetsHeavy() {

  // If there are no extra jets, then accept
  if (jetMomenta.empty()) return 0;

  // Number of hard partons
  int nParton = typeIdx[1].size();

  // Remove jets that are close to heavy quarks
  set < int > removeJets;

  // Jet matching procedure: (1) deltaR between partons and jets
  if (jetMatch == 1) {

    // Loop over heavy hard partons and get 4-momentum
    for (int i = 0; i < nParton; i++) {
      Vec4 p1 = eventProcess[typeIdx[1][i]].p();

      // Loop over all jets, find dR and mark for removal if match
      for (int j = 0; j < int(jetMomenta.size()); j++) {
        double dR = (jetAlgorithm == 1) ?
            deltaReta(p1, jetMomenta[j]) : deltaRy(p1, jetMomenta[j]);
        if (dR < coneRadiusHeavy * coneMatchHeavy)
          removeJets.insert(j);

      } // for (j)
    } // for (i)

  // Jet matching procedure: (2) ghost particles in SlowJet
  } else {

    // Loop over added 'ghost' particles and if assigned to a jet
    // then mark this jet for removal
    for (int i = workEventJet.size() - nParton;
        i < workEventJet.size(); i++) {
      int jMin = slowJet->jetAssignment(i);
      if (jMin >= 0) removeJets.insert(jMin);
    }
      
  }

  // Remove jets (backwards order to not disturb indices)
  for (set < int >::reverse_iterator it  = removeJets.rbegin();
                                     it != removeJets.rend(); it++)
    jetMomenta.erase(jetMomenta.begin() + *it);

  // Handle case if there are still extra jets
  if (!jetMomenta.empty()) {

    // Exclusive mode, so immediate veto
    if (exclusive) return 1;

    // Inclusive mode; extra jets must be softer than any matched light jet
    else if (eTpTlightMin >= 0.)
      for (size_t j = 0; j < jetMomenta.size(); j++) {
        // CellJet uses eT, SlowJet uses pT
        if ( (jetAlgorithm == 1 && jetMomenta[j].eT() > eTpTlightMin) ||
             (jetAlgorithm == 2 && jetMomenta[j].pT() > eTpTlightMin) )
          return 2;
      }

  } // if (!jetMomenta.empty())

  // No extra jets were present so no veto
  return 0;
}

//==========================================================================

#endif // MLMHOOKS_H
