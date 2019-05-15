// MergingHooks.h is a part of the PYTHIA event generator.
// Copyright (C) 2012 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file is written by Stefan Prestel.
// Header file to allow user access to program at different stages.
// HardProcess: Container class for the hard process to be merged. Holds the
//              bookkeeping of particles not be be reclustered
// MergingHooks: Steering class for matrix element merging. Some functions can
//               be redefined in a derived class to have access to the merging

#ifndef Pythia8_MergingHooks_H
#define Pythia8_MergingHooks_H

#include "Basics.h"
#include "BeamParticle.h"
#include "Event.h"
#include "Info.h"
#include "ParticleData.h"
#include "PartonSystems.h"
#include "PythiaStdlib.h"
#include "Settings.h"

namespace Pythia8 {

//==========================================================================

// Declaration of hard process class
// This class holds information on the desired hard 2->2 process 
// for the merging.
// This class is a container class for History class use.

class HardProcess {

public:

   // Flavour of the first incoming particle
  int hardIncoming1;
  // Flavour of the second incoming particle
  int hardIncoming2;
  // Flavours of the outgoing particles
  vector<int> hardOutgoing1;
  vector<int> hardOutgoing2;
  // Flavour of intermediate bosons in the hard 2->2
  vector<int> hardIntermediate;

  // Current reference event
  Event state;
  // Potential positions of outgoing particles in reference event
  vector<int> PosOutgoing1;
  vector<int> PosOutgoing2;
  // Potential positions of intermediate bosons in reference event
  vector<int> PosIntermediate;

  // Information on merging scale read from LHE file
  double tms;
  // Type of ME generator
  int meGenType;

  // Default constructor
  HardProcess(){}
  // Default destructor
  ~HardProcess(){}

  // Copy constructor
  HardProcess( const HardProcess& hardProcessIn )
    : state(hardProcessIn.state),
      tms(hardProcessIn.tms), meGenType(hardProcessIn.meGenType) {
    hardIncoming1 = hardProcessIn.hardIncoming1;
    hardIncoming2 = hardProcessIn.hardIncoming2;
    for(int i =0; i < int(hardProcessIn.hardOutgoing1.size());++i)
      hardOutgoing1.push_back( hardProcessIn.hardOutgoing1[i]);
    for(int i =0; i < int(hardProcessIn.hardOutgoing2.size());++i)
      hardOutgoing2.push_back( hardProcessIn.hardOutgoing2[i]);
    for(int i =0; i < int(hardProcessIn.hardIntermediate.size());++i)
      hardIntermediate.push_back( hardProcessIn.hardIntermediate[i]);
    for(int i =0; i < int(hardProcessIn.PosOutgoing1.size());++i)
      PosOutgoing1.push_back( hardProcessIn.PosOutgoing1[i]);
    for(int i =0; i < int(hardProcessIn.PosOutgoing2.size());++i)
      PosOutgoing2.push_back( hardProcessIn.PosOutgoing2[i]);
    for(int i =0; i < int(hardProcessIn.PosIntermediate.size());++i)
      PosIntermediate.push_back( hardProcessIn.PosIntermediate[i]);
  }

  // Constructor with path to LHE file
  HardProcess( string LHEfile, ParticleData* particleData) {
    state = Event();
    state.init("(hard process)", particleData);
    translateLHEFString(LHEfile);
  }

  // Constructor with core process input
  void initOnProcess( string process, ParticleData* particleData);

  // Constructor with path to LHE file input
  void initOnLHEF( string LHEfile, ParticleData* particleData);

  // Function to access the LHE file and read relevant information
  void translateLHEFString( string LHEpath);

  // Function to translate the process string (in MG/ME notation)
  void translateProcessString( string process);

  // Function to clear hard process information
  void clear();

  // Function to check whether the sets of candidates Pos1, Pos2, together
  // with the proposed candidate iPos give an allowed hard process state
  bool allowCandidates(int iPos, vector<int> Pos1, vector<int> Pos2,
    const Event& event);
  // Function to identify the hard subprocess in the current event
  void storeCandidates( const Event& event, string process);
  // Function to check if the particle event[iPos] matches any of
  // the stored outgoing particles of the hard subprocess
  bool matchesAnyOutgoing(int iPos, const Event& event);

  // Function to return the type of the ME generator
  int MEgenType();
  // Function to get the number of coloured final state partons in the
  // hard process
  int nQuarksOut();
  // Function to get the number of uncoloured final state particles in the 
  // hard process
  int nLeptonOut();
  // Function to get the number of electroweak final state bosons in the 
  // hard process
  int nBosonsOut();

  // Function to get the number of coloured initial state partons in the 
  // hard process
  int nQuarksIn();
  // Function to get the number of uncoloured initial state particles in the 
  // hard process
  int nLeptonIn();
  // Function to report if a resonace decay was found in the 2->2 sub-process 
  // of the  current state
  int hasResInCurrent();
  // Function to report the number of resonace decays in the 2->2 sub-process 
  // of the  current state
  int nResInCurrent();
  // Function to report if a resonace decay was found in the 2->2 hard process
  bool hasResInProc();
  // Function to print the hard process (for debug)
  void list() const;
  // Function to print the hard process candidates in the
  // Matrix element state (for debug)
  void listCandidates() const;

};

//==========================================================================

// MergingHooks is base class for user input to the merging procedure.

class MergingHooks {

public:

  // Constructor.
  MergingHooks() : 
    doUserMergingSave(false),
    doMGMergingSave(false),
    doKTMergingSave(false),
    doPTLundMergingSave(false),
    doCutBasedMergingSave(false),
    forceOrderedSave(false),
    forceUnorderedSave(false),
    doOrderHistoriesSave(true),
    doForceUnorderedHistoriesSave(false),
    doCutOnRecStateSave(false),
    doWClusteringSave(false) {}

  // Destructor.
  virtual ~MergingHooks() {}
  // Function encoding the functional definition of the merging scale
  virtual double tmsDefinition( const Event& event){ return event[0].e();}

  // Function returning the value of the merging scale.
  double tms() {
    if(doCutBasedMergingSave) return 0.;
    else return tmsValueSave;
  }
  // Function returning the value of the Delta R_{ij} cut for
  // cut based merging scale definition.
  double dRijMS() {
    return ((tmsListSave.size() == 3) ? tmsListSave[0] : 0.);
  }
  // Function returning the value of the pT_{i} cut for
  // cut based merging scale definition.
  double pTiMS() {
    return ((tmsListSave.size() == 3) ? tmsListSave[1] : 0.);
  }
  // Function returning the value of the pT_{i} cut for
  // cut based merging scale definition.
  double QijMS() {
    return ((tmsListSave.size() == 3) ? tmsListSave[2] : 0.);
  }
  // Function returning the value of the maximal number of merged jets.
  int nMaxJets() { return nJetMaxSave;}

  // Function to return hard process string.
  string getProcessString() { return processSave;}

  // Function to return the number of outgoing partons in the core process
  int nHardOutPartons(){ return hardProcess.nQuarksOut();}
  // Function to return the number of outgoing leptons in the core process
  int nHardOutLeptons(){ return hardProcess.nLeptonOut();}
  // Function to return the number of outgoing electroweak bosons in the core
  // process.
  int nHardOutBosons(){ return hardProcess.nBosonsOut();}

  // Function to return the number of incoming partons (hadrons) in the core
  // process.
  int nHardInPartons(){ return hardProcess.nQuarksIn();}
  // Function to return the number of incoming leptons in the core process.
  int nHardInLeptons(){ return hardProcess.nLeptonIn();}
  // Function to report the number of resonace decays in the 2->2 sub-process 
  // of the  current state.
  int nResInCurrent(){ return hardProcess.nResInCurrent();}
  // Function to check if event contains an emission not present in the hard
  // process.
  bool isFirstEmission(const Event& event);
  // Function to allow effective gg -> EW boson couplings.
  bool hasEffectiveG2EW() {
    if ( getProcessString().compare("pp>h") == 0 ) return true;
    return false; }

  // Function to return the number of clustering steps for the current event
  int getNumberOfClusteringSteps(const Event& event);

  // Function to determine if user defined merging should be applied.
  bool doUserMerging(){ return doUserMergingSave;}
  // Function to determine if automated MG/ME merging should be applied.
  bool doMGMerging() { return doMGMergingSave;}
  // Function to determine if KT merging should be applied.
  bool doKTMerging() { return doKTMergingSave;}
  // Function to determine if PTLund merging should be applied.
  bool doPTLundMerging() { return doPTLundMergingSave;}
  // Function to determine if cut based merging should be applied.
  bool doCutBasedMerging() { return doCutBasedMergingSave;}

  // Function to dampen weights calculated from histories with lowest 
  // multiplicity reclustered events that do not pass the ME cuts
  virtual double dampenIfFailCuts( const Event& inEvent ) {
    // Dummy statement to avoid compiler warnings
    if(false) cout << inEvent[0].e();
    return 1.;
  }

  // Hooks to disallow states in the construction of all histories, e.g.
  // because jets are below the merging scale or fail the matrix element cuts
  // Function to allow interference in the construction of histories 
  virtual bool canCutOnRecState() { return doCutOnRecStateSave; }
  // Function to check reclustered state while generating all possible
  // histories
  // Function implementing check of reclustered events while constructing
  // all possible histories
  virtual bool doCutOnRecState( const Event& event ) {
    // Dummy statement to avoid compiler warnings.
    if(false) cout << event[0].e();
    // Count number of final state partons.
    int nPartons = 0;
    for( int i=0; i < int(event.size()); ++i)
      if(  event[i].isFinal()
        && (event[i].isGluon() || event[i].isQuark()) )
        nPartons++;
    // For gg -> h, allow only histories with gluons in initial state
    if( hasEffectiveG2EW() && nPartons < 2){
      if(event[3].id() != 21 || event[4].id() != 21)
        return true;
    }
    return false;
  }

  // Function to allow not counting a trial emission.
  virtual bool canVetoTrialEmission() { return false;}
  // Function to check if trial emission should be rejected.
  virtual bool doVetoTrialEmission( const Event&, const Event& ) {
    return false; }

  // Make History class friend to allow access to advanced switches
  friend class History;
  // Make Pythia class friend
  friend class Pythia;
  // Make PartonLevel class friend
  friend class PartonLevel;
  // Make SpaceShower class friend
  friend class SpaceShower;
  // Make TimeShower class friend
  friend class TimeShower;

  // Function to force preferred picking of ordered histories. By default,
  // unordered histories will only be considered if no ordered histories
  // were found. 
  void orderHistories( bool doOrderHistoriesIn) {
    doOrderHistoriesSave = doOrderHistoriesIn; }
  // Function to force preferred picking of unordered histories.
  void forceUnorderedHistories( bool doForceUnorderedHistoriesIn) {
    doForceUnorderedHistoriesSave = doForceUnorderedHistoriesIn; }
  // Function to force cut on reconstructed states internally, as needed
  // for gg -> Higgs to ensure that e.g. uu~ -> Higgs is not constructed.
  void allowCutOnRecState( bool doCutOnRecStateIn) {
    doCutOnRecStateSave = doCutOnRecStateIn; }

  // Function to allow final state clusterings of W-bosons
  void doWClustering( bool doWClusteringIn ) {
    doWClusteringSave = doWClusteringIn; }

protected:

  // Functions for internal use inside Pythia source code
  // Initialize.
  void init(  Settings& settings, Info* infoPtrIn, 
    ParticleData* particleDataPtrIn, ostream& os = cout);
  // Function storing candidates for the hard process in the current event
  // Needed in order not to cluster members of the core process
  void storeHardProcessCandidates(const Event& event){
    hardProcess.storeCandidates(event,getProcessString());
  }
  // Function to set the path to the LHE file, so that more automated merging
  // can be used.
  // Remove "_1.lhe" suffix from LHE file name.
  // This is done so that the HarsProcess class can access both the +0 and +1
  // LHE files to find both the merging scale and the core process string
  // Store.
  void setLHEInputFile( string lheFile) {
    lheInputFile = lheFile.substr(0,lheFile.size()-6); }

  // Function to save the current CKKW-L weight
  void setWeight(double wgt){ weightSave = wgt;}

  // Save CKKW-L weight / O(\alpha_s) weight.
  double weightCKKWLSave;
  // Return CKKW-L weight.
  double getWeightCKKWL() { return weightCKKWLSave; }
  // Set CKKW-L weight.
  void setWeightCKKWL( double weightIn) { weightCKKWLSave = weightIn; }

  // Function to calculate the minimal kT in the event
  double kTms(const Event & event);

  // Find the minimal Lund pT between coloured partons in the event
  double rhoms( const Event& event, bool withColour);

  // Function to check if the properties of the input particle should be
  // checked against the cut-based merging scale defintion.
  bool checkAgainstCut( const Particle& particle);

  // Find the if the event passes the Delta R_{ij}, pT_{i} and Q_{ij} cuts on
  // the matrix element, as needed for cut-based merging scale definition
  double cutbasedms( const Event& event );

  // Function to compute Delta R separation from 4-vector input
  double deltaRij(Vec4 jet1, Vec4 jet2);

  // Function to decide if (too) many events are significantly above the
  // merging scale.
  bool stats() { double ALLOWEDFRACTION = 0.75;
    double fraction = double(infoPtr->getCounter(41))
                    / max(1., double(infoPtr->nAccepted()));
    return ( fraction > ALLOWEDFRACTION ); } 

  // Function to get the CKKW-L weight for the current event
  double getWeight() { return weightSave;}

  // Helper class doing all the core process book-keeping
  HardProcess hardProcess;

  // Pointer to various information on the generation.
  Info*          infoPtr;

  // Pointer to the particle data table.
  ParticleData*  particleDataPtr;

  // AlphaS objects for alphaS reweighting use
  AlphaStrong AlphaS_FSRSave;
  AlphaStrong AlphaS_ISRSave;
  AlphaEM AlphaEM_FSRSave;

  // Return AlphaStrong objects
  AlphaStrong* AlphaS_FSR() { return &AlphaS_FSRSave;}
  AlphaStrong* AlphaS_ISR() { return &AlphaS_ISRSave;}
  AlphaEM* AlphaEM_FSR() { return &AlphaEM_FSRSave;}

  // Saved path to LHE file for more automated merging
  string lheInputFile;

  bool   doUserMergingSave, doMGMergingSave, doKTMergingSave,
         doPTLundMergingSave, doCutBasedMergingSave,
         includeMassiveSave, enforceStrongOrderingSave, orderInRapiditySave,
         pickByFullPSave, pickByPoPT2Save, includeRedundantSave,
         pickBySumPTSave, allowColourShufflingSave, resetHardQRenSave,
         resetHardQFacSave;
  int    unorderedScalePrescipSave, unorderedASscalePrescipSave,
         unorderedPDFscalePrescipSave, incompleteScalePrescipSave,
         ktTypeSave;
  double scaleSeparationFactorSave, nonJoinedNormSave,
         fsrInRecNormSave, herwigAcollFSRSave, herwigAcollISRSave,
         pT0ISRSave, pTcutSave;
  bool   forceOrderedSave, forceUnorderedSave;

  // Functions to return advanced merging switches
  // Include masses in definition of evolution pT and splitting kernels
  bool includeMassive() { return includeMassiveSave;}
  // Prefer strongly ordered histories
  bool enforceStrongOrdering() { return enforceStrongOrderingSave;}
  // Prefer histories ordered in rapidity and evolution pT
  bool orderInRapidity() { return orderInRapiditySave;}
  // Pick history probabilistically by full (correct) splitting probabilities
  bool pickByFull() { return pickByFullPSave;}
  // Pick history probabilistically, with easier form of probabilities
  bool pickByPoPT2() { return pickByPoPT2Save;}
  // Include redundant terms (e.g. PDF ratios) in the splitting probabilities
  bool includeRedundant(){ return includeRedundantSave;}
  // Pick by winner-takes-it-all, with minimum sum of scalar evolution pT
  bool pickBySumPT(){ return pickBySumPTSave;}

  bool forceOrdered() { return forceOrderedSave;}
  bool forceUnordered() { return forceUnorderedSave;}

  // Prescription for combined scale of unordered emissions
  // 0 : use larger scale
  // 1 : use smaller scale
  int unorderedScalePrescip() { return unorderedScalePrescipSave;}
  // Prescription for combined scale used in alpha_s for unordered emissions
  // 0 : use combined emission scale in alpha_s weight for both (!) splittings
  // 1 : use original reclustered scales of each emission in alpha_s weight
  int unorderedASscalePrescip() { return unorderedASscalePrescipSave;}
  // Prescription for combined scale used in PDF ratios for unordered
  // emissions
  // 0 : use combined emission scale in PDFs for both (!) splittings
  // 1 : use original reclustered scales of each emission in PDF ratiost
  int unorderedPDFscalePrescip() { return unorderedPDFscalePrescipSave;}
  // Prescription for starting scale of incomplete histories
  // 0: use factorization scale
  // 1: use sHat
  // 2: use s
  int incompleteScalePrescip() { return incompleteScalePrescipSave;}

  // Allow swapping one colour index while reclustering
  bool allowColourShuffling() { return allowColourShufflingSave;}

  // Allow use of dynamical renormalisation scale of the core 2-> 2 process.
  bool resetHardQRen() { return resetHardQRenSave; }
  // Allow use of dynamical factorisation scale of the core 2-> 2 process.
  bool resetHardQFac() { return resetHardQFacSave; }

  // Factor by which two scales should differ to be classified strongly
  // ordered.
  double scaleSeparationFactor() { return scaleSeparationFactorSave;}
  // Absolute normalization of splitting probability for non-joined
  // evolution.
  double nonJoinedNorm() { return nonJoinedNormSave;}
  // Absolute normalization of splitting probability for final state
  // splittings with initial state recoiler
  double fsrInRecNorm() { return fsrInRecNormSave;}
  // Factor multiplying scalar evolution pT for FSR splitting, when picking 
  // history by minimum scalar pT (see Jonathan Tully's thesis)
  double herwigAcollFSR() { return herwigAcollFSRSave;}
  // Factor multiplying scalar evolution pT for ISR splitting, when picking 
  // history by minimum scalar pT (see Jonathan Tully's thesis)
  double herwigAcollISR() { return herwigAcollISRSave;}
  // ISR regularisation scale
  double pT0ISR() { return pT0ISRSave;}
  // Shower cut-off scale
  double pTcut() { return pTcutSave;}

  // Function to calculate the kT separation between two particles
  double kTdurham(const Particle& RadAfterBranch,
    const Particle& EmtAfterBranch, int Type, double D );
  // Function to compute "pythia pT separation" from Particle input
  double rhoPythia(const Particle& RadAfterBranch,
    const Particle& EmtAfterBranch, const Particle& RecAfterBranch, 
    int ShowerType);
  // Function to find a colour (anticolour) index in the input event,
  // used to find colour-connected recoilers
  int findColour(int col, int iExclude1, int iExclude2,
    const Event& event, int type, bool isHardIn);

  // Saved members.
  double tmsValueSave;
  int nJetMaxSave;
  string processSave;
  double weightSave;

  // List of cut values to used to define a merging scale. Ordering:
  // 0: DeltaR_{jet_i,jet_j,min}
  // 1: p_{T,jet_i,min}
  // 2: Q_{jet_i,jet_j,min}
  vector<double> tmsListSave;

  // INTERNAL Hooks to allow construction of all histories,
  // including un-ordered ones
  bool doOrderHistoriesSave, doForceUnorderedHistoriesSave;
  bool orderHistories() { return doOrderHistoriesSave;}

  // INTERNAL Hooks to force construction of only unordered histories,
  bool forceUnorderedHistories() { return doForceUnorderedHistoriesSave;}

  // INTERNAL Hooks to disallow states in the construction of all histories,
  // e.g. because jets are below the merging scale, of to avoid the
  // construction of uu~ -> Higgs histories.
  bool doCutOnRecStateSave;
  bool allowCutOnRecState() { return doCutOnRecStateSave;}

  // INTERNAL Hooks to allow clustering W bosons.
  bool doWClusteringSave;
  bool doWClustering() { return doWClusteringSave;}

};

//==========================================================================

} // end namespace Pythia8

#endif // Pythia8_MergingHooks_H
