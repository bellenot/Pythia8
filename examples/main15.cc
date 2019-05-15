// main15.cc is a part of the PYTHIA event generator.
// Copyright (C) 2007 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. 
// It illustrates how B decays could be repeated a number of times 
// for each event to improve statistics when this could be a problem.

#include "Pythia.h"
using namespace Pythia8;
 
int main() {

  // Number of events, generated and listed ones.
  int nEvent = 1000;
  int nList = 1;

  // Number of times B decays should be redone for each event.
  int nRepeat = 10; 

  // List of weakly decaying B hadrons.
  // Note: this list is overkill; some will never be produced.
  int bCodes[28] = {511, 521, 531, 541, 5122, 5132, 5142, 5232, 5242,
    5332, 5342, 5412, 5414, 5422, 5424, 5432, 5434, 5442, 5444, 5512,
    5514, 5522, 5524, 5532, 5534, 5542, 5544, 5544 }; 
  int nCodes = 28;

  // Generator. Shorthand for event.
  Pythia pythia;
  Event& event = pythia.event;

  // Simulate b production above given pTmin scale.
  // Warning: these processes do not catch all possible production modes.
  // You would need to use HardQCD:all or even SoftQCD:minBias for that.
  pythia.readString("HardQCD:gg2bbbar = on");    
  pythia.readString("HardQCD:qqbar2bbbar = on");    
  pythia.readString("PhaseSpace:pTHatMin = 50.");  

  // Initialize for LHC energies.    
  pythia.init( 2212, 2212, 14000.);

  // Histogram invariant mass of muon pairs.
  Hist nBperEvent("number of B hadrons in an event", 10, -0.5, 9.5); 
  Hist nSameEvent("number of times same event is used", 10, -0.5, 9.5); 
  Hist oppSignMass("mass of opposite-sign muon pair", 100, 0.0, 100.0);
  Hist sameSignMass("mass of same-sign muon pair", 100, 0.0, 100.0);

  // Begin event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Switch off decays of weakly decaying B hadrons.
    // (More compact solution than repeated readString(..).)
    for (int iC = 0; iC < nCodes; ++iC) 
      ParticleDataTable::mayDecay( bCodes[iC], false);

    // Generate event. Skip it if error.
    if (!pythia.next()) continue;

    // List first few events.
    if (iEvent < nList) {
      pythia.info.list(); 
      pythia.process.list();
      pythia.event.list();
    }

    // Find all locations where B hadrons are stored.
    vector<int> iBLoc;
    for (int i = 0; i < event.size(); ++i) {
      int idAbs = event[i].idAbs();
      for (int iC = 0; iC < 28; ++iC) 
      if (idAbs == bCodes[iC]) {       
        iBLoc.push_back(i);
        break;
      }
    }

    // Histogram number of B's in event.
    int nB = iBLoc.size();
    nBperEvent.fill( nB );

    // Store size of current event (also as local copy).
    event.saveSize();
    int saveSize = event.size();

    // Switch back on decays of weakly decaying B hadrons.
    // (More compact solution than repeated readString(..).)
    for (int iC = 0; iC < nCodes; ++iC) 
      ParticleDataTable::mayDecay( bCodes[iC], true);

    // Begin loop over rounds of B decay tries for same event.
    int nWithPair = 0;
    for (int iRepeat = 0; iRepeat < nRepeat; ++iRepeat) {

      // Remove B decay products from previous round. (Added at end!)
      if (iRepeat > 0) {
        event.restoreSize();

        // Mark decayed B hadrons as undecayed.
        for (int iB = 0; iB < nB; ++iB) event[ iBLoc[iB] ].statusPos(); 
      } 
  
      // Do decays of B hadrons, sequentially for products. 
      // Note: this routine does not work for bottomonium (or heavier) states,
      // since there decays like Upsilon -> g g g also need hadronization.
      // Also, there is no provision for Bose-Einstein effects.
      pythia.moreDecays();
  
      // Look for muons among B decay products (also from charm/tau).
      vector<int> iMuNeg, iMuPos;
      for (int i = saveSize; i < event.size(); ++i) {
        int id = event[i].id();  
        if (id ==  13) iMuNeg.push_back(i);
        if (id == -13) iMuPos.push_back(i);
      }
      
      // Check whether pair(s) present.
      int nMuNeg = iMuNeg.size();
      int nMuPos = iMuPos.size();
      if (nMuNeg + nMuPos > 1) {
        ++nWithPair;

        // Fill masses of opposite-sign pairs.
        for (int iN = 0; iN < nMuNeg; ++iN)
        for (int iP = 0; iP < nMuPos; ++iP) 
          oppSignMass.fill(
            (event[iMuNeg[iN]].p() + event[iMuPos[iP]].p()).mCalc() );

        // Fill masses of same-sign pairs.
        for (int i1 = 0; i1 < nMuNeg - 1; ++i1)
        for (int i2 = i1 + 1; i2 < nMuNeg; ++i2) 
          sameSignMass.fill(
            (event[iMuNeg[i1]].p() + event[iMuNeg[i2]].p()).mCalc() );
        for (int i1 = 0; i1 < nMuPos - 1; ++i1)
        for (int i2 = i1 + 1; i2 < nMuPos; ++i2) 
          sameSignMass.fill(
            (event[iMuPos[i1]].p() + event[iMuPos[i2]].p()).mCalc() );

      // Finished analysis of current round. 
      }

    // End of loop over many rounds. fill number of rounds with pairs.
    }
    nSameEvent.fill( nWithPair );

  // End of event loop.
  }

  // Statistics. Histograms. 
  pythia.statistics();
  cout << nBperEvent << nSameEvent << oppSignMass << sameSignMass << endl;

  // Done. 
  return 0;
}
