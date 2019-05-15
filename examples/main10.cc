// File: main10.cc
// This is a simple test program. 
// It studies jet production at the LHC, using CellJet.
// Copyright C 2006 Torbjorn Sjostrand

#include "Pythia.h"
using namespace Pythia8;
 
int main() {

  // Generator. Process selection. LHC initialization.
  Pythia pythia;
  pythia.readString("Pythia6:ckin(3) = 200.");    
  pythia.init( 2212, 2212, 14000.);

  // Jet finder.
  CellJet cellJet;
  // CellJet cellJet(20., 0.7, 2, 5., 100, 100);

  // Histograms.
  Hist nJets("number of jets",20,-0.5,19.5);
  Hist eTjets("eT for jets",100,0.,500.);
  Hist etaJets("eta for jets",100,-5.,5.);
  Hist phiJets("phi for jets",100,-M_PI,M_PI);  
  Hist distJets("R distance between jets",100,0.,10.);
  Hist eTdiff("eT difference",100,-100.,400.);

  // Begin event loop. Generate event. Skip if error. list first few. 
  for (int iEvent = 0; iEvent < 100; ++iEvent) {
    if (!pythia.next()) continue;
    if (iEvent < 1) pythia.event.list();

    // Analyze jet properties. 
    cellJet. analyze( pythia.event );
    if (iEvent < 3) cellJet.list();
    nJets.fill( cellJet.size() );
    for (int i = 0; i < cellJet.size(); ++i) {
      eTjets.fill( cellJet.eT(i) );
      etaJets.fill( cellJet.etaWeighted(i) );
      phiJets.fill( cellJet.phiWeighted(i) );
    }
    for (int i = 0; i < cellJet.size() - 1; ++i)
    for (int j = i +1; j < cellJet.size(); ++j) {
      double dEta = cellJet.etaWeighted(i) 
        - cellJet.etaWeighted(j);
      double dPhi = abs( cellJet.phiWeighted(i) 
        - cellJet.phiWeighted(j) );
      if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
      double dR = sqrt( pow2(dEta) + pow2(dPhi) );
      distJets.fill( dR );
    }
    for (int i = 1; i < cellJet.size(); ++i) 
      eTdiff.fill( cellJet.eT(i-1)- cellJet.eT(i) );

  // End of event loop. Statistics. Histogram. 
  }
  pythia.statistics();
  cout << nJets << eTjets << etaJets << phiJets 
       << distJets << eTdiff;

  // Done. 
  return 0;
}
