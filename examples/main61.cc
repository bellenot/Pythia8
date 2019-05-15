// main61.cc is a part of the PYTHIA event generator.
// Copyright (C) 2008 Richard Corke.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

/*
 * Simple example of fastjet analysis. Roughly follows analysis of:
 * T. Aaltonen et al. [CDF Collaboration],
 * Measurement of the cross section for W-boson production in association
 * with jets in ppbar collisions at sqrt(s)=1.96$ TeV
 * Phys. Rev. D 77 (2008) 011108
 * arXiv:0711.4044 [hep-ex]
 *
 * Cuts:
 *   ET(elec)     > 20GeV
 *   |eta(elec)|  < 1.1
 *   ET(missing)  > 30GeV
 *   ET(jet)      > 20GeV
 *   |eta(jet)|   < 2.0
 *   deltaR(elec, jet) > 0.52
 * Not used:
 *   mT(W)        > 20GeV
 */

#include "Pythia.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

using namespace Pythia8;

// sigma(W -> ev + >= n-jet; ET(n'th-jet) > 25GeV), n = 0, 1, 2, 3, 4
const double expCrossSec[] = { 798.0, 53.5, 6.8, 0.84, 0.074 };

int main() {
  // Settings
  int  nEvent = 1000;
  bool doMI = true;
    
  // Generator
  Pythia pythia;

  // Single W production
  pythia.readString("WeakSingleBoson:ffbar2W = on");
  // Force decay W->ev
  pythia.readString("24:onMode = off");
  pythia.readString("24:onIfAny = 11 12");
  // Multiple Interactions
  if (doMI == false)
    pythia.readString("PartonLevel:MI = off");

  // Initialisation
  pythia.init( 2212, -2212, 1960.);

  // Histograms
  Hist dSigma1("1-jet cross-section (E_jet1 > 20 GeV)", 66, 20.0, 350.0);
  Hist dSigma2("2-jet cross-section (E_jet2 > 20 GeV)", 34, 20.0, 190.0);
  Hist dSigma3("3-jet cross-section (E_jet3 > 20 GeV)", 12, 20.0, 80.0);
  Hist dSigma4("4-jet cross-section (E_jet4 > 20 GeV)",  3, 20.0, 35.0);
  Hist *dSigmaHist[5] = { NULL, &dSigma1, &dSigma2, &dSigma3, &dSigma4 };

  // Fastjet analysis - select algorithm and parameters
  double Rparam = 0.4;
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recomb_scheme = fastjet::E_scheme;
  fastjet::JetDefinition         *jet_def = NULL;
  jet_def = new fastjet::JetDefinition(fastjet::kt_algorithm, Rparam,
                                       recomb_scheme, strategy);

  // Fastjet input
  std::vector <fastjet::PseudoJet> fjInputs;

  int nEventAccept25[5] = { 0, 0, 0, 0, 0 };
  int vetoCount[4] = { 0, 0, 0, 0 };
  const char *vetoStr[] = { "ET(elec)", "|eta(elec)|", "ET(missing)", 
    "deltaR(elec, jet)" };
  bool firstEvent = true;

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    if (iEvent < 1) {pythia.info.list(); pythia.event.list();} 
    if (iEvent % 10000 == 0) printf("Event %d\n", iEvent);

    // Reset Fastjet input
    fjInputs.resize(0);

    // Cuts
    double missingET = 0.0;
    Vec4   missingETvec = 0.0;
    bool   vetoEvent = false;
    int    haveElectron = 0, haveNeutrino = 0, haveW = 0;

    for (int i = 0; i < pythia.event.size(); ++i) {
      // Keep track of the W
      if (pythia.event[i].idAbs() == 24)
        haveW = i;

      // Final state only
      if (!pythia.event[i].isFinal()) continue;

      // Check for W decay products
      if (pythia.event[i].status() == 23 || pythia.event[i].status() == 51
       || pythia.event[i].status() == 52) {

        // Electron cuts
        if (pythia.event[i].idAbs() == 11) {
          if (pythia.event[i].pT() <= 20.0) {
            vetoCount[0]++;
            vetoEvent = true;
            break;
          }
          if (fabs(pythia.event[i].eta()) >= 1.1) {
            vetoCount[1]++;
            vetoEvent = true;
            break;
          }
          missingETvec += pythia.event[i].p();
          haveElectron = i;

        // Neutrino cuts
        } else if (pythia.event[i].idAbs() == 12) {
          haveNeutrino = i;

        }

        // Don't include decay products for jet algorithm
        continue;
      }

      // Only |eta| < 3.6
      if (fabs(pythia.event[i].eta()) >= 3.6) continue;

      // No neutrinos
      if (pythia.event[i].id() == 12 || pythia.event[i].id() == 14 ||
          pythia.event[i].id() == 16) continue;

      // Missing ET
      missingETvec += pythia.event[i].p();

      // Store as input to Fastjet
      fjInputs.push_back( fastjet::PseudoJet (pythia.event[i].px(),
						pythia.event[i].py(), pythia.event[i].pz(),
						pythia.event[i].e() ) );
    }

    // Check if event meets cuts
    if (vetoEvent) continue;
    if (haveElectron == 0 || haveNeutrino == 0) {
      pythia.event.list();
      cout << "Error: Couldn't find electron and neutrino decay products. "
           << "Exiting.." << endl;
      return 1;
    }
    missingET = missingETvec.pT();
    if (missingET <= 30.0) { vetoCount[2]++; continue; }

    // Run Fastjet algorithm
    vector <fastjet::PseudoJet> inclusive_jets, sorted_jets;
    if (fjInputs.size() != 0) {
      fastjet::ClusterSequence clust_seq(fjInputs, *jet_def);

      if (firstEvent) {
        cout << "Ran " << jet_def->description() << endl;
        cout << "Strategy adopted by FastJet was "
             << clust_seq.strategy_string() << endl << endl;
        firstEvent = false;
      }

      // Extract inclusive jets sorted by pT
      inclusive_jets = clust_seq.inclusive_jets();
      sorted_jets    = sorted_by_pt(inclusive_jets);  
    }

    int jetCount20 = 0, jetCount25 = 0;
    for (unsigned int i = 0; i < sorted_jets.size(); i++) {
      // Check deltaR between W decay electron and jets
      fastjet::PseudoJet fjElec(pythia.event[haveElectron].px(),
                                pythia.event[haveElectron].py(),
                                pythia.event[haveElectron].pz(),
                                pythia.event[haveElectron].e());
      double deltaPhi = fjElec.phi() - sorted_jets[i].phi();
      double deltaEta = fjElec.eta() - sorted_jets[i].eta();
      double deltaR = sqrt(deltaPhi * deltaPhi + deltaEta * deltaEta);
      if (deltaR <= 0.52) { vetoEvent = true; break; }

      // Jet cut of |eta| < 2.0
      if (fabs(sorted_jets[i].rap()) >= 2.0) continue;

      // Fill dSigma histograms and count jets with ET > 25.0
      if (sorted_jets[i].perp() > 20.0) {
        if (sorted_jets[i].perp() > 25.0)
          jetCount25++;

        if (jetCount20 <= 3)
          dSigmaHist[++jetCount20]->fill(sorted_jets[i].perp());
      } else
        break;

    }
    if (vetoEvent) { vetoCount[3]++; continue; }

    if (jetCount25 > 4) jetCount25 = 4;
    for (int i = jetCount25; i >= 0; i--)
      nEventAccept25[i]++;

  // End of event loop.
  }

  // Statistics
  pythia.statistics();

  // Output histograms
  double sigmapb = pythia.info.sigmaGen() * 1.0E9;

  for (int i = 1; i <= 4; i++)
    (*dSigmaHist[i]) = ((*dSigmaHist[i]) * sigmapb) / nEvent;
  cout << dSigma1 << dSigma2 << dSigma3 << dSigma4 << endl;

  // Output cross-sections
  cout << "Jet algorithm is kT" << endl;
  cout << "Multiple interactions are switched "
       << ( (doMI) ? "on" : "off" ) << endl;
  cout << endl << nEvent << " events generated. " << nEventAccept25[0]
       << " events passed cuts." << endl;
  cout << "Vetos:" << endl;
  for (int i = 0; i < 4; i++)
    cout << "  " << vetoStr[i] << " = " << vetoCount[i] << endl;

  cout << endl << "Inclusive cross-sections (pb):" << endl;
  for (int i = 0; i < 5; i++) {
    cout << scientific << setprecision(3)
         << "  " << i << "-jet - Pythia = "
         << ((double) nEventAccept25[i] / (double) nEvent) * sigmapb;
    cout << ", Experimental = " << expCrossSec[i];
    if (i != 0) {
      cout << scientific << setprecision(3)
           << ", Pythia ratio to " << i - 1 << "-jet = "
           << ((double) nEventAccept25[i] / (double) nEventAccept25[i - 1]);
      cout << scientific << setprecision(3)
           << ", Experimental ratio to " << i - 1 << "-jet = "
           << expCrossSec[i] / expCrossSec[i - 1];
    }
    cout << endl;
  }

  return 0;
}

