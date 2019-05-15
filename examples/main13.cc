// main13.cc is a part of the PYTHIA event generator.
// Copyright (C) 2007 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Comparison with Pythia6.408 cross sections process by process.
// Note: processes g g / q qbar -> H0 t tbar and g g / q qbar -> H0 b bbar 
// are so slow that they have been left out to keep reasonable execution time.
// Subruns  0 -  5 : QCD jets
//          6 - 10 : prompt photons.
//         11 - 12 : t-channel gamma/Z/W exchange.
//         13 - 23 : gamma*/Z^0/W^+-, singly, in pairs or with parton
//         24 - 25 : onia.
//         26 - 30 : top.
//         31 - 40 : Standard Model Higgs.

#include "Pythia.h"

using namespace Pythia8;
 
int main() {

  // First and last process to test: can run from 0 through 40.
  int iFirst = 0;
  int iLast  = 40;   

  // Statistics. Pythia6 run was with 10000, so no point to use more.
  int nEvent = 10000; 

  // Normally one subprocess word per subrun, but exceptions exist. 
  int nSub[41];
  for (int i = 0; i < 41; ++i) nSub[i] = 1; 
  nSub[1] = 3;  
  nSub[5] = 3;  
  // Starting positions in subprocess words list, recursively defined.
  int iBeg[42] = { 0 };
  for (int i = 0; i < 41; ++i) iBeg[i + 1] = iBeg[i] + nSub[i]; 
  // List of subprocess words. 
  string processes[45] = { "HardQCD:gg2gg", "HardQCD:gg2qqbar",
    "HardQCD:gg2ccbar",  "HardQCD:gg2bbbar","HardQCD:qg2qg" , 
    "HardQCD:qq2qq", "HardQCD:qqbar2gg", "HardQCD:qqbar2qqbarNew", 
    "HardQCD:qqbar2ccbar", "HardQCD:qqbar2bbbar", "PromptPhoton:qg2qgamma",
    "PromptPhoton:qqbar2ggamma", "PromptPhoton:gg2ggamma", 
    "PromptPhoton:ffbar2gammagamma", "PromptPhoton:gg2gammagamma",
    "WeakBosonExchange:ff2ff(t:gmZ)", "WeakBosonExchange:ff2ff(t:W)",
    "WeakSingleBoson:ffbar2gmZ", "WeakSingleBoson:ffbar2W", 
    "WeakDoubleBoson:ffbar2gmZgmZ", "WeakDoubleBoson:ffbar2ZW",
    "WeakDoubleBoson:ffbar2WW", "WeakBosonAndParton:qqbar2gmZg",
    "WeakBosonAndParton:qg2gmZq", "WeakBosonAndParton:ffbar2gmZgm",
    "WeakBosonAndParton:qqbar2Wg", "WeakBosonAndParton:qg2Wq",
    "WeakBosonAndParton:ffbar2Wgm", "Charmonium:all", "Bottomonium:all",
    "Top:gg2ttbar",  "Top:qqbar2ttbar",  "Top:qq2tq(t:W)",    
    "Top:ffbar2ttbar(s:gmZ)", "Top:ffbar2tqbar(s:W)",
    "HiggsSM:ffbar2H", "HiggsSM:gg2H", "HiggsSM:ffbar2HZ", 
    "HiggsSM:ffbar2HW", "HiggsSM:ff2Hff(t:ZZ)", "HiggsSM:ff2Hff(t:WW)",
    "HiggsSM:qg2Hq", "HiggsSM:gg2Hg(l:t)", "HiggsSM:qg2Hq(l:t)", 
    "HiggsSM:qqbar2Hg(l:t)" }; 

  // List of cross sections from Pythia6.
  double sigma6[41] = {   4.960e-01, 1.627e-02, 2.790e-01, 2.800e-02, 
    3.310e-04, 3.653e-04, 1.697e-04, 1.163e-05, 1.065e-07, 8.259e-08, 
    8.237e-08, 2.544e-05, 5.321e-06, 5.571e-05, 1.621e-04, 9.039e-09,
    2.247e-08, 5.893e-08, 3.781e-06, 1.078e-05, 4.551e-08, 1.025e-05,
    3.208e-05, 5.435e-08, 1.038e-04, 3.929e-05, 4.155e-07, 6.685e-08,
    1.898e-07, 4.240e-10, 7.142e-09, 1.547e-10, 7.064e-09, 1.316e-10,
    2.332e-10, 5.105e-10, 1.316e-09, 4.462e-11, 5.557e-09, 1.966e-09,
    8.725e-12};

  // Generator. 
  Pythia pythia;

  // Standard set of masses for comparison with Fortran code.
  pythia.readString("5:m0  = 4.2");
  pythia.readString("6:m0  = 175.");
  pythia.readString("23:m0 = 91.2");
  pythia.readString("24:m0 = 80.");
  pythia.readString("25:m0 = 200.");

  // Same kinematics cuts as Fortran code.
  pythia.readString("PhaseSpace:pTHatMin = 20.");
  pythia.readString("6:mMin = 20.");
  pythia.readString("23:mMin = 20.");
  pythia.readString("24:mMin = 20.");
  pythia.readString("25:mMin = 20.");

  // Also same renormalization and factorization scale.
  pythia.readString("SigmaProcess:renormScale2 = 3");
  pythia.readString("SigmaProcess:factorScale2 = 3");

  // Switch off unnecessary parts.
  pythia.readString("PartonLevel:all = off");

  // Debug.
  //pythia.readString("PhaseSpace:showSearch = on");
  //pythia.readString("PhaseSpace:showViolation = on");

  // Loop over processes.
  for (int iProc = iFirst; iProc <= iLast; ++iProc) {
    cout << "\n Begin subrun number " << iProc << endl;

    // Switch off previous process(es) and switch on new one(s).
    if (iProc > iFirst) for (int i = iBeg[iProc - 1]; i < iBeg[iProc]; ++i)
      pythia.readString( processes[i] + " = off" );
    for (int i = iBeg[iProc]; i < iBeg[iProc + 1]; ++i)
      pythia.readString( processes[i] + " = on" );
 
    // Initialize for LHC.
    pythia.init( 2212, 2212, 14000.);

    // Generate events to get cross section statistics.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) pythia.next();
 
    // Show statistics.
    //pythia.statistics();
    double sigma = pythia.info.sigmaGen();
    cout << " Cross section is " << scientific << setprecision(3)  
         << sigma << " and in Pythia6 was " << sigma6[iProc] 
         << ",\n i.e. now is factor >>> " << fixed 
         << sigma / sigma6[iProc] << " <<< different" <<endl;

  // End of loop over processes.
  }

  // Done. 
  return 0;
}
