// File: test32.cc
// Studies of Prob(n_ch) and <pT>(n_ch),
// comparing internal default PDF with one from LHAPDF.
// Major differences indicate need for major retuning. 
// Copyright C 2007 Torbjorn Sjostrand.

#include "Pythia.h"
using namespace Pythia8; 

int main() {

  // Machine: 1 = Tevatron, 2 = LHC. Statistics.
  int machine = 2;
  int nEvent  = 1000;

  // Histograms.
  double nMax = (machine == 1) ? 199.5 : 399.5;
  Hist nChargedOld("n_charged old PDF", 100, -0.5, nMax);
  Hist nChargedNew("n_charged new PDF", 100, -0.5, nMax);
  Hist avgPTnChOld("<pT>(n_charged) old PDF", 100, -0.5, nMax);  
  Hist avgPTnChNew("<pT>(n_charged) new PDF", 100, -0.5, nMax);  
  Hist xDistOld("log(x) distribution old PDF", 100, -8., 0.); 
  Hist xDistNew("log(x) distribution new PDF", 100, -8., 0.); 
  Hist pTDistOld("pT (=Q) distribution old PDF", 100, 0., 20.); 
  Hist pTDistNew("pT (=Q) distribution new PDF", 100, 0., 20.); 

  // Loop over default run and one with new PDF.
  for (int iRun = 0; iRun < 2; ++iRun) {

    // Generator.
    Pythia pythia;
    Event& event = pythia.event;

    // Generate minimum-bias events 
    pythia.readString("SoftQCD:minBias = on");  
    //pythia.readString("SoftQCD:doubleDiffractive = on"); 

    // In second run pick new PDF set.
    if (iRun == 1) {
      pythia.readString("Pythia:useLHAPDF = on");
      //pythia.readString("Pythia:LHAPDFset = cteq5l.LHgrid");
      pythia.readString("Pythia:LHAPDFset = cteq61.LHpdf");
      //pythia.readString("Pythia:LHAPDFset = cteq61.LHgrid");
      //pythia.readString("Pythia:LHAPDFset = MRST2004nlo.LHgrid");
    }

    // Tevatron/LHC initialization. 
    double eCM = (machine == 1) ? 1960. : 14000.;
    if (machine == 1) pythia.init( 2212, -2212, eCM);
    else              pythia.init( 2212,  2212, eCM);
    pythia.settings.listChanged();
   
    // Begin event loop.
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      // Generate events.  Skip if error.
      if (!pythia.next()) continue;

      // Statistics on multiplicity and pT.
      int    nCh   = 0;
      double pTsum = 0.;
      for (int i = 0; i < event.size(); ++i) 
      if (event[i].isFinal() && event[i].isCharged()) {
        ++nCh;
        pTsum += event[i].pT();
      }
      
      // Fill histograms.
      if (iRun == 0) {
        nChargedOld.fill( nCh );
        avgPTnChOld.fill( nCh, pTsum / max(1, nCh) );
      } else {
        nChargedNew.fill( nCh );
        avgPTnChNew.fill( nCh, pTsum / max(1, nCh) );
      }

      // Loop through event record and fill x of all incoming partons.
      for (int i = 1; i < event.size(); ++i) 
      if (event[i].status() == -21 || event[i].status() == -31) {
        double x = 2. * event[i].e() / eCM;
        if (iRun == 0) xDistOld.fill( log10(x) );
        else           xDistNew.fill( log10(x) );
      }

      // Loop through multiple interactions list and fill pT of all MI's.
      for (int i = 0; i < pythia.info.nMI(); ++i) {
        double pT = pythia.info.pTMI(i);
        if (iRun == 0) pTDistOld.fill( pT );
        else           pTDistNew.fill( pT );
      }

    // End of event loop.
    }

  // Statistics. End of loop over two runs.
  pythia.statistics( true );
  }

  // Histograms.
  avgPTnChOld /= nChargedOld;
  avgPTnChNew /= nChargedNew;
  cout << nChargedOld << nChargedNew << avgPTnChOld << avgPTnChNew
       << xDistOld << xDistNew << pTDistOld << pTDistNew;

  // Done.
  return 0;
}
