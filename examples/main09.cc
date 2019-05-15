// File: main09.cc
// This is a simple test program. 
// It studies event properties of LEP1 events.
// Copyright C 2006 Torbjorn Sjostrand

#include "Pythia.h"

using namespace Pythia8; 

int main() {

  // Generator.
  Pythia pythia;

  // Process selection.
  pythia.readString("Pythia6:msel = 1");    

  // No ISR in e+e- - beam remnants not yet implemented.
  pythia.readString("Pythia6:mstp(11) = 0"); 
  pythia.readString("PartonLevel:ISR = off"); 

  // Switch off Z0 decay to charged leptons and neutrinos. 
  //pythia.readString("Pythia6:mdme(182,1) = 0"); 
  //pythia.readString("Pythia6:mdme(183,1) = 0"); 
  //pythia.readString("Pythia6:mdme(184,1) = 0"); 
  //pythia.readString("Pythia6:mdme(185,1) = 0"); 
  //pythia.readString("Pythia6:mdme(186,1) = 0"); 
  //pythia.readString("Pythia6:mdme(187,1) = 0"); 

  // LEP1 initialization.   
  pythia.init( 11, -11, 91.2);

  // Histograms.
  Hist nCharge("charged multiplicity", 100, -0.5, 99.5);
  Hist spheri("Sphericity", 100, 0., 1.);
  Hist linea("Linearity", 100, 0., 1.);
  Hist saxis("cos(theta_Sphericity)", 100, -1., 1.);
  Hist laxis("cos(theta_Linearity)", 100, -1., 1.);

  // Set up Sphericity and "Linearity" analyses.
  Sphericity sph;  
  Sphericity lin(1.);  

  // Begin event loop. Generate event. Skip if error. List first few.
  for (int iEvent = 0; iEvent < 10000; ++iEvent) {
    if (!pythia.next()) continue;
    if (iEvent < 1) pythia.event.list();

    // Find and histogram charged multiplicity. 
    int nCh = 0;
    for (int i = 0; i < pythia.event.size(); ++i) 
      if (pythia.event[i].remains() && pythia.event[i].isCharged()) ++nCh;
    nCharge.fill( nCh );

    // Find and histogram sphericity. 
    if (sph.analyze( pythia.event )) { 
      // if (iEvent < 3) sph.list();
      spheri.fill( sph.sph() ); 
      saxis.fill( sph.eigenVector(1).pz() );
      double e1 = sph.eigenValue(1);
      double e2 = sph.eigenValue(2);
      double e3 = sph.eigenValue(3);
      if (e2 > e1 || e3 > e2) cout << "eigenvalues out of order: "
      << e1 << "  " << e2 << "  " << e3 << endl;
    }
    if (lin.analyze( pythia.event )) {
      // if (iEvent < 3) lin.list();
      linea.fill( lin.sph() ); 
      laxis.fill( lin.eigenVector(1).pz() );
      double e1 = lin.eigenValue(1);
      double e2 = lin.eigenValue(2);
      double e3 = lin.eigenValue(3);
      if (e2 > e1 || e3 > e2) cout << "eigenvalues out of order: "
      << e1 << "  " << e2 << "  " << e3 << endl;
    }

  // End of event loop. Statistics. Output histograms. 
  }
  pythia.statistics();
  cout << nCharge << spheri << linea << saxis << laxis; 

  // Done.
  return 0;
}
