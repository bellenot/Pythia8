// main52.cc is a part of the PYTHIA event generator.
// Copyright (C) 2009 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a more extensive program linking to Pythia6 for hard processes.
// It illustrates how input can be provided from external files both for
// the Fortran and the C++ part, see main52.fcmnd and main52.ccmnd.

#include "Pythia.h"
#include "LHAFortran.h"
#include "Pythia6Interface.h"

using namespace Pythia8; 

//**************************************************************************

// Implement initialization fillHepRup method for Pythia6 example.

bool LHAupFortran::fillHepRup() { 

  // Open file with commands to Pythia6. Check that it worked.
  ifstream is("main52.fcmnd");  
  if (!is) {
    infoPtr->errorMsg("Error in LHAupFortran::fillHepRup: "
      "did not find file");
    return false;
  } 

  // Read in one line at a time and send it on to pygive and Pythia6.
  string line;
  while ( getline(is, line) ) Pythia6Interface::pygive( line);
    
  // Initialize for 14 TeV pp collider.
  Pythia6Interface::pyinit("cms","p","p",14000.);   

  // Fill initialization information in HEPRUP.
  Pythia6Interface::pyupin();

  // Done.
  return true;

}

//**************************************************************************

// Implement event generation fillHepEup method for Pythia6 example.

bool LHAupFortran::fillHepEup() { 

  // Generate and fill the next Pythia6 event in HEPEUP.
  Pythia6Interface::pyupev();

  // Done.
  return true;

}

//**************************************************************************

int main() {

  // Generator. Shorthand for the event.
  Pythia pythia;
  Event& event = pythia.event;

  // Read in Pythia8 commands from external file.
  pythia.readFile("main52.ccmnd");    

  // Initialize to access Pythia6 generator by Les Houches interface.
  LHAupFortran pythia6;
  pythia.init(&pythia6);    

  // Extract settings to be used in the main program.
  int  nEvent  = pythia.mode("Main:numberOfEvents");
  int  nList   = pythia.mode("Main:numberToList");
  int  nShow   = pythia.mode("Main:timesToShow");
  int  nAbort  = pythia.mode("Main:timesAllowErrors");
  bool showCS  = pythia.flag("Main:showChangedSettings");
  bool showAS  = pythia.flag("Main:showAllSettings");
  bool showCPD = pythia.flag("Main:showChangedParticleData");
  bool showAPD = pythia.flag("Main:showAllParticleData");

  // List changed data.
  if (showCS) pythia.settings.listChanged();
  if (showAS) pythia.settings.listAll();

  // List particle data.  
  if (showCPD) pythia.particleData.listChanged();
  if (showAPD) pythia.particleData.listAll();

  // Histograms.
  double eCM = 14000.;
  double epTol = 1e-7 * eCM;
  Hist epCons("deviation from energy-momentum conservation",100,0.,epTol);
  Hist nFinal("final particle multiplicity",100,-0.5,1599.5);
  Hist nChg("final charged multiplicity",100,-0.5,799.5);
  Hist nISR("number of ISR emissions for hard system",40,-0.5,39.5);
  Hist nMI("number of MI (excluding hard system)",100,-0.5,99.5);
  Hist nISRMI("number of ISR emissions per MI",40,-0.5,39.5);
  Hist nFSR("total number of FSR emissions",100,-0.5,299.5);
  Hist nJUN("number of junctions",10,-0.5,9.5);
  Hist pThard("ISR pT kick to hard system",100,0.,400.);
  Hist sumETparticle("summed ET of particles",100,0.,2000.);
  Hist dnCHparticleDy("dn_charged/dy for particles",100,-10.,10.);
  Hist dETparticleDy("dET/dy for particles",100,-10.,10.);

  // Begin event loop.
  int nPace  = max(1, nEvent / max(1, nShow) ); 
  int iAbort = 0; 
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (nShow > 0 && iEvent%nPace == 0) 
      cout << " Now begin event " << iEvent << endl;

    // Generate events. Quit if too many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }
 
    // List first few events, both hard process and complete events.
    if (iEvent < nList) { 
      pythia.info.list();
      // This call to Pythia6 is superfluous, but shows it can be done.
      Pythia6Interface::pylist(1);
      pythia.process.list();
      event.list();
    }

    // Number of ISR for hard subprocess.
    int nisr = -1;
    int iNow = 3;
    do { iNow = event[iNow].mother1(); ++nisr;}
    while (iNow > 1 && abs(event[iNow].status()) < 50 ); 
    nISR.fill(nisr);

    // Total pT kick of hard subsystem.
    Vec4 pHard;
    for (int i = 0; i < event.size(); ++i) {
      if (abs(event[i].status()) > 30) break;
      if (event[i].status() == -22 || event[i].status() == -23) {      
        iNow = i;
        while (event[iNow].daughter2() == event[iNow].daughter1() &&
          event[iNow].daughter1() > iNow) iNow = event[iNow].daughter1();
        pHard += event[iNow].p();
      }
    }
    pThard.fill(pHard.pT());

    // Reset quantities to be summed over event.
    int nfin = 0;
    int nch = 0;
    int nmi = 0;
    int nfsr = 0;
    int nfinqg = 0;
    Vec4 pSum = - (event[1].p() + event[2].p());
    double eTsum = 0.;

    // Loop over particles in the event. 
    for (int i = 0; i < event.size(); ++i) {

      // Number of MI and of ISR per MI.
      if (i < event.size() - 1 && event[i].status() == -31 
        && event[i+1].status() == -31) {
        ++nmi;
        int inow = i;
        int nisrmi = -1;
        do { inow = event[inow].mother1(); ++ nisrmi;}
        while (inow > 1 && abs(event[inow].status()) < 50) ; 
        nISRMI.fill(nisrmi);
      }
    
      // Number of FSR branchings.
      if (event[i].status() == -52) ++nfsr; 

      // Specialize to final particles. Total multiplicity and momentum.
      if (event[i].isFinal()) {
        ++nfin;
        if (event[i].isQuark() || event[i].isGluon()) ++nfinqg;
        if (event[i].isCharged()) ++nch;
        pSum += event[i].p();

        // Final-state particle spectra.
        double eTnow = event[i].pT();
        double ynow = event[i].y();
        eTsum += eTnow;
        if (event[i].isCharged()) dnCHparticleDy.fill(ynow);
        dETparticleDy.fill(ynow,eTnow);

      // End of loop over (final/all) particles.
      }
    }

    // Energy-momentum deviation.
    double epDev = abs(pSum.e()) + abs(pSum.px()) + abs(pSum.py())
      + abs(pSum.pz());
    epCons.fill(epDev);
      
    // Fill summed quantities.
    nFinal.fill(nfin);
    nChg.fill(nch);
    nMI.fill(nmi);
    nFSR.fill(nfsr);
    nJUN.fill( event.sizeJunction() );
    sumETparticle.fill(eTsum);

    // Debug.
    if (nfinqg > 0) {
      cout << " Error: number of unframented q/qbar/g = " << nfinqg << "\n";
      event.list();
    }

  // End of event loop.
  }

  // Final statistics. Must do call to Pythia6 explicitly.
  pythia.statistics();
  Pythia6Interface::pystat(1);  

  // Histogram normalization.
  double factor = 5. / (nEvent - nAbort);  
  dnCHparticleDy *= factor;
  dETparticleDy *= factor;

  // Histogram output.
  cout << epCons << nFinal<< nChg << nISR << nMI << nISRMI << nFSR 
       << nJUN << pThard << sumETparticle << dnCHparticleDy 
       << dETparticleDy; 

  // Done.
  return 0;
}
