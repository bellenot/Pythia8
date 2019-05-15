// main54.cc is a part of the PYTHIA event generator.
// Copyright (C) 2009 Mikhail Kirsanov, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. 
// It illustrates the chain Pythia6 -> Pythia8 -> HepMC.
// All input is specified in the main54.cmnd file.
// HepMC events are output to the hepmcout54.dat file.
// Written by Mikhail Kirsanov based on main52.cc.

#include "Pythia.h"
#include "LHAFortran.h"
#include "Pythia6Interface.h"
#include "HepMCInterface.h"

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "HepMC/Units.h"
// Following two lines are deprecated alternative.
//#include "HepMC/IO_Ascii.h"
//#include "HepMC/IO_AsciiParticles.h"

using namespace Pythia8; 

//**************************************************************************

// Implement initialization fillHepRup method for Pythia6 example.

bool LHAupFortran::fillHepRup() { 

  // Set process to generate.
  // Example: Z+jet production, must set pTmin, canset mMin.
  Pythia6Interface::pygive("msel = 13"); 
  Pythia6Interface::pygive("ckin(3) = 20."); 
  Pythia6Interface::pygive("ckin(1) = 50."); 

  // Switch off everything but Z0 -> leptons. 
  // Warning: only works with version Pythia 6.411 onwards.
  Pythia6Interface::pygive("23:alloff"); 
  Pythia6Interface::pygive("23:onifany = 11 13 15"); 

  // Speed up initialization: multiple interactions only in C++ code.
  Pythia6Interface::pygive("mstp(81)=0");
    
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

  // Interface for conversion from Pythia8::Event to HepMC one. 
  HepMC::I_Pythia8 ToHepMC;
  //  ToHepMC.set_crash_on_problem();

  // Specify file where HepMC events will be stored.
  HepMC::IO_GenEvent ascii_io("hepmcout54.dat", std::ios::out);
  // Following two lines are deprecated alternative.
  // HepMC::IO_Ascii ascii_io("hepmcout54.dat", std::ios::out);
  // HepMC::IO_AsciiParticles ascii_io("hepmcout54.dat", std::ios::out);

  // Generator. Shorthand for the event.
  Pythia8::Pythia pythia;
  Event& event = pythia.event;

  // Read in commands from external file.
  pythia.readFile("main54.cmnd");

  // Extract settings to be used in the main program.
  int  nEvent  = pythia.mode("Main:numberOfEvents");
  int  nList   = pythia.mode("Main:numberToList");
  int  nShow   = pythia.mode("Main:timesToShow");
  int  nAbort  = pythia.mode("Main:timesAllowErrors");
  bool showCS  = pythia.flag("Main:showChangedSettings");
  bool showAS  = pythia.flag("Main:showAllSettings");
  bool showCPD = pythia.flag("Main:showChangedParticleData");
  bool showAPD = pythia.flag("Main:showAllParticleData");

  // Initialize to access Pythia6 generator by Les Houches interface.
  LHAupFortran pythia6;
  pythia.init(&pythia6);    

  // List changed data.
  if (showCS) pythia.settings.listChanged();
  if (showAS) pythia.settings.listAll();

  // List particle data.  
  if (showCPD) pythia.particleData.listChanged();
  if (showAPD) pythia.particleData.listAll();

  // Histograms.
  double eCM   = 14000.;
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
  bool generated;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (nShow > 0 && iEvent%nPace == 0) 
      cout << " Now begin event " << iEvent << endl;

    // Generate events. Quit if too many failures.
    generated = pythia.next();
    if (!generated) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }
    cout << " successfully generated = " << generated << endl;
 
    // List first few events, both hard process and complete events.
    if (iEvent < nList) { 
      pythia.process.list();
      event.list();
    }

    // Construct new empty HepMC event. Arguments superfluous 
    // if HepMC was built with GeV and mm as units. 
    HepMC::GenEvent* hepmcevt = new HepMC::GenEvent(
      HepMC::Units::GEV, HepMC::Units::MM); 

    // Fill HepMC event, including PDF info.
    ToHepMC.fill_next_event( pythia, hepmcevt );
    // This alternative older method fills event, without PDF info.
    // ToHepMC.fill_next_event( pythia.event, hepmcevt );

    // Write the HepMC event to file. Done with it.
    ascii_io << hepmcevt;
    delete hepmcevt;

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
        int iNow = i;
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
      if (event[i].status() > 0) {
        ++nfin;
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

  // End of event loop.
  }

  // Final statistics.
  pythia.statistics();

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
