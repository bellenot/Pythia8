// File: main01.cc
// This is a simple test program. 
// It illustrates how Pythia6 processes can be used in Pythia8.
// All input is specified in the main01.cmnd file.
// Copyright C 2006 Torbjorn Sjostrand

#include "Pythia.h"

using namespace Pythia8; 

//**************************************************************************

int main() {

  // Generator. Shorthand for the event and for settings.
  Pythia pythia;
  Event& event = pythia.event;
  Settings& settings = pythia.settings;

  // Read in commands from external file.
  pythia.readFile("main01.cmnd");    

  // Extract settings to be used in the main program.
  int idBeamA = settings.mode("Main:idBeamA");
  int idBeamB = settings.mode("Main:idBeamB");
  bool inCMframe = settings.flag("Main:inCMframe");
  double eCM = settings.parameter("Main:eCM");
  double eBeamA = settings.parameter("Main:eBeamA");
  double eBeamB = settings.parameter("Main:eBeamB");
  int nEvent = settings.mode("Main:numberOfEvents");
  int nList = settings.mode("Main:numberToList");
  int nShow = settings.mode("Main:timesToShow");
  int nAbort = settings.mode("Main:timesAllowErrors");
  bool showChangedSettings = settings.flag("Main:showChangedSettings");
  bool showAllSettings = settings.flag("Main:showAllSettings");
  bool showParticleData = settings.flag("Main:showParticleData");

  // Initialization for Pythia6 event input.
  if (inCMframe) pythia.init( idBeamA, idBeamB, eCM);
  else pythia.init( idBeamA, idBeamB, eBeamA, eBeamB);

  // List changed data.
  if (showChangedSettings) settings.listChanged();
  if (showAllSettings) settings.listAll();

  // List particle data.  
  if (showParticleData) ParticleDataTable::list();

  // Histograms.
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
  int nShowPace = max(1,nEvent/nShow); 
  int iAbort = 0; 
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (iEvent%nShowPace == 0) cout << " Now begin event " 
      << iEvent << endl;

    // Generate events. Quit if too many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }
 
    // List first few events, both hard process and complete events.
    if (iEvent < nList) { 
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
      if (event[i].remains()) {
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
