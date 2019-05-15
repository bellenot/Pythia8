// File: main15.cc
// This is a simple test program. 
// It illustrates how Les Houches Accord input can be used 
// to feed in toy parton-level configurations.
// Copyright C 2007 Torbjorn Sjostrand

#include "Pythia.h"
using namespace Pythia8; 

//**************************************************************************

// Derived classes for external feed-in of parton configurations
// to be hadronized, without any parton-level evolution.

//*********

class LHAinitToy : public LHAinit {

public:

  // Constructor.
  LHAinitToy() {}

  // Routine for doing the job of setting initialization info.  
  // The strategy = 10 is special option that does not require 
  // beams to have been specified, and does not add showers.
  bool set() {strategy(10); return true;} 

};

//*********

class LHAevntToy: public LHAevnt {

public:

  // Constructor.
  LHAevntToy(int typeIn = 0, double energy = 10.) : type(typeIn),
    ee(energy) {}

  // Routine for doing the job of setting info on next event.  
  bool set(); 

  // Store which kind of toy events to generate.
  int type;
  double ee;

};

//*********

bool LHAevntToy::set() {

  // Information on process: needs to be called but OK to use defaults.
  process();

  // Information on a q qbar or g g system, to be hadronized.
  if (type == 1 || type == 2) {
    int id1 = (type == 1) ?  2 : 21;
    int id2 = (type == 1) ? -2 : 21;
    particle( id1, 1, 0, 0, 101, 0, 0., 0.,  ee, ee); 
    particle( id2, 1, 0, 0, 0, 101, 0., 0., -ee, ee);

  // Information on a g g g system, to be hadronized.
  } else if (type == 3) {  
    particle( 21, 1, 0, 0, 101, 102, 0., 0., ee, ee); 
    particle( 21, 1, 0, 0, 102, 103,  0.8 * ee, 0., -0.6 * ee, 
      ee); 
    particle( 21, 1, 0, 0, 103, 101, -0.8 * ee, 0., -0.6 * ee, 
      ee); 

  // Information on a q q q junction system, to be hadronized.
  } else if (type == 4 || type == 5) { 

    // The three endpoint q q q; the minimal system. 
    double rt75 = sqrt(0.75);  
    particle( 2, 1, 0, 0, 101, 0, 0., 0., 1.01 * ee, 1.01 * ee); 
    particle( 2, 1, 0, 0, 102, 0,  rt75 * ee, 0., -0.5 * ee, ee ); 
    particle( 1, 1, 0, 0, 103, 0, -rt75 * ee, 0., -0.5 * ee, ee );

    // Define the qqq configuration as starting point for adding gluons.
    if (type == 5) {
      int colNow[4] = {0, 101, 102, 103}; 
      Vec4 pQ[4];
      pQ[1] = Vec4(0., 0., 1., 0.);
      pQ[2] = Vec4( rt75, 0., -0.5, 0.);
      pQ[3] = Vec4(-rt75, 0., -0.5, 0.);

      // Minimal cos(q-g opening angle), allows more or less nasty events.
      double cosThetaMin =0.;
      
      // Add a few gluons (almost) at random. 
      for (int nglu = 0; nglu < 5; ++nglu) { 
        int iq = 1 + int( 2.99999 * Rndm::flat() );
        double px, py, pz, e, prod; 
        do {
          e =  ee * Rndm::flat();
          double cThe = 2. * Rndm::flat() - 1.;
          double phi = 2. * M_PI * Rndm::flat(); 
          px = e * sqrt(1. - cThe*cThe) * cos(phi);
          py = e * sqrt(1. - cThe*cThe) * sin(phi);
          pz = e * cThe;
          prod = ( px * pQ[iq].px() + py * pQ[iq].py() + pz * pQ[iq].pz() ) 
            / e;
        } while (prod < cosThetaMin); 
        int colNew = 104 + nglu;
        particle( 21, 1, 0, 0, colNew, colNow[iq], px, py, pz, e, 0.);
        colNow[iq] = colNew;   
      }
    }

  // Information on a q q qbar qbar dijunction system, to be hadronized.
  } else if (type >= 6) {

    // The two fictitious beam remnant particles; needed for documentation.
    particle( 2212, 2, 0, 0, 0, 0, 0., 0., ee, ee, 0.);
    particle(-2212, 2, 0, 0, 0, 0, 0., 0., ee, ee, 0.);

    // Opening angle between "diquark" legs.
    double theta = 0.2;
    double cThe = cos(theta);
    double sThe = sin(theta);  

    // Set one colour depending on whether more gluons or not.
    int acol = (type == 6) ? 103 : 106;

    // The four endpoint q q qbar qbar; the minimal system. 
    // Two additional fictitious partons to make up original beams.
    particle(  2, 1, 1, 0, 101, 0,  ee * sThe, 0.,  ee * cThe, ee, 0.); 
    particle(  1, 1, 1, 0, 102, 0, -ee * sThe, 0.,  ee * cThe, ee, 0.); 
    particle(  2, 2, 1, 0, 103, 0,         0., 0.,  ee       , ee, 0.); 
    particle( -2, 1, 2, 0, 0, 104,  ee * sThe, 0., -ee * cThe, ee, 0.); 
    particle( -1, 1, 2, 0, 0, 105, -ee * sThe, 0., -ee * cThe, ee, 0.); 
    particle( -2, 2, 2, 0, 0, acol,        0., 0., -ee       , ee, 0.); 

    // Add extra gluons on string between junctions.
    if (type == 7) {
      particle( 21, 1, 5, 8, 103, 106, 0., ee, 0., ee, 0.); 
    } else if (type == 8) {
      particle( 21, 1, 5, 8, 103, 108, 0., ee, 0., ee, 0.); 
      particle( 21, 1, 5, 8, 108, 106, 0.,-ee, 0., ee, 0.); 
    } else if (type == 9) {
      particle( 21, 1, 5, 8, 103, 107, 0., ee, 0., ee, 0.); 
      particle( 21, 1, 5, 8, 107, 108, ee, 0., 0., ee, 0.); 
      particle( 21, 1, 5, 8, 108, 106, 0.,-ee, 0., ee, 0.); 
    } else if (type == 10) {
      particle( 21, 1, 5, 8, 103, 107, 0., ee, 0., ee, 0.); 
      particle( 21, 1, 5, 8, 107, 108, ee, 0., 0., ee, 0.); 
      particle( 21, 1, 5, 8, 108, 109, 0.,-ee, 0., ee, 0.); 
      particle( 21, 1, 5, 8, 109, 106,-ee, 0., 0., ee, 0.); 
    }

  // No more cases: done.
  } 
  return true; 
}

//**************************************************************************

int main() {

  // Pick kind of events to generate:
  // 1 = q qbar.
  // 2 = g g. 
  // 3 = g g g.
  // 4 = minimal q q q junction topology.
  // 5 = q q q junction topology with gluons on the strings.
  // 6 = q q qbar qbar dijunction topology, no gluons.
  // 7 - 10 : ditto, but with 1 - 4 gluons on string between junctions.
  int eventType = 1;

  // Set typical energy per parton.
  double energy = 20.0;

  // Set number of events to generate and to list.
  int nEvent = 10000;
  int nList = 1;

  // Generator; shorthand for event.                           
  Pythia pythia;  
  Event& event = pythia.event;

  // Standard checks not meaningful (assume incoming beams in lines 1 and 2).
  pythia.readString("Pythia:checkEvent = off");

  // Expand event listing with junctions.                           
  pythia.readString("Event:listJunctions = on");

  // Optionally switch off decays.
  pythia.readString("HadronLevel:Decay = off");
                             
  // Initialize generation with LHA input.                           
  LHAinitToy lhaInit;       
  LHAevntToy lhaEvnt(eventType, energy);      
  pythia.init(&lhaInit, &lhaEvnt);    

  // Provide printout of initial information.        
  pythia.settings.listChanged();

  // Book histograms.                          
  Hist epCons("deviation from energy-momentum conservation",100,0.,1e-4);
  Hist chgCons("deviation from charge conservation",57,-9.5,9.5);
  Hist nFinal("final particle multiplicity",100,-0.5,99.5);   
  Hist dnparticledp("dn/dp for particles",100,0.,energy);
  Hist status85("multiplicity status code 85",50,-0.5,49.5);
  Hist status86("multiplicity status code 86",50,-0.5,49.5);
  Hist status83("multiplicity status code 83",50,-0.5,49.5);
  Hist status84("multiplicity status code 84",50,-0.5,49.5);
  Hist dndtheta("particle flow in event plane",100,-M_PI,M_PI);
  Hist dedtheta("energy flow in event plane",100,-M_PI,M_PI);
  Hist dpartondtheta("parton flow in event plane",100,-M_PI,M_PI);
  Hist dndyAnti("dn/dy primaries antijunction",100, -10., 10.);
  Hist dndyJun("dn/dy primaries junction",100, -10., 10.);
  Hist dndySum("dn/dy all primaries",100, -10., 10.);
  
  // Begin of event loop.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (iEvent%(max(1,nEvent/20)) == 0) cout << " Now begin event " 
      << iEvent << endl;

    // Generate events. Quit if failure.
    if (!pythia.next()) {
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }
 
    // List first few events, both hard process and complete events.
    if (iEvent < nList) { 
      pythia.LHAevntList();
      pythia.process.list();
      pythia.event.list();
    }

    // Initialize statistics. 
    Vec4 pSum = - event[0].p();
    double chargeSum = 0.;
    if (eventType == 4 || eventType == 5) chargeSum = -1;
    int nFin = 0;  
    int n85 = 0;
    int n86 = 0;
    int n83 = 0;
    int n84 = 0;
                          
    // Loop over all particles.
    for (int i = 0; i < event.size(); ++i) {
      int status = event[i].statusAbs();

      // Find any unrecognized particle codes.
      int id = event[i].id();
      if (id == 0 || !ParticleDataTable::isParticle(id))
        cout << " Error! Unknown code id = " << id << "\n"; 

      // Find particles with E-p mismatch.
      double eCalc = event[i].eCalc();
      if (abs(eCalc/event[i].e() - 1.) > 1e-6) cout << " e mismatch, i = "
        << i << " e_nominal = " << event[i].e() << " e-from-p = " 
        << eCalc << " m-from-e " << event[i].mCalc() << "\n";

      // Parton flow in event plane.
      if (status == 71 || status == 72) { 
        double thetaXZ = event[i].thetaXZ();
        dpartondtheta.fill(thetaXZ);
      }

      // Origin of primary hadrons.
      if (status == 85) ++n85;
      if (status == 86) ++n86;
      if (status == 83) ++n83;
      if (status == 84) ++n84;

      // Flow of primary hadrons in the event plane.
      if (status > 80 && status < 90) {
        double eAbs = event[i].e();
        if (eAbs < 0.) {cout << " e < 0 line " << i; event.list();}
        double thetaXZ = event[i].thetaXZ();
        dndtheta.fill(thetaXZ);
        dedtheta.fill(thetaXZ, eAbs);
 
        // Rapidity distribution of primary hadrons.
        double y = event[i].y();
        dndySum.fill(y);
        if (eventType >= 6) {
          int motherId = event[event[i].mother1()].id();
          if (motherId > 0 ) dndyJun.fill(event[i].y()); 
          else dndyAnti.fill(event[i].y());
	}
      }

      // Study final-state particles.
      if (event[i].isFinal()) {
        pSum += event[i].p();
        chargeSum += event[i].charge();
        nFin++;
        double pAbs = event[i].pAbs();
        dnparticledp.fill(pAbs);
      }
    }
 
    // Fill histograms once for each event.
    double epDev = abs(pSum.e()) + abs(pSum.px()) + abs(pSum.py())
      + abs(pSum.pz());
    epCons.fill(epDev);
    chgCons.fill(chargeSum);
    nFinal.fill(nFin); 
    status85.fill(n85);
    status86.fill(n86);
    status83.fill(n83);
    status84.fill(n84);
    if (epDev > 1e-3  || abs(chargeSum) > 0.1) event.list(); 
                       
  // End of event loop.
  }                                           

  // Print histogram and done.
  cout << epCons << chgCons << nFinal << dnparticledp
       << dndtheta << dedtheta << dpartondtheta << dndySum;
  if (eventType >= 4) cout << status85 << status86 << status83 
       << status84; 
  if (eventType >= 6) cout << dndyJun << dndyAnti;

  // Done.
  return 0;
}
