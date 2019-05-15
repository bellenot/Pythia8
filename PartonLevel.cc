// Function definitions (not found in the header) for the PartonLevel class.
// Copyright C 2006 Torbjorn Sjostrand

#include "PartonLevel.h"

namespace Pythia8 {
 
//**************************************************************************

// The PartonLevel class.

//*********

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Maximum number of tries to produce parton level from given input..
const int PartonLevel::NTRY = 10; 

//*********

// Main routine to initialize the parton-level generation process.

bool PartonLevel::init( BeamParticle& beamA, BeamParticle& beamB, 
  int strategyIn) {

  // Save input.
  strategyLHA = strategyIn;

  // Main flags.
  ISR = Settings::flag("PartonLevel:ISR");
  FSRinProcess = Settings::flag("PartonLevel:FSRinProcess");
  FSRinResonances = Settings::flag("PartonLevel:FSRinResonances");
  MI = Settings::flag("PartonLevel:MI");

  // Set info in the respective program elements.
  if (strategyLHA < 10) {
    times.init();
    if (ISR) space.init( beamA, beamB);
    if (MI) MI = multi.init( beamA, beamB);
  }

  // Succeeded. (Check return values from other classes??)
  return true;
}

//*********

// Main routine to do the parton-level evolution.

bool PartonLevel::next( BeamParticle& beamA, BeamParticle& beamB, 
  Event& process, Event& event) {

  // Special case if all partons already given.
  if (strategyLHA >= 10) return setupSimpleSys( process, event);

  // Allow up to ten tries; failure possible for beam remnants.
  // Currently caused by lack of junction handling?? (Or colours!)
  bool physical = true;
  for (int iTry = 0; iTry < NTRY; ++ iTry) {
    physical = true;

    // Identify hard interaction system for showers.
    setupHardSys( beamA, beamB, process, event);

    // Set hard scale, maximum for showers and multiple interactions,
    double pTmax = process.scale();

    // Prepare the classes to begin the generation.
    // Need to redo everything??
    if  (MI) multi.prepare( pTmax);
    if (ISR) space.prepare( event);
    //if (FSRinProcess) times.prepare( event);

    // Begin evolution down in pT from hard pT scale.  
    do {

      // Find next pT value for FSR, MI and ISR.
      // Order calls to minimize time expenditure.
      int sizeOld = event.size();
      double pTgen = 0.;
      // Currently timelike showers are not interleaved??
      //double pTtimes = (FSRinProcess) ? times.pTnext( event, pTmax, pTgen) 
      //  : -1.;
      double pTtimes = -1.;
      pTgen = max( pTgen, pTtimes);
      double pTmulti =  (MI) ? multi.pTnext( beamA, beamB, pTmax, pTgen) : -1.;
      pTgen = max( pTgen, pTmulti);
      double pTspace = (ISR) ? space.pTnext( beamA, beamB, pTmax, pTgen) : -1.;

      // One further possibility is that all pT have fallen
      // below a pTveto scale, and that a special veto routine 
      // is called to decide whether to keep the event. ??

      // Do a multiple interaction (if allowed).
      if (pTmulti > 0. && pTmulti > pTspace && pTmulti > pTtimes) {
        if (multi.scatter( beamA, beamB, event)) { 
          if (ISR) space.prepare( event, sizeOld);
          // times.update( event);
	}
        pTmax = pTmulti;
      }
   
      // Do an initial-state emission (if allowed).
      else if (pTspace > 0. && pTspace > pTtimes) { 
      //   if (space.branch()) times.update( event); 
        space.branch( beamA, beamB, event);
        pTmax = pTspace;
      }

      // Do a final-state emission (if allowed).
      // Yet to be done here; currently moved to end??
      //else if (pTtimes > 0) {
      //  times.branch( event); 
      //  pTmax = pTtimes;
      //}
   
      // Keep on evolving until nothing is left to be done.
      else pTmax = 0.;

    } while (pTmax > 0.);   

    // Now add beam remnants. Includes primordial kT kick and colour tracing.
    if (!remnants.add( beamA, beamB, event)) physical = false;

    // Temporary position for final-state emissions??
    if (FSRinProcess && physical) {
      // Find largest scale for final partons.
      double pTmaxFSR = 0.;
      for (int i = 0; i < event.size(); ++i) 
        if (event[i].remains() && event[i].scale() > pTmaxFSR)
          pTmaxFSR = event[i].scale();     
      // Let all partons shower up to their individual scale.
      times.shower( event, 0, 0, pTmaxFSR);
    }

    // If no problems then done, else restore and loop.
    if (physical) break;
    event.clear();
    beamA.clear();
    beamB.clear();

  // End loop over ten tries. Hopefully it worked
  }
  if (!physical) return false;

  // Perform showers in resonance decay chains.
  resonanceShowers( process, event); 

  // Done.
  return true;
}

//*********

// Set up the hard process, excluding subsequent resonance decays.

void PartonLevel::setupHardSys( BeamParticle& beamA, BeamParticle& beamB, 
  Event& process, Event& event) {

  // Incoming partons to hard process are stored in slots 3 and 4. Scale.
  int inP = 3;
  int inM = 4;
  double scale = process.scale();

  // If incoming partons are massive then recalculate to put them massless.
  if (process[inP].m() != 0. || process[inM].m() != 0.) { 
    double pPlus  = process[inP].pPlus() + process[inM].pPlus();
    double pMinus = process[inP].pMinus() + process[inM].pMinus(); 
    process[inP].pz( 0.5 * pPlus);
    process[inP].e(  0.5 * pPlus);
    process[inP].m( 0.);
    process[inM].pz(-0.5 * pMinus);
    process[inM].e(  0.5 * pMinus);
    process[inM].m( 0.);
  }

  // Add incoming hard-scattering partons to list in beam remnants.
  // Not valid if not in rest frame??
  double x1 = process[inP].pPlus() / process[0].e();
  beamA.append( inP, process[inP].id(), x1);
  double x2 = process[inM].pMinus() / process[0].e();
  beamB.append( inM, process[inM].id(), x2);

  // Find whether incoming partons are valence or sea.
  beamA.xfISR( 0, process[inP].id(), x1, scale*scale);
  beamA.pickValSeaComp(); 
  beamB.xfISR( 0, process[inM].id(), x2, scale*scale);
  beamB.pickValSeaComp(); 

  // Initialize info needed for subsequent sequential decays + showers.
  nHardDone = 0;
  iPosBefShow.resize( process.size() );

  // Add the beam and hard subprocess partons to the event record.
  for (int i = 0; i < process.size(); ++ i) { 
    if (process[i].mother1() > inM) break;
    event.append(process[i]);
    iPosBefShow[i] = i;

    // Currently outgoing ones should not count as decayed.
    if (event[i].status() == -22) { 
      event[i].statusPos(); 
      event[i].daughters(0, 0);
    }

    // Complete task of copying hard subsystem into event record.
    ++nHardDone;
  }

  // Update event colour tag to maximum in whole process.
  int maxColTag = 0;
  for (int i = 0; i < process.size(); ++ i) { 
    if (process[i].col() > maxColTag) maxColTag = process[i].col();
    if (process[i].acol() > maxColTag) maxColTag = process[i].acol();
  }
  event.initColTag(maxColTag); 

  // Copy junctions from process to event.
  for (int i = 0; i < process.sizeJunction(); ++i) 
    event.appendJunction( process.getJunction(i));

  // Done. 
}

//*********

// Set up the hard process, special case if all partons already given.

bool PartonLevel::setupSimpleSys( Event& process, Event& event) {

  // Copy particles from process to event.
  for (int i = 0; i < process.size(); ++ i) event.append( process[i]);

  // Copy junctions from process to event.
  for (int i = 0; i < process.sizeJunction(); ++i) 
    event.appendJunction( process.getJunction(i));

  // Done.
  return true;
}

//*********

// Handle showers in successive resonance decays.

void PartonLevel::resonanceShowers( Event& process, Event& event) {

  // Isolate next system to be processed, if anything remains.
  while (nHardDone < process.size()) {
    int iBegin = nHardDone;

    // Mother in hard process and in complete event (after shower).
    int iHardMother = process[iBegin].mother1();
    Particle& hardMother = process[iHardMother];
    int iBefMother = iPosBefShow[iHardMother];
    int iAftMother = event.iBotCopyId(iBefMother);
    Particle& aftMother = event[iAftMother];

    // From now mother counts as decayed.
    aftMother.statusNeg();

    // Mother can have been moved by showering (in any of previous steps), 
    // so prepare to update colour and momentum information for system.
    int colBef = hardMother.col();
    int acolBef = hardMother.acol();
    int colAft = aftMother.col();
    int acolAft = aftMother.acol();
    RotBstMatrix M;
    M.bst( hardMother.p(), aftMother.p());

    // Extract next partons from hard event into normal event record.
    for (int i = iBegin; i < process.size(); ++ i) { 
      if (process[i].mother1() != iHardMother) break;
      int iNow = event.append( process[i] );
      iPosBefShow[i] = iNow;
      Particle& now = event.back();

      // Currently outgoing ones should not count as decayed.
      if (now.status() == -22) { 
        now.statusPos(); 
        now.daughters(0, 0);
      }
      
      // Update daughter and mother information.
      if (i == iBegin) aftMother.daughter1( iNow);
      else aftMother.daughter2( iNow); 
      now.mother1(iAftMother); 

      // Update colour and momentum information.
      if (now.col() == colBef) now.col( colAft);
      if (now.acol() == acolBef) now.acol( acolAft);
      now.rotbst( M);   

      // Complete task of copying next subsystem into event record.
      ++nHardDone;
    }
    int iEnd = nHardDone - 1;

    // Do parton showers inside subsystem.
    if (FSRinResonances) {
      double pTmax = 0.5 * hardMother.m();
      times.shower( event, iPosBefShow[iBegin], iPosBefShow[iEnd], pTmax);
    }    

  // No more systems to be processed.
  }

}

//*********

// Print statistics, if any.

void PartonLevel::statistics() {

  // Preliminary list, to expand.
  if (MI) multi.statistics();

}
 
//**************************************************************************

} // end namespace Pythia8
