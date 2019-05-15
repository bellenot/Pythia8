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

bool PartonLevel::init( Info* infoPtrIn, BeamParticle* beamAPtrIn, 
  BeamParticle* beamBPtrIn, int strategyIn) {

  // Store input pointers and modes for future use. 
  infoPtr = infoPtrIn;
  beamAPtr = beamAPtrIn;
  beamBPtr = beamBPtrIn;
  strategyLHA = strategyIn;

  // Main flags.
  ISR = Settings::flag("PartonLevel:ISR");
  FSRinProcess = Settings::flag("PartonLevel:FSRinProcess");
  FSRinResonances = Settings::flag("PartonLevel:FSRinResonances");
  MI = Settings::flag("PartonLevel:MI");
  if (!Settings::flag("Pythia:afterProcessLevel")) MI = false;

  // Set info in the respective program elements.
  if (strategyLHA < 10) {
    times.init();
    if (ISR) space.init( beamAPtr, beamBPtr);
    if (MI) MI = multi.init( beamAPtr, beamBPtr);
  }

  // Succeeded. (Check return values from other classes??)
  return true;
}

//*********

// Main routine to do the parton-level evolution.

bool PartonLevel::next( Event& process, Event& event) {

  // Special case if all partons already given.
  if (strategyLHA >= 10) return setupSimpleSys( process, event);

  // Special case if unresolved = elastic/diffractive event.
  if (!infoPtr->isResolved()) return setupUnresolvedSys( process, event);

  // Special case if minimum bias: do hardest interaction.
  multi.clear();
  if (infoPtr->isMinBias()) {
    multi.pTfirst();
    multi.setupFirstSys( infoPtr, process);
  }

  // Allow up to ten tries; failure possible for beam remnants.
  // Currently caused by lack of junction handling?? (Or colours!)
  bool physical = true;
  for (int iTry = 0; iTry < NTRY; ++ iTry) {

    // Reset counters and flag.
    nMI = 1;
    nISR = 0;
    nFSRinProc = 0;
    nFSRinRes = 0;
    physical = true;

    // Identify hard interaction system for showers.
    setupHardSys( process, event);

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
      double pTmulti =  (MI) ? multi.pTnext( pTmax, pTgen) 
        : -1.;
      pTgen = max( pTgen, pTmulti);
      double pTspace = (ISR) ? space.pTnext( pTmax, pTgen) 
        : -1.;

      // One further possibility is that all pT have fallen
      // below a pTveto scale, and that a special veto routine 
      // is called to decide whether to keep the event. ??

      // Do a multiple interaction (if allowed).
      if (pTmulti > 0. && pTmulti > pTspace && pTmulti > pTtimes) {
        multi.scatter( event);  
        ++nMI;
        if (ISR) space.prepare( event, sizeOld);
        // times.update( event);
        pTmax = pTmulti;
      }
   
      // Do an initial-state emission (if allowed).
      else if (pTspace > 0. && pTspace > pTtimes) { 
      //   if (space.branch()) times.update( event); 
        space.branch( event);
        ++nISR;
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
    if (!remnants.add( *beamAPtr, *beamBPtr, event)) physical = false;

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
    beamAPtr->clear();
    beamBPtr->clear();

  // End loop over ten tries. Hopefully it worked
  }
  if (!physical) return false;

  // Perform showers in resonance decay chains.
  resonanceShowers( process, event); 

  // Store statistics.
  infoPtr->setImpact( multi.bMI(), multi.enhanceMI());
  infoPtr->setCounters( nMI, nISR, nFSRinProc, nFSRinRes);
 
  // Done.
  return true;
}

//*********

// Set up the hard process, excluding subsequent resonance decays.

void PartonLevel::setupHardSys( Event& process, Event& event) {

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
  beamAPtr->append( inP, process[inP].id(), x1);
  double x2 = process[inM].pMinus() / process[0].e();
  beamBPtr->append( inM, process[inM].id(), x2);

  // Find whether incoming partons are valence or sea.
  beamAPtr->xfISR( 0, process[inP].id(), x1, scale*scale);
  beamAPtr->pickValSeaComp(); 
  beamBPtr->xfISR( 0, process[inM].id(), x2, scale*scale);
  beamBPtr->pickValSeaComp(); 

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

// Set up an unresolved process, i.e. elastic or diffractive.

bool PartonLevel::setupUnresolvedSys( Event& process, Event& event) {

  // Copy particles from process to event.
  for (int i = 0; i < process.size(); ++ i) event.append( process[i]);

  // Loop to find diffractively excited beams.
  for (int i = 0; i < 2; ++i)  
  if ( (i == 0 && infoPtr->isDiffractiveA()) 
    || (i == 1 && infoPtr->isDiffractiveB()) ) {
    int iBeam = i + 3;
    BeamParticle* beamPtr = (i == 0) ? beamAPtr : beamBPtr;

    // Diffractive mass. Reconstruct boost and rotation to event cm frame.
    double mDiff = process[iBeam].m();  
    double m2Diff = mDiff*mDiff;  
    double beta = process[iBeam].pAbs() / process[iBeam].e();
    double theta = process[iBeam].theta();
    double phi = process[iBeam].phi();
  
    // Pick quark or gluon kicked out and flavour subdivision.
    bool gluonIsKicked = beamPtr->pickGluon(mDiff);
    int id1 = beamPtr->pickValence();
    int id2 = beamPtr->pickRemnant();

    // Find flavour masses. Scale them down if too big.
    double m1 = ParticleDataTable::constituentMass(id1);
    double m2 = ParticleDataTable::constituentMass(id2);
    if (m1 + m2 > 0.5 * mDiff) { 
      double reduce = 0.5 * mDiff / (m1 + m2);
      m1 *= reduce;
      m2 *= reduce;
    }

    // If quark is kicked out, then trivial kinematics in rest frame.
    if (!gluonIsKicked) { 
      double pAbs = sqrt( pow2(m2Diff - m1*m1 - m2*m2) 
        - pow2(2. * m1 * m2) ) / (2. * mDiff);
      double e1 = (m2Diff + m1*m1 - m2*m2) / (2. * mDiff);
      double e2 = (m2Diff + m2*m2 - m1*m1) / (2. * mDiff);
      Vec4 p1(0.,0., -pAbs, e1);
      Vec4 p2(0.,0., pAbs, e2);

      // Boost and rotate to event cm frame.
      p1.bst(0., 0., beta); p1.rot(theta, phi);   
      p2.bst(0., 0., beta); p2.rot(theta, phi);   

      // Set colours.
      int col1, acol1, col2, acol2;
      if (ParticleDataTable::colType(id1) == 1) {
        col1 = event.nextColTag(); acol1 = 0;
        col2 = 0; acol2 = col1;
      } else {  
        col1 = 0; acol1 = event.nextColTag();
        col2 = acol1; acol2 = 0;
      }    
    
      // Store partons of diffractive system and mark system decayed.
      int iDauBeg = event.append( id1, 23, iBeam, 0, 0, 0, col1, acol1, 
        p1, m1);
      int iDauEnd = event.append( id2, 63, iBeam, 0, 0, 0, col2, acol2, 
        p2, m2);
      event[iBeam].statusNeg();
      event[iBeam].daughters(iDauBeg, iDauEnd);   


    // If gluon is kicked out: share momentum between two remnants.
    } else {
      double m2Sys, zSys, pxSys, pySys, mTS1, mTS2;
      zSys = beamPtr->zShare(mDiff, m1, m2);

      // Provide relative pT kick in remnant. Construct (transverse) masses.
      pxSys = beamPtr->pxShare(); 
      pySys = beamPtr->pyShare(); 
      mTS1 = m1*m1 + pxSys*pxSys + pySys*pySys;
      mTS2 = m2*m2 + pxSys*pxSys + pySys*pySys;
      m2Sys = mTS1 / zSys + mTS2 / (1. - zSys);

      // Momentum of kicked-out massless gluon in diffractive rest frame.
      double pAbs = (m2Diff - m2Sys) / (2. * mDiff);
      Vec4 pG(0., 0., -pAbs, pAbs);
      Vec4 pRem(0., 0., pAbs, mDiff - pAbs);

      // Momenta of the two beam remnant flavours. (Lightcone p+ = m_diff!)
      double e1 = 0.5 * (zSys * mDiff + mTS1 / (zSys * mDiff));    
      double pL1 = 0.5 * (zSys * mDiff - mTS1 / (zSys * mDiff));  
      Vec4 p1(pxSys, pySys, pL1, e1);
      Vec4 p2 = pRem - p1;
  
      // Boost and rotate to event cm frame.
      pG.bst(0., 0., beta); pG.rot(theta, phi);   
      p1.bst(0., 0., beta); p1.rot(theta, phi);   
      p2.bst(0., 0., beta); p2.rot(theta, phi); 

      // Set colours.
      int colG, acolG, col1, acol1, col2, acol2;
      if (ParticleDataTable::colType(id1) == 1) {
        col1 = event.nextColTag(); acol1 = 0;
        colG = event.nextColTag(); acolG = col1;
        col2 = 0; acol2 = colG;
      } else {  
        col1 = 0; acol1 = event.nextColTag();
        colG = acol1; acolG = event.nextColTag();
        col2 = acolG; acol2 = 0;
      } 
       
      // Store partons of diffractive system and mark system decayed.
      int iDauBeg = event.append( 21, 23, iBeam, 0, 0, 0, colG, acolG, 
        pG, m1);
      event.append( id1, 63, iBeam, 0, 0, 0, col1, acol1, p1, m1);
      int iDauEnd = event.append( id2, 63, iBeam, 0, 0, 0, col2, acol2, 
        p2, m2);
      event[iBeam].statusNeg();
      event[iBeam].daughters(iDauBeg, iDauEnd);   
    }

  // End loop over beams. Done.
  }
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
