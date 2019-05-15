// Function definitions (not found in the header) for 
// the ResonanceDecays class.
// Copyright C 2007 Torbjorn Sjostrand

#include "ResonanceDecays.h"

namespace Pythia8 {
 
//**************************************************************************

// The ResonanceDecays class.
// Do all resonance decays sequentially.
// Remains: three-body; Breit-Wigner masses * phase space ??

//*********

// Number of tries to pick a decay channel.
const int    ResonanceDecays::NTRYDECAY = 100;

// Mass above threshold for allowed decays.
const double ResonanceDecays::MSAFETY   = 1.; 

//*********
  
bool ResonanceDecays::next( Event& process) {

  // Loop over all entries to find resonances that should decay.
  int iDec = 0;
  do {
    Particle& decayer = process[iDec];
    if (decayer.isFinal() && decayer.canDecay() && decayer.mayDecay() 
    && decayer.isResonance() ) {

      // Particle data for decaying particle.
      id0    = decayer.id();
      id0Abs = abs(id0);
      m0     = decayer.m();

      // Pointer to dynamically defined resonances.
      ResonanceProperties* resonancePtr = 0;
      if (id0Abs == 23)    resonancePtr = new ResonanceGmZ;
      if (id0Abs == 24)    resonancePtr = new ResonanceW;
      if (id0Abs == 25)    resonancePtr = new ResonanceH; 

      // Prepare decay selection for other resonances.
      if (resonancePtr == 0) decayer.particleData().decay.preparePick(id0);

      // Mother flavour.
      int idIn = process[decayer.mother1()].id();

      // Default properties of a decay channel.
      bool physical         = false;
      DecayChannel* channel = 0;

      // Pick a decay channel; allow up to ten tries.
      for (int iTryChannel = 0; iTryChannel < NTRYDECAY; ++iTryChannel) {

        // Correct dynamic treatment so far only for gamma*/Z0, W+- and H0.
        if (resonancePtr > 0) 
             channel = &resonancePtr->dynamicDecay(m0, idIn);
        else channel = &decayer.particleData().decay.pickChannel();  

        // Consider for now only two-body decay.
        int mult = channel->multiplicity();
        if (mult != 2) continue;

        // Read out flavours.
        id1 = channel->product(0);
        if (id0 < 0 && ParticleDataTable::hasAnti(id1)) id1 = -id1;
        id2 = channel->product(1);
        if (id0 < 0 && ParticleDataTable::hasAnti(id2)) id2 = -id2;

        // Read out masses. Check phase space. Breit-Wigner phase-space??
        m1 = ParticleDataTable::mass(id1);          
        m2 = ParticleDataTable::mass(id2); 
        if (m1 + m2 + MSAFETY > m0) continue;

        // Working decay channel found.
        physical = true;
        break;
      }

      // Release pointer to resonance.
      delete resonancePtr;

      // Failed to find acceptable decays.
      if (physical == false) {
        ErrorMsg::message("Error in ResonanceDecays::next:"
          " failed to find workable decay channel");         
        return false;
      }

      // Select colours in decay.
      if (!pickColours(decayer, process)) return false;

      // Select four-momenta in mother rest frame and boost them.
      pickKinematics();
      p1.bst( decayer.p() );
      p2.bst( decayer.p() );

      // Append decay products to the event record.
      int i1 = process.append( id1, 23, iDec, 0, 0, 0, col1, acol1, 
        p1, m1, m0);
      int i2 = process.append( id2, 23, iDec, 0, 0, 0, col2, acol2, 
        p2, m2, m0);

      // Modify mother status and daughters.
      decayer.status(-22);
      decayer.daughter1(i1); 
      decayer.daughter2(i2); 
                 
    // End of loop over all entries.
    }
  } while (++iDec < process.size());

  // Done.
  return true;
}

//*********

// Select colours in decay.
  
bool ResonanceDecays::pickColours(Particle& decayer, Event& process) {

  // Setup to find colours.
  col0     = decayer.col();
  acol0    = decayer.acol();
  col1     = 0;
  acol1    = 0;
  col2     = 0;
  acol2    = 0;
  colType1 = ParticleDataTable::colType(id1);
  colType2 = ParticleDataTable::colType(id2);

  // Colour singlet decay.
  if (col0 == 0 && acol0 == 0) {
    if (colType1 == 0 && colType2 == 0) ;
    else if (colType1 == 1 && colType2 == -1) {
      col1 = process.nextColTag(); 
      acol2 = col1;
    } else if (colType1 == -1 && colType2 == 1) {
      acol1 = process.nextColTag(); 
      col2 = acol1;
    } else if (colType1 == 2 && colType2 == 2) {
      col1 = process.nextColTag();
      acol1 = process.nextColTag(); 
      col2 = acol1;
      acol2 = col1;
    } else {
      ErrorMsg::message("Error in ResonanceDecays::pickColours:"
        " inconsistent colour tags");
      return false;
    }

  // Colour triplet decay.
  } else if (col0 > 0 && acol0 == 0) {
    if (colType1 == 1 && colType2 == 0) col1 = col0;
    else if (colType1 == 0 && colType2 == 1) col2 = col0;
    else if (colType1 == 1 && colType2 == 2) {
      col2 = col0;
      col1 = process.nextColTag(); 
      acol2 = col1;
    } else if (colType1 == 2 && colType2 == 1) {
      col1 = col0;
      acol1 = process.nextColTag(); 
      col2 = acol1;
   } else {
      ErrorMsg::message("Error in ResonanceDecays::pickColours:"
        " inconsistent colour tags");
      return false;
    }

  // Colour antitriplet decay.
  } else if (col0 == 0 && acol0 > 0) {
    if (colType1 == -1 && colType2 == 0) acol1 = acol0;
    else if (colType1 == 0 && colType2 == -1) acol2 = acol0;
    else if (colType1 == -1 && colType2 == 2) {
      acol2 = acol0;
      acol1 = process.nextColTag(); 
      col2 = acol1;
    } else if (colType1 == 2 && colType2 == -1) {
      acol1 = acol0;
      col1 = process.nextColTag(); 
      acol2 = col1;
    } else {
      ErrorMsg::message("Error in ResonanceDecays::pickColours:"
        " inconsistent colour tags");
      return false;
    }

  // Colour octet decay.
  } else {
    if (colType1 == 1 && colType2 == -1) {
      col1 = col0;
      acol2 = acol0;
    } else if (colType1 == -1 && colType2 == 1) {
      acol1 = acol0;
      col2 = col0;
    } else if (colType1 == 2 && colType2 == 0) {
      col1 = col0;
      acol1 = acol0;
    } else if (colType1 == 0 && colType2 == 2) {
      col2 = col0;
      acol2 = acol0;
    } else if (colType1 == 2 && colType2 == 2 && Rndm::flat() > 0.5) {
      col1 = col0;
      acol1 = process.nextColTag(); 
      col2 = acol1;
      acol2 = acol0;
    } else if (colType1 == 2 && colType2 == 2) {
      col1 = process.nextColTag(); 
      acol1 = acol0;
      col2 = col0;
      acol2 = col1;
    } else {
      ErrorMsg::message("Error in ResonanceDecays::pickColours:"
        " inconsistent colour tags");
      return false;
    }
  }

  // Done.
  return true;

}

//*********

// Select momenta isotropically in rest frame of decaying particle.
// Later on, process-dependent angular distributions may be imposed.
  
void ResonanceDecays::pickKinematics() {

  // Energies and absolute momentum in the rest frame.
  double e1   = 0.5 * (m0*m0 + m1*m1 - m2*m2) / m0;
  double e2   = 0.5 * (m0*m0 + m2*m2 - m1*m1) / m0;
  double pAbs = 0.5 * sqrtpos( (m0 - m1 - m2) * (m0 + m1 + m2)
    * (m0 + m1 - m2) * (m0 - m1 + m2) ) / m0;  

  // Pick isotropic angles to give three-momentum. 
  double cosTheta = 2. * Rndm::flat() - 1.;
  double sinTheta = sqrt(1. - cosTheta*cosTheta);
  double phi      = 2. * M_PI * Rndm::flat();
  double pX       = pAbs * sinTheta * cos(phi);  
  double pY       = pAbs * sinTheta * sin(phi);  
  double pZ       = pAbs * cosTheta;  

  // Fill four-momenta in mother rest frame. 
  p1 = Vec4(  pX,  pY,  pZ, e1);
  p2 = Vec4( -pX, -pY, -pZ, e2);

}

//**************************************************************************

} // end namespace Pythia8
