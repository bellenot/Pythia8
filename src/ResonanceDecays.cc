// Function definitions (not found in the header) for 
// the ResonanceDecays class.
// Copyright C 2007 Torbjorn Sjostrand

#include "ResonanceDecays.h"

namespace Pythia8 {
 
//**************************************************************************

// Do all resonance decays. First draft??
  
bool ResonanceDecays::next( Event& process) {

  // Loop over all entries to find resonances that should decay.
  int iDec = 0;
  do {
    Particle& decayer = process[iDec];
    if (decayer.isFinal() && decayer.canDecay() && decayer.mayDecay() 
    && decayer.isResonance() ) {

      // Particle data for decaying particle.
      int id0 = decayer.id();
      double m0 = decayer.m();

      // Special case for Z0: need dynamic decay channel selection.
      ResonanceGmZ gammaZRes;
      ResonanceW WRes;
      int idIn = process[decayer.mother1()].id();

      // Pick a decay channel; allow up to ten tries.
      int NTRYDECAY = 10;
      double mSafety = 1.; 
      int id1 = 0;
      int id2 = 0;
      double m1 = 0.;
      double m2 = 0.;
      bool physical = false;
      if (id0 != 23 && abs(id0) != 24) 
        decayer.particleData().decay.preparePick(id0);
      for (int iTryChannel = 0; iTryChannel < NTRYDECAY; ++iTryChannel) {
        // Correct dynamic treatment so far only for gamma*/Z0 and W+-.
        DecayChannel& channel = (id0 == 23 || abs(id0) == 24) 
          ? ( (id0 == 23) ? gammaZRes.dynamicDecay(m0, idIn) 
            : WRes.dynamicDecay(m0, idIn) ) 
          : decayer.particleData().decay.pickChannel();
        // int mode = channel.modeME();
        int mult = channel.multiplicity();

        // Consider for now only two-body decay. Check phase space. 
        if (mult != 2) continue;
        id1 = channel.product(0);
        if (id0 < 0 && ParticleDataTable::hasAnti(id1)) id1 = -id1;
        id2 = channel.product(1);
        if (id0 < 0 && ParticleDataTable::hasAnti(id2)) id2 = -id2;
        m1 = ParticleDataTable::mass(id1);          
        m2 = ParticleDataTable::mass(id2); 
        if (m1 + m2 + mSafety > m0) continue;

        // End of loop over tries.
        physical = true;
        break;
      }
      if (physical == false) return false;


      // Energies and absolute momentum in the rest frame.
      double e1 = 0.5 * (m0*m0 + m1*m1 - m2*m2) / m0;
      double e2 = 0.5 * (m0*m0 + m2*m2 - m1*m1) / m0;
      double pAbs = 0.5 * sqrtpos( (m0 - m1 - m2) * (m0 + m1 + m2)
        * (m0 + m1 - m2) * (m0 - m1 + m2) ) / m0;  

      // Isotropic angles give three-momentum. 
      double cosTheta = 2. * Rndm::flat() - 1.;
      // Correct treatment so far only for gamma*/Z0 and W+-.
      if (id0 == 23) cosTheta = gammaZRes.cosTheta(m0, idIn, id1);
      if (abs(id0) == 24) cosTheta = WRes.cosTheta(m0, idIn, id1, id2);        
      double sinTheta = sqrt(1. - cosTheta*cosTheta);
      double phi = 2. * M_PI * Rndm::flat();
      double pX = pAbs * sinTheta * cos(phi);  
      double pY = pAbs * sinTheta * sin(phi);  
      double pZ = pAbs * cosTheta;  

      // Fill four-momenta and boost them away from mother rest frame.
      Vec4 p1( pX, pY, pZ, e1);
      Vec4 p2( -pX, -pY, -pZ, e2);
      p1.bst( decayer.p() );
      p2.bst( decayer.p() );

      // Setup to find colours.
      int col0 = decayer.col();
      int acol0 = decayer.acol();
      int col1 = 0;
      int acol1 = 0;
      int col2 = 0;
      int acol2 = 0;
      int colType1 = ParticleDataTable::colType(id1);
      int colType2 = ParticleDataTable::colType(id2);

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
          ErrorMsg::message("Error in ResonanceDecays::next:"
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
          ErrorMsg::message("Error in ResonanceDecays::next:"
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
          ErrorMsg::message("Error in ResonanceDecays::next:"
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
          ErrorMsg::message("Error in ResonanceDecays::next:"
            " inconsistent colour tags");
          return false;
        }
      }

      // Append decay products to the event record.
      process.append( id1, 23, iDec, 0, 0, 0, col1, acol1, p1, m1, m0);
      process.append( id2, 23, iDec, 0, 0, 0, col2, acol2, p2, m2, m0);

      // Modify mother status to show it is a decayed resonance.
      decayer.status(-22);
                 
    // End of loop over all entries.
    }
  } while (++iDec < process.size());

  // Done.
  return true;
}

//**************************************************************************

} // end namespace Pythia8
