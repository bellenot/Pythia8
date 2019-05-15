// ResonanceProperties.cc is a part of the PYTHIA event generator.
// Copyright (C) 2007 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for 
// the ResonanceProperties class and classes derived from it.

#include "ResonanceProperties.h"
#include "PythiaComplex.h"

namespace Pythia8 {

//**************************************************************************

// The ResonanceProperties class.
// Base class for the various resonances.

//*********
 
// Definitions of static variables and functions.
// (Values will be overwritten in initStatic call, so are purely dummy.)

int          ResonanceProperties::alphaSorder    = 1;
int          ResonanceProperties::alphaEMorder   = 1;
double       ResonanceProperties::alphaSvalue    = 0.1265;
AlphaStrong  ResonanceProperties::alphaS;
AlphaEM      ResonanceProperties::alphaEM;

// Number of points in integration direction for numInt routines.
const int    ResonanceProperties::NPOINT         = 100;

// The sum of product masses must not be too close to the resonance mass.
const double ResonanceProperties::MASSMARGIN     = 0.1;

//*********

// Initialize static data members.

void ResonanceProperties::initStatic() {

  // Parameters of alphaStrong generation .
  alphaSvalue  = Settings::parm("SigmaProcess:alphaSvalue");
  alphaSorder  = Settings::mode("SigmaProcess:alphaSorder");

  // Initialize alphaStrong generation.
  alphaS.init( alphaSvalue, alphaSorder); 

  // Parameters of alphaEM generation.
  alphaEMorder = Settings::mode("SigmaProcess:alphaEMorder");

  // Initialize alphaEM generation.
  alphaEM.init( alphaEMorder); 

}

//*********

// Numerical integration of matrix-element in two-body decay,
// where one particle is described by a Breit-Wigner mass distribution.
// Normalization to unit integral if matrix element is unity 
// and there are no phase-space restrictions. 
 
double ResonanceProperties::numInt1BW(double mHat, double m1, double Gamma1, 
  double mMin1, double m2, int meMode) {

  // Check that phase space is open for integration.
  if (mMin1 + m2 > mHat) return 0.;

  // Precalculate coefficients for Breit-Wigner selection.
  double s1       = m1 * m1;
  double mG1      = m1 * Gamma1;
  double mMax1    = mHat - m2; 
  double atanMin1 = atan( (mMin1 * mMin1 - s1) / mG1 );
  double atanMax1 = atan( (mMax1 * mMax1 - s1) / mG1 );
  double atanDif1 = atanMax1 - atanMin1;
  double wtDif1   = atanDif1 / (M_PI * NPOINT); 

  // Step size in atan-mapped variable.
  double xStep   = 1. / NPOINT;

  // Variables used in loop over integration points.
  double sum = 0.;
  double xNow1, sNow1, mNow1, ratNow1, phasespace, value; 
  double rat2 = pow2(m2 / mHat);

  // Loop with first-particle mass selection.
  for (int ip1 = 0; ip1 < NPOINT; ++ip1) {
    xNow1   = xStep * (ip1 + 0.5);
    sNow1   = s1 + mG1 * tan(atanMin1 + xNow1 * atanDif1);
    mNow1   = min( mMax1, max( mMin1, sqrtpos(sNow1) ) );
    ratNow1 = pow2(mNow1 / mHat);

    // Evaluate value and add to sum. Different matrix elements.
    phasespace = sqrtpos( pow2(1. - ratNow1 - rat2) - 4. * ratNow1 * rat2);
    value      = phasespace;
    if (meMode == 1) value *= (pow2(1. - ratNow1 - rat2) 
       + 8. * ratNow1 * rat2);
    sum += value;

  // End of  loop over integration points. Overall normalization.
  }  
  sum *= wtDif1;
 
  // Done.
  return sum;
}

//*********

// Numerical integration of matrix-element in two-body decay,
// where both particles are described by Breit-Wigner mass distributions.
// Normalization to unit integral if matrix element is unity
// and there are no phase-space restrictions. 
 
double ResonanceProperties::numInt2BW(double mHat, double m1, double Gamma1, 
  double mMin1, double m2, double Gamma2, double mMin2, int meMode) {

  // Check that phase space is open for integration.
  if (mMin1 + mMin2 > mHat) return 0.;

  // Precalculate coefficients for Breit-Wigner selection.
  double s1       = m1 * m1;
  double mG1      = m1 * Gamma1;
  double mMax1    = mHat - mMin2; 
  double atanMin1 = atan( (mMin1 * mMin1 - s1) / mG1 );
  double atanMax1 = atan( (mMax1 * mMax1 - s1) / mG1 );
  double atanDif1 = atanMax1 - atanMin1;
  double wtDif1   = atanDif1 / (M_PI * NPOINT); 
  double s2       = m2 * m2;
  double mG2      = m2 * Gamma2;
  double mMax2    = mHat - mMin1; 
  double atanMin2 = atan( (mMin2 * mMin2 - s2) / mG2 );
  double atanMax2 = atan( (mMax2 * mMax2 - s2) / mG2 );
  double atanDif2 = atanMax2 - atanMin2;
  double wtDif2   = atanDif2 / (M_PI * NPOINT); 

  // If on-shell decay forbidden then split integration range
  // to ensure that low-mass region is not forgotten.
  bool mustDiv    = false;
  double mDiv1    = 0.;
  double atanDiv1 = 0.;
  double atanDLo1 = 0.; 
  double atanDHi1 = 0.;
  double wtDLo1   = 0.;
  double wtDHi1   = 0.;
  double mDiv2    = 0.;
  double atanDiv2 = 0.;
  double atanDLo2 = 0.;
  double atanDHi2 = 0.;
  double wtDLo2   = 0.;
  double wtDHi2   = 0.;
  if (m1 + m2 > mHat) {
    mustDiv = true;
    double tmpDiv = (mHat - m1 - m2) / (Gamma1 + Gamma2);
    mDiv1         = m1 + Gamma1 * tmpDiv;
    atanDiv1      = atan( (mDiv1 * mDiv1 - s1) / mG1 );
    atanDLo1      = atanDiv1 - atanMin1;
    atanDHi1      = atanMax1 - atanDiv1;
    wtDLo1        = atanDLo1 / (M_PI * NPOINT); 
    wtDHi1        = atanDHi1 / (M_PI * NPOINT); 
    mDiv2         = m2 + Gamma2 * tmpDiv;
    atanDiv2      = atan( (mDiv2 * mDiv2 - s2) / mG2 ); 
    atanDLo2      = atanDiv2 - atanMin2;
    atanDHi2      = atanMax2 - atanDiv2;
    wtDLo2        = atanDLo2 / (M_PI * NPOINT); 
    wtDHi2        = atanDHi2 / (M_PI * NPOINT); 
  }

  // Step size in atan-mapped variable.
  double xStep   = 1. / NPOINT;
  int nIter      = (mustDiv) ? 2 * NPOINT : NPOINT;

  // Variables used in loop over integration points.
  double sum = 0.;
  double xNow1, sNow1, mNow1, ratNow1,
         xNow2, sNow2, mNow2, ratNow2, 
         phasespace, value; 
  double wtNow1 = wtDif1; 
  double wtNow2 = wtDif2; 

  // Outer loop with first-particle mass selection.
  for (int ip1 = 0; ip1 < nIter; ++ip1) {
    if (!mustDiv) {
      xNow1  = xStep * (ip1 + 0.5);
      sNow1  = s1 + mG1 * tan(atanMin1 + xNow1 * atanDif1);
    } else if (ip1 < NPOINT) {
      xNow1  = xStep * (ip1 + 0.5);
      sNow1  = s1 + mG1 * tan(atanMin1 + xNow1 * atanDLo1);
      wtNow1 = wtDLo1;
    } else {
      xNow1  = xStep * (ip1 - NPOINT + 0.5);
      sNow1  = s1 + mG1 * tan(atanDiv1 + xNow1 * atanDHi1);
      wtNow1 = wtDHi1;
    }
    mNow1    = min( mMax1, max( mMin1, sqrtpos(sNow1) ) );
    ratNow1  = pow2(mNow1 / mHat);

    // Inner loop with second-particle mass selection.
    for (int ip2 = 0; ip2 < nIter; ++ip2) {
      if (!mustDiv) {
        xNow2  = xStep * (ip2 + 0.5);
        sNow2  = s2 + mG2 * tan(atanMin2 + xNow2 * atanDif2); 
      } else if (ip2 < NPOINT) {
        xNow2  = xStep * (ip2 + 0.5);
        sNow2  = s2 + mG2 * tan(atanMin2 + xNow2 * atanDLo2);
        wtNow2 = wtDLo2;
      } else {
        xNow2  = xStep * (ip2 - NPOINT + 0.5);
        sNow2  = s2 + mG2 * tan(atanDiv2 + xNow2 * atanDHi2);
        wtNow2 = wtDHi2;
      }
      mNow2 = min( mMax2, max( mMin2, sqrtpos(sNow2) ) );
      ratNow2 = pow2(mNow2 / mHat);

      // Check that point is inside phase space.
      if (mNow1 + mNow2 > mHat) break;

      // Evaluate value and add to sum. Different matrix elements.
      phasespace = sqrtpos( pow2(1. - ratNow1 - ratNow2) 
        - 4. * ratNow1 * ratNow2);
      value      = phasespace;
      // This meMode is dummy so far.
      if (meMode == 1) value *= (pow2(1. - ratNow1 - ratNow2) 
        + 8. * ratNow1 * ratNow2);
      sum += value * wtNow1 * wtNow2;

    // End of second and first loop over integration points.
    }
  }  

  // Done.
  return sum;
}
 
//**************************************************************************

// The ResonanceGmZ class.
// Derived class for gamma*/Z0 properties.

//*********
 
// Definitions of static variables and functions.
int    ResonanceGmZ::idRes   = 23;
int    ResonanceGmZ::gmZmode = 0; 
double ResonanceGmZ::mRes;
double ResonanceGmZ::GammaRes;
double ResonanceGmZ::m2Res; 
double ResonanceGmZ::GamMRat; 
double ResonanceGmZ::thetaWRat; 
double ResonanceGmZ::openPos = 1.; 
double ResonanceGmZ::openNeg = 0.; 
ParticleDataEntry* ResonanceGmZ::particlePtr;

//*********

// Initialize static data members.

void ResonanceGmZ::initStatic() {

  // Set pointer to particle properties and decay table.
  particlePtr = ParticleDataTable::particleDataPtr(idRes);

  // Z0 properties.
  mRes      = ParticleDataTable::m0(idRes);
  GammaRes  = ParticleDataTable::mWidth(idRes);
  m2Res     = mRes*mRes;
  GamMRat   = GammaRes / mRes;
  thetaWRat = 1. / (16. * CoupEW::sin2thetaW() * CoupEW::cos2thetaW());

  // Allow to pick only gamma* or Z0 part of full gamma*/Z0 expression.
  gmZmode   = Settings::mode("SigmaProcess:gmZmode");

}

//*********

// Calculate and store partial and total widths at the nominal mass. 
// Only consider the Z0 part, as appropriate for Breit-Wigner width.

void ResonanceGmZ::widthInit() {

  // Common coupling factors.
  double mHat   = mRes;
  double sH     = mHat*mHat;
  double alpEM  = alphaEM.alphaEM(sH);
  double alpS   = alphaS.alphaS(sH);
  double colQ   = 3. * (1. + alpS / M_PI);
  double preFac = alpEM * thetaWRat * mHat / 3.;

  // Reset quantities to sum. Declare variables inside loop.
  double widTot = 0.; 
  double widPos = 0.;
  int    idAbs, onMode;
  double widNow, mf, mr, betaf, psvec, psaxi;

  // Loop over all decay channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow = 0.;
    idAbs  = abs( particlePtr->decay[i].product(0) );

    // Only contributions from three fermion generations, except top.
    if ( (idAbs > 0 && idAbs < 6) || ( idAbs > 10 && idAbs < 17)) {
      mf = ParticleDataTable::m0(idAbs);

      // Check that above threshold. Phase space.
      if (mHat > 2. * mf + MASSMARGIN) {
        mr    = pow2(mf / mHat);
        betaf = sqrtpos(1. - 4. * mr); 
        psvec = betaf * (1. + 2. * mr);
        psaxi = pow3(betaf);

        // Combine phase space with colour factor and couplings.
        widNow = preFac * (CoupEW::vf2(idAbs) * psvec 
          + CoupEW::af2(idAbs) * psaxi); 
        if (idAbs < 6) widNow *= colQ;
      }
    }

    // Store partial widths and update sum, also for open channels only.
    particlePtr->decay[i].onShellWidth(widNow); 
    widTot += widNow;    
    onMode = particlePtr->decay[i].onMode();
    if (onMode == 1 || onMode == 2) widPos += widNow;
  }

  // Update total width and branching ratios.
  particlePtr->setMWidth(widTot, false);
  double bRatio;
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    bRatio = particlePtr->decay[i].onShellWidth() / widTot;
    particlePtr->decay[i].bRatio( bRatio, false);
  }
  particlePtr->setIsInitResonance();

  // Update the width in this object itself. Secondary widths for open.
  GammaRes = widTot;
  GamMRat  = GammaRes / mRes;  
  openPos  = widPos / widTot;
    

  }

//*********

// Return fraction of width open for one/two-resonance combinations.
  
double ResonanceGmZ::openFrac( int idResSgn1, int idResSgn2) {
  
  // Only one resonance.
  if (idResSgn2 == 0) return (idResSgn1 > 0) ? openPos : openNeg;

  // Pair of resonances.
  if (idResSgn1 * idResSgn2 < 0) return openPos * openNeg;
  return (idResSgn1 > 0) ? pow2(openPos) : pow2(openNeg);  

}

//*********

// Calculate the total width and store phase-space-weighted coupling sums.
// For idIn > 0 provides relative weight of specific incoming flavour.
// Use full gamma*/Z0, as appropriate for final-state flavour choice.

double ResonanceGmZ::width(double mHat, int , int idIn, 
  bool openOnly, bool setBR) {

  // Common coupling factors.
  double sH     = mHat*mHat;
  double alpEM  = alphaEM.alphaEM(sH); 
  double alpS   = alphaS.alphaS(sH);
  double colQ   = 3. * (1. + alpS / M_PI);
  double preFac = alpEM * thetaWRat * mHat / 3.;

  // Couplings when an incoming fermion is specified.
  int idInAbs   = abs(idIn);
  double ei2    = (idInAbs > 0) ? CoupEW::ef2(idInAbs) : 0;
  double eivi   = (idInAbs > 0) ? CoupEW::efvf(idInAbs) : 0;
  double vi2ai2 = (idInAbs > 0) ? CoupEW::vf2af2(idInAbs) : 0; 

  // Calculate prefactors for gamma/interference/Z0 terms.
  double gamNorm = ei2;
  double intNorm = 2. * eivi * thetaWRat * sH * (sH - m2Res)
    / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  double resNorm = vi2ai2 * pow2(thetaWRat * sH) 
    / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );

  // Optionally only keep gamma* or Z0 term.
  if (gmZmode == 1) {intNorm = 0.; resNorm = 0.;}
  if (gmZmode == 2) {gamNorm = 0.; intNorm = 0.;}

  // Reset quantities to sum. Declare variables inside loop.
  if (setBR) particlePtr->decay.resetDynamicBR(); 
  double widSum = 0.; 
  int    idAbs, onMode;
  double widNow, mf, mr, betaf, psvec, psaxi, vf2af2;

  // Loop over all decay channels; optionally only open channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow = 0.;
    onMode = particlePtr->decay[i].onMode();
    if (openOnly && onMode !=1 && onMode != 2) continue;
    idAbs  = abs( particlePtr->decay[i].product(0) );

    // Only contributions from three fermion generations, except top.
    if ( (idAbs > 0 && idAbs < 6) || ( idAbs > 10 && idAbs < 17)) {
      mf = ParticleDataTable::m0(idAbs);

      // Check that above threshold. Phase space.
      if (mHat > 2. * mf + MASSMARGIN) {
        mr     = pow2(mf / mHat);
        betaf  = sqrtpos(1. - 4. * mr); 
        psvec  = betaf * (1. + 2. * mr);
        psaxi  = pow3(betaf);
        vf2af2 = CoupEW::vf2(idAbs) * psvec + CoupEW::af2(idAbs) * psaxi; 

        // Z0 width: combine phase space with couplings.
        if (idIn == 0) widNow = preFac * vf2af2;

        // Relative outwidths: combine instate, propagator and outstate.
        else widNow = gamNorm * CoupEW::ef2(idAbs) * psvec 
          + intNorm * CoupEW::efvf(idAbs) * psvec + resNorm * vf2af2;

        // Colour factor.
        if (idAbs < 6) widNow *= colQ;
      }
    }

    // Sum back up.
    widSum += widNow;

    // Optionally store partial widths for later decay channel choice.
    if (setBR) particlePtr->decay[i].dynamicBR(widNow); 
  }

  // Done.
  return widSum;
  
}
 
//**************************************************************************

// The ResonanceW class.
// Derived class for W+- properties.

//*********
 
// Definitions of static variables and functions.
int    ResonanceW::idRes = 24;
double ResonanceW::mRes;
double ResonanceW::GammaRes;
double ResonanceW::m2Res;
double ResonanceW::GamMRat;
double ResonanceW::thetaWRat;
double ResonanceW::openPos = 1.; 
double ResonanceW::openNeg = 1.; 
ParticleDataEntry* ResonanceW::particlePtr;

//*********

// Initialize static data members.

void ResonanceW::initStatic() {

  // Set pointer to particle properties and decay table.
  particlePtr = ParticleDataTable::particleDataPtr(idRes);

  // W properties.
  mRes      = ParticleDataTable::m0(idRes);
  GammaRes  = ParticleDataTable::mWidth(idRes);
  m2Res     = mRes*mRes;
  GamMRat   = GammaRes / mRes;
  thetaWRat = 1. / (12. * CoupEW::sin2thetaW());

}

//*********

// Calculate and store partial and total widths at the nominal mass. 

void ResonanceW::widthInit() {

  // Common coupling factors.
  double mHat   = mRes;
  double sH     = mHat*mHat;
  double alpEM  = alphaEM.alphaEM(sH);
  double alpS   = alphaS.alphaS(sH);
  double colQ   = 3. * (1. + alpS / M_PI);
  double preFac = alpEM * thetaWRat * mHat;

  // Reset quantities to sum. Declare others.
  double widTot = 0.; 
  double widPos = 0.;
  double widNeg = 0.;
  int    id1Abs, id2Abs, onMode;
  double widNow, mf1, mf2, mr1, mr2, ps;

  // Loop over all decay channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow = 0.;
    id1Abs = abs( particlePtr->decay[i].product(0) );
    id2Abs = abs( particlePtr->decay[i].product(1) );

    // Only contributions from three fermion generations, except top.
    if ( (id1Abs > 0 && id1Abs < 6 && id2Abs > 0 && id2Abs < 6) 
      || (id1Abs > 10 && id1Abs < 17 && id2Abs > 10 && id2Abs < 17)) {
      mf1 = ParticleDataTable::m0(id1Abs);
      mf2 = ParticleDataTable::m0(id2Abs);

      // Check that above threshold. Phase space.
      if (mHat > mf1 + mf2 + MASSMARGIN) {
        mr1 = pow2(mf1 / mHat);
        mr2 = pow2(mf2 / mHat);
        ps = (1. - 0.5 * (mr1 + mr2) - 0.5 * pow2(mr1 - mr2))
          * sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 ); 

        // Combine phase space with colour factor and CKM couplings.
        widNow = preFac * ps;
        if (id1Abs < 6) widNow *= colQ * VCKM::V2id(id1Abs, id2Abs);
      }
    }

    // Store partial widths and update sum, also for open channels only.
    particlePtr->decay[i].onShellWidth(widNow); 
    widTot += widNow;    
    onMode = particlePtr->decay[i].onMode();
    if (onMode == 1 || onMode == 2) widPos += widNow;
    if (onMode == 1 || onMode == 3) widNeg += widNow;
  }

  // Update total width and branching ratios.
  particlePtr->setMWidth(widTot, false);
  double bRatio;
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    bRatio = particlePtr->decay[i].onShellWidth() / widTot;
    particlePtr->decay[i].bRatio( bRatio, false);
  }
  particlePtr->setIsInitResonance();

  // Update the width in this object itself. Secondary widths for open.
  GammaRes = widTot;
  GamMRat  = GammaRes / mRes;  
  openPos  = widPos / widTot;
  openNeg  = widNeg / widTot;

}

//*********

// Return fraction of width open for one/two-resonance combinations.
  
double ResonanceW::openFrac( int idResSgn1, int idResSgn2) {
  
  // Only one resonance.
  if (idResSgn2 == 0) return (idResSgn1 > 0) ? openPos : openNeg;

  // Pair of resonances.
  if (idResSgn1 * idResSgn2 < 0) return openPos * openNeg;
  return (idResSgn1 > 0) ? pow2(openPos) : pow2(openNeg);  

}

//*********

// Calculate the total width and store phases-space-weighted coupling sums.

double ResonanceW::width(double mHat, int idResSgn, int , 
  bool openOnly, bool setBR) {

  // Common coupling factors.
  double sH     = mHat*mHat;
  double alpEM  = alphaEM.alphaEM(sH); 
  double alpS   = alphaS.alphaS(sH);
  double colQ   = 3. * (1. + alpS / M_PI);
  double preFac = alpEM * thetaWRat * mHat;

  // Reset quantities to sum. Declare variables inside loop.
  if (setBR) particlePtr->decay.resetDynamicBR(); 
  double widSum = 0.; 
  int    id1Abs, id2Abs, onMode;
  double widNow, mf1, mf2, mr1, mr2, ps;

  // Loop over all decay channels; optionally only open channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow = 0.;
    onMode = particlePtr->decay[i].onMode();
    if (openOnly) {
      if (idResSgn > 0 && onMode !=1 && onMode != 2) continue;
      if (idResSgn < 0 && onMode !=1 && onMode != 3) continue;
    }
    id1Abs = abs( particlePtr->decay[i].product(0) );
    id2Abs = abs( particlePtr->decay[i].product(1) );

    // Only contributions from three fermion generations, except top.
    if ( (id1Abs > 0 && id1Abs < 6 && id2Abs > 0 && id2Abs < 6) 
      || (id1Abs > 10 && id1Abs < 17 && id2Abs > 10 && id2Abs < 17)) {
      mf1 = ParticleDataTable::m0(id1Abs);
      mf2 = ParticleDataTable::m0(id2Abs);

      // Check that above threshold. Phase space.
      if (mHat > mf1 + mf2 + MASSMARGIN) {
        mr1 = pow2(mf1 / mHat);
        mr2 = pow2(mf2 / mHat);
        ps = (1. - 0.5 * (mr1 + mr2) - 0.5 * pow2(mr1 - mr2))
          * sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 ); 

        // Combine phase space with colour factor and CKM couplings.
        widNow = preFac * ps;
        if (id1Abs < 6) widNow *= colQ * VCKM::V2id(id1Abs, id2Abs);
      }
    }

    // Sum back up.
    widSum += widNow;

    // Optionally store partial widths for later decay channel choice.
    if (setBR) particlePtr->decay[i].dynamicBR(widNow); 
  }

  // Done.
  return widSum;
  
}
 
//**************************************************************************

// The ResonanceTop class.
// Derived class for top/antitop properties.

//*********
 
// Definitions of static variables and functions.
int    ResonanceTop::idRes = 6;
double ResonanceTop::mRes;
double ResonanceTop::GammaRes;
double ResonanceTop::m2Res;
double ResonanceTop::GamMRat;
double ResonanceTop::thetaWRat;
double ResonanceTop::m2W;
double ResonanceTop::openPos = 1.; 
double ResonanceTop::openNeg = 1.; 
ParticleDataEntry* ResonanceTop::particlePtr;

//*********

// Initialize static data members.

void ResonanceTop::initStatic() {

  // Set pointer to particle properties and decay table.
  particlePtr = ParticleDataTable::particleDataPtr(idRes);

  // Top properties and other constants.
  mRes      = ParticleDataTable::m0(idRes);
  GammaRes  = ParticleDataTable::mWidth(idRes);
  m2Res     = mRes*mRes;
  GamMRat   = GammaRes / mRes;
  thetaWRat = 1. / (16. * CoupEW::sin2thetaW());
  m2W       = pow2(ParticleDataTable::m0(24));

}

//*********

// Calculate and store partial and total widths at the nominal mass. 

void ResonanceTop::widthInit() {

  // Common coupling factors.
  double mHat   = mRes;
  double sH     = mHat*mHat;
  double alpEM  = alphaEM.alphaEM(sH);
  double alpS   = alphaS.alphaS(sH);
  double colQ   = 1. - 2.5 * alpS / M_PI;
  double preFac = alpEM * thetaWRat * (sH / m2W) * mHat;

  // Reset quantities to sum. Declare others.
  double widTot = 0.; 
  double widPos = 0.;
  double widNeg = 0.;
  int    id1Abs, id2Abs, onMode;
  double widNow, widSecPos, widSecNeg, mf1, mf2, mr1, mr2, ps;

  // Loop over all decay channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow    = 0.;
    widSecPos = 1.;
    widSecNeg = 1.;
    id1Abs = abs( particlePtr->decay[i].product(0) );
    id2Abs = abs( particlePtr->decay[i].product(1) );
    if (id2Abs > id1Abs) swap( id1Abs, id2Abs);

    // Only contributions from W + d/s/b.
    if ( id1Abs == 24 && (id2Abs == 1 || id2Abs == 3 || id2Abs ==5) ) { 
      mf1 = ParticleDataTable::m0(id1Abs);
      mf2 = ParticleDataTable::m0(id2Abs);

      // Check that above threshold. Phase space.
      if (mHat > mf1 + mf2 + MASSMARGIN) {
        mr1 = pow2(mf1 / mHat);
        mr2 = pow2(mf2 / mHat);
        ps = ( pow2(1. - mr2) + (1. + mr2) * mr1 - 2. * mr1 * mr1 )
           * sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 ); 

        // Combine phase space with colour factor and CKM couplings.
        widNow = preFac * ps * colQ * VCKM::V2id(6, id2Abs);
 
        // Secondary width from W+ and W- decay.
        widSecPos = ResonanceW::openFrac(24); 
        widSecNeg = ResonanceW::openFrac(-24); 
      }
    }

    // Store partial widths and update sum
    particlePtr->decay[i].onShellWidth(widNow); 
    widTot += widNow;    
    onMode = particlePtr->decay[i].onMode();
    if (onMode == 1 || onMode == 2) widPos += widNow * widSecPos;
    if (onMode == 1 || onMode == 3) widNeg += widNow * widSecNeg;
  }

  // Update total width and branching ratios.
  particlePtr->setMWidth(widTot, false);
  double bRatio;
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    bRatio = particlePtr->decay[i].onShellWidth() / widTot;
    particlePtr->decay[i].bRatio( bRatio, false);
  }
  particlePtr->setIsInitResonance();

  // Update the width in this object itself. Secondary widths for open.
  GammaRes = widTot;
  GamMRat  = GammaRes / mRes;  
  openPos  = (widPos / widTot);
  openNeg  = (widNeg / widTot);

}

//*********

// Return fraction of width open for one/two-resonance combinations.
  
double ResonanceTop::openFrac( int idResSgn1, int idResSgn2) {
  
  // Only one resonance.
  if (idResSgn2 == 0) return (idResSgn1 > 0) ? openPos : openNeg;

  // Pair of resonances.
  if (idResSgn1 * idResSgn2 < 0) return openPos * openNeg;
  return (idResSgn1 > 0) ? pow2(openPos) : pow2(openNeg);  

}

//*********

// Calculate the total width and store phases-space-weighted coupling sums.

double ResonanceTop::width(double mHat, int idResSgn, int , 
  bool openOnly, bool setBR) {

  // Common coupling factors.
  double sH     = mHat*mHat;
  double alpEM  = alphaEM.alphaEM(sH); 
  double alpS   = alphaS.alphaS(sH);
  double colQ   = 3. * (1. + alpS / M_PI);
  double preFac = alpEM * thetaWRat * (sH / m2W) * mHat;
 
  // Reset quantities to sum. Declare variables inside loop.
  if (setBR) particlePtr->decay.resetDynamicBR(); 
  double widSum = 0.; 
  int    id1Abs, id2Abs, onMode;
  double widNow, mf1, mf2, mr1, mr2, ps;

  // Loop over all decay channels; optionally only open channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow = 0.;
    onMode = particlePtr->decay[i].onMode();
    if (openOnly) {
      if (idResSgn > 0 && onMode !=1 && onMode != 2) continue;
      if (idResSgn < 0 && onMode !=1 && onMode != 3) continue;
    }
    id1Abs = abs( particlePtr->decay[i].product(0) );
    id2Abs = abs( particlePtr->decay[i].product(1) );
    if (id2Abs > id1Abs) swap( id1Abs, id2Abs);

    // Only contributions from W + d/s/b.
    if ( id1Abs == 24 && (id2Abs == 1 || id2Abs == 3 || id2Abs ==5) ) { 
      mf1 = ParticleDataTable::m0(id1Abs);
      mf2 = ParticleDataTable::m0(id2Abs);

      // Check that above threshold. Phase space.
      if (mHat > mf1 + mf2 + MASSMARGIN) {
        mr1 = pow2(mf1 / mHat);
        mr2 = pow2(mf2 / mHat);
        ps = ( pow2(1. - mr2) + (1. + mr2) * mr1 - 2. * mr1 * mr1 )
           * sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 ); 

        // Combine phase space with colour factor and CKM couplings.
        widNow = preFac * ps * colQ * VCKM::V2id(6, id2Abs);

        // Secondary width from W decay.
        if (openOnly) widNow *= ResonanceW::openFrac(idResSgn); 
      }
    }

    // Sum back up.
    widSum += widNow;

    // Optionally store partial widths for later decay channel choice.
    if (setBR) particlePtr->decay[i].dynamicBR(widNow); 
  }

  // Done.
  return widSum;
  
}

//*********

// Evaluate weight for W decay distribution in t -> W b -> f fbar b.

double ResonanceTop::weightDecayAngles( Event& process, int iResBeg, 
  int iResEnd) {

  // If not pair W d/s/b and mother t then return unit weight.
  if (iResEnd - iResBeg != 2) return 1.;
  int iW = iResBeg;
  int iB = iResBeg + 1;
  int idW = process[iW].idAbs();
  int idB = process[iB].idAbs();
  if (idW != 24) {
    swap(iW, iB); 
    swap(idW, idB);
  } 
  if (idW != 24 || (idB != 1 && idB != 3 && idB != 5)) return 1.;
  int iT = process[iW].mother1(); 
  if (iT <= 0 || process[iT].idAbs() != 6) return 1.;

  // Find sign-matched order of W decay products. 
  int iF    = process[iW].daughter1(); 
  int iFbar = process[iW].daughter2();
  if (iFbar - iF != 1) return 1.; 
  if (process[iT].id() * process[iF].id() < 0) swap(iF, iFbar);

  // Weight and maximum weight.
  double wt    = (process[iT].p() * process[iFbar].p()) 
               * (process[iF].p() * process[iB].p());
  double wtMax = ( pow4(process[iT].m()) - pow4(process[iW].m()) ) / 8.;  

  // Done.
  return wt / wtMax;

}
 
//**************************************************************************

// The ResonanceSMH class.
// Derived class for SM Higgs properties.

//*********
 
// Definitions of static variables and functions.
int    ResonanceSMH::idRes = 25;
bool   ResonanceSMH::linearWidthWWZZ;
double ResonanceSMH::mRes;
double ResonanceSMH::GammaRes;
double ResonanceSMH::m2Res;
double ResonanceSMH::GamMRat;
double ResonanceSMH::sin2tW;
double ResonanceSMH::cos2tW;
double ResonanceSMH::mZ;
double ResonanceSMH::mW;
double ResonanceSMH::GammaZ;
double ResonanceSMH::GammaW;
double ResonanceSMH::GammaT;
double ResonanceSMH::widTable[25]     = {0.};
int    ResonanceSMH::SMHiggsParity = 1;
double ResonanceSMH::SMHiggsEta    = 0.;
double ResonanceSMH::openPos = 1.; 
double ResonanceSMH::openNeg = 0.; 
ParticleDataEntry* ResonanceSMH::particlePtr;

// Minimal mass for W, Z, top in integration over respective Breit-Wigner.
const double ResonanceSMH::MASSMIN = 10.;

// Number of widths above threshold where B-W integration not needed.
const double ResonanceSMH::GAMMAMARGIN = 10.;

//*********

// Initialize static data members.

void ResonanceSMH::initStatic() {

  // Set pointer to particle properties and decay table.
  particlePtr = ParticleDataTable::particleDataPtr(idRes);

  // Settings.
  linearWidthWWZZ = Settings::flag("ResonanceSMH:linearWidthWWZZ");

  // H properties, or other properties used to derive those.
  mRes      = ParticleDataTable::m0(idRes);
  GammaRes  = ParticleDataTable::mWidth(idRes);
  m2Res     = mRes*mRes;
  GamMRat   = GammaRes / mRes;
  sin2tW    = CoupEW::sin2thetaW();
  cos2tW    = 1. - sin2tW;
  mZ        = ParticleDataTable::m0(23);
  mW        = ParticleDataTable::m0(24);
  GammaZ    = ParticleDataTable::mWidth(23);
  GammaW    = ParticleDataTable::mWidth(24);
  GammaT    = ParticleDataTable::mWidth(6);

  // Higgs parity assumption.
  SMHiggsParity  = Settings::mode("ResonanceSMH:parity");
  SMHiggsEta     = Settings::parm("ResonanceSMH:etaParity");

}

//*********

// Calculate and store partial and total widths at the nominal mass. 

void ResonanceSMH::widthInit() {

  // Common coupling factors.
  double mHat   = mRes;
  double sH     = mHat*mHat;
  double alpEM  = alphaEM.alphaEM(sH);
  double alpS   = alphaS.alphaS(sH);
  double colQ   = 3. * (1. + alpS / M_PI);
  double preFac = (alpEM / (8. * sin2tW)) * pow3(mHat) / pow2(mW); 

  // Reset quantities to sum. Declare others.
  double widTot = 0.; 
  double widPos = 0.;
  int    id1Abs, id2Abs, onMode;
  double widNow, widSec, mf1, mr1, ps;

  // Loop over all decay channels.
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow = 0.;
    widSec = 1.;
    id1Abs = abs( particlePtr->decay[i].product(0) );
    id2Abs = abs( particlePtr->decay[i].product(1) );
    mf1    = ParticleDataTable::m0(id1Abs);
    mr1    = pow2(mf1 / mHat);

    // Widths of decays H0 -> f + fbar.
    if ( id2Abs == id1Abs && ( (id1Abs > 0 && id1Abs < 7) 
      || (id1Abs > 10 && id1Abs < 17) ) ) {
      ps  = 0.;

      // Check that above threshold. Phase space.
      if (id1Abs != 6 && mHat > 2. * mf1 + MASSMARGIN)
        ps = sqrtpos(1. - 4. * pow2(mf1 / mHat)); 

      // For top near threshold use numerical integration.     
      else if (id1Abs == 6 && mHat > 2. * (mf1 - GAMMAMARGIN * GammaT)) 
        ps = ( mHat > 2. * (mf1 + GAMMAMARGIN * GammaT) ) 
          ? sqrtpos(1. - 4. * mr1) 
          : numInt2BW( mHat, mf1, GammaT, MASSMIN, mf1, GammaT, MASSMIN, 0);

      // Combine phase space with colour factor and running masses.
      double ratioRun = pow2(ParticleDataTable::mRun(id1Abs, mHat) / mHat);
      if (ps > 0.) widNow = preFac * ps * ratioRun; 
      if (id1Abs < 7) widNow *= colQ;

      // Secondary width from top decay.
      if (id1Abs == 6) widSec = ResonanceTop::openFrac(6, -6); 

      // Store widths for in-state: no phase space and no colour factors.
      widTable[id1Abs] = preFac * ratioRun;
    }

    // Widths of decays H0 -> g + g. No colour factor for instate.
    else if (id1Abs == 21 && id2Abs == 21) { 
      widNow    = preFac * pow2(alpS / M_PI) * eta2gg(mHat); 
      widTable[21] = widNow / 8.;
    }

    // Widths of decays H0 -> gamma + gamma.
    else if (id1Abs == 22 && id2Abs == 22) { 
      widNow    = preFac * pow2(alpEM / M_PI) * 0.5 * eta2gaga(mHat); 
      widTable[22] = widNow;
    }
 
    // Widths of decays H0 -> gamma + Z0.
    else if (id1Abs == 22 && id2Abs == 23 && mHat > mZ) {
      widNow = preFac * pow2(alpEM / M_PI) * eta2gaZ(mHat)
        * pow3(1. - pow2(mZ / mHat)); 

      // Secondary width from Z0 decay.
      widSec = ResonanceGmZ::openFrac(23); 
    }
 
    // Widths of decays H0 -> Z0 + Z0.
    else if (id1Abs == 23 && id2Abs == 23) {
      // If H0 heavy use on-shell expression, else numerical integration.
      widNow = ( mHat > 2. * (mZ + GAMMAMARGIN * GammaZ) ) 
        ? (1. - 4. * mr1 + 12. * pow2(mr1)) * sqrtpos(1. - 4. * mr1) 
        : numInt2BW( mHat, mZ, GammaZ, MASSMIN, mZ, GammaZ, MASSMIN, 1);
      widNow   *= 0.25 * preFac;
      widTable[23] = 0.25 * preFac;

      // Secondary width from Z0 Z0 decay.
      widSec = ResonanceGmZ::openFrac(23, 23); 
    }
 
    // Widths of decays H0 -> W+ + W-.
    else if (id1Abs == 24 && id2Abs == 24) {
      // If H0 heavy use on-shell expression, else numerical integration.
      widNow = (mHat > 2. * ( mW + GAMMAMARGIN * GammaW)) 
        ? (1. - 4. * mr1 + 12. * pow2(mr1)) * sqrtpos(1. - 4. * mr1) 
        : numInt2BW( mHat, mW, GammaW, MASSMIN, mW, GammaW, MASSMIN, 1);
      widNow   *= 0.5 * preFac;
      widTable[24] = 0.5 * preFac;

      // Secondary width from W+ W- decay.
      widSec = ResonanceW::openFrac(24, -24); 
    }

    // Store partial widths and update sum, also for open channels only.
    particlePtr->decay[i].onShellWidth(widNow); 
    widTot += widNow;    
    onMode = particlePtr->decay[i].onMode();
    if (onMode == 1 || onMode == 2) widPos += widNow * widSec;
  }

  // Update total width and branching ratios.
  particlePtr->setMWidth(widTot, false);
  double bRatio;
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    bRatio = particlePtr->decay[i].onShellWidth() / widTot;
    particlePtr->decay[i].bRatio( bRatio, false);
  }
  particlePtr->setIsInitResonance();

  // Update the width in this object itself. Secondary widths for open.
  GammaRes = widTot;
  GamMRat  = GammaRes / mRes;  
  openPos  = widPos / widTot;

}

//*********
 
// Calculate the partial width for a specific mass and in/out-state.
// This width has no phase-space suppression and is given without
// colour factors. Uses pretabulation from widthInit() to speed up.

double ResonanceSMH::widthChan(double mHat, int id1Abs, int ) { 

  // Read out from stored values.
  double widNow = (id1Abs > 0 && id1Abs < 25) ? widTable[id1Abs] : 0.;

  // Simple rescaling with third power of actual mass.
  widNow *= pow3(mHat / mRes);

  // Normally replace mHat^3 by mHat * mHiggs^2 for WW/ZZ decays.
  if ( linearWidthWWZZ && (id1Abs == 23 || id1Abs == 24) )
    widNow *= pow2(mRes / mHat);   

  // Done.
  return widNow;

}

//*********

// Return fraction of width open for one/two-resonance combinations.
  
double ResonanceSMH::openFrac( int idResSgn1, int idResSgn2) {
  
  // Only one resonance.
  if (idResSgn2 == 0) return (idResSgn1 > 0) ? openPos : openNeg;

  // Pair of resonances.
  if (idResSgn1 * idResSgn2 < 0) return openPos * openNeg;
  return (idResSgn1 > 0) ? pow2(openPos) : pow2(openNeg);  

}

//*********
 
// Calculate the total or open width at a given energy.
// Has approximate rescaling since the correct width calculation 
// (in widthInit()) is slow owing to numerical integrations.

double ResonanceSMH::width(double mHat, int , int , 
  bool openOnly, bool setBR) { 

  // Reset quantities to sum. Declare variables inside loop.
  if (setBR) particlePtr->decay.resetDynamicBR(); 
  double widSum = 0.; 
  int    id1Abs, id2Abs, onMode;
  double widNow, widSec;

  // Loop over all decay channels; optionally only open channels.
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow = 0.;
    widSec = 1.;
    onMode = particlePtr->decay[i].onMode();
    if (openOnly && onMode !=1 && onMode != 2) continue;
    
    // Start out from on-shell partial width.     
    widNow = particlePtr->decay[i].onShellWidth();

    // Simple rescaling with third power of actual mass. (Primitive??)
    widNow *= pow3(mHat / mRes);

    // Normally replace mHat^3 by mHat * mHiggs^2 for WW/ZZ decays.
    id1Abs = abs( particlePtr->decay[i].product(0) );
    id2Abs = abs( particlePtr->decay[i].product(1) );
    if (linearWidthWWZZ && (id1Abs == 23 || id1Abs == 24) 
      && id2Abs == id1Abs) widNow *= pow2(mRes / mHat);   

    // Secondary width from t tbar, gamma Z0, Z0 Z0 and W+ W- decays.
    if (openOnly) {
      if (id1Abs  == 6 && id2Abs ==  6) widSec = ResonanceTop::openFrac(6, -6);
      if (id1Abs == 22 && id2Abs == 23) widSec = ResonanceGmZ::openFrac(23);
      if (id1Abs == 23 && id2Abs == 23) widSec = ResonanceGmZ::openFrac(23, 23); 
      if (id1Abs == 24 && id2Abs == 24) widSec = ResonanceW::openFrac(24, -24);
    }
    widNow *= widSec;

    // Sum back up.
    widSum += widNow;

    // Optionally store partial widths for later decay channel choice.
    if (setBR) particlePtr->decay[i].dynamicBR(widNow); 
  }

  // Done.
  return widSum;

} 

//*********

// Evaluate weight for Z0/W+- decay distributions in H -> Z0/W+ Z0/W- -> 4f.

double ResonanceSMH::weightDecayAngles( Event& process, int iResBeg, 
  int iResEnd) {

  // If not pair Z0 Z0 or W+ W- or not mother Higgs then return unit weight.
  if (iResEnd - iResBeg != 2) return 1.;
  int iZW1  = iResBeg;
  int iZW2  = iResBeg + 1;
  int idZW1 = process[iZW1].id();
  int idZW2 = process[iZW2].id();
  if (idZW1 < 0) {
    swap(iZW1, iZW2); 
    swap(idZW1, idZW2);
  } 
  if ( (idZW1 != 23 || idZW2 != 23) && (idZW1 != 24 || idZW2 != -24) )
    return 1.;
  int iH = process[iZW1].mother1(); 
  if (iH <= 0 || process[iH].id() != 25) return 1.;

  // Option with isotropic decays.
  if (SMHiggsParity == 0) return 1.;

  // Maximum and initial weight. 
  double wtMax = pow4(process[iH].m());
  double wt    = wtMax; 

  // Find sign-matched order of Z0/W+- decay products. 
  int i3 = process[iZW1].daughter1();
  int i4 = process[iZW1].daughter2();
  if (process[i3].id() < 0) swap( i3, i4); 
  int i5 = process[iZW2].daughter1();
  int i6 = process[iZW2].daughter2();
  if (process[i5].id() < 0) swap( i5, i6); 

  // Evaluate four-vector products and find masses..
  double p35  = 2. * process[i3].p() * process[i5].p(); 
  double p36  = 2. * process[i3].p() * process[i6].p(); 
  double p45  = 2. * process[i4].p() * process[i5].p(); 
  double p46  = 2. * process[i4].p() * process[i6].p(); 
  double p34  = 2. * process[i3].p() * process[i4].p(); 
  double p56  = 2. * process[i5].p() * process[i6].p(); 
  double mZW1 = process[iZW1].m();
  double mZW2 = process[iZW2].m();

  // For mixed CP states need epsilon product and gauge boson masses.
  double epsilonProd = 0.;
  if (SMHiggsParity == 3) {
    double p[4][4];
    for (int i = 0; i < 4; ++i) {
      int         ii = i3;
      if (i == 1) ii = i4;
      if (i == 2) ii = i5;
      if (i == 3) ii = i6;
      p[i][0] = process[ii].e();
      p[i][1] = process[ii].px();
      p[i][2] = process[ii].py();
      p[i][3] = process[ii].pz();
    }     
    epsilonProd 
      = p[0][0]*p[1][1]*p[2][2]*p[3][3] - p[0][0]*p[1][1]*p[2][3]*p[3][2] 
      - p[0][0]*p[1][2]*p[2][1]*p[3][3] + p[0][0]*p[1][2]*p[2][3]*p[3][1]
      + p[0][0]*p[1][3]*p[2][1]*p[3][2] - p[0][0]*p[1][3]*p[2][2]*p[3][1]
      - p[0][1]*p[1][0]*p[2][2]*p[3][3] + p[0][1]*p[1][0]*p[2][3]*p[3][2]
      + p[0][1]*p[1][2]*p[2][0]*p[3][3] - p[0][1]*p[1][2]*p[2][3]*p[3][0]
      - p[0][1]*p[1][3]*p[2][0]*p[3][2] + p[0][1]*p[1][3]*p[2][2]*p[3][0]
      + p[0][2]*p[1][0]*p[2][1]*p[3][3] - p[0][2]*p[1][0]*p[2][3]*p[3][1]
      - p[0][2]*p[1][1]*p[2][0]*p[3][3] + p[0][2]*p[1][1]*p[2][3]*p[3][0] 
      + p[0][2]*p[1][3]*p[2][0]*p[3][1] - p[0][2]*p[1][3]*p[2][1]*p[3][0]
      - p[0][3]*p[1][0]*p[2][1]*p[3][2] + p[0][3]*p[1][0]*p[2][2]*p[3][1] 
      + p[0][3]*p[1][1]*p[2][0]*p[3][2] - p[0][3]*p[1][1]*p[2][2]*p[3][0] 
      - p[0][3]*p[1][2]*p[2][0]*p[3][1] + p[0][3]*p[1][2]*p[2][1]*p[3][0];
  }

  // Z0 Z0 decay: vector and axial couplings of two fermion pairs.
  if (idZW1 == 23) {
    double vf1 = CoupEW::vf(process[i3].idAbs());
    double af1 = CoupEW::af(process[i3].idAbs());
    double vf2 = CoupEW::vf(process[i5].idAbs());
    double af2 = CoupEW::af(process[i5].idAbs());
    double va12asym = 4. * vf1 * af1 * vf2 * af2 
      / ( (vf1*vf1 + af1*af1) * (vf2*vf2 + af2*af2) );
    double etaMod = SMHiggsEta / pow2( ParticleDataTable::m0(23) );
    
    // Normal CP-even decay.
    if (SMHiggsParity == 1) wt = 8. * (1. + va12asym) * p35 * p46 
      + 8. * (1. - va12asym) * p36 * p45;

    // CP-odd decay.
    else if (SMHiggsParity == 2) wt = ( pow2(p35 + p46) 
      + pow2(p36 + p45) - 2. * p34 * p56
      - 2. * pow2(p35 * p46 - p36 * p45) / (p34 * p56) 
      + va12asym * (p35 + p36 - p45 - p46) * (p35 + p45 - p36 - p46) )
      / (1. +  va12asym);

    // Mixed CP states. 
    else wt = 32. * ( 0.25 * ( (1. + va12asym) * p35 * p46 
      + (1. - va12asym) * p36 * p45 ) - 0.5 * etaMod * epsilonProd
      * ( (1. + va12asym) * (p35 + p46) - (1. - va12asym) * (p36 + p45) )
      + 0.0625 * etaMod * etaMod * (-2. * pow2(p34 * p56) 
      - 2. * pow2(p35 * p46 - p36 * p45) 
      + p34 * p56 * (pow2(p35 + p46) + pow2(p36 + p45)) 
      + va12asym * p34 * p56 * (p35 + p36 - p45 - p46) 
      * (p35 + p45 - p36 - p46) ) ) / ( 1. * 2. * etaMod * mZW1 * mZW2 
      + 2. * pow2(etaMod * mZW1 * mZW2) * (1. + va12asym) );

  // W+ W- decay.
  } else if (idZW1 == 24) {
    double etaMod = SMHiggsEta / pow2( ParticleDataTable::m0(24) );
    
    // Normal CP-even decay.
    if (SMHiggsParity == 1) wt = 16. * p35 * p46; 

    // CP-odd decay.
    else if (SMHiggsParity == 2) wt = 0.5 * ( pow2(p35 + p46) 
      + pow2(p36 + p45) - 2. * p34 * p56  
      - 2. * pow2(p35 * p46 - p36 * p45) / (p34 * p56) 
      + (p35 + p36 - p45 - p46) * (p35 + p45 - p36 - p46) );

    // Mixed CP states. 
    else wt = 32. * ( 0.25 * 2. * p35 * p46 
      - 0.5 * etaMod * epsilonProd * 2. * (p35 + p46)
      + 0.0625 * etaMod * etaMod * (-2. * pow2(p34 * p56) 
      - 2. * pow2(p35 * p46 - p36 * p45) 
      + p34 * p56 * (pow2(p35 + p46) + pow2(p36 + p45)) 
      + p34 * p56 * (p35 + p36 - p45 - p46) * (p35 + p45 - p36 - p46) ) ) 
      / ( 1. * 2. * etaMod * mZW1 * mZW2 + 2. * pow2(etaMod * mZW1 * mZW2) );
  }

  // Done.
  return wt / wtMax;

}

//*********

// Sum up quark loop contributions in H0 -> g + g.

double ResonanceSMH::eta2gg(double mHat) {

  // Initial values.
  complex eta = complex(0., 0.);
  double epsilon, root, rootLog;
  complex phi, etaNow;

  // Loop over s, c, b, t quark flavours.
  for (int id = 3; id < 7; ++id) {
    epsilon = pow2(2. * ParticleDataTable::mRun(id, mHat) / mHat);

    // Value of loop integral.
    if (epsilon <= 1.) {
      if (epsilon < 1e-4) rootLog = log(4. / epsilon - 2.);
      else { 
        root    = sqrt(1. - epsilon);
        rootLog = log( (1. + root) / (1. - root) );
      }
      phi = complex( -0.25 * (pow2(rootLog) - pow2(M_PI)), 
                     0.5 * M_PI * rootLog );
    } 
    else phi = complex( pow2( asin(1. / sqrt(epsilon)) ), 0.);
    etaNow = -0.5 * epsilon * (complex(1., 0.) + (1. - epsilon) * phi);
    
    // Sum up contribution and return square of absolute value.
    eta += etaNow;
  }
  return (pow2(eta.real()) + pow2(eta.imag()));

}

//*********

// Sum up quark, lepton and W+- loop contributions in H0 -> gamma + gamma.

double ResonanceSMH::eta2gaga(double mHat) {

  // Initial values.
  complex eta = complex(0., 0.);
  int id;
  double ef, epsilon, root, rootLog;
  complex phi, etaNow;

  // Loop over s, c, b, t, mu , tau, W flavours.
  for (int idLoop = 0; idLoop < 7; ++idLoop) {
    if      (idLoop < 4) id = idLoop + 3;
    else if (idLoop < 6) id = 2 * idLoop + 5;
    else                 id = 24;
 
    // Charge and loop integral parameter.
    ef      = (id < 20) ? CoupEW::ef(id) : 1.;
    epsilon = pow2(2. * ParticleDataTable::mRun(id, mHat) / mHat);

    // Value of loop integral.
    if (epsilon <= 1.) {
      if (epsilon < 1e-4) rootLog = log(4. / epsilon - 2.);
      else { 
        root    = sqrt(1. - epsilon);
        rootLog = log( (1. + root) / (1. - root) );
      }
      phi = complex( -0.25 * (pow2(rootLog) - pow2(M_PI)), 
                     0.5 * M_PI * rootLog );
    } 
    else phi = complex( pow2( asin(1. / sqrt(epsilon)) ), 0.);

    // Separate expressions for quarks, leptons and W.
    if (idLoop < 6) { 
      etaNow = -0.5 * epsilon * (complex(1., 0.) + (1. - epsilon) * phi);
      if (idLoop < 4) etaNow *= 3. * pow2(ef);
      else            etaNow *=      pow2(ef);
    } 
    else etaNow = complex(0.5 + 0.75 * epsilon, 0.)
                + 0.75 * epsilon * (2. - epsilon) * phi;       
    
    // Sum up contribution and return square of absolute value.
    eta += etaNow;
  }
  return (pow2(eta.real()) + pow2(eta.imag()));

}

//*********

// Sum up quark, lepton and W+- loop contributions in H0 -> gamma + Z0.

double ResonanceSMH::eta2gaZ(double mHat) {

  // Initial values.
  complex eta = complex(0., 0.);
  int id;
  double ef, vf, mRun, epsilon, epsPrime, root, rootLog, asinEps;
  complex phi, psi, phiPrime, psiPrime, fXY, f1, etaNow;

  // Loop over s, c, b, t, mu , tau, W flavours.
  for (int idLoop = 0; idLoop < 7; ++idLoop) {
    if      (idLoop < 4) id = idLoop + 3;
    else if (idLoop < 6) id = 2 * idLoop + 5;
    else                 id = 24;

    // Electroweak charges and loop integral parameters.
    ef       = (id < 20) ? CoupEW::ef(id) : 1.;
    vf       = (id < 20) ? CoupEW::vf(id) : 0.;
    mRun     = ParticleDataTable::mRun(id, mHat);
    epsilon  = pow2(2. * mRun / mHat);
    epsPrime = pow2(2. * mRun / mZ);

    // Value of loop integral for epsilon = 4 m^2 / sHat.
    if (epsilon <= 1.) {
      root    = sqrt(1. - epsilon);
      rootLog = (epsilon < 1e-4) ? log(4. / epsilon - 2.)
                : log( (1. + root) / (1. - root) );
      phi = complex( -0.25 * (pow2(rootLog) - pow2(M_PI)), 
                     0.5 * M_PI * rootLog );
      psi = 0.5 * root * complex( rootLog, -M_PI); 
    } else {
      asinEps = asin(1. / sqrt(epsilon));
      phi = complex( pow2(asinEps), 0.);
      psi = complex( sqrt(epsilon - 1.) * asinEps, 0.);
    }

    // Value of loop integral for epsilonPrime = 4 m^2 / m_Z^2.
    if (epsPrime <= 1.) {
      root    = sqrt(1. - epsPrime);
      rootLog = (epsPrime < 1e-4) ? log(4. / epsPrime - 2.)
              : log( (1. + root) / (1. - root) );
      phiPrime = complex( -0.25 * (pow2(rootLog) - pow2(M_PI)), 
                          0.5 * M_PI * rootLog );
      psiPrime = 0.5 * root * complex( rootLog, -M_PI); 
    } else {
      asinEps = asin(1. / sqrt(epsPrime));
      phiPrime = complex( pow2(asinEps), 0.);
      psiPrime = complex( sqrt(epsPrime - 1.) * asinEps, 0.);
    }

    // Combine the two loop integrals.
    fXY = (epsilon * epsPrime / (8. * pow2(epsilon - epsPrime)))
      * ( complex(epsilon - epsPrime, 0) 
      + epsilon * epsPrime * (phi - phiPrime)
      + 2. * epsilon * (psi - psiPrime) );
    f1 = - (epsilon * epsPrime / (2. * (epsilon - epsPrime)))
      * (phi - phiPrime);    

    // Separate expressions for quarks, leptons and W.
    if (idLoop < 6) { 
      etaNow = -fXY + 0.25 * f1;
      if (idLoop < 4) etaNow *= 3. * ef * vf;
      else            etaNow *=      ef * vf;
    } else {
      double coef1  = 3. - sin2tW / cos2tW;
      double coefXY = (1. + 2. / epsilon) * sin2tW / cos2tW 
        - (5. + 2. / epsilon);
      etaNow = -cos2tW * (coef1 * f1 + coefXY * fXY); 
    }
    
    // Sum up contribution and return square of absolute value.
    eta += etaNow;
  }
  return ( (pow2(eta.real()) + pow2(eta.imag())) / (sin2tW * cos2tW) );

}

//**************************************************************************

} // end namespace Pythia8
