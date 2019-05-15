// ResonanceWidths.cc is a part of the PYTHIA event generator.
// Copyright (C) 2007 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for 
// the ResonanceWidths class and classes derived from it.

#include "ResonanceWidths.h"
#include "PythiaComplex.h"

namespace Pythia8 {

//**************************************************************************

// The ResonanceWidths class.
// Base class for the various resonances.

//*********
 
// Definitions of static variables and functions.
// (Values will be overwritten in initStatic call, so are purely dummy.)

int          ResonanceWidths::alphaSorder    = 1;
int          ResonanceWidths::alphaEMorder   = 1;
double       ResonanceWidths::alphaSvalue    = 0.1265;
AlphaStrong  ResonanceWidths::alphaS;
AlphaEM      ResonanceWidths::alphaEM;

// Number of points in integration direction for numInt routines.
const int    ResonanceWidths::NPOINT         = 100;

// The sum of product masses must not be too close to the resonance mass.
const double ResonanceWidths::MASSMARGIN     = 0.1;

// If the total width of a particle is too low (or 0) then make it stable.
const double ResonanceWidths::MINWIDTH       = 1e-20;

// When meMode = 103 the beta at nominal mass should not be taken too tiny. 
const double ResonanceWidths::BETAMIN        = 0.1;

//*********

// Initialize static data members.

void ResonanceWidths::initStatic() {

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

// Initialize particle properties always present.

void ResonanceWidths::initBasic(ParticleDataEntry* particlePtrIn) {

  // Save pointer to particle species.
  if (particlePtrIn != 0) particlePtr  = particlePtrIn; 

  // Resonance propertiesL identity code, mass, width
  idRes        = particlePtr->id();
  mRes         = particlePtr->m0();
  GammaRes     = particlePtr->mWidth();
  m2Res        = mRes*mRes;

  // For very narrow resonances assign fictitious small width.
  if (GammaRes < MINWIDTH) GammaRes = 0.1 * MINWIDTH;  
  GamMRat      = GammaRes / mRes;

  // Secondary widths by default all on.
  openPos      = 1.;
  openNeg      = 1.;

  // Allow option where on-shell width is forced to current value.
  doForceWidth = particlePtr->doForceWidth();
  forceFactor  = 1.;
  
}

//*********

// Numerical integration of matrix-element in two-body decay,
// where one particle is described by a Breit-Wigner mass distribution.
// Normalization to unit integral if matrix element is unity 
// and there are no phase-space restrictions. 
 
double ResonanceWidths::numInt1BW(double mHat, double m1, double Gamma1, 
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
  double xStep    = 1. / NPOINT;

  // Variables used in loop over integration points.
  double sum      = 0.;
  double mr2      = pow2(m2 / mHat);
  double xNow1, sNow1, mNow1, mrNow1, ps, value; 

  // Loop with first-particle mass selection.
  for (int ip1 = 0; ip1 < NPOINT; ++ip1) {
    xNow1         = xStep * (ip1 + 0.5);
    sNow1         = s1 + mG1 * tan(atanMin1 + xNow1 * atanDif1);
    mNow1         = min( mMax1, max( mMin1, sqrtpos(sNow1) ) );
    mrNow1        = pow2(mNow1 / mHat);

    // Evaluate value and add to sum. Different matrix elements.
    ps            = sqrtpos( pow2(1. - mrNow1 - mr2) - 4. * mrNow1 * mr2);
    value         = 1.;
    if (meMode == 1) value = ps;
    if (meMode == 2) value = ps * ps;
    if (meMode == 3) value = pow3(ps);
    if (meMode == 5) value = ps * 
      (pow2(1. - mrNow1 - mr2) + 8. * mrNow1 * mr2);
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
 
double ResonanceWidths::numInt2BW(double mHat, double m1, double Gamma1, 
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
  double xStep    = 1. / NPOINT;
  int nIter       = (mustDiv) ? 2 * NPOINT : NPOINT;

  // Variables used in loop over integration points.
  double sum      = 0.;
  double xNow1, sNow1, mNow1, mrNow1, xNow2, sNow2, mNow2, mrNow2, ps, 
         value; 
  double wtNow1   = wtDif1; 
  double wtNow2   = wtDif2; 

  // Outer loop with first-particle mass selection.
  for (int ip1 = 0; ip1 < nIter; ++ip1) {
    if (!mustDiv) {
      xNow1       = xStep * (ip1 + 0.5);
      sNow1       = s1 + mG1 * tan(atanMin1 + xNow1 * atanDif1);
    } else if (ip1 < NPOINT) {
      xNow1       = xStep * (ip1 + 0.5);
      sNow1       = s1 + mG1 * tan(atanMin1 + xNow1 * atanDLo1);
      wtNow1      = wtDLo1;
    } else {
      xNow1       = xStep * (ip1 - NPOINT + 0.5);
      sNow1       = s1 + mG1 * tan(atanDiv1 + xNow1 * atanDHi1);
      wtNow1      = wtDHi1;
    }
    mNow1         = min( mMax1, max( mMin1, sqrtpos(sNow1) ) );
    mrNow1        = pow2(mNow1 / mHat);

    // Inner loop with second-particle mass selection.
    for (int ip2 = 0; ip2 < nIter; ++ip2) {
      if (!mustDiv) {
        xNow2     = xStep * (ip2 + 0.5);
        sNow2     = s2 + mG2 * tan(atanMin2 + xNow2 * atanDif2); 
      } else if (ip2 < NPOINT) {
        xNow2     = xStep * (ip2 + 0.5);
        sNow2     = s2 + mG2 * tan(atanMin2 + xNow2 * atanDLo2);
        wtNow2    = wtDLo2;
      } else {
        xNow2     = xStep * (ip2 - NPOINT + 0.5);
        sNow2     = s2 + mG2 * tan(atanDiv2 + xNow2 * atanDHi2);
        wtNow2    = wtDHi2;
      }
      mNow2       = min( mMax2, max( mMin2, sqrtpos(sNow2) ) );
      mrNow2      = pow2(mNow2 / mHat);

      // Check that point is inside phase space.
      if (mNow1 + mNow2 > mHat) break;

      // Evaluate value and add to sum. Different matrix elements.
      ps          = sqrtpos( pow2(1. - mrNow1 - mrNow2) 
                  - 4. * mrNow1 * mrNow2);
      value       = 1.;
      if      (meMode == 1) value = ps;
      else if (meMode == 2) value = ps * ps;
      else if (meMode == 3) value = pow3(ps);
      else if (meMode == 5) value = ps 
        * (pow2(1. - mrNow1 - mrNow2) + 8. * mrNow1 * mrNow2);
      sum        += value * wtNow1 * wtNow2;

    // End of second and first loop over integration points.
    }
  }  

  // Done.
  return sum;
}
 
//**************************************************************************

// The ResonanceNeutral class.
// Derived class for a generic resonance its own antiparticle, i.e. neutral.

//*********

// Initialize data members.
// Calculate and store partial and total widths at the nominal mass. 

void ResonanceNeutral::init() {

  // Do or redo initialization of basic quantities.
  initBasic();

  // Reset quantities to sum. Declare variables inside loop.
  double widTot = 0.; 
  double widPos = 0.;
  int    idNow, onMode;
  double widNow, widSec;

  // Loop over all decay channels. Find partial width.  
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow      = GammaRes * particlePtr->decay[i].bRatio();
    widSec      = 1.;

    // Find secondary widths.
    for (int j = 0; j < particlePtr->decay[i].multiplicity(); ++j) {
      idNow     = particlePtr->decay[i].product(j);
      widSec   *= ParticleDataTable::resOpenFrac(idNow); 
    }

    // Store partial widths and update sum, also for open channels only.
    particlePtr->decay[i].onShellWidth(widNow); 
    widTot += widNow;    
    onMode  = particlePtr->decay[i].onMode();
    if (onMode == 1 || onMode == 2) widPos += widNow * widSec;
  }

  // If no decay channels are open then set particle stable.
  bool mayDecay = (widTot > MINWIDTH); 
  if (!mayDecay) particlePtr->setMayDecay(false, false);

  // Secondary widths for open.
  openPos       = (mayDecay) ? (widPos / widTot) : 1.;
    
}

//*********

// Calculate the total width and store phase-space-weighted coupling sums.

double ResonanceNeutral::width(int , double mHat, int , 
  bool openOnly, bool setBR) {

  // Reset quantities to sum. Declare variables inside loop.
  double widSum  = 0.; 
  int    onMode, meMode, mult, idNow;
  double widNow, mfSum, mf1, mf2, mr1, mr2, betaf, betaN;

  // Loop over all decay channels; optionally only open channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow   = 0.;
    onMode   = particlePtr->decay[i].onMode();
    if (openOnly && onMode !=1 && onMode != 2) continue;
  
    // Matrix element mode and channel multiplicity.
    meMode   = particlePtr->decay[i].meMode();
    mult     = particlePtr->decay[i].multiplicity();

    // Correction by step at threshold.
    if (meMode == 101) {
      mfSum  = 0.;
      for (int j = 0; j < mult; ++j) mfSum 
            += ParticleDataTable::m0( particlePtr->decay[i].product(j) );
      if (mfSum < mHat) widNow = GammaRes * particlePtr->decay[i].bRatio(); 
    }  

    // Correction by a phase space factor for two-body decays.
    else if ( (meMode == 102 || meMode == 103) && mult == 2) { 
      mf1    = ParticleDataTable::m0( particlePtr->decay[i].product(0) );
      mf2    = ParticleDataTable::m0( particlePtr->decay[i].product(1) );
      mr1    = pow2(mf1 / mHat);    
      mr2    = pow2(mf2 / mHat);    
      betaf  = sqrtpos( pow2(1.- mr1 -mr2) - 4. * mr1 * mr2);
      mr1    = pow2(mf1 / mRes);   
      mr2    = pow2(mf2 / mRes);   
      betaN  = (meMode == 102) ? 1. : max( BETAMIN, 
               sqrtpos( pow2(1.- mr1 -mr2) - 4. * mr1 * mr2) );
      widNow = GammaRes * particlePtr->decay[i].bRatio() * betaf / betaN;
    }
    
    // Correction by simple threshold factor for multibody decay.
    else if (meMode == 102 || meMode == 103) {
      mfSum  = 0.;
      for (int j = 0; j < mult; ++j) mfSum 
            += ParticleDataTable::m0( particlePtr->decay[i].product(j) );
      betaf  = sqrtpos(1. - mfSum / mHat);
      betaN  = (meMode == 102) ? 1. : max( BETAMIN, 
               sqrtpos(1. - mfSum / mRes) );
      widNow = GammaRes * particlePtr->decay[i].bRatio() * betaf / betaN;
    } 

    // Default: no correction. 
    else widNow = GammaRes * particlePtr->decay[i].bRatio();

    // Find secondary widths.
    if (openOnly) {
      for (int j = 0; j < particlePtr->decay[i].multiplicity(); ++j) {
        idNow     = particlePtr->decay[i].product(j);
        widNow   *= ParticleDataTable::resOpenFrac(idNow); 
      }
    }

    // Sum back up.
    widSum += widNow;

    // Optionally store partial widths for later decay channel choice.
    if (setBR) particlePtr->decay[i].currentBR(widNow); 
  }

  // Done.
  return widSum;
  
}
 
//**************************************************************************

// The ResonanceCharged class.
// Derived class for a generic resonance with nonidentical antiparticle.

//*********

// Initialize data members.
// Calculate and store partial and total widths at the nominal mass. 

void ResonanceCharged::init() {

  // Do or redo initialization of basic quantities.
  initBasic();

  // Reset quantities to sum. Declare variables inside loop.
  double widTot = 0.; 
  double widPos = 0.;
  double widNeg = 0.;
  int    idNow, idAnti, onMode;
  double widNow, widSecPos, widSecNeg;

  // Loop over all decay channels. Find partial width.  
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow      = GammaRes * particlePtr->decay[i].bRatio();
    widSecPos   = 1.;
    widSecNeg   = 1.;

    // Find secondary widths.
    for (int j = 0; j < particlePtr->decay[i].multiplicity(); ++j) {
      idNow     = particlePtr->decay[i].product(j);
      idAnti    = (ParticleDataTable::hasAnti(idNow)) ? -idNow : idNow;
      widSecPos *= ParticleDataTable::resOpenFrac(idNow); 
      widSecNeg *= ParticleDataTable::resOpenFrac(idAnti); 
    }

    // Store partial widths and update sum, also for open channels only.
    particlePtr->decay[i].onShellWidth(widNow); 
    widTot += widNow;    
    onMode  = particlePtr->decay[i].onMode();
    if (onMode == 1 || onMode == 2) widPos += widNow * widSecPos;
    if (onMode == 1 || onMode == 3) widNeg += widNow * widSecNeg;
  }

  // If no decay channels are open then set particle stable.
  bool mayDecay = (widTot > MINWIDTH); 
  if (!mayDecay) particlePtr->setMayDecay(false, false);

  // Secondary widths for open.
  openPos  = (mayDecay) ? (widPos / widTot) : 1.;
  openNeg  = (mayDecay) ? (widNeg / widTot) : 1.;
    
}

//*********

// Calculate the total width and store phase-space-weighted coupling sums.

double ResonanceCharged::width(int idSgn, double mHat, int , 
  bool openOnly, bool setBR) {

  // Reset quantities to sum. Declare variables inside loop.
  double widSum  = 0.; 
  int    idNow, onMode, meMode, mult;
  double widNow, mfSum, mf1, mf2, mr1, mr2, betaf, betaN;

  // Loop over all decay channels; optionally only open channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow   = 0.;
    onMode   = particlePtr->decay[i].onMode();
    if (openOnly) {
      if (idSgn > 0 && onMode !=1 && onMode != 2) continue;
      if (idSgn < 0 && onMode !=1 && onMode != 3) continue;
    }
  
    // Matrix element mode and channel multiplicity.
    meMode   = particlePtr->decay[i].meMode();
    mult     = particlePtr->decay[i].multiplicity();

    // Correction by step at threshold.
    if (meMode == 101) {
      mfSum  = 0.;
      for (int j = 0; j < mult; ++j) mfSum 
            += ParticleDataTable::m0( particlePtr->decay[i].product(j) );
      if (mfSum < mHat) widNow = GammaRes * particlePtr->decay[i].bRatio(); 
    }  

    // Correction by a phase space factor for two-body decays.
    else if ( (meMode == 102 || meMode == 103) && mult == 2) { 
      mf1    = ParticleDataTable::m0( particlePtr->decay[i].product(0) );
      mf2    = ParticleDataTable::m0( particlePtr->decay[i].product(1) );
      mr1    = pow2(mf1 / mHat);    
      mr2    = pow2(mf2 / mHat);    
      betaf  = sqrtpos( pow2(1.- mr1 -mr2) - 4. * mr1 * mr2);
      mr1    = pow2(mf1 / mRes);   
      mr2    = pow2(mf2 / mRes);   
      betaN  = (meMode == 102) ? 1. : max( BETAMIN, 
               sqrtpos( pow2(1.- mr1 -mr2) - 4. * mr1 * mr2) );
      widNow = GammaRes * particlePtr->decay[i].bRatio() * betaf / betaN;
    }
    
    // Correction by simple threshold factor for multibody decay.
    else if (meMode == 102 || meMode == 103) {
      mfSum  = 0.;
      for (int j = 0; j < mult; ++j) mfSum 
            += ParticleDataTable::m0( particlePtr->decay[i].product(j) );
      betaf  = sqrtpos(1. - mfSum / mHat);
      betaN  = (meMode == 102) ? 1. : max( BETAMIN, 
               sqrtpos(1. - mfSum / mRes) );
      widNow = GammaRes * particlePtr->decay[i].bRatio() * betaf / betaN;
    } 

    // Default: no correction. 
    else widNow = GammaRes * particlePtr->decay[i].bRatio();

    // Find secondary widths.
    if (openOnly) {
      for (int j = 0; j < particlePtr->decay[i].multiplicity(); ++j) {
        idNow   = particlePtr->decay[i].product(j);
        if (ParticleDataTable::hasAnti(idNow)) idNow = -idNow;
        widNow *= ParticleDataTable::resOpenFrac(idNow); 
      }
    }

    // Sum back up.
    widSum += widNow;

    // Optionally store partial widths for later decay channel choice.
    if (setBR) particlePtr->decay[i].currentBR(widNow); 
  }

  // Done.
  return widSum;
  
}
 
//**************************************************************************

// The ResonanceGmZ class.
// Derived class for gamma*/Z0 properties.

//*********

// Initialize data members.
// Calculate and store partial and total widths at the nominal mass. 
// Only consider the Z0 part, as appropriate for Breit-Wigner width.

void ResonanceGmZ::init() {

  // Do or redo initialization of basic quantities.
  initBasic();

  // Locally stored properties and couplings.
  gmZmode       = Settings::mode("SigmaProcess:gmZmode");
  thetaWRat     = 1. / (16. * CoupEW::sin2thetaW() * CoupEW::cos2thetaW());

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
  double widNow, mf, mr, betaf, kinFacVec, kinFacAxi;

  // Loop over all decay channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow      = 0.;
    idAbs       = abs( particlePtr->decay[i].product(0) );

    // Only contributions from three fermion generations, except top.
    if ( (idAbs > 0 && idAbs < 6) || ( idAbs > 10 && idAbs < 17)) {
      mf        = ParticleDataTable::m0(idAbs);

      // Check that above threshold. Kinematical factors.
      if (mHat > 2. * mf + MASSMARGIN) {
        mr        = pow2(mf / mHat);
        betaf     = sqrtpos(1. - 4. * mr); 
        kinFacVec = betaf * (1. + 2. * mr);
        kinFacAxi = pow3(betaf);

        // Combine kinematics with colour factor and couplings.
        widNow    = preFac * (CoupEW::vf2(idAbs) * kinFacVec 
                  + CoupEW::af2(idAbs) * kinFacAxi); 
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
  if (!doForceWidth) particlePtr->setMWidth(widTot, false);
  double bRatio;
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    bRatio = particlePtr->decay[i].onShellWidth() / widTot;
    particlePtr->decay[i].bRatio( bRatio, false);
  }

  // Update the width in this object itself. Secondary widths for open.
  forceFactor   = GammaRes / widTot;
  if (doForceWidth) for (int i = 0; i < particlePtr->decay.size(); ++i)
    particlePtr->decay[i].onShellWidthFactor( forceFactor); 
  else GammaRes = widTot;
  GamMRat       = GammaRes / mRes;  
  openPos       = widPos / widTot;

}

//*********

// Calculate the total width and store phase-space-weighted coupling sums.
// For idIn fermion provides relative weight of specific incoming flavour.
// Use full gamma*/Z0, as appropriate for final-state flavour choice.

double ResonanceGmZ::width(int , double mHat, int idIn, 
  bool openOnly, bool setBR) {

  // Common coupling factors.
  double sH      = mHat*mHat;
  double alpEM   = alphaEM.alphaEM(sH); 
  double alpS    = alphaS.alphaS(sH);
  double colQ    = 3. * (1. + alpS / M_PI);
  double preFac  = alpEM * thetaWRat * mHat / 3.;
  if (doForceWidth) preFac *= forceFactor;

  // Couplings when an incoming fermion is specified; elso only pure Z0.
  double ei2     = 0.;
  double eivi    = 0.;  
  double vi2ai2  = 1.;
  int idInAbs    = abs(idIn);
  if (idInAbs > 0 && idInAbs < 19) {
    ei2          = CoupEW::ef2(idInAbs);
    eivi         = CoupEW::efvf(idInAbs);
    vi2ai2       = CoupEW::vf2af2(idInAbs); 
  }

  // Calculate prefactors for gamma/interference/Z0 terms.
  double gamNorm = ei2;
  double intNorm = 2. * eivi * thetaWRat * sH * (sH - m2Res)
    / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  double resNorm = vi2ai2 * pow2(thetaWRat * sH) 
    / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );

  // Rescale Z0 height normalization to compensate for a width one.
  if (doForceWidth) {
    intNorm     *= forceFactor;
    resNorm     *= forceFactor;
  }

  // Optionally only keep gamma* or Z0 term.
  if (gmZmode == 1) {intNorm = 0.; resNorm = 0.;}
  if (gmZmode == 2) {gamNorm = 0.; intNorm = 0.;}

  // Reset quantities to sum. Declare variables inside loop.
  double widSum  = 0.; 
  int    idAbs, onMode;
  double widNow, mf, mr, betaf, kinFacVec, kinFacAxi, vf2af2;

  // Loop over all decay channels; optionally only open channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow       = 0.;
    onMode       = particlePtr->decay[i].onMode();
    if (openOnly && onMode !=1 && onMode != 2) continue;
    idAbs        = abs( particlePtr->decay[i].product(0) );

    // Only contributions from three fermion generations, except top.
    if ( (idAbs > 0 && idAbs < 6) || ( idAbs > 10 && idAbs < 17)) {
      mf = ParticleDataTable::m0(idAbs);

      // Check that above threshold. Kinematical factors.
      if (mHat > 2. * mf + MASSMARGIN) {
        mr        = pow2(mf / mHat);
        betaf     = sqrtpos(1. - 4. * mr); 
        kinFacVec = betaf * (1. + 2. * mr);
        kinFacAxi = pow3(betaf);
        vf2af2    = CoupEW::vf2(idAbs) * kinFacVec 
                  + CoupEW::af2(idAbs) * kinFacAxi; 

        // Z0 width: combine kinematics with couplings.
        if (idIn == 0) widNow = preFac * vf2af2;

        // Relative outwidths: combine instate, propagator and outstate.
        else widNow = gamNorm * CoupEW::ef2(idAbs) * kinFacVec 
                    + intNorm * CoupEW::efvf(idAbs) * kinFacVec 
                    + resNorm * vf2af2;

        // Colour factor.
        if (idAbs < 6) widNow *= colQ;
      }
    }

    // Sum back up.
    widSum += widNow;

    // Optionally store partial widths for later decay channel choice.
    if (setBR) particlePtr->decay[i].currentBR(widNow); 
  }

  // Done.
  return widSum;
  
}
 
//**************************************************************************

// The ResonanceW class.
// Derived class for W+- properties.

//*********

// Initialize static data members.
// Calculate and store partial and total widths at the nominal mass. 

void ResonanceW::init() {

  // Do or redo initialization of basic quantities.
  initBasic();

  // Locally stored properties and couplings.
  thetaWRat     = 1. / (12. * CoupEW::sin2thetaW());

  // Common coupling factors.
  double mHat   = mRes;
  double sH     = mHat*mHat;
  double alpEM  = alphaEM.alphaEM(sH);
  double alpS   = alphaS.alphaS(sH);
  double colQ   = 3. * (1. + alpS / M_PI);
  double preFac = alpEM * thetaWRat * mHat;

  // Reset quantities to sum. Declare variables inside loop.
  double widTot = 0.; 
  double widPos = 0.;
  double widNeg = 0.;
  int    id1Abs, id2Abs, onMode;
  double widNow, mf1, mf2, mr1, mr2, kinFac;

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

      // Check that above threshold. Kinematical factor.
      if (mHat > mf1 + mf2 + MASSMARGIN) {
        mr1    = pow2(mf1 / mHat);
        mr2    = pow2(mf2 / mHat);
        kinFac = (1. - 0.5 * (mr1 + mr2) - 0.5 * pow2(mr1 - mr2))
               * sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 ); 

        // Combine kinematics with colour factor and CKM couplings.
        widNow = preFac * kinFac;
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
  if (!doForceWidth) particlePtr->setMWidth(widTot, false);
  double bRatio;
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    bRatio = particlePtr->decay[i].onShellWidth() / widTot;
    particlePtr->decay[i].bRatio( bRatio, false);
  }

  // Update the width in this object itself. Secondary widths for open.
  forceFactor   = GammaRes / widTot;
  if (doForceWidth) for (int i = 0; i < particlePtr->decay.size(); ++i)
    particlePtr->decay[i].onShellWidthFactor( forceFactor); 
  else GammaRes = widTot;
  GamMRat       = GammaRes / mRes;  
  openPos       = widPos / widTot;
  openNeg       = widNeg / widTot;

}

//*********

// Calculate the total width and store phases-space-weighted coupling sums.

double ResonanceW::width(int idSgn, double mHat, int , 
  bool openOnly, bool setBR) {

  // Common coupling factors.
  double sH     = mHat*mHat;
  double alpEM  = alphaEM.alphaEM(sH); 
  double alpS   = alphaS.alphaS(sH);
  double colQ   = 3. * (1. + alpS / M_PI);
  double preFac = alpEM * thetaWRat * mHat;
  if (doForceWidth) preFac *= forceFactor;

  // Reset quantities to sum. Declare variables inside loop.
  double widSum = 0.; 
  int    id1Abs, id2Abs, onMode;
  double widNow, mf1, mf2, mr1, mr2, ps;

  // Loop over all decay channels; optionally only open channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow = 0.;
    onMode = particlePtr->decay[i].onMode();
    if (openOnly) {
      if (idSgn > 0 && onMode !=1 && onMode != 2) continue;
      if (idSgn < 0 && onMode !=1 && onMode != 3) continue;
    }
    id1Abs = abs( particlePtr->decay[i].product(0) );
    id2Abs = abs( particlePtr->decay[i].product(1) );

    // Only contributions from three fermion generations, except top.
    if ( (id1Abs > 0 && id1Abs < 6 && id2Abs > 0 && id2Abs < 6) 
      || (id1Abs > 10 && id1Abs < 17 && id2Abs > 10 && id2Abs < 17)) {
      mf1 = ParticleDataTable::m0(id1Abs);
      mf2 = ParticleDataTable::m0(id2Abs);

      // Check that above threshold. Kinematical factor.
      if (mHat > mf1 + mf2 + MASSMARGIN) {
        mr1 = pow2(mf1 / mHat);
        mr2 = pow2(mf2 / mHat);
        ps = (1. - 0.5 * (mr1 + mr2) - 0.5 * pow2(mr1 - mr2))
          * sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 ); 

        // Combine kinematics with colour factor and CKM couplings.
        widNow = preFac * ps;
        if (id1Abs < 6) widNow *= colQ * VCKM::V2id(id1Abs, id2Abs);
      }
    }

    // Sum back up.
    widSum += widNow;

    // Optionally store partial widths for later decay channel choice.
    if (setBR) particlePtr->decay[i].currentBR(widNow); 
  }

  // Done.
  return widSum;
  
}
 
//**************************************************************************

// The ResonanceTop class.
// Derived class for top/antitop properties.

//*********

// Initialize static data members.
// Calculate and store partial and total widths at the nominal mass. 

void ResonanceTop::init() {

  // Do or redo initialization of basic quantities.
  initBasic();

  // Locally stored properties and couplings.
  thetaWRat     = 1. / (16. * CoupEW::sin2thetaW());
  m2W           = pow2(ParticleDataTable::m0(24));

  // Common coupling factors.
  double mHat   = mRes;
  double sH     = mHat*mHat;
  double alpEM  = alphaEM.alphaEM(sH);
  double alpS   = alphaS.alphaS(sH);
  double colQ   = 1. - 2.5 * alpS / M_PI;
  double preFac = alpEM * thetaWRat * (sH / m2W) * mHat;

  // Reset quantities to sum. Declare variables inside loop.
  double widTot = 0.; 
  double widPos = 0.;
  double widNeg = 0.;
  int    id1Abs, id2Abs, onMode;
  double widNow, widSecPos, widSecNeg, mf1, mf2, mr1, mr2, kinFac;

  // Loop over all decay channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow      = 0.;
    widSecPos   = 1.;
    widSecNeg   = 1.;
    id1Abs      = abs( particlePtr->decay[i].product(0) );
    id2Abs      = abs( particlePtr->decay[i].product(1) );
    if (id2Abs > id1Abs) swap( id1Abs, id2Abs);

    // Only contributions from W + d/s/b.
    if ( id1Abs == 24 && (id2Abs == 1 || id2Abs == 3 || id2Abs ==5) ) { 
      mf1 = ParticleDataTable::m0(id1Abs);
      mf2 = ParticleDataTable::m0(id2Abs);

      // Check that above threshold. Kinematical factor.
      if (mHat > mf1 + mf2 + MASSMARGIN) {
        mr1    = pow2(mf1 / mHat);
        mr2    = pow2(mf2 / mHat);
        kinFac = ( pow2(1. - mr2) + (1. + mr2) * mr1 - 2. * mr1 * mr1 )
           * sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 ); 

        // Combine kinematics with colour factor and CKM couplings.
        widNow = preFac * kinFac * colQ * VCKM::V2id(6, id2Abs);
 
        // Secondary width from W+ and W- decay.
        widSecPos = ParticleDataTable::resOpenFrac(24); 
        widSecNeg = ParticleDataTable::resOpenFrac(-24); 
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
  if (!doForceWidth) particlePtr->setMWidth(widTot, false);
  double bRatio;
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    bRatio = particlePtr->decay[i].onShellWidth() / widTot;
    particlePtr->decay[i].bRatio( bRatio, false);
  }

  // Update the width in this object itself. Secondary widths for open.
  forceFactor   = GammaRes / widTot;
  if (doForceWidth) for (int i = 0; i < particlePtr->decay.size(); ++i)
    particlePtr->decay[i].onShellWidthFactor( forceFactor); 
  else GammaRes = widTot;
  GamMRat       = GammaRes / mRes;  
  openPos       = widPos / widTot;
  openNeg       = (widNeg / widTot);

}

//*********

// Calculate the total width and store phases-space-weighted coupling sums.

double ResonanceTop::width(int idSgn, double mHat, int , 
  bool openOnly, bool setBR) {

  // Common coupling factors.
  double sH     = mHat*mHat;
  double alpEM  = alphaEM.alphaEM(sH); 
  double alpS   = alphaS.alphaS(sH);
  double colQ   = 3. * (1. + alpS / M_PI);
  double preFac = alpEM * thetaWRat * (sH / m2W) * mHat;
  if (doForceWidth) preFac *= forceFactor;
 
  // Reset quantities to sum. Declare variables inside loop.
  double widSum = 0.; 
  int    id1Abs, id2Abs, onMode;
  double widNow, mf1, mf2, mr1, mr2, kinFac;

  // Loop over all decay channels; optionally only open channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow = 0.;
    onMode = particlePtr->decay[i].onMode();
    if (openOnly) {
      if (idSgn > 0 && onMode !=1 && onMode != 2) continue;
      if (idSgn < 0 && onMode !=1 && onMode != 3) continue;
    }
    id1Abs = abs( particlePtr->decay[i].product(0) );
    id2Abs = abs( particlePtr->decay[i].product(1) );
    if (id2Abs > id1Abs) swap( id1Abs, id2Abs);

    // Only contributions from W + d/s/b.
    if ( id1Abs == 24 && (id2Abs == 1 || id2Abs == 3 || id2Abs ==5) ) { 
      mf1 = ParticleDataTable::m0(id1Abs);
      mf2 = ParticleDataTable::m0(id2Abs);

      // Check that above threshold. Kinematical factor.
      if (mHat > mf1 + mf2 + MASSMARGIN) {
        mr1    = pow2(mf1 / mHat);
        mr2    = pow2(mf2 / mHat);
        kinFac = ( pow2(1. - mr2) + (1. + mr2) * mr1 - 2. * mr1 * mr1 )
           * sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 ); 

        // Combine kinematics with colour factor and CKM couplings.
        widNow = preFac * kinFac * colQ * VCKM::V2id(6, id2Abs);

        // Secondary width from W decay.
        int idSgnW = (idSgn > 0) ? 24 : -24;
        if (openOnly) widNow *= ParticleDataTable::resOpenFrac(idSgnW); 
      }
    }

    // Sum back up.
    widSum += widNow;

    // Optionally store partial widths for later decay channel choice.
    if (setBR) particlePtr->decay[i].currentBR(widNow); 
  }

  // Done.
  return widSum;
  
}
 
//**************************************************************************

// The ResonanceFour class.
// Derived class for fourth-generation-fermion properties.

//*********

// Initialize static data members.
// Calculate and store partial and total widths at the nominal mass. 

void ResonanceFour::init() {

  // Do or redo initialization of basic quantities.
  initBasic();

  // Locally stored properties and couplings.
  thetaWRat     = 1. / (16. * CoupEW::sin2thetaW());
  m2W           = pow2(ParticleDataTable::m0(24));

  // Common coupling factors.
  double mHat   = mRes;
  double sH     = mHat*mHat;
  double alpEM  = alphaEM.alphaEM(sH);
  double alpS   = alphaS.alphaS(sH);
  double colQ   = (idRes < 9) ? 1. - 2.5 * alpS / M_PI : 1.;
  double preFac = alpEM * thetaWRat * (sH / m2W) * mHat;

  // Reset quantities to sum. Declare variables inside loop.
  double widTot = 0.; 
  double widPos = 0.;
  double widNeg = 0.;
  int    id1Abs, id2Abs, onMode;
  double widNow, widSecPos, widSecNeg, mf1, mf2, mr1, mr2, kinFac;

  // Loop over all decay channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow      = 0.;
    widSecPos   = 1.;
    widSecNeg   = 1.;
    id1Abs      = abs( particlePtr->decay[i].product(0) );
    id2Abs      = abs( particlePtr->decay[i].product(1) );
    if (id2Abs > id1Abs) swap( id1Abs, id2Abs);

    // Contributions from W + fermion.
    if ( id1Abs == 24 && id2Abs < 19) { 
      mf1 = ParticleDataTable::m0(id1Abs);
      mf2 = ParticleDataTable::m0(id2Abs);

      // Check that above threshold. Kinematical factor.
      if (mHat > mf1 + mf2 + MASSMARGIN) {
        mr1    = pow2(mf1 / mHat);
        mr2    = pow2(mf2 / mHat);
        kinFac = ( pow2(1. - mr2) + (1. + mr2) * mr1 - 2. * mr1 * mr1 )
           * sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 ); 

        // Combine kinematics with colour factor and CKM couplings.
        widNow = preFac * kinFac * colQ;
        if (idRes < 9) widNow *= VCKM::V2id(idRes, id2Abs);
 
        // Secondary width from W and fermion decay.
        int idSgnW = (idRes%2 == 0) ? 24 : -24;
        widSecPos = ParticleDataTable::resOpenFrac( idSgnW,  id2Abs); 
        widSecNeg = ParticleDataTable::resOpenFrac(-idSgnW, -id2Abs); 
      }
    }

    // Store partial widths and update sum
    particlePtr->decay[i].onShellWidth(widNow); 
    widTot += widNow;    
    onMode = particlePtr->decay[i].onMode();
    if (onMode == 1 || onMode == 2) widPos += widNow * widSecPos;
    if (onMode == 1 || onMode == 3) widNeg += widNow * widSecNeg;
  }

  // If no decay channels are open then set particle stable.
  bool mayDecay = (widTot > MINWIDTH); 
  if (!mayDecay) particlePtr->setMayDecay(false, false);

  // Update total width and branching ratios.
  if (!doForceWidth) particlePtr->setMWidth(widTot, false);
  double bRatio;
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    bRatio = (mayDecay) ? (particlePtr->decay[i].onShellWidth() / widTot) 
      : 0.;
    particlePtr->decay[i].bRatio( bRatio, false);
  }

  // Update the width in this object itself. Secondary widths for open.
  forceFactor   = GammaRes / widTot;
  if (doForceWidth) for (int i = 0; i < particlePtr->decay.size(); ++i)
    particlePtr->decay[i].onShellWidthFactor( forceFactor); 
  else GammaRes = widTot;
  GamMRat       = GammaRes / mRes;  
  openPos       = (mayDecay) ? (widPos / widTot) : 1.;
  openNeg       = (mayDecay) ? (widNeg / widTot) : 1.;

}

//*********

// Calculate the total width and store phases-space-weighted coupling sums.

double ResonanceFour::width(int idSgn, double mHat, int , 
  bool openOnly, bool setBR) {

  // Common coupling factors.
  double sH     = mHat*mHat;
  double alpEM  = alphaEM.alphaEM(sH); 
  double alpS   = alphaS.alphaS(sH);
  double colQ   = (idRes < 9) ? 1. - 2.5 * alpS / M_PI : 1.;
  double preFac = alpEM * thetaWRat * (sH / m2W) * mHat;
  if (doForceWidth) preFac *= forceFactor;
 
  // Reset quantities to sum. Declare variables inside loop.
  double widSum = 0.; 
  int    id1Abs, id2Abs, onMode;
  double widNow, mf1, mf2, mr1, mr2, kinFac;

  // Loop over all decay channels; optionally only open channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow = 0.;
    onMode = particlePtr->decay[i].onMode();
    if (openOnly) {
      if (idSgn > 0 && onMode !=1 && onMode != 2) continue;
      if (idSgn < 0 && onMode !=1 && onMode != 3) continue;
    }
    id1Abs = abs( particlePtr->decay[i].product(0) );
    id2Abs = abs( particlePtr->decay[i].product(1) );
    if (id2Abs > id1Abs) swap( id1Abs, id2Abs);

    // Contributions from W + fermion.
    if ( id1Abs == 24 && id2Abs < 19 ) { 
      mf1 = ParticleDataTable::m0(id1Abs);
      mf2 = ParticleDataTable::m0(id2Abs);

      // Check that above threshold. Kinematical factor.
      if (mHat > mf1 + mf2 + MASSMARGIN) {
        mr1    = pow2(mf1 / mHat);
        mr2    = pow2(mf2 / mHat);
        kinFac = ( pow2(1. - mr2) + (1. + mr2) * mr1 - 2. * mr1 * mr1 )
           * sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 ); 

        // Combine kinematics with colour factor and CKM couplings.
        widNow = preFac * kinFac * colQ;
        if (idRes < 9) widNow *= VCKM::V2id(idRes, id2Abs);

        // Secondary width from W and fermion decay.
        int idSgnW = ((idSgn > 0 && idRes%2 == 0) 
          || (idSgn < 0 && idRes%2 == 1)) ? 24 : -24;
        int idSgnF = (idSgn > 0) ? id2Abs : -id2Abs;
        if (openOnly) widNow 
          *= ParticleDataTable::resOpenFrac(idSgnW, idSgnF); 
      }
    }

    // Sum back up.
    widSum += widNow;

    // Optionally store partial widths for later decay channel choice.
    if (setBR) particlePtr->decay[i].currentBR(widNow); 
  }

  // Done.
  return widSum;
  
}
 
//**************************************************************************

// The ResonanceH class.
// Derived class for SM and BSM Higgs properties.

//*********

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Minimal mass for W, Z, top in integration over respective Breit-Wigner.
const double ResonanceH::MASSMIN = 10.;

// Number of widths above threshold where B-W integration not needed.
const double ResonanceH::GAMMAMARGIN = 10.;

//*********

// Initialize data members.
// Calculate and store partial and total widths at the nominal mass. 

void ResonanceH::init() {

  // Do or redo initialization of basic quantities.
  initBasic();

  // Locally stored properties and couplings.
  useCubicWidth = Settings::flag("Higgs:cubicWidth");
  sin2tW        = CoupEW::sin2thetaW();
  cos2tW        = 1. - sin2tW;
  mZ            = ParticleDataTable::m0(23);
  mW            = ParticleDataTable::m0(24);
  mHchg         = ParticleDataTable::m0(37);
  GammaZ        = ParticleDataTable::mWidth(23);
  GammaW        = ParticleDataTable::mWidth(24);
  GammaT        = ParticleDataTable::mWidth(6);

  // Couplings to fermions, Z and W, depending on Higgs type.
  coup2d        = 1.;
  coup2u        = 1.;
  coup2l        = 1.;
  coup2Z        = 1.;
  coup2W        = 1.;
  coup2Hchg     = 0.;
  coup2H1H1     = 0.;
  coup2A3A3     = 0.;
  coup2H1Z      = 0.;
  coup2A3Z      = 0.;
  coup2A3H1     = 0.;
  coup2HchgW    = 0.;
  if (higgsType == 1) {
    coup2d      = Settings::parm("HiggsH1:coup2d");
    coup2u      = Settings::parm("HiggsH1:coup2u");
    coup2l      = Settings::parm("HiggsH1:coup2l");
    coup2Z      = Settings::parm("HiggsH1:coup2Z");
    coup2W      = Settings::parm("HiggsH1:coup2W");
    coup2Hchg   = Settings::parm("HiggsH1:coup2Hchg");
  } else if (higgsType == 2) {
    coup2d      = Settings::parm("HiggsH2:coup2d");
    coup2u      = Settings::parm("HiggsH2:coup2u");
    coup2l      = Settings::parm("HiggsH2:coup2l");
    coup2Z      = Settings::parm("HiggsH2:coup2Z");
    coup2W      = Settings::parm("HiggsH2:coup2W");
    coup2Hchg   = Settings::parm("HiggsH2:coup2Hchg");
    coup2H1H1   = Settings::parm("HiggsH2:coup2H1H1");
    coup2A3A3   = Settings::parm("HiggsH2:coup2A3A3");
    coup2H1Z    = Settings::parm("HiggsH2:coup2H1Z");
    coup2A3Z    = Settings::parm("HiggsA3:coup2H2Z");
    coup2A3H1   = Settings::parm("HiggsH2:coup2A3H1");
    coup2HchgW  = Settings::parm("HiggsH2:coup2HchgW");
  } else if (higgsType == 3) {
    coup2d      = Settings::parm("HiggsA3:coup2d");
    coup2u      = Settings::parm("HiggsA3:coup2u");
    coup2l      = Settings::parm("HiggsA3:coup2l");
    coup2Z      = Settings::parm("HiggsA3:coup2Z");
    coup2W      = Settings::parm("HiggsA3:coup2W");
    coup2Hchg   = Settings::parm("HiggsA3:coup2Hchg");
    coup2H1H1   = Settings::parm("HiggsA3:coup2H1H1");
    coup2H1Z    = Settings::parm("HiggsA3:coup2H1Z");
    coup2HchgW  = Settings::parm("HiggsA3:coup2Hchg");
  }

  // Common coupling factors.
  double mHat   = mRes;
  double sH     = mHat*mHat;
  double alpEM  = alphaEM.alphaEM(sH);
  double alpS   = alphaS.alphaS(sH);
  double colQ   = 3. * (1. + alpS / M_PI);
  double preFac = (alpEM / (8. * sin2tW)) * pow3(mHat) / pow2(mW); 

  // Reset quantities to sum. Declare variables inside loop.
  double widTot = 0.; 
  double widPos = 0.;
  int    id1, id2, id1Abs, id2Abs, meMode, onMode;
  double widNow, widSec, kinFac, coupFac, mf1, mf2, mr1, mr2, ps;

  // Loop over all decay channels.
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow      = 0.;
    widSec      = 1.;
    kinFac      = 0.;
    id1         = particlePtr->decay[i].product(0);
    id2         = particlePtr->decay[i].product(1);
    id1Abs      = abs(id1);
    id2Abs      = abs(id2);
    mf1         = ParticleDataTable::m0(id1Abs);
    mf2         = ParticleDataTable::m0(id2Abs);
    mr1         = pow2(mf1 / mHat);
    mr2         = pow2(mf2 / mHat);
    ps          = (mf1 + mf2 + MASSMARGIN > mHat) ? 0.   
                : sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 );

    // Widths of decays Higgs -> f + fbar.
    if ( id2Abs == id1Abs && ( (id1Abs > 0 && id1Abs < 7) 
      || (id1Abs > 10 && id1Abs < 17) ) ) {

      // Check that above threshold. Kinematical factor.
      if ( (id1Abs != 6 && mHat > 2. * mf1 + MASSMARGIN)
        || (id1Abs == 6 && mHat > 2. * (mf1 + GAMMAMARGIN * GammaT) ) ) {
        // A0 behaves like beta, h0 and H0 like beta**3.
        kinFac = (higgsType < 3) ? pow3(ps) : ps;
      }

      // For top near threshold use numerical integration.   
      else if (id1Abs == 6 && mHat > 2. * (mf1 - GAMMAMARGIN * GammaT)) {
        meMode = (higgsType < 3) ? 3 : 1;       
        kinFac = numInt2BW( mHat, mf1, GammaT, MASSMIN, mf1, GammaT, 
        MASSMIN, meMode);
      }

      // Coupling from mass and from BSM deviation from SM.
      coupFac = pow2(ParticleDataTable::mRun(id1Abs, mHat) / mHat);
      if (id1Abs < 7 && id1Abs%2 == 1) coupFac *= coup2d * coup2d;
      else if (id1Abs < 7)             coupFac *= coup2u * coup2u;   
      else                             coupFac *= coup2l * coup2l;

      // Combine couplings and phase space with colour factor.
      widNow                  = preFac * coupFac * kinFac; 
      if (id1Abs < 7) widNow *= colQ;

      // Secondary width from top decay.
      if (id1Abs == 6) widSec = ParticleDataTable::resOpenFrac(6, -6); 

      // Store widths for in-state: no phase space and no colour factors.
      widTable[id1Abs]        = preFac * coupFac;
    }

    // Widths of decays Higgs -> g + g. No colour factor for instate.
    else if (id1Abs == 21 && id2Abs == 21) { 
      widNow       = preFac * pow2(alpS / M_PI) * eta2gg(mHat); 
      widTable[21] = widNow / 8.;
    }

    // Widths of decays Higgs -> gamma + gamma.
    else if (id1Abs == 22 && id2Abs == 22) { 
      widNow       = preFac * pow2(alpEM / M_PI) * 0.5 * eta2gaga(mHat); 
      widTable[22] = widNow;
    }
 
    // Widths of decays Higgs -> Z0 + gamma0. Secondary width.
    else if (max(id1Abs, id2Abs) == 23 && min(id1Abs, id2Abs) == 22) {
      widNow = preFac * pow2(alpEM / M_PI) * pow3(ps) * eta2gaZ(mHat); 
      widSec = ParticleDataTable::resOpenFrac(23); 
    }
 
    // Widths of decays Higgs (h0, H0) -> Z0 + Z0. Secondary width.
    else if (id1Abs == 23 && id2Abs == 23) {
      // If Higgs heavy use on-shell expression, else numerical integration.
      widNow       = ( mHat > 2. * (mZ + GAMMAMARGIN * GammaZ) ) 
        ? (1. - 4. * mr1 + 12. * mr1 * mr1) * ps
        : numInt2BW( mHat, mZ, GammaZ, MASSMIN, mZ, GammaZ, MASSMIN, 5);
      widNow      *= 0.25 * preFac * pow2(coup2Z);
      widTable[23] = 0.25 * preFac * pow2(coup2Z);
      widSec       = ParticleDataTable::resOpenFrac(23, 23); 
    }
 
    // Widths of decays Higgs (h0, H0) -> W+ + W-. Secondary width.
    else if (id1Abs == 24 && id2Abs == 24) {
      // If Higgs heavy use on-shell expression, else numerical integration.
      widNow       = (mHat > 2. * ( mW + GAMMAMARGIN * GammaW)) 
        ? (1. - 4. * mr1 + 12. * mr1 * mr1) * ps
        : numInt2BW( mHat, mW, GammaW, MASSMIN, mW, GammaW, MASSMIN, 5);
      widNow      *= 0.5 * preFac * pow2(coup2W);
      widTable[24] = 0.5 * preFac * pow2(coup2W);
      widSec       = ParticleDataTable::resOpenFrac(24, -24); 
    }
 
    // Widths of decays Higgs (H0) -> h0 + h0. Secondary width.
    else if (id1Abs == 25 && id2Abs == 25) {
      widNow = 0.25 * preFac * pow4(mZ / mHat) * ps * pow2(coup2H1H1);
      widSec = ParticleDataTable::resOpenFrac(25, 25); 
    }
     
    // Widths of decays Higgs (H0) -> A0 + A0. Secondary width.
    else if (id1Abs == 36 && id2Abs == 36) {
      widNow = 0.5 * preFac * pow4(mZ / mHat) * ps * pow2(coup2A3A3);
      widSec = ParticleDataTable::resOpenFrac(36, 36); 
    }
 
    // Widths of decays Higgs (A0) -> h0 + Z0. Secondary width.
    else if (max(id1Abs, id2Abs) == 25 && min(id1Abs, id2Abs) == 23) {
      widNow = 0.5 * preFac * pow3(ps) * pow2(coup2H1Z);
      widSec = ParticleDataTable::resOpenFrac(25, 23); 
    }
 
    // Widths of decays Higgs (H0) -> A0 + Z0. Secondary width.
    else if (max(id1Abs, id2Abs) == 36 && min(id1Abs, id2Abs) == 23) {
      widNow = 0.5 * preFac * pow3(ps) * pow2(coup2A3Z);
      widSec = ParticleDataTable::resOpenFrac(36, 23); 
    }
     
    // Widths of decays Higgs (H0) -> A0 + h0. Secondary width.
    else if (max(id1Abs, id2Abs) == 36 && min(id1Abs, id2Abs) == 25) {
      widNow = 0.25 * preFac * pow4(mZ / mHat) * ps * pow2(coup2A3H1);
      widSec = ParticleDataTable::resOpenFrac(36, 25); 
    }
 
    // Widths of decays Higgs -> H+- + W-+. Secondary width.
    else if (max(id1Abs, id2Abs) == 37 && min(id1Abs, id2Abs) == 24) {
      widNow = 0.5 * preFac * pow3(ps) * pow2(coup2HchgW);
      widSec = (id1 == 37 || id2 == 37)
             ? ParticleDataTable::resOpenFrac( 37, -24)  
             : ParticleDataTable::resOpenFrac(-37,  24); 
    }
    
    // Store partial widths and update sum, also for open channels only.
    particlePtr->decay[i].onShellWidth(widNow); 
    widTot += widNow;    
    onMode = particlePtr->decay[i].onMode();
    if (onMode == 1 || onMode == 2) widPos += widNow * widSec;
  }

  // Update total width and branching ratios.
  if (!doForceWidth) particlePtr->setMWidth(widTot, false);
  double bRatio;
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    bRatio = particlePtr->decay[i].onShellWidth() / widTot;
    particlePtr->decay[i].bRatio( bRatio, false);
  }

  // Update the width in this object itself. Secondary widths for open.
  forceFactor   = GammaRes / widTot;
  if (doForceWidth) for (int i = 0; i < particlePtr->decay.size(); ++i)
    particlePtr->decay[i].onShellWidthFactor( forceFactor); 
  else GammaRes = widTot;
  GamMRat       = GammaRes / mRes;  
  openPos       = widPos / widTot;

}

//*********
 
// Calculate the total or open width at a given energy.
// Has approximate rescaling since the correct width calculation 
// (in widthInit()) is slow owing to numerical integrations.

double ResonanceH::width(int , double mHat, int , 
  bool openOnly, bool setBR) { 

  // Reset quantities to sum. Declare variables inside loop.
  double widSum = 0.; 
  int    id1, id2, id1Abs, id2Abs, onMode;
  double widNow;

  // Loop over all decay channels; optionally only open channels.
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow = 0.;
    onMode = particlePtr->decay[i].onMode();
    if (openOnly && onMode !=1 && onMode != 2) continue;
    id1    = particlePtr->decay[i].product(0);
    id2    = particlePtr->decay[i].product(1);
    id1Abs = abs(id1);
    id2Abs = abs(id2);
    if (id2Abs > id1Abs) swap( id1Abs, id2Abs);
    
    // Start out from on-shell partial width.     
    widNow = particlePtr->decay[i].onShellWidth();

    // Simple rescaling with first power of actual mass.
    widNow *= mHat / mRes;

    // Optionally two more powers for decays to gauge bosons and Higgses.
    if (useCubicWidth && id1Abs > 22 && id2Abs > 22) 
      widNow *= pow2(mHat / mRes);   

    // Approximate empirical threshold for H -> Z0 Z0 or W+ W-.
    if ( (id1Abs == 23 || id1Abs == 24) && id2Abs == id1Abs) {
      double mZW     = (id1Abs == 23) ? mZ : mW;  
      double GammaZW = (id1Abs == 23) ? GammaZ : GammaW;  
      double mRatio5 = pow5(0.5 * mHat / mZW);
      double thrNow  = (atan((mHat - 2. * mZW) / GammaZW) / M_PI + 0.5)
                     * mRatio5 / (1. + mRatio5);
      mRatio5        = pow5(0.5 * mRes / mZW);
      double thrRes  = (atan((mRes - 2. * mZW) / GammaZW) / M_PI + 0.5)
                     * mRatio5 / (1. + mRatio5);
      widNow        *= thrNow / thrRes;  
    }  

    // Approximate phase-space threshold for decays to other Higgses.
    if (id1Abs >= 25 && widNow > 0.) {
      double mf1     = ParticleDataTable::m0(id1Abs);
      double mf2     = ParticleDataTable::m0(id2Abs);
      double mr1     = pow2(mf1 / mHat);
      double mr2     = pow2(mf2 / mHat);
      double ps      = (mf1 + mf2 + MASSMARGIN > mHat) ? 0.   
                     : sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 );
      mr1            = pow2(mf1 / mRes);
      mr2            = pow2(mf2 / mRes);
      double psRes   = sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 );
      widNow *= (id2Abs == 23 || id2Abs == 24) 
              ? pow3(ps / psRes) : ps / psRes;        
    }

    // Secondary widths.
    if (openOnly) widNow *= ParticleDataTable::resOpenFrac(id1, id2);

    // Sum back up.
    widSum += widNow;

    // Optionally store partial widths for later decay channel choice.
    if (setBR) particlePtr->decay[i].currentBR(widNow); 
  }

  // Done.
  return widSum;

} 

//*********
 
// Calculate the partial width for a specific mass and in/out-state.
// This width has no phase-space suppression and is given without
// colour factors. Uses pretabulation from widthInit() to speed up.

double ResonanceH::widthChan(double mHat, int id1Abs, int ) { 

  // Read out from stored values.
  double widNow = (id1Abs > 0 && id1Abs < 25) ? widTable[id1Abs] : 0.;

  // Simple rescaling with first power of actual mass.
  widNow *= mHat / mRes;

  // Optionally two more powers for decays to gauge bosons and Higgses.
  if (useCubicWidth && id1Abs > 22) widNow *= pow2(mHat / mRes);   

  // Done.
  return widNow;

}

//*********

// Sum up quark loop contributions in Higgs -> g + g.
// Note: running quark masses are used, unlike Pythia6 (not negligible shift). 

double ResonanceH::eta2gg(double mHat) {

  // Initial values.
  complex eta = complex(0., 0.);
  double  epsilon, root, rootLog;
  complex phi, etaNow;

  // Loop over s, c, b, t quark flavours.
  for (int id = 3; id < 7; ++id) {
    epsilon = pow2(2. * ParticleDataTable::mRun(id, mHat) / mHat);

    // Value of loop integral.
    if (epsilon <= 1.) {
      root    = sqrt(1. - epsilon);
      rootLog = (epsilon < 1e-4) ? log(4. / epsilon - 2.)
                : log( (1. + root) / (1. - root) );
      phi = complex( -0.25 * (pow2(rootLog) - pow2(M_PI)), 
                     0.5 * M_PI * rootLog );
    } 
    else phi = complex( pow2( asin(1. / sqrt(epsilon)) ), 0.);
  
    // Factors that depend on Higgs and flavour type.
    if (higgsType < 3) etaNow = -0.5 * epsilon 
      * (complex(1., 0.) + (1. - epsilon) * phi);
    else etaNow = -0.5 * epsilon * phi;    
    if (id%2 == 1) etaNow *= coup2d;
    else           etaNow *= coup2u;   
    
    // Sum up contribution and return square of absolute value.
    eta += etaNow;
  }
  return (pow2(eta.real()) + pow2(eta.imag()));

}

//*********

// Sum up quark, lepton, W+- and (for BSM) H+- loop contributions 
// in Higgs -> gamma + gamma.

double ResonanceH::eta2gaga(double mHat) {

  // Initial values.
  complex eta = complex(0., 0.);
  int     id;
  double  ef, epsilon, root, rootLog;
  complex phi, etaNow;

  // Loop over s, c, b, t, mu, tau, W+-, H+- flavours.
  for (int idLoop = 0; idLoop < 8; ++idLoop) {
    if      (idLoop < 4) id = idLoop + 3;
    else if (idLoop < 6) id = 2 * idLoop + 5;
    else if (idLoop < 7) id = 24;
    else                 id = 37;
    if (id == 37 && higgsType == 0) continue;
 
    // Charge and loop integral parameter.
    ef      = (id < 20) ? CoupEW::ef(id) : 1.;
    epsilon = pow2(2. * ParticleDataTable::mRun(id, mHat) / mHat);

    // Value of loop integral.
    if (epsilon <= 1.) {
      root    = sqrt(1. - epsilon);
      rootLog = (epsilon < 1e-4) ? log(4. / epsilon - 2.)
                : log( (1. + root) / (1. - root) );
      phi = complex( -0.25 * (pow2(rootLog) - pow2(M_PI)), 
                     0.5 * M_PI * rootLog );
    } 
    else phi = complex( pow2( asin(1. / sqrt(epsilon)) ), 0.);

    // Expressions for quarks and leptons that depend on Higgs type.
    if (id < 17) { 
      if (higgsType < 3) etaNow = -0.5 * epsilon 
        * (complex(1., 0.) + (1. - epsilon) * phi);
      else etaNow = -0.5 * epsilon * phi;    
      if (id < 7 && id%2 == 1) etaNow *= 3. * pow2(ef) * coup2d;
      else if (id < 7 )        etaNow *= 3. * pow2(ef) * coup2u;
      else                     etaNow *=      pow2(ef) * coup2l;
    } 

    // Expression for W+-.
    else if (id == 24) etaNow = (complex(0.5 + 0.75 * epsilon, 0.)
      + 0.75 * epsilon * (2. - epsilon) * phi) * coup2W;  
 
    // Expression for H+-.
   else etaNow = (complex(epsilon, 0.) - epsilon * epsilon * phi)
     * pow2(mW / mHchg) * coup2Hchg;      
    
    // Sum up contribution and return square of absolute value.
    eta += etaNow;
  }
  return (pow2(eta.real()) + pow2(eta.imag()));

}

//*********

// Sum up quark, lepton, W+- and (for BSM) H+- loop contributions 
// in Higgs -> gamma + Z0.

double ResonanceH::eta2gaZ(double mHat) {

  // Initial values.
  complex eta = complex(0., 0.);
  int     id;
  double  ef, vf, mRun, epsilon, epsPrime, root, rootLog, asinEps;
  complex phi, psi, phiPrime, psiPrime, fXY, f1, etaNow;

  // Loop over s, c, b, t, mu , tau, W+-, H+- flavours.
  for (int idLoop = 0; idLoop < 7; ++idLoop) {
    if      (idLoop < 4) id = idLoop + 3;
    else if (idLoop < 6) id = 2 * idLoop + 5;
    else if (idLoop < 7) id = 24;
    else                 id = 37;

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

    // Expressions for quarks and leptons that depend on Higgs type.
    if (id < 17) { 
      etaNow = (higgsType < 3) ? -fXY + 0.25 * f1 : 0.25 * f1;
      if (id < 7 && id%2 == 1) etaNow *= 3. * ef * vf * coup2d;
      else if (id < 7)         etaNow *= 3. * ef * vf * coup2u;
      else                     etaNow *=      ef * vf * coup2l;

    // Expression for W+-.
    } else if (id == 24) {
      double coef1  = 3. - sin2tW / cos2tW;
      double coefXY = (1. + 2. / epsilon) * sin2tW / cos2tW 
        - (5. + 2. / epsilon);
      etaNow = -cos2tW * (coef1 * f1 + coefXY * fXY) * coup2W; 

    // Expression for H+-.
    } else etaNow = (1. - 2. * sin2tW) * fXY * pow2(mW / mHchg) 
      * coup2Hchg;     
    
    // Sum up contribution and return square of absolute value.
    eta += etaNow;
  }
  return ( (pow2(eta.real()) + pow2(eta.imag())) / (sin2tW * cos2tW) );

}
 
//**************************************************************************

// The ResonanceHchg class.
// Derived class for H+- properties.

//*********

// Initialize static data members.
// Calculate and store partial and total widths at the nominal mass. 

void ResonanceHchg::init() {

  // Do or redo initialization of basic quantities.
  initBasic();

  // Locally stored properties and couplings.
  useCubicWidth = Settings::flag("Higgs:cubicWidth");
  thetaWRat     = 1. / (8. * CoupEW::sin2thetaW());
  mW            = ParticleDataTable::m0(24);
  tanBeta       = Settings::parm("HiggsHchg:tanBeta");
  tan2Beta      = tanBeta * tanBeta;
  coup2H1W      = Settings::parm("HiggsHchg:coup2H1W");
  
  // Common coupling factors.
  double mHat   = mRes;
  double sH     = mHat*mHat;
  double alpEM  = alphaEM.alphaEM(sH);
  double alpS   = alphaS.alphaS(sH);
  double colQ   = 3. * (1. + alpS / M_PI);
  double preFac = alpEM * thetaWRat * pow3(mHat) / pow2(mW); 

  // Reset quantities to sum. Declare variables inside loop.
  double widTot = 0.; 
  double widPos = 0.;
  double widNeg = 0.;
  int    id1Abs, id2Abs, onMode;
  double widNow, widSecPos, widSecNeg, mf1, mf2, mr1, mr2, ps, 
         mRun1, mRun2, mrRunDn, mrRunUp;

  // Loop over all decay channels. Check that above threshold.
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow      = 0.;
    widSecPos   = 1.;
    widSecNeg   = 1.;
    id1Abs      = abs( particlePtr->decay[i].product(0) );
    id2Abs      = abs( particlePtr->decay[i].product(1) );
    mf1         = ParticleDataTable::m0(id1Abs);
    mf2         = ParticleDataTable::m0(id2Abs);
    if (mHat < mf1 + mf2 + MASSMARGIN) continue;
    mr1         = pow2(mf1 / mHat);
    mr2         = pow2(mf2 / mHat);
    ps          = sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 );

    // H+- decay to fermions involves running masses.
    if (id1Abs < 17 && id2Abs < 17 && (id1Abs < 7 || id1Abs > 10)) {
      mRun1     = ParticleDataTable::mRun(id1Abs, mHat);
      mRun2     = ParticleDataTable::mRun(id2Abs, mHat);
      mrRunDn   = pow2(mRun1 / mHat);
      mrRunUp   = pow2(mRun2 / mHat);
      if (id1Abs%2 == 0) swap( mrRunDn, mrRunUp);

      // Width to fermions. 
      widNow    = preFac * max( 0., (mrRunDn * tan2Beta + mrRunUp / tan2Beta) 
                * (1. - mrRunDn - mrRunUp) - 4. *mrRunDn * mrRunUp ) * ps;
      if (id1Abs < 7) widNow *= colQ;

      // Secondary width to top/antitop.
      if (id1Abs == 6 || id2Abs == 6) {
        widSecPos = ParticleDataTable::resOpenFrac( 6);
        widSecNeg = ParticleDataTable::resOpenFrac(-6);
      }
    }

    // H+- decay to h0 + W+-.
    else if (max(id1Abs, id2Abs) == 25 && min(id1Abs, id2Abs) == 24) {
      widNow    = 0.5 * preFac * pow3(ps) * pow2(coup2H1W);
      widSecPos = ParticleDataTable::resOpenFrac(25,  24); 
      widSecNeg = ParticleDataTable::resOpenFrac(25, -24); 
    }

    // Store partial widths and update sum, also for open channels only.
    particlePtr->decay[i].onShellWidth(widNow); 
    widTot += widNow;    
    onMode = particlePtr->decay[i].onMode();
    if (onMode == 1 || onMode == 2) widPos += widNow * widSecPos;
    if (onMode == 1 || onMode == 3) widNeg += widNow * widSecNeg;
  }

  // Update total width and branching ratios.
  if (!doForceWidth) particlePtr->setMWidth(widTot, false);
  double bRatio;
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    bRatio = particlePtr->decay[i].onShellWidth() / widTot;
    particlePtr->decay[i].bRatio( bRatio, false);
  }

  // Update the width in this object itself. Secondary widths for open.
  forceFactor   = GammaRes / widTot;
  if (doForceWidth) for (int i = 0; i < particlePtr->decay.size(); ++i)
    particlePtr->decay[i].onShellWidthFactor( forceFactor); 
  else GammaRes = widTot;
  GamMRat       = GammaRes / mRes;  
  openPos       = widPos / widTot;
  openNeg       = widNeg / widTot;

}

//*********

// Calculate the total width and store phases-space-weighted coupling sums.

double ResonanceHchg::width(int idSgn, double mHat, int , 
  bool openOnly, bool setBR) {

  // Common coupling factors.
  double sH     = mHat*mHat;
  double alpEM  = alphaEM.alphaEM(sH); 
  double alpS   = alphaS.alphaS(sH);
  double colQ   = 3. * (1. + alpS / M_PI);
  double preFac = alpEM * thetaWRat * pow3(mHat) / pow2(mW); 
  if (doForceWidth) preFac *= forceFactor;

  // Reset quantities to sum. Declare variables inside loop.
  double widSum = 0.; 
  int    id1Abs, id2Abs, onMode;
  double widNow, mf1, mf2, mr1, mr2, ps, mRun1, mRun2, mrRunDn, mrRunUp;

  // Loop over all decay channels; optionally only open channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow = 0.;
    onMode = particlePtr->decay[i].onMode();
    if (openOnly) {
      if (idSgn > 0 && onMode !=1 && onMode != 2) continue;
      if (idSgn < 0 && onMode !=1 && onMode != 3) continue;
    }
    id1Abs = abs( particlePtr->decay[i].product(0) );
    id2Abs = abs( particlePtr->decay[i].product(1) );

    // Check that above threshold.
    mf1         = ParticleDataTable::m0(id1Abs);
    mf2         = ParticleDataTable::m0(id2Abs);
    if (mHat < mf1 + mf2 + MASSMARGIN) continue;
    mr1         = pow2(mf1 / mHat);
    mr2         = pow2(mf2 / mHat);
    ps          = sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 );

    // H+- decay to fermions involves running masses.
    if (id1Abs < 17 && id2Abs < 17 && (id1Abs < 7 || id1Abs > 10)) {
      mRun1     = ParticleDataTable::mRun(id1Abs, mHat);
      mRun2     = ParticleDataTable::mRun(id2Abs, mHat);
      mrRunDn   = pow2(mRun1 / mHat);
      mrRunUp   = pow2(mRun2 / mHat);
      if (id1Abs%2 == 0) swap( mrRunDn, mrRunUp);

      // Width to fermions. 
      widNow    = preFac * max( 0., (mrRunDn * tan2Beta + mrRunUp / tan2Beta) 
                * (1. - mrRunDn - mrRunUp) - 4. *mrRunDn * mrRunUp ) * ps;
      if (id1Abs < 7) widNow *= colQ;

      // Secondary width to top/antitop.
      if (openOnly && (id1Abs == 6 || id2Abs == 6)) {
        int idTop = (idSgn > 0) ? 6 : -6;
        widNow   *= ParticleDataTable::resOpenFrac(idTop);
      }
    }

    // H+- decay to h0 + W+-.
    else if (max(id1Abs, id2Abs) == 25 && min(id1Abs, id2Abs) == 24) {
      widNow    = 0.5 * preFac * pow3(ps) * pow2(coup2H1W);
      if (!useCubicWidth) widNow *= pow2(mRes / mHat);
      int idW = (idSgn > 0) ? 24 : -24;
      widNow   *= ParticleDataTable::resOpenFrac(25, idW);
    }

    // Sum back up.
    widSum += widNow;

    // Optionally store partial widths for later decay channel choice.
    if (setBR) particlePtr->decay[i].currentBR(widNow); 
  }

  // Done.
  return widSum;
  
}
 
//**************************************************************************

// The ResonanceGraviton class.
// Derived class for excited Graviton properties, for extra dimensions.

//*********

// Initialize data members.
// Calculate and store partial and total widths at the nominal mass. 

void ResonanceGraviton::init() {

  // Do or redo initialization of basic quantities.
  initBasic();

  // Locally stored properties and couplings: kappa * m_G*.
  kappaMG       = Settings::parm("ExtraDimensionsG*:kappaMG");

  // Common coupling factors.
  double mHat   = mRes;
  double sH     = mHat*mHat;
  double alpS   = alphaS.alphaS(sH);
  double colQ   = 3. * (1. + alpS / M_PI);
  double preFac = pow2(kappaMG) * mHat / M_PI;

  // Reset quantities to sum. Declare variables inside loop.
  double widTot = 0.; 
  double widPos = 0.;
  int    idAbs, onMode;
  double widNow, widSec, mf, mr;

  // Loop over all decay channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow      = 0.;
    widSec      = 1.;
    idAbs       = abs( particlePtr->decay[i].product(0) );

    // Check that above threshold.
    mf          = ParticleDataTable::m0(idAbs);
    if (mHat > 2. * mf + MASSMARGIN) {
      mr        = pow2(mf / mHat);

      // Widths to fermion pairs.
      if (idAbs < 19) {
        widNow  = preFac * pow3( sqrtpos(1. - 4. * mr) ) 
          * (1. + 8. * mr / 3.) / 320.;        
        if (idAbs < 9) widNow *= colQ;
        widSec = ParticleDataTable::resOpenFrac(idAbs, -idAbs);
      }

      // Widths to gluon and photon pair.
      else if (idAbs == 21) widNow = preFac / 20.;
      else if (idAbs == 22) widNow = preFac / 160.;
     
      // Widths to Z0 Z0 and W+ W- pair.
      else if (idAbs == 23 || idAbs == 24) {
        widNow  = preFac * sqrtpos(1. - 4. * mr)
          * (13. / 12. + 14. * mr / 3. + 4. * mr * mr) / 80.;
        if (idAbs == 23) widNow *= 0.5;
        widSec = (idAbs == 23) ? ParticleDataTable::resOpenFrac(23,  23)
                               : ParticleDataTable::resOpenFrac(24, -24);
      }
    }

    // Store partial widths and update sum, also for open channels only.
    particlePtr->decay[i].onShellWidth(widNow); 
    widTot     += widNow;    
    onMode      = particlePtr->decay[i].onMode();
    if (onMode == 1 || onMode == 2) widPos += widNow * widSec;
  }

  // Update total width and branching ratios.
  if (!doForceWidth) particlePtr->setMWidth(widTot, false);
  double bRatio;
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    bRatio = particlePtr->decay[i].onShellWidth() / widTot;
    particlePtr->decay[i].bRatio( bRatio, false);
  }

  // Update the width in this object itself. Secondary widths for open.
  forceFactor   = GammaRes / widTot;
  if (doForceWidth) for (int i = 0; i < particlePtr->decay.size(); ++i)
    particlePtr->decay[i].onShellWidthFactor( forceFactor); 
  else GammaRes = widTot;
  GamMRat       = GammaRes / mRes;  
  openPos       = widPos / widTot;
    
}

//*********

// Calculate the total width and store phase-space-weighted coupling sums.

double ResonanceGraviton::width(int , double mHat, int , 
  bool openOnly, bool setBR) {

  // Common coupling factors.
  double sH      = mHat*mHat;
  double alpS    = alphaS.alphaS(sH);
  double colQ    = 3. * (1. + alpS / M_PI);
  double preFac  = pow2(kappaMG / M_PI) * mHat;
  if (doForceWidth) preFac *= forceFactor;

  // Reset quantities to sum. Declare variables inside loop.
  double widSum  = 0.; 
  int    idAbs, onMode;
  double widNow, mf, mr;

  // Loop over all decay channels; optionally only open channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow       = 0.;
    onMode       = particlePtr->decay[i].onMode();
    if (openOnly && onMode !=1 && onMode != 2) continue;
    idAbs        = abs( particlePtr->decay[i].product(0) );

    // Check that above threshold.
    mf          = ParticleDataTable::m0(idAbs);
    if (mHat > 2. * mf + MASSMARGIN) {
      mr        = pow2(mf / mHat);

      // Widths to fermion pairs.
      if (idAbs < 19) {
        widNow  = preFac * pow3( sqrtpos(1. - 4. * mr) ) 
          * (1. + 8. * mr / 3.) / 320.;        
        if (idAbs < 9) widNow *= colQ;
        if (openOnly) widNow 
                       *= ParticleDataTable::resOpenFrac(idAbs, -idAbs);
      }

      // Widths to gluon and photon pair.
      else if (idAbs == 21) widNow = preFac / 20.;
      else if (idAbs == 22) widNow = preFac / 160;
     
      // Widths to Z0 Z0 and W+ W- pair.
      else if (idAbs == 23 || idAbs == 24) {
        widNow  = preFac * sqrtpos(1. - 4. * mr) 
          * (13. / 12. + 14. + mr / 3. + 4. * mr * mr) / 80.;
        if (idAbs == 23) widNow *= 0.5;
        if (openOnly) widNow *= (idAbs == 23) 
                              ? ParticleDataTable::resOpenFrac(23,  23)
                              : ParticleDataTable::resOpenFrac(24, -24);
      }
    }

    // Sum back up.
    widSum += widNow;

    // Optionally store partial widths for later decay channel choice.
    if (setBR) particlePtr->decay[i].currentBR(widNow); 
  }

  // Done.
  return widSum;
  
}
 
//**************************************************************************

// The ResonanceLeptoquark class.
// Derived class for H+- properties.

//*********

// Initialize static data members.
// Calculate and store partial and total widths at the nominal mass. 

void ResonanceLeptoquark::init() {

  // Do or redo initialization of basic quantities.
  initBasic();

  // Locally stored properties and couplings.
  kCoup         = Settings::parm("LeptoQuark:kCoup");

  // Check that flavour info in decay channel is correctly set.
  int id1 = particlePtr->decay[0].product(0);
  int id2 = particlePtr->decay[0].product(1);
  if (id1 < 1 || id1 > 5) {
    ErrorMsg::message("Error in ResonanceLeptoquark::init:"
      " unallowed input quark flavour reset to u"); 
    id1   = 2;
    particlePtr->decay[0].product(0, id1);
  }
  if (abs(id2) < 11 || abs(id2) > 16) {
    ErrorMsg::message("Error in ResonanceLeptoquark::init:"
      " unallowed input lepton flavour reset to e-"); 
    id2   = 11;
    particlePtr->decay[0].product(1, id2);
  }

  // Set/overwrite charge and name of particle.
  bool changed  = particlePtr->hasChanged();
  int chargeLQ  = ParticleDataTable::chargeType(id1) 
                + ParticleDataTable::chargeType(id2);
  particlePtr->setChargeType(chargeLQ); 
  string nameLQ = "LQ_" + ParticleDataTable::name(id1) + ","
                + ParticleDataTable::name(id2);
  particlePtr->setNames(nameLQ, nameLQ + "bar"); 
  if (!changed) particlePtr->setHasChanged(false);
  
  // Common coupling factors.
  double mHat   = mRes;
  double sH     = mHat*mHat;
  double alpEM  = alphaEM.alphaEM(sH);
  double preFac = 0.25 * alpEM * kCoup * mHat; 

  // Reset quantities to sum. Declare variables inside loop.
  double widTot = 0.; 
  double widPos = 0.;
  double widNeg = 0.;
  int    id1Abs, id2Abs, onMode;
  double widNow, mf1, mf2, mr1, mr2, ps;

  // Loop over all decay channels. Check that above threshold.
  // Note: only one decay channel for now, but don't exclude more later.
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow      = 0.;
    id1Abs      = abs( particlePtr->decay[i].product(0) );
    id2Abs      = abs( particlePtr->decay[i].product(1) );
    mf1         = ParticleDataTable::m0(id1Abs);
    mf2         = ParticleDataTable::m0(id2Abs);
    if (mHat < mf1 + mf2 + MASSMARGIN) continue;
    mr1         = pow2(mf1 / mHat);
    mr2         = pow2(mf2 / mHat);
    ps          = sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 );

    // Width. Only allow first to be quark and second lepton.
    if (id1Abs < 6 && id2Abs > 10 && id2Abs < 17) 
      widNow    = preFac * pow3(ps);

    // Store partial widths and update sum, also for open channels only.
    particlePtr->decay[i].onShellWidth(widNow); 
    widTot += widNow;    
    onMode = particlePtr->decay[i].onMode();
    if (onMode == 1 || onMode == 2) widPos += widNow;
    if (onMode == 1 || onMode == 3) widNeg += widNow;
  }

  // Update total width and branching ratios.
  if (!doForceWidth) particlePtr->setMWidth(widTot, false);
  double bRatio;
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    bRatio = particlePtr->decay[i].onShellWidth() / widTot;
    particlePtr->decay[i].bRatio( bRatio, false);
  }

  // Update the width in this object itself. Secondary widths for open.
  forceFactor   = GammaRes / widTot;
  if (doForceWidth) for (int i = 0; i < particlePtr->decay.size(); ++i)
    particlePtr->decay[i].onShellWidthFactor( forceFactor); 
  else GammaRes = widTot;
  GamMRat       = GammaRes / mRes;  
  openPos       = widPos / widTot;
  openNeg       = widNeg / widTot;

}

//*********

// Calculate the total width and store phases-space-weighted coupling sums.

double ResonanceLeptoquark::width(int idSgn, double mHat, int , 
  bool openOnly, bool setBR) {

  // Common coupling factors.
  double sH     = mHat*mHat;
  double alpEM  = alphaEM.alphaEM(sH); 
  double preFac = 0.25 * alpEM * kCoup * mHat; 
  if (doForceWidth) preFac *= forceFactor;

  // Reset quantities to sum. Declare variables inside loop.
  double widSum = 0.; 
  int    id1Abs, id2Abs, onMode;
  double widNow, mf1, mf2, mr1, mr2, ps;

  // Loop over all decay channels; optionally only open channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    widNow = 0.;
    onMode = particlePtr->decay[i].onMode();
    if (openOnly) {
      if (idSgn > 0 && onMode !=1 && onMode != 2) continue;
      if (idSgn < 0 && onMode !=1 && onMode != 3) continue;
    }
    id1Abs = abs( particlePtr->decay[i].product(0) );
    id2Abs = abs( particlePtr->decay[i].product(1) );

    // Check that above threshold.
    mf1         = ParticleDataTable::m0(id1Abs);
    mf2         = ParticleDataTable::m0(id2Abs);
    if (mHat < mf1 + mf2 + MASSMARGIN) continue;
    mr1         = pow2(mf1 / mHat);
    mr2         = pow2(mf2 / mHat);
    ps          = sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 );

    // Width. Only allow first to be quark and second lepton.
    if (id1Abs < 6 && id2Abs > 10 && id2Abs < 17) 
      widNow    = preFac * pow3(ps);

    // Sum back up.
    widSum += widNow;

    // Optionally store partial widths for later decay channel choice.
    if (setBR) particlePtr->decay[i].currentBR(widNow); 
  }

  // Done.
  return widSum;
  
}

//**************************************************************************

} // end namespace Pythia8
