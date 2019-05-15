// Function definitions (not found in the header) for the AlphaStrong class.
// Copyright © 2005 Torbjörn Sjöstrand

#include "StandardModel.h"

namespace Pythia8 {

//**************************************************************************

// The AlphaStrong class.

//*********

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of iterations to determine Lambda from given alpha_s.
const int AlphaStrong::NITER = 10;

//*********

// Initialize alpha_strong calculation by finding Lambda values etc.

void AlphaStrong::init( double valueIn, int orderIn) {

  // Order of alpha_s evaluation. Charm, bottom and Z0 masses. 
  // Pick defaults if ParticleDataTable not properly initialized.
  valueRef = valueIn;
  order = max( 0, min( 2, orderIn ) );
  mc = ParticleDataTable::m0(4);
  if (mc < 1. || mc > 2.) mc = 1.5;
  mb = ParticleDataTable::m0(5);
  if (mb < 3. || mb > 6.) mb = 4.8;
  mZ = ParticleDataTable::m0(23);
  if (mZ < 90. || mZ > 92.) mZ = 91.188;

  // Fix alpha_s.
  if (order == 0) {
    Lambda3Save = Lambda4Save = Lambda5Save = 0.;

  // First order alpha_s: match at flavour thresholds.
  } else if (order == 1) {
    Lambda5Save = mZ * exp( -6. * M_PI / (23. * valueRef) );
    Lambda4Save = Lambda5Save * pow(mb/Lambda5Save, 2./25.); 
    Lambda3Save = Lambda4Save * pow(mc/Lambda4Save, 2./27.); 

  // Second order alpha_s: iterative match at flavour thresholds.
  } else {

    // Find Lambda_5 at m_Z.
    Lambda5Save = mZ * exp( -6. * M_PI / (23. * valueRef) );
    for (int iter = 0; iter < NITER; ++iter) {
      double logScale = 2. * log(mZ/Lambda5Save);
      double correction = 1. - (348./529.) * log(logScale) / logScale;
      double valueIter = valueRef / correction; 
      Lambda5Save = mZ * exp( -6. * M_PI / (23. * valueIter) );
    }

    // Find Lambda_4 at m_b.
    double logScaleB = 2. * log(mb/Lambda5Save);
    double valueB = 12. * M_PI / (23. * logScaleB) 
        * (1. - (348./529.) * log(logScaleB) / logScaleB);
    Lambda4Save = Lambda5Save;
    for (int iter = 0; iter < NITER; ++iter) {
      double logScale = 2. * log(mb/Lambda4Save);
      double correction = 1. - (462./625.) * log(logScale) / logScale;
      double valueIter = valueB / correction; 
      Lambda4Save = mb * exp( -6. * M_PI / (25. * valueIter) );
    }

    // Find Lambda_3 at m_c.
    double logScaleC = 2. * log(mc/Lambda4Save);
    double valueC = 12. * M_PI / (25. * logScaleC) 
        * (1. - (462./625.) * log(logScaleC) / logScaleC);
    Lambda3Save = Lambda4Save;
    for (int iter = 0; iter < NITER; ++iter) {
      double logScale = 2. * log(mc/Lambda3Save);
      double correction = 1. - (64./81.) * log(logScale) / logScale;
      double valueIter = valueC / correction; 
      Lambda3Save = mc * exp( -6. * M_PI / (27. * valueIter) );
    }
  }

  // Save squares of mass and Lambda values as well.
  mc2 = pow2(mc);
  mb2 = pow2(mb);
  Lambda3Save2 = pow2(Lambda3Save);
  Lambda4Save2 = pow2(Lambda4Save);
  Lambda5Save2 = pow2(Lambda5Save);
  valueNow = valueIn;
  scale2Now = mZ*mZ;
  isInit = true;

};

//*********

// Calculate alpha_s value    

double AlphaStrong::alphaS( double scale2) {

  // If equal to old scale then same answer.
  if (!isInit) return 0.;
  if (scale2 == scale2Now && (order < 2 || lastCallToFull)) return valueNow;
  scale2Now = scale2;
  lastCallToFull = true;

  // Fix alpha_s.
  if (order == 0) {
    valueNow = valueRef;        
  
  // First order alpha_s: differs by mass region.  
  } else if (order == 1) {
    if (scale2 > mb2) valueNow 
      = 12. * M_PI / (23. * log(scale2/Lambda5Save2));
    else if (scale2 > mc2) valueNow 
      = 12. * M_PI / (25. * log(scale2/Lambda4Save2));
    else valueNow = 12. * M_PI / (27. * log(scale2/Lambda3Save2));
  
  // Second order alpha_s: differs by mass region.  
  } else {
    if (scale2 > mb2) { double logScale = log(scale2/Lambda5Save2);
      valueNow = 12. * M_PI / (23. * logScale) 
        * (1. - (348./529.) * log(logScale) / logScale); }
    else if (scale2 > mc2) { double logScale = log(scale2/Lambda4Save2); 
      valueNow = 12. * M_PI / (25. * logScale)
        * (1. - (462./625.) * log(logScale) / logScale); }
    else { double logScale = log(scale2/Lambda3Save2); 
      valueNow = 12. * M_PI / (27. * logScale)
        * (1. - (64./81.) * log(logScale) / logScale); }
  }

  // Done.
  return valueNow;
} ;

//*********

// Calculate alpha_s value, but only use up to first-order piece.
// (To be combined with alphaS2OrdCorr.)

double  AlphaStrong::alphaS1Ord( double scale2) {

  // If equal to old scale then same answer.
  if (!isInit) return 0.;
  if (scale2 == scale2Now && (order < 2 || !lastCallToFull)) return valueNow;
  scale2Now = scale2;
  lastCallToFull = false;

  // Fix alpha_S.
  if (order == 0) {
    valueNow = valueRef;        
  
  // First/second order alpha_s: differs by mass region.  
  } else {
    if (scale2 > mb2) valueNow 
      = 12. * M_PI / (23. * log(scale2/Lambda5Save2));
    else if (scale2 > mc2) valueNow 
      = 12. * M_PI / (25. * log(scale2/Lambda4Save2));
    else valueNow = 12. * M_PI / (27. * log(scale2/Lambda3Save2));
  }

  // Done.
  return valueNow;
} ;

//*********

// Calculates the second-order extra factor in alpha_s.
// (To be combined with alphaS1Ord.)

double AlphaStrong::alphaS2OrdCorr( double scale2) {

  // Only meaningful for second order calculations.
  if (!isInit) return 1.;
  if (order < 2) return 1.; 
  
  // Second order correction term: differs by mass region.  
  if (scale2 > mb2) { double logScale = log(scale2/Lambda5Save2);
    return 1. - (348./529.) * log(logScale) / logScale; }
  if (scale2 > mc2) { double logScale = log(scale2/Lambda4Save2); 
    return 1. - (462./625.) * log(logScale) / logScale; }
  double logScale = log(scale2/Lambda3Save2); 
    return 1. - (64./81.) * log(logScale) / logScale; 
} ;

//**************************************************************************

} // end namespace Pythia8
