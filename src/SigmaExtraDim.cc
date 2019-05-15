// SigmaExtraDim.cc is a part of the PYTHIA event generator.
// Copyright (C) 2009 Torbjorn Sjostrand.
// Copyright (C) 2009 Stefan Ask for the *LED* routines.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the 
// extra-dimensional simulation classes. 

#include "SigmaExtraDim.h"

namespace Pythia8 {

//**************************************************************************

// Sigma1gg2GravitonStar class.
// Cross section for g g -> G* (excited graviton state). 

//*********

// Initialize process. 
  
void Sigma1gg2GravitonStar::initProc() {

  // Store G* mass and width for propagator. 
  idGstar  = 5100039;
  mRes     = ParticleDataTable::m0(idGstar);
  GammaRes = ParticleDataTable::mWidth(idGstar);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

  // Overall coupling strength kappa * m_G*.
  kappaMG  = Settings::parm("ExtraDimensionsG*:kappaMG");

  // Set pointer to particle properties and decay table.
  gStarPtr = ParticleDataTable::particleDataPtr(idGstar);

} 

//*********

// Evaluate sigmaHat(sHat), part independent of incoming flavour. 

void Sigma1gg2GravitonStar::sigmaKin() { 

  // Incoming width for gluons.
  double widthIn  = pow2(kappaMG) * mH / (160. * M_PI);

  // Set up Breit-Wigner. Width out only includes open channels. 
  double sigBW    = 5. * M_PI/ ( pow2(sH - m2Res) + pow2(sH * GamMRat) );    
  double widthOut = gStarPtr->resWidthOpen(idGstar, mH);

  // Modify cross section in wings of peak. Done.
  sigma           = widthIn * sigBW * widthOut * pow2(sH / m2Res);    

}

//*********

// Select identity, colour and anticolour.

void Sigma1gg2GravitonStar::setIdColAcol() {

  // Flavours trivial.
  setId( 21, 21, idGstar);

  // Colour flow topology.
  setColAcol( 1, 2, 2, 1, 0, 0);

}

//*********

// Evaluate weight for G* decay angle.
  
double Sigma1gg2GravitonStar::weightDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For top decay hand over to standard routine.
  if (idMother == 6) 
    return weightTopDecay( process, iResBeg, iResEnd);

  // G* should sit in entry 5.
  if (iResBeg != 5 || iResEnd != 5) return 1.;

  // Phase space factors. Reconstruct decay angle.
  double mr1    = pow2(process[6].m()) / sH;
  double mr2    = pow2(process[7].m()) / sH;
  double betaf  = sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2); 
  double cosThe = (process[3].p() - process[4].p()) 
    * (process[7].p() - process[6].p()) / (sH * betaf);

  // Default is isotropic decay.
  double wt     = 1.;

  // Angular weight for g + g -> G* -> f + fbar.
  if (process[6].idAbs() < 19) wt = 1. - pow4(cosThe);

  // Angular weight for g + g -> G* -> g + g or gamma + gamma.
  else if (process[6].id() == 21 || process[6].id() == 22)
    wt = (1. + 6. * pow2(cosThe) + pow4(cosThe)) / 8.;
 
  // Done.
  return wt;

}

//**************************************************************************

// Sigma1ffbar2GravitonStar class.
// Cross section for f fbar -> G* (excited graviton state). 

//*********

// Initialize process. 
  
void Sigma1ffbar2GravitonStar::initProc() {

  // Store G* mass and width for propagator. 
  idGstar  = 5100039;
  mRes     = ParticleDataTable::m0(idGstar);
  GammaRes = ParticleDataTable::mWidth(idGstar);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

  // Overall coupling strength kappa * m_G*.
  kappaMG  = Settings::parm("ExtraDimensionsG*:kappaMG");

  // Set pointer to particle properties and decay table.
  gStarPtr = ParticleDataTable::particleDataPtr(idGstar);

} 

//*********

// Evaluate sigmaHat(sHat), part independent of incoming flavour. 

void Sigma1ffbar2GravitonStar::sigmaKin() { 

  // Incoming width for fermions, disregarding colour factor.
  double widthIn  = pow2(kappaMG) * mH / (80. * M_PI);

  // Set up Breit-Wigner. Width out only includes open channels. 
  double sigBW    = 5. * M_PI/ ( pow2(sH - m2Res) + pow2(sH * GamMRat) );    
  double widthOut = gStarPtr->resWidthOpen(idGstar, mH);

  // Modify cross section in wings of peak. Done.
  sigma0          = widthIn * sigBW * widthOut * pow2(sH / m2Res);    

}

//*********

// Select identity, colour and anticolour.

void Sigma1ffbar2GravitonStar::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, idGstar);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//*********

// Evaluate weight for G* decay angle.
  
double Sigma1ffbar2GravitonStar::weightDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For top decay hand over to standard routine.
  if (idMother == 6) 
    return weightTopDecay( process, iResBeg, iResEnd);

  // G* should sit in entry 5.
  if (iResBeg != 5 || iResEnd != 5) return 1.;

  // Phase space factors. Reconstruct decay angle.
  double mr1    = pow2(process[6].m()) / sH;
  double mr2    = pow2(process[7].m()) / sH;
  double betaf  = sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2); 
  double cosThe = (process[3].p() - process[4].p()) 
    * (process[7].p() - process[6].p()) / (sH * betaf);

  // Default is isotropic decay.
  double wt     = 1.;

  // Angular weight for f + fbar -> G* -> f + fbar.
  if (process[6].idAbs() < 19)
    wt = (1. - 3. * pow2(cosThe) + 4. * pow4(cosThe)) / 2.;

  // Angular weight for f + fbar -> G* -> g + g or gamma + gamma.
  else if (process[6].id() == 21 || process[6].id() == 22)
    wt = 1. - pow4(cosThe);
 
  // Done.
  return wt;

}

//**************************************************************************

// Sigma2gg2GravitonStarg class.
// Cross section for g g -> G* g (excited graviton state). 

//*********

// Initialize process. 
  
void Sigma2gg2GravitonStarg::initProc() {

  // Store G* mass and width for propagator. 
  idGstar  = 5100039;
  mRes     = ParticleDataTable::m0(idGstar);
  GammaRes = ParticleDataTable::mWidth(idGstar);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

  // Overall coupling strength kappa * m_G*.
  kappaMG  = Settings::parm("ExtraDimensionsG*:kappaMG");

   // Secondary open width fraction.
  openFrac = ParticleDataTable::resOpenFrac(idGstar);

} 

//*********

// Evaluate sigmaHat(sHat), part independent of incoming flavour. 

void Sigma2gg2GravitonStarg::sigmaKin() { 

  //  Evaluate cross section. Secondary width for G*.
  sigma = (3. * pow2(kappaMG) * alpS) / (32. * sH * s3)
    * ( pow2(tH2 + tH * uH + uH2) / (sH2 * tH * uH) 
    + 2. * (tH2 / uH + uH2 / tH) / sH + 3. * (tH / uH + uH / tH)
    + 2. * (sH / uH + sH/tH) + sH2 / (tH * uH) );
  sigma *= openFrac;

}

//*********

// Select identity, colour and anticolour.

void Sigma2gg2GravitonStarg::setIdColAcol() {

  // Flavours trivial.
  setId( 21, 21, idGstar, 21);

  // Colour flow topologies: random choice between two mirrors.
  if (Rndm::flat() < 0.5) setColAcol( 1, 2, 2, 3, 0, 0, 1, 3);
  else                    setColAcol( 1, 2, 3, 1, 0, 0, 3, 2);

}

//*********

// Evaluate weight for decay angles: currently G* assumed isotropic.
  
double Sigma2gg2GravitonStarg::weightDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For top decay hand over to standard routine.
  if (idMother == 6) 
    return weightTopDecay( process, iResBeg, iResEnd);

  // No equations for G* decay so assume isotropic.
  return 1.;

}

//**************************************************************************

// Sigma2qg2GravitonStarq class.
// Cross section for q g -> G* q (excited graviton state). 

//*********

// Initialize process. 
  
void Sigma2qg2GravitonStarq::initProc() {

  // Store G* mass and width for propagator. 
  idGstar  = 5100039;
  mRes     = ParticleDataTable::m0(idGstar);
  GammaRes = ParticleDataTable::mWidth(idGstar);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

  // Overall coupling strength kappa * m_G*.
  kappaMG  = Settings::parm("ExtraDimensionsG*:kappaMG");

   // Secondary open width fraction.
  openFrac = ParticleDataTable::resOpenFrac(idGstar);

} 

//*********

// Evaluate sigmaHat(sHat), part independent of incoming flavour. 

void Sigma2qg2GravitonStarq::sigmaKin() { 

  //  Evaluate cross section. Secondary width for G*.
  sigma = -(pow2(kappaMG) * alpS) / (192. * sH * s3)
    * ( 4. * (sH2 + uH2) / (tH * sH) + 9. * (sH + uH) / sH + sH / uH
    + uH2 / sH2 + 3. * tH * (4. + sH / uH + uH / sH) / sH
    + 4. * tH2 * (1. / uH + 1. / sH) / sH + 2. * tH2 * tH / (uH * sH2) );
  sigma *= openFrac;

}

//*********

// Select identity, colour and anticolour.

void Sigma2qg2GravitonStarq::setIdColAcol() {

  // Flavour set up for q g -> H q.
  int idq = (id2 == 21) ? id1 : id2;
  setId( id1, id2, idGstar, idq);

  // tH defined between f and f': must swap tHat <-> uHat if q g in.
  swapTU = (id2 == 21); 

  // Colour flow topologies. Swap when antiquarks.
  if (id2 == 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else           setColAcol( 2, 1, 1, 0, 0, 0, 2, 0);
  if (idq < 0) swapColAcol();

}

//*********

// Evaluate weight for decay angles: currently G* assumed isotropic.
  
double Sigma2qg2GravitonStarq::weightDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For top decay hand over to standard routine.
  if (idMother == 6) 
    return weightTopDecay( process, iResBeg, iResEnd);

  // No equations for G* decay so assume isotropic.
  return 1.;

}

//**************************************************************************

// Sigma2qqbar2GravitonStarg class.
// Cross section for q qbar -> G* g (excited graviton state). 

//*********

// Initialize process. 
  
void Sigma2qqbar2GravitonStarg::initProc() {

  // Store G* mass and width for propagator. 
  idGstar  = 5100039;
  mRes     = ParticleDataTable::m0(idGstar);
  GammaRes = ParticleDataTable::mWidth(idGstar);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

  // Overall coupling strength kappa * m_G*.
  kappaMG  = Settings::parm("ExtraDimensionsG*:kappaMG");

   // Secondary open width fraction.
  openFrac = ParticleDataTable::resOpenFrac(idGstar);

} 

//*********

// Evaluate sigmaHat(sHat), part independent of incoming flavour. 

void Sigma2qqbar2GravitonStarg::sigmaKin() { 

  // Evaluate cross section. Secondary width for G*.
  sigma = (pow2(kappaMG) * alpS) / (72. * sH * s3)
    * ( 4. * (tH2 + uH2) / sH2 + 9. * (tH + uH) / sH 
    + (tH2 / uH + uH2 / tH) / sH + 3. * (4. + tH / uH + uH/ tH)
    + 4. * (sH / uH + sH / tH) + 2. * sH2 / (tH * uH) );
  sigma *= openFrac;

}

//*********

// Select identity, colour and anticolour.

void Sigma2qqbar2GravitonStarg::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, idGstar, 21);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 2, 0, 0, 1, 2);
  if (id1 < 0) swapColAcol();

}

//*********

// Evaluate weight for decay angles: currently G* assumed isotropic.
  
double Sigma2qqbar2GravitonStarg::weightDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // Identity of mother of decaying reseonance(s).
  int idMother = process[process[iResBeg].mother1()].idAbs();

  // For top decay hand over to standard routine.
  if (idMother == 6) 
    return weightTopDecay( process, iResBeg, iResEnd);

  // No equations for G* decay so assume isotropic.
  return 1.;

}

//**************************************************************************

// Sigma2gg2LEDGravitong class.
// Cross section for g g -> G g (real graviton emission in large extra 
// dimensions). 

//*********

void Sigma2gg2LEDGravitong::initProc() {
  
  //+++ Init model parameters.
  m_idG    = 5000039;
  m_trunc  = Settings::flag("ExtraDimensionsLED:Trunc"); 
  m_nGrav  = Settings::mode("ExtraDimensionsLED:n");
  m_MD     = Settings::parm("ExtraDimensionsLED:MD");

  //+++ Torus Surface
  double tmpTerm1 = sqrt( pow(M_PI, double(m_nGrav)) );
  double tmpTerm2 = funcGammaIntHalf(m_nGrav);
  double tmpExp   = m_nGrav + 2.0;
  m_constantTerm  =  tmpTerm1 / ( tmpTerm2 * pow(m_MD, tmpExp) );

} 

//*********

void Sigma2gg2LEDGravitong::sigmaKin() { 

  //+++ Set graviton mass and mandelstam variables
  mG        = m3;
  mGS       = mG*mG;

  double A0 = 1/sH;    
  double xH = tH/sH;
  double yH = mGS/sH;
  double xHS = pow2(xH);
  double yHS = pow2(yH);
  double xHC = pow(xH,3);
  double yHC = pow(yH,3);
  double xHQ = pow(xH,4);
  double yHQ = pow(yH,4);
  
  double T0 = 1/(xH*(yH-1-xH));
  double T1 = 1 + 2*xH + 3*xHS + 2*xHC + xHQ;
  double T2 = -2*yH*(1 + xHC);
  double T3 = 3*yHS*(1 + xHS);
  double T4 = -2*yHC*(1 + xH);
  double T5 = yHQ;
  
  m_sigma0 = A0 * T0 *( T1 + T2 + T3 + T4 + T5 );

  //+++ Couplings..
  m_sigma0 *= alpS * 3 / 16;

  //+++ Mass Spectrum, (m)^(n-2)
  double tmpExp = m_nGrav - 2;
  m_sigma0 *= pow(mG, tmpExp);

  //+++ Constants
  m_sigma0 *= m_constantTerm;

}

//*********

double Sigma2gg2LEDGravitong::sigmaHat() { 

  //+++ Mass spactrum weighting.
  double sigma = m_sigma0 /runBW3;      

  //+++ Truncation, to test perturbative region
  if (m_trunc) {
    if (sH > pow2(m_MD) ) { sigma *= pow(m_MD,4)/pow2(sH); }
  }
  
  return sigma;  
}

//*********

void Sigma2gg2LEDGravitong::setIdColAcol() {

 //+++ Flavours trivial.
  setId( 21, 21, m_idG, 21);

  //+++ Colour flow topologies: random choice between two mirrors.
  if (Rndm::flat() < 0.5) setColAcol( 1, 2, 2, 3, 0, 0, 1, 3);
  else                    setColAcol( 1, 2, 3, 1, 0, 0, 3, 2);

}

//*********

double Sigma2gg2LEDGravitong::funcGammaIntHalf( int n ) {

  double gamma = 1;
  if (n%2 == 0) {
    int k = int(n/2);
    for (int i = 1; i <= k-1; ++i)  {
      gamma *= i;
    }  
  } else {
    int k = int((n-1)/2);
    for (int i = 0; i <= k-1; ++i)  {
      gamma *= i+0.5;
    } 
    gamma *= sqrt(M_PI);
  }  

  return gamma;
}

//**************************************************************************

// Sigma2qg2LEDGravitonq class.
// Cross section for q g -> G q (real graviton emission in large extra 
// dimensions). 

//*********

void Sigma2qg2LEDGravitonq::initProc() {
  
  //+++ Init model parameters.
  m_idG    = 5000039;
  m_trunc  = Settings::flag("ExtraDimensionsLED:Trunc"); 
  m_nGrav  = Settings::mode("ExtraDimensionsLED:n");
  m_MD     = Settings::parm("ExtraDimensionsLED:MD");

  //+++ Torus Surface
  double tmpTerm1 = sqrt( pow(M_PI, double(m_nGrav)) );
  double tmpTerm2 = funcGammaIntHalf(m_nGrav);
  double tmpExp   = m_nGrav + 2.0;
  m_constantTerm  =  tmpTerm1 / ( tmpTerm2 * pow(m_MD, tmpExp) );

} 

//*********

void Sigma2qg2LEDGravitonq::sigmaKin() { 

  //+++ Set graviton mass and mandelstam variables
  mG        = m3;
  mGS       = mG*mG;

  double A0 = 1/sH;    
  double xH = tH/sH;
  double yH = mGS/sH;
  double x2H = xH/(yH - 1 - xH);
  double y2H = yH/(yH - 1 - xH);
  double x2HS = pow2(x2H);
  double y2HS = pow2(y2H);
  double x2HC = pow(x2H,3);
  double y2HC = pow(y2H,3);

  double T0 = -(yH - 1 - xH);
  double T20 = 1/(x2H*(y2H-1-x2H));
  double T21 = -4*x2H*(1 + x2H)*(1 + 2*x2H + 2*x2HS);
  double T22 = y2H*(1 + 6*x2H + 18*x2HS + 16*x2HC);
  double T23 = -6*y2HS*x2H*(1+2*x2H);
  double T24 = y2HC*(1 + 4*x2H);

  m_sigma0 = A0 * T0 * T20 * ( T21 + T22 + T23 + T24 );

  //+++ Couplings..
  m_sigma0 *= alpS / 96;

  //+++ Mass Spectrum, (m)^(n-2)
  double tmpExp = m_nGrav - 2;
  m_sigma0 *= pow(mG, tmpExp);

  //+++ Constants
  m_sigma0 *= m_constantTerm;

}

//*********

double Sigma2qg2LEDGravitonq::sigmaHat() { 

  //+++ Mass spactrum weighting.
  double sigma = m_sigma0 /runBW3;      

  //+++ Truncation, to test perturbative region
  if (m_trunc) {
    if (sH > pow2(m_MD) ) { sigma *= pow(m_MD,4)/pow2(sH); }
  }
  
  return sigma;  
}

//*********

void Sigma2qg2LEDGravitonq::setIdColAcol() {

  //+++ Flavour set up for q g -> G* q.
  int idq = (id2 == 21) ? id1 : id2;
  setId( id1, id2, m_idG, idq);

  //+++ tH defined between f and f': must swap tHat <-> uHat if q g in.
  swapTU = (id2 == 21); 

  //+++ Colour flow topologies. Swap when antiquarks.
  if (id2 == 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else           setColAcol( 2, 1, 1, 0, 0, 0, 2, 0);
  if (idq < 0) swapColAcol();

}

//*********

double Sigma2qg2LEDGravitonq::funcGammaIntHalf( int n ) {

  double gamma = 1;
  if (n%2 == 0) {
    int k = int(n/2);
    for (int i = 1; i <= k-1; ++i)  {
      gamma *= i;
    }  
  } else {
    int k = int((n-1)/2);
    for (int i = 0; i <= k-1; ++i)  {
      gamma *= i+0.5;
    } 
    gamma *= sqrt(M_PI);
  }  

  return gamma;
}

//**************************************************************************

// Sigma2qqbar2LEDGravitong class.
// Cross section for q qbar -> G g (real graviton emission in large extra 
// dimensions). 

//*********

void Sigma2qqbar2LEDGravitong::initProc() {
  
  //+++ Init model parameters.
  m_idG    = 5000039;
  m_trunc  = Settings::flag("ExtraDimensionsLED:Trunc"); 
  m_nGrav  = Settings::mode("ExtraDimensionsLED:n");
  m_MD     = Settings::parm("ExtraDimensionsLED:MD");

  //+++ Torus Surface
  double tmpTerm1 = sqrt( pow(M_PI, double(m_nGrav)) );
  double tmpTerm2 = funcGammaIntHalf(m_nGrav);
  double tmpExp   = m_nGrav + 2.0;
  m_constantTerm  =  tmpTerm1 / ( tmpTerm2 * pow(m_MD, tmpExp) );

} 

//*********

void Sigma2qqbar2LEDGravitong::sigmaKin() { 

  //+++ Set graviton mass and mandelstam variables
  mG        = m3;
  mGS       = mG*mG;

  double A0 = 1/sH;    
  double xH = tH/sH;
  double yH = mGS/sH;
  double xHS = pow2(xH);
  double yHS = pow2(yH);
  double xHC = pow(xH,3);
  double yHC = pow(yH,3);
  
  double T0 = 1/(xH*(yH-1-xH));
  double T1 = -4*xH*(1 + xH)*(1 + 2*xH + 2*xHS);
  double T2 = yH*(1 + 6*xH + 18*xHS + 16*xHC);
  double T3 = -6*yHS*xH*(1+2*xH);
  double T4 = yHC*(1 + 4*xH);
  
  m_sigma0 = A0 * T0 *( T1 + T2 + T3 + T4 );

  //+++ Couplings..
  m_sigma0 *= alpS / 36;

  //+++ Mass Spectrum, (m)^(n-2)
  double tmpExp = m_nGrav - 2;
  m_sigma0 *= pow(mG, tmpExp);

  //+++ Constants
  m_sigma0 *= m_constantTerm;

}

//*********

double Sigma2qqbar2LEDGravitong::sigmaHat() { 

  //+++ Mass spactrum weighting.
  double sigma = m_sigma0 /runBW3;      

  //+++ Truncation, to test perturbative region
  if (m_trunc) {
    if (sH > pow2(m_MD) ) { sigma *= pow(m_MD,4)/pow2(sH); }
  }
  
  return sigma;  
}

//*********

void Sigma2qqbar2LEDGravitong::setIdColAcol() {

  //+++ Flavours trivial.
  setId( id1, id2, m_idG, 21);

  //+++ Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 2, 0, 0, 1, 2);
  if (id1 < 0) swapColAcol();

}

//*********

double Sigma2qqbar2LEDGravitong::funcGammaIntHalf( int n ) {

  double gamma = 1;
  if (n%2 == 0) {
    int k = int(n/2);
    for (int i = 1; i <= k-1; ++i)  {
      gamma *= i;
    }  
  } else {
    int k = int((n-1)/2);
    for (int i = 0; i <= k-1; ++i)  {
      gamma *= i+0.5;
    } 
    gamma *= sqrt(M_PI);
  }  

  return gamma;
}

//**************************************************************************

// Sigma2ffbar2LEDUnparticleZ class.
// Cross section for f fbar -> U/G Z (real LED graviton or unparticle 
// emission).

//*********

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// FIXRATIO:
// Ratio between the two possible coupling constants of the spin-2 ME. 
// A value different from one give rise to an IR divergence which makes 
// the event generation very slow, so this values is fixed to 1 until
// investigated further.
const double Sigma2ffbar2LEDUnparticleZ::FIXRATIO = 1.;

//*********

Sigma2ffbar2LEDUnparticleZ::Sigma2ffbar2LEDUnparticleZ( bool Graviton ):
  m_graviton(Graviton) {

}

//*********

void Sigma2ffbar2LEDUnparticleZ::initProc() {

  //+++ Init model parameters.
  m_idG        = 5000039;
  if (m_graviton) {
    m_trunc    = Settings::flag("ExtraDimensionsLED:Trunc"); 
    m_nGrav    = Settings::mode("ExtraDimensionsLED:n");
    m_LambdaU  = Settings::parm("ExtraDimensionsLED:MD");
    m_dU       = m_nGrav/2 + 1; 
    m_spin     = 2;
  } else {
    m_trunc    = Settings::flag("ExtraDimensionsUnpart:Trunc"); 
    m_spin     = Settings::mode("ExtraDimensionsUnpart:spinU");
    m_dU       = Settings::parm("ExtraDimensionsUnpart:dU");
    m_LambdaU  = Settings::parm("ExtraDimensionsUnpart:LambdaU");
    m_lambda   = Settings::parm("ExtraDimensionsUnpart:lambda");
    m_ratio    = FIXRATIO; 
    //         = Settings::parm("ExtraDimensionsUnpart:ratio");
  }

  //+++ Store Z0 mass and width for propagator.
  mZ        = ParticleDataTable::m0(23);
  widZ      = ParticleDataTable::mWidth(23);
  mZS       = mZ*mZ;
  mwZS      = pow2(mZ * widZ);

  //+++ Init spin-2 parameters
  if ( m_spin != 2 ){
    m_graviton = false;
    m_lambdaPrime = 0;
  } else if (m_graviton) {
    m_lambda = 1;
    m_ratio = 1;
    m_lambdaPrime = m_lambda;
  } else {
    m_lambdaPrime = m_ratio * m_lambda;
  }

  //+++ The A(dU) or S'(n) value
  double tmpArg1   = m_dU + 0.5;
  double tmpGamma1 = funcGammaReal( tmpArg1, 100000);
  double tmpArg2   = m_dU - 1;
  double tmpGamma2 = funcGammaReal( tmpArg2, 100000);  
  double tmpArg3   = 2 * m_dU;
  double tmpGamma3 = funcGammaReal( tmpArg3, 100000);
  double tmpBase   = 2 * M_PI;
  double tmpExp1   = 2 * m_dU;
  double tmpAdU    = 16 * pow2(M_PI) * sqrt(M_PI) / pow(tmpBase, tmpExp1)
                   * tmpGamma1 / (tmpGamma2 * tmpGamma3);
  if (m_graviton) { 
    m_nGrav = 2*(m_dU - 1);
    tmpAdU  = 2 * M_PI * sqrt( pow(M_PI, m_nGrav) ) 
            / funcGammaIntHalf( int(m_nGrav) ); 
  } 

  //+++ Standard 2 to 2 cross section related constants
  double tmpTerm1 = 1/(2 * 16 * pow2(M_PI));
  double tmpLS    = pow2(m_LambdaU);

  //+++ Spin dependent constants from ME.
  //+++ Spin-0 is only place holder.
  double tmpTerm2 = 0;
  if ( m_spin == 0 ) { 

  } else if (m_spin == 1) {
    tmpTerm2 = 4 * pow2(m_lambda);
  } else if (m_spin == 2) {
    tmpTerm2 = pow2(m_lambda)/(4 * 3 * tmpLS);
  }

  //+++ Unparticle phase space related
  double tmpExp2 = m_dU - 2;
  double tmpTerm3 = tmpAdU / (tmpLS * pow(tmpLS, tmpExp2));

  //+++ All in total
  m_constantTerm = tmpTerm1 * tmpTerm2 * tmpTerm3;

} 

//*********

void Sigma2ffbar2LEDUnparticleZ::sigmaKin() { 

  //+++ Set graviton mass and some powers of mandelstam variables
  mU        = m3;
  mUS       = mU*mU;

  sHS = pow2(sH);
  tHS = pow2(tH);
  uHS = pow2(uH);
  tHC = pow(tH,3);
  uHC = pow(uH,3);
  tHQ = pow(tH,4);
  uHQ = pow(uH,4);
  tHuH = tH+uH;

  //+++ Evaluate (m**2, t, u) part of differential cross section.
  //+++ Extra 1/sHS comes from standard 2 to 2 cross section 
  //+++ phase space factors.

  if ( m_spin == 1 ) {
    
    double A0 = 1/sHS;
    double T1 = 0.5 * (tH/uH + uH/tH); 
    double T2 =  pow2(mZS + mUS)/(tH * uH); 
    double T3 = - 0.5 * mUS * (mZS/tHS + mZS/uHS) ;
    double T4 = - (mZS+mUS)*(1/tH + 1/uH);
    
    m_sigma0 = A0 * ( T1 + T2 + T3 + T4 );

  } else if ( m_spin == 2 ) {

    double A0   = 1 / ( sHS * uHS * tHS * pow2(sH-mZS) ); 
    double F0 = 2*tHS*uHS*( 16*pow(mZS,3) +  mUS*(7*tHS + 12*tH*uH + 7*uHS)
              - 3*(3*tHC + 11*tHS*uH + 11*tH*uHS + 3*uHC)
              + 6*pow(mZS,2)*(7*mUS - 2*tHuH) + mZS*(14*pow(mUS,2) 
              - 15*tHS - 44*tH*uH - 15*uHS + 2*mUS*tHuH) );
    double F2 = 2*tHS*uHS*tHuH*( -8*pow(mZS,2)*tHuH 
              + 4*mZS*(tHS + 3*tH*uH + uHS) 
              + 3*(tHC + 5*tHS*uH + 5*tH*uHS + uHC) );
    double F4 = -2*tHS*uHS*pow(tHuH,3)*(tHS + uHS - mZS*tHuH);

    double G0 = 4*tH*uH*( 6*pow(mZS,3)*(mUS - tH - uH)*tHuH
	      + pow(mZS,2)*( 9*tHC + 7*tHS*uH + 7*tH*uHS + 9*uHC 
              + 15*pow2(mUS)*tHuH - 2*mUS*(12*tHS + 19*tH*uH + 12*uHS) ) 
	      + tH*uH*( 6*pow(mUS,3) - 9*pow(mUS,2)*tHuH - mUS*(tHS 
              + 12*tH*uH + uHS) + 6*(tHC + 6*tHS*uH + 6*tH*uHS + uHC) ) 
	      + mZS*(-3*tHQ + 25*tHC*uH + 58*tHS*uHS + 25*tH*uHC 
              - 3*uHQ + 6*pow(mUS,3)*tHuH 
	      - pow(mUS,2)*(15*tHS + 2*tH*uH + 15*uHS) + 2*mUS*(6*tHC 
              - 11*tHS*uH - 11*tH*uHS + 6*uHC)) );
    double G2 = -4*tHS*uHS*tHuH*( -10*pow2(mZS)*tHuH + 2*mZS*(3*tHS 
              + 7*tH*uH + 3*uHS) + 3*(tHC + 5*tHS*uH + 5*tH*uHS + uHC) ); 
    double G4 = -2*F4;

    double H0 = 24*pow(mZS,3)*tH*uH*pow2(-mUS + tHuH) 
              - 6*pow(mZS,2)*tH*uH*( -9*pow(mUS,3) + 24*pow(mUS,2)*tHuH 
              - mUS*(21*tHS + 38*tH*uH + 21*uHS) 
              + 2*(3*tHC + 5*tHS*uH + 5*tH*uHS + 3*uHC) )
              - mZS*( 3*pow(mUS,4)*(tHS - 12*tH*uH + uHS) 
              - 2*tH*uH*pow2(tHuH)*(6*tHS - 29*tH*uH + 6*uHS) 
	      - 6*pow(mUS,3)*(tHC - 16*tHS*uH - 16*tH*uHS + uHC) 
              + 54*mUS*tH*uH*(tHC + tHS*uH + tH*uHS + uHC) 
	      + pow2(mUS)*(3*tHQ - 102*tHC*uH - 166*tHS*uHS 
              - 102*tH*uHC + 3*uHQ) )
              + tH*uH*( 6*pow(mUS,5) - 18*pow(mUS,4)*tHuH 
              - 12*pow(mUS,2)*pow(tHuH,3) 
              + 3*pow(mUS,3)*(7*tHS + 12*tH*uH + 7*uHS) 
	      - 18*tH*uH*(tHC + 5*tHS*uH + 5*tH*uHS + uHC) 
              + mUS*(3*tHQ + 32*tHC*uH + 78*tHS*uHS + 32*tH*uHC + 3*uHQ) );
    double H2 = 2*tHS*uHS*pow2(tHuH)*( -12*pow2(mZS) + 8*mZS*tHuH 
              + 3*(tHS + 4*tH*uH + uHS) );
    double H4 = F4;

    m_sigma0 = A0*( F0 + 1/mUS*F2 + 1/pow2(mUS)*F4 
	     + m_ratio*(G0 + 1/mUS*G2 + 1/pow2(mUS)*G4) 
	     + pow2(m_ratio)*(H0 + 1/mUS*H2 + 1/pow2(mUS)*H4) );

  } else {
    
    m_sigma0 = 0;
  
  }

}

//*********

double Sigma2ffbar2LEDUnparticleZ::sigmaHat() { 

  //+++ Electroweak couplings.
  int idAbs    = abs(id1);
  //+++ Note: 1/2 * (g_L^2 + g_R^2) = (g_v^2 + g_a^2) 
  double facEWS  = 4 * M_PI * alpEM  
                   / (CoupEW::sin2thetaW() * CoupEW::cos2thetaW()) 
                   * ( 0.25 * 0.25 * CoupEW::vf2af2(idAbs) );   

  //+++ Mass Spectrum, (m^2)^(d-2)
  double tmpExp = m_dU - 2;
  double facSpect = pow(mUS, tmpExp);

  //+++ Total cross section
  double sigma = m_constantTerm * facEWS * facSpect * m_sigma0;  
  
  //+++ If f fbar are quarks (1/N_c)
  if (idAbs < 9) sigma /= 3.;

  //+++ Related to mass spactrum weighting.
  sigma /= runBW3;   

  //+++ Truncation, to test perturbative region
  if(m_trunc) {
    if (sH > pow2(m_LambdaU) ) { sigma *= pow(m_LambdaU,4)/pow2(sH); }
  }

  return sigma;  

}

//*********

void Sigma2ffbar2LEDUnparticleZ::setIdColAcol() {

  //+++ Flavours trivial.
  setId( id1, id2, m_idG, 23);

  //+++ Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}
  
//*********

double Sigma2ffbar2LEDUnparticleZ::funcGammaIntHalf( int n ) {

  double gamma = 1;
  if (n%2 == 0) {
    int k = int(n/2);
    for (int i = 1; i <= k-1; ++i)  {
      gamma *= i;
    }  
  } else {
    int k = int((n-1)/2);
    for (int i = 0; i <= k-1; ++i)  {
      gamma *= i+0.5;
    } 
    gamma *= sqrt(M_PI);
  }  

  return gamma;
}

//*********

double Sigma2ffbar2LEDUnparticleZ::funcGammaReal( double x , int nmax ){

  if (x<0) {
    return 0;
  } else if (x == 0) {
    x = 0.000001;
    infoPtr->errorMsg("Warning in Sigma2ffbar2LEDUnparticleZ::funcGammaReal: "
		      "Zero argument, Gamma( x = 0.000001) is used");
  }

  double gamma = 1;
  double eulerGamma = 0.5772156649;

  for (int n=1; n < nmax+1; ++n){
    gamma *= exp(x/n)/(1 + x/n);
  }
  gamma *= exp(-eulerGamma * x)/x;

  return gamma;
}

//**************************************************************************

// Sigma2ffbar2LEDUnparticlegamma class.
// Cross section for f fbar -> U/G gamma (real LED graviton or unparticle 
// emission).

//*********

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// FIXRATIO:
// Ratio between the two possible coupling constants of the spin-2 ME. 
// A value different from one give rise to an IR divergence which makes 
// the event generation very slow, so this values is fixed to 1 until
// investigated further.
const double Sigma2ffbar2LEDUnparticlegamma::FIXRATIO = 1.;

//*********

Sigma2ffbar2LEDUnparticlegamma::Sigma2ffbar2LEDUnparticlegamma( 
  bool Graviton ): m_graviton(Graviton) {

}

//*********

void Sigma2ffbar2LEDUnparticlegamma::initProc() {

  //+++ WARNING: Keep in mind that this class uses the photon limit 
  //+++          of the Z+G/U ME code. This might give rise to some 
  //+++          confusing things, e.g. mZ = ParticleDataTable::m0(22);          
  
  //+++ Init model parameters.
  m_idG        = 5000039;
  if (m_graviton) {
    m_trunc    = Settings::flag("ExtraDimensionsLED:Trunc"); 
    m_nGrav    = Settings::mode("ExtraDimensionsLED:n");
    m_LambdaU  = Settings::parm("ExtraDimensionsLED:MD");
    m_dU       = m_nGrav/2 + 1; 
    m_spin     = 2;
  } else {
    m_trunc    = Settings::flag("ExtraDimensionsUnpart:Trunc"); 
    m_spin     = Settings::mode("ExtraDimensionsUnpart:spinU");
    m_dU       = Settings::parm("ExtraDimensionsUnpart:dU");
    m_LambdaU  = Settings::parm("ExtraDimensionsUnpart:LambdaU");
    m_lambda   = Settings::parm("ExtraDimensionsUnpart:lambda");
    m_ratio    = FIXRATIO; 
    //         = Settings::parm("ExtraDimensionsUnpart:ratio");
  }

  //+++ Store Z0 mass.
  mZ        = ParticleDataTable::m0(22);
  mZS       = mZ*mZ;  

  //+++ Init spin-2 parameters
  if ( m_spin != 2 ){
    m_graviton = false;
    m_lambdaPrime = 0;
  } else if (m_graviton) {
    m_lambda = 1;
    m_ratio = 1;
    m_lambdaPrime = m_lambda;
  } else {
    m_lambdaPrime = m_ratio * m_lambda;
  }

  //+++ The A(dU) or S'(n) value
  double tmpArg1   = m_dU + 0.5;
  double tmpGamma1 = funcGammaReal( tmpArg1, 100000);
  double tmpArg2   = m_dU - 1;
  double tmpGamma2 = funcGammaReal( tmpArg2, 100000);  
  double tmpArg3   = 2 * m_dU;
  double tmpGamma3 = funcGammaReal( tmpArg3, 100000);
  double tmpBase   = 2 * M_PI;
  double tmpExp1   = 2 * m_dU;
  double tmpAdU    = 16 * pow2(M_PI) * sqrt(M_PI) / pow(tmpBase, tmpExp1)
                   * tmpGamma1 / (tmpGamma2 * tmpGamma3);
  if (m_graviton) { 
    m_nGrav = 2*(m_dU - 1);
    tmpAdU  = 2 * M_PI * sqrt( pow(M_PI, m_nGrav) ) 
            / funcGammaIntHalf( int(m_nGrav) ); 
  } 

  //+++ Standard 2 to 2 cross section related constants
  double tmpTerm1 = 1/(2 * 16 * pow2(M_PI));
  double tmpLS    = pow2(m_LambdaU);

  //+++ Spin dependent constants from ME.
  //+++ Spin-0 is only place holder.
  double tmpTerm2 = 0;
  if ( m_spin == 0 ) {
    
  } else if (m_spin == 1) {
    tmpTerm2 = 4 * pow2(m_lambda);
  } else if (m_spin == 2) {
    tmpTerm2 = pow2(m_lambda)/(4 * 3 * tmpLS);
  } 

  //+++ Unparticle phase space related
  double tmpExp2 = m_dU - 2;
  double tmpTerm3 = tmpAdU / (tmpLS * pow(tmpLS, tmpExp2));

  //+++ All in total
  m_constantTerm = tmpTerm1 * tmpTerm2 * tmpTerm3;

} 

//*********

void Sigma2ffbar2LEDUnparticlegamma::sigmaKin() { 

  //+++ Set graviton mass and some powers of mandelstam variables
  mU        = m3;
  mUS       = mU*mU;

  sHS = pow2(sH);
  tHS = pow2(tH);
  uHS = pow2(uH);
  tHC = pow(tH,3);
  uHC = pow(uH,3);
  tHQ = pow(tH,4);
  uHQ = pow(uH,4);
  tHuH = tH+uH;

  //+++ Evaluate (m**2, t, u) part of differential cross section.
  //+++ Extra 1/sHS comes from standard 2 to 2 cross section 
  //+++ phase space factors.

  if ( m_spin == 1 ) {
    
    double A0 = 1/sHS;
    double T1 = 0.5 * (tH/uH + uH/tH); 
    double T2 =  pow2(mZS + mUS)/(tH * uH); 
    double T3 = - 0.5 * mUS * (mZS/tHS + mZS/uHS) ;
    double T4 = - (mZS+mUS)*(1/tH + 1/uH);
    
    m_sigma0 = A0 * ( T1 + T2 + T3 + T4 );

  } else if ( m_spin == 2 ) {

    double A0 = 1 / ( sHS * uHS * tHS * pow2(sH-mZS) ); 
    double F0 = 2*tHS*uHS*( 16*pow(mZS,3) +  mUS*(7*tHS + 12*tH*uH + 7*uHS) 
              - 3*(3*tHC + 11*tHS*uH + 11*tH*uHS + 3*uHC)
              + 6*pow(mZS,2)*(7*mUS - 2*tHuH) + mZS*(14*pow(mUS,2) 
              - 15*tHS - 44*tH*uH - 15*uHS + 2*mUS*tHuH) );
    double F2 = 2*tHS*uHS*tHuH*( -8*pow(mZS,2)*tHuH 
              + 4*mZS*(tHS + 3*tH*uH + uHS) 
              + 3*(tHC + 5*tHS*uH + 5*tH*uHS + uHC) );
    double F4 = -2*tHS*uHS*pow(tHuH,3)*(tHS + uHS - mZS*tHuH);

    double G0 = 4*tH*uH*( 6*pow(mZS,3)*(mUS - tH - uH)*tHuH
	      + pow(mZS,2)*( 9*tHC + 7*tHS*uH + 7*tH*uHS + 9*uHC 
              + 15*pow2(mUS)*tHuH - 2*mUS*(12*tHS + 19*tH*uH + 12*uHS) ) 
	      + tH*uH*( 6*pow(mUS,3) - 9*pow(mUS,2)*tHuH 
              - mUS*(tHS + 12*tH*uH + uHS) 
              + 6*(tHC + 6*tHS*uH + 6*tH*uHS + uHC) ) 
	      + mZS*(-3*tHQ + 25*tHC*uH + 58*tHS*uHS + 25*tH*uHC 
              - 3*uHQ + 6*pow(mUS,3)*tHuH 
	      - pow(mUS,2)*(15*tHS + 2*tH*uH + 15*uHS) 
              + 2*mUS*(6*tHC - 11*tHS*uH - 11*tH*uHS + 6*uHC)) );
    double G2 = -4*tHS*uHS*tHuH*( -10*pow2(mZS)*tHuH 
              + 2*mZS*(3*tHS + 7*tH*uH + 3*uHS) 
              + 3*(tHC + 5*tHS*uH + 5*tH*uHS + uHC) ); 
    double G4 = -2*F4;

    double H0 = 24*pow(mZS,3)*tH*uH*pow2(-mUS + tHuH) 
              - 6*pow(mZS,2)*tH*uH*( -9*pow(mUS,3) + 24*pow(mUS,2)*tHuH 
              - mUS*(21*tHS + 38*tH*uH + 21*uHS) 
              + 2*(3*tHC + 5*tHS*uH + 5*tH*uHS + 3*uHC) )
              - mZS*( 3*pow(mUS,4)*(tHS - 12*tH*uH + uHS) 
              - 2*tH*uH*pow2(tHuH)*(6*tHS - 29*tH*uH + 6*uHS) 
	      - 6*pow(mUS,3)*(tHC - 16*tHS*uH - 16*tH*uHS + uHC) 
              + 54*mUS*tH*uH*(tHC + tHS*uH + tH*uHS + uHC) 
	      + pow2(mUS)*(3*tHQ - 102*tHC*uH - 166*tHS*uHS 
              - 102*tH*uHC + 3*uHQ) )
              + tH*uH*( 6*pow(mUS,5) - 18*pow(mUS,4)*tHuH 
              - 12*pow(mUS,2)*pow(tHuH,3) 
              + 3*pow(mUS,3)*(7*tHS + 12*tH*uH + 7*uHS) 
	      - 18*tH*uH*(tHC + 5*tHS*uH + 5*tH*uHS + uHC) 
              + mUS*(3*tHQ + 32*tHC*uH + 78*tHS*uHS + 32*tH*uHC + 3*uHQ) );
    double H2 = 2*tHS*uHS*pow2(tHuH)*( -12*pow2(mZS) + 8*mZS*tHuH 
              + 3*(tHS + 4*tH*uH + uHS) );
    double H4 = F4;

    m_sigma0 = A0*( F0 + 1/mUS*F2 + 1/pow2(mUS)*F4 
	     + m_ratio*(G0 + 1/mUS*G2 + 1/pow2(mUS)*G4) 
	     + pow2(m_ratio)*(H0 + 1/mUS*H2 + 1/pow2(mUS)*H4) );

  } else {
    
    m_sigma0 = 0;
  
  }

}

//*********

double Sigma2ffbar2LEDUnparticlegamma::sigmaHat() { 

  //+++ Electroweak couplings..
  int idAbs    = abs(id1);
  double facEWS = 4 * M_PI * alpEM * CoupEW::ef2(idAbs);

  //+++ Mass Spectrum, (m^2)^(d-2)
  double tmpExp = m_dU - 2;
  double facSpect = pow(mUS, tmpExp);

  //+++ Total cross section
  double sigma = m_constantTerm * facEWS * facSpect * m_sigma0;  

  //+++ If f fbar are quarks
  if (idAbs < 9) sigma /= 3.;

  //+++ Related to mass spactrum weighting.
  sigma /= runBW3;      

  //+++ Truncation, to test perturbative region
  if (m_trunc) {
    if (sH > pow2(m_LambdaU) ) { sigma *= pow(m_LambdaU,4)/pow2(sH); }
  }
  return sigma;  

}

//*********

void Sigma2ffbar2LEDUnparticlegamma::setIdColAcol() {

  //+++ Flavours trivial.
  setId( id1, id2, m_idG, 22);

  //+++ Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//*********

double Sigma2ffbar2LEDUnparticlegamma::funcGammaIntHalf( int n ) {

  double gamma = 1;
  if (n%2 == 0) {
    int k = int(n/2);
    for (int i = 1; i <= k-1; ++i)  {
      gamma *= i;
    }  
  } else {
    int k = int((n-1)/2);
    for (int i = 0; i <= k-1; ++i)  {
      gamma *= i+0.5;
    } 
    gamma *= sqrt(M_PI);
  }  

  return gamma;
}

//*********

double Sigma2ffbar2LEDUnparticlegamma::funcGammaReal( double x , int nmax ){

  if (x<0) {
    return 0;
  } else if (x == 0) {
    x = 0.000001;
    infoPtr->errorMsg("Warning in Sigma2ffbar2LEDUnparticlegamma::"
      "funcGammaReal: Zero argument, Gamma( x = 0.000001) is used");
  }

  double gamma = 1;
  double eulerGamma = 0.5772156649;

  for (int n=1; n < nmax+1; ++n){
    gamma *= exp(x/n)/(1 + x/n);
  }
  gamma *= exp(-eulerGamma * x)/x;

  return gamma;
}

//**************************************************************************

// Sigma2ffbar2LEDgammagamma class.
// Cross section for f fbar -> (LED G*/U*) -> gamma gamma 
// (virtual graviton/unparticle exchange).

//*********

Sigma2ffbar2LEDgammagamma::Sigma2ffbar2LEDgammagamma( bool Graviton ) 
  : m_graviton(Graviton) {
  
}

//*********

void Sigma2ffbar2LEDgammagamma::initProc() {
  
  //+++ Init model parameters.
  if (m_graviton) {
    m_spin     = 2;
    m_dU       = 2;
    m_LambdaU  = Settings::parm("ExtraDimensionsLED:LambdaT");
    m_lambda   = 1;
  } else {
    m_spin     = Settings::mode("ExtraDimensionsUnpart:spinU");
    m_dU       = Settings::parm("ExtraDimensionsUnpart:dU");
    m_LambdaU  = Settings::parm("ExtraDimensionsUnpart:LambdaU");
    m_lambda   = Settings::parm("ExtraDimensionsUnpart:lambda");
  }

  //+++ Model dependent constants.
  if (m_graviton) {
    m_lambda2chi = 4*M_PI;
  } else {
    double tmp_arg1   = m_dU + 0.5;
    double tmp_Gamma1 = funcGammaReal( tmp_arg1, 100000);
    double tmp_arg2   = m_dU - 1;
    double tmp_Gamma2 = funcGammaReal( tmp_arg2, 100000);  
    double tmp_arg3   = 2 * m_dU;
    double tmp_Gamma3 = funcGammaReal( tmp_arg3, 100000);
    double tmp_base   = 2 * M_PI;
    double tmp_exp1   = 2 * m_dU;
    double tmp_AdU    = 16 * pow2(M_PI) * sqrt(M_PI) 
                      / pow(tmp_base, tmp_exp1)
                      * tmp_Gamma1 / (tmp_Gamma2 * tmp_Gamma3);
    double tmp_dUpi = m_dU * M_PI;
    m_lambda2chi = pow2(m_lambda) * tmp_AdU / (2 * sin(tmp_dUpi));
  }

  //+++ Model parameter check (if not applicable, sigma = 0).
  //+++ Note: SM contribution still generated.
  if ( !(m_spin==0 || m_spin==2) ) {
    m_lambda2chi = 0;
    infoPtr->errorMsg("Error in Sigma2ffbar2LEDgammagamma::initProc: "
		      "Incorrect spin value (turn process off)!");
  } else if ( !m_graviton && (m_dU >= 2)) {
    m_lambda2chi = 0;
    infoPtr->errorMsg("Error in Sigma2ffbar2LEDgammagamma::initProc: "
		      "This process requires dU < 2 (turn process off)!");
  }

} 

//*********

void Sigma2ffbar2LEDgammagamma::sigmaKin() { 

  //+++ Mandelstam variables.
  double sHS = pow2(sH);
  double sHQ = pow(sH, 4);
  double tHS = pow2(tH);
  double uHS = pow2(uH);

  //+++ ME from spin-0 and spin-2 unparticles
  //+++ including extra 1/sHS from 2-to-2 phase space.
  if (m_spin == 0) {
    double tmp_sLambda2 = sH / pow2(m_LambdaU);
    double tmp_exp = 2 * m_dU - 1;
    m_term1 = pow(tmp_sLambda2,tmp_exp);
    m_term1 /= sHS;
  } else {
    m_term1 = (uH / tH + tH / uH);
    m_term1 /= sHS;
    double tmp_sLambda2 = sH / pow2(m_LambdaU);
    double tmp_exp = m_dU;
    m_term2 = pow(tmp_sLambda2,tmp_exp) * (uHS + tHS) / sHS;
    m_term2 /= sHS;
    tmp_exp = 2 * m_dU;
    m_term3 = pow(tmp_sLambda2,tmp_exp) * tH * uH * (uHS + tHS) / sHQ;
    m_term3 /= sHS;
  }

}

//*********

double Sigma2ffbar2LEDgammagamma::sigmaHat() { 

  //+++ Incoming fermion flavor.
  int idAbs      = abs(id1);

  //+++ Couplings and constants.
  //+++ Note: ME already contain 1/2 for identical 
  //+++       particles in the final state.
  double sigma = 0;
  if (m_spin == 0) {
    sigma = pow2(m_lambda2chi) * m_term1 / 8;
  } else {
    double tmp_e2Q2 = 4 * M_PI * alpEM * CoupEW::ef2(idAbs);
    double tmp_dUpi = m_dU * M_PI;
    sigma = pow2(tmp_e2Q2) * m_term1
          - tmp_e2Q2 * m_lambda2chi * cos(tmp_dUpi) * m_term2
          + pow2(m_lambda2chi) * m_term3 / 4;
  } 

  //+++ dsigma/dt, 2-to-2 phase space factors.
  sigma /= 16 * M_PI;

  //+++ If f fbar are quarks.
  if (idAbs < 9) sigma /= 3.;

  return sigma;  
}

//*********

void Sigma2ffbar2LEDgammagamma::setIdColAcol() {

  //+++ Flavours trivial.
  setId( id1, id2, 22, 22);

  //+++ Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//*********

double Sigma2ffbar2LEDgammagamma::funcGammaReal( double x , int nmax ){

  if (x<0) {
    return 0;
  } else if (x == 0) {
    x = 0.000001;
    infoPtr->errorMsg("Warning in Sigma2ffbar2LEDgammagamma::"
      "funcGammaReal: Zero argument, Gamma( x = 0.000001) is used");
  }

  double gamma = 1;
  double eulerGamma = 0.5772156649;

  for (int n=1; n < nmax+1; ++n){
    gamma *= exp(x/n)/(1 + x/n);
  }
  gamma *= exp(-eulerGamma * x)/x;

  return gamma;
}


//**************************************************************************

// Sigma2gg2LEDgammagamma class.
// Cross section for g g -> (LED G*/U*) -> gamma gamma 
// (virtual graviton/unparticle exchange).

//*********

Sigma2gg2LEDgammagamma::Sigma2gg2LEDgammagamma( bool Graviton ) 
  : m_graviton(Graviton) {

}

//*********

void Sigma2gg2LEDgammagamma::initProc() {

  //+++ Init model parameters.
  if (m_graviton) {
    m_spin     = 2;
    m_dU       = 2;
    m_LambdaU  = Settings::parm("ExtraDimensionsLED:LambdaT");
    m_lambda   = 1;
  } else {
    m_spin     = Settings::mode("ExtraDimensionsUnpart:spinU");
    m_dU       = Settings::parm("ExtraDimensionsUnpart:dU");
    m_LambdaU  = Settings::parm("ExtraDimensionsUnpart:LambdaU");
    m_lambda   = Settings::parm("ExtraDimensionsUnpart:lambda");
  }

  //+++ Model dependent constants.
  if (m_graviton) {
    m_lambda2chi = 4 * M_PI;

  } else {
    double tmp_arg1   = m_dU + 0.5;
    double tmp_Gamma1 = funcGammaReal( tmp_arg1, 100000);
    double tmp_arg2   = m_dU - 1;
    double tmp_Gamma2 = funcGammaReal( tmp_arg2, 100000);  
    double tmp_arg3   = 2 * m_dU;
    double tmp_Gamma3 = funcGammaReal( tmp_arg3, 100000);
    double tmp_base   = 2 * M_PI;
    double tmp_exp1   = 2 * m_dU;
    double tmp_AdU    = 16 * pow2(M_PI) * sqrt(M_PI) 
                      / pow(tmp_base, tmp_exp1)
                      * tmp_Gamma1 / (tmp_Gamma2 * tmp_Gamma3);
    double tmp_dUpi = m_dU * M_PI;
    m_lambda2chi = pow2(m_lambda) * tmp_AdU / (2 * sin(tmp_dUpi));
  }

  //+++ Model parameter check (if not applicable, sigma = 0).
  if ( !(m_spin==0 || m_spin==2) ) {
    m_lambda2chi = 0;
    infoPtr->errorMsg("Error in Sigma2gg2LEDgammagamma::initProc: "
		      "Incorrect spin value (turn process off)!");
  } else if ( !m_graviton && (m_dU >= 2)) {
    m_lambda2chi = 0;
    infoPtr->errorMsg("Error in Sigma2gg2LEDgammagamma::initProc: "
		      "This process requires dU < 2 (turn process off)!");
  }

} 

//*********

void Sigma2gg2LEDgammagamma::sigmaKin() { 
  
  //+++ Mandelstam variables.
  double sHS = pow2(sH);
  double sHQ = pow(sH, 4);
  double tHQ = pow(tH, 4);
  double uHQ = pow(uH, 4);

  //+++ ME from spin-0 and spin-2 unparticles.
  if (m_spin == 0) {
    double tmp_sLambda2 = sH / pow2(m_LambdaU);
    double tmp_exp = 2 * m_dU;
    m_sigma0 = pow(tmp_sLambda2,tmp_exp);
  } else {
    double tmp_sLambda2 = sH / pow2(m_LambdaU);
    double tmp_exp = 2 * m_dU;
    m_sigma0 = pow(tmp_sLambda2,tmp_exp) * (uHQ + tHQ) / sHQ;
  }

  //+++ extra 1/sHS from 2-to-2 phase space.
  m_sigma0 /= sHS;

}

//*********

double Sigma2gg2LEDgammagamma::sigmaHat() { 

  //+++ Couplings and constants.
  //+++ Note: ME already contain 1/2 for identical 
  //+++       particles in the final state.
  double sigma = m_sigma0;
  if (m_spin == 0) {
    sigma *= pow2(m_lambda2chi) / 256;
  } else {
    sigma *= pow2(m_lambda2chi) / 32;
  } 

  //+++ dsigma/dt, 2-to-2 phase space factors.
  sigma /= 16 * M_PI;

  return sigma;  
}

//*********

void Sigma2gg2LEDgammagamma::setIdColAcol() {

  //+++ Flavours trivial.
  setId( 21, 21, 22, 22);

  //+++ Colour flow topologies. 
  setColAcol( 1, 2, 2, 1, 0, 0, 0, 0);

}

//*********

double Sigma2gg2LEDgammagamma::funcGammaReal( double x , int nmax ){

  if (x<0) {
    return 0;
  } else if (x == 0) {
    x = 0.000001;
    infoPtr->errorMsg("Warning in Sigma2gg2LEDgammagamma::"
      "funcGammaReal: Zero argument, Gamma( x = 0.000001) is used");
  }

  double gamma = 1;
  double eulerGamma = 0.5772156649;

  for (int n=1; n < nmax+1; ++n){
    gamma *= exp(x/n)/(1 + x/n);
  }
  gamma *= exp(-eulerGamma * x)/x;

  return gamma;
}

//**************************************************************************

} // end namespace Pythia8
