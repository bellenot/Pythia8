// Function definitions (not found in the header) for the 
// Higgs simulation classes. 
// Copyright C 2007 Torbjorn Sjostrand

#include "SigmaHiggs.h"

namespace Pythia8 {

//**************************************************************************

// Sigma1ffbar2H class.
// Cross section for f fbar -> H0 (f is quark or lepton, H0 SM Higgs). 

//*********

// Initialize process. 
  
void Sigma1ffbar2H::initProc() {

  // Store H0 mass and width for propagator. 
  mRes     = ParticleDataTable::m0(25);
  GammaRes = ParticleDataTable::mWidth(25);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

} 

//*********

// Initialize parton-flux object. 
  
void Sigma1ffbar2H::initFlux() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxffbarSame();

  // Multiply by colour factor 1/3.
  inFluxPtr->weightInvCol();

} 

//*********

// Evaluate sigmaHat(sHat). 
// Note: owing to the fermion-mass dependence,
// the incoming partial width is put in inFlux. 

double Sigma1ffbar2H::sigmaHat() { 

  // Initial values.
  double mHat  = sqrt(sH);
  double widthIn;

  // Loop over incoming fermion flavour species.
  for (int i = 0; i < 5; ++i) {
    int idAbs = i + 1;
    if (hasLeptonBeams && i == 0) idAbs = abs(idA);
    else if (hasLeptonBeams) continue;

    // Calculate and store incoming width for each fermion species.
    widthIn = HRes.widthIn( mHat, idAbs);
    inFluxPtr->weightInState( idAbs, -idAbs, widthIn, true, true, false);
  }

  // Set up Breit-Wigner. Width out only includes open channels. 
  double sigBW    = 4. * M_PI / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );    
  double widthOut = HRes.width( mHat, true);

  // Done.
  return CONVERT2MB * sigBW * widthOut;    

}

//*********

// Select identity, colour and anticolour.

void Sigma1ffbar2H::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, 25);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 10) setColAcol( 1, 0, 0, 1, 0, 0);
  else               setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//*********

// Evaluate weight for Z0 Z0 or W+W- decay angles in Higgs decay.

double Sigma1ffbar2H::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // For Higgs decay hand over to standard routine, else done.
  if (process[process[iResBeg].mother1()].idAbs() == 25) 
       return weightHiggsDecay( process, iResBeg, iResEnd);
  else return 1.; 

}

//**************************************************************************

// Sigma1gg2H class.
// Cross section for g g -> H0 (H0 SM Higgs). 

//*********

// Initialize process. 
  
void Sigma1gg2H::initProc() {

  // Store H0 mass and width for propagator. 
  mRes     = ParticleDataTable::m0(25);
  GammaRes = ParticleDataTable::mWidth(25);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

} 

//*********

// Initialize parton-flux object. 
  
void Sigma1gg2H::initFlux() {

  // Set up for g g initial state.
  inFluxPtr = new InFluxgg();

  // Multiply by colour factor 1/8.
  inFluxPtr->weightInvCol();

} 

//*********

// Evaluate sigmaHat(sHat). 

double Sigma1gg2H::sigmaHat() { 

  // Initial values. Incoming width for gluons excludes factor of 8.
  double mHat     = sqrt(sH);
  double widthIn  = HRes.widthIn( mHat, 21);

  // Set up Breit-Wigner. Width out only includes open channels. 
  double sigBW    = 8. * M_PI/ ( pow2(sH - m2Res) + pow2(sH * GamMRat) );    
  double widthOut = HRes.width( mHat, true);

  // Done.
  return CONVERT2MB * widthIn * sigBW * widthOut;    

}

//*********

// Select identity, colour and anticolour.

void Sigma1gg2H::setIdColAcol() {

  // Flavours trivial.
  setId( 21, 21, 25);

  // Colour flow topology.
  setColAcol( 1, 2, 2, 1, 0, 0);

}

//*********

// Evaluate weight for Z0 Z0 or W+W- decay angles in Higgs decay.

double Sigma1gg2H::weightDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // For Higgs decay hand over to standard routine, else done.
  if (process[process[iResBeg].mother1()].idAbs() == 25) 
       return weightHiggsDecay( process, iResBeg, iResEnd);
  else return 1.; 

}

//**************************************************************************

// Sigma1gmgm2H class.
// Cross section for gamma gamma -> H0 (H0 SM Higgs). 

//*********

// Initialize process. 
  
void Sigma1gmgm2H::initProc() {

  // Store H0 mass and width for propagator. 
  mRes     = ParticleDataTable::m0(25);
  GammaRes = ParticleDataTable::mWidth(25);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

} 

//*********

// Initialize parton-flux object. 
  
void Sigma1gmgm2H::initFlux() {

  // Set up for gamma gamma initial state.
  inFluxPtr = new InFluxgmgm();

} 

//*********

// Evaluate sigmaHat(sHat). 

double Sigma1gmgm2H::sigmaHat() { 

  // Initial values. Incoming width for gluons excludes factor of 8.
  double mHat     = sqrt(sH);
  double widthIn  = HRes.widthIn( mHat, 22);

  // Set up Breit-Wigner. Width out only includes open channels. 
  double sigBW    = 8. * M_PI/ ( pow2(sH - m2Res) + pow2(sH * GamMRat) );    
  double widthOut = HRes.width( mHat, true);

  // Done.
  return CONVERT2MB * widthIn * sigBW * widthOut;    

}

//*********

// Select identity, colour and anticolour.

void Sigma1gmgm2H::setIdColAcol() {

  // Flavours trivial.
  setId( 22, 22, 25);

  // Colour flow trivial.
  setColAcol( 0, 0, 0, 0, 0, 0);

}

//*********

// Evaluate weight for Z0 Z0 or W+W- decay angles in Higgs decay.

double Sigma1gmgm2H::weightDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // For Higgs decay hand over to standard routine, else done.
  if (process[process[iResBeg].mother1()].idAbs() == 25) 
       return weightHiggsDecay( process, iResBeg, iResEnd);
  else return 1.; 

}

//**************************************************************************

// Sigma2ffbar2HZ class.
// Cross section for f fbar -> H0 Z0 (H0 SM Higgs). 

//*********

// Initialize process. 
  
void Sigma2ffbar2HZ::initProc() {

  // Store Z0 mass and width for propagator. Common coupling factor.
  mZ        = ParticleDataTable::m0(23);
  widZ      = ParticleDataTable::mWidth(23);
  mZS       = mZ*mZ;
  mwZS      = pow2(mZ * widZ);
  thetaWRat = 1. / (16. * CoupEW::sin2thetaW() * CoupEW::cos2thetaW());  

} 

//*********

// Initialize parton-flux object. 
  
void Sigma2ffbar2HZ::initFlux() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxffbarSame();

  // Multiply by coupling a_f^2 + v_f^2 to s-channel Z0.
  if (hasLeptonBeams) {
    inFluxPtr->weightFixed( idA, -idA, CoupEW::vf2af2(abs(idA)) );
  } else {    
    inFluxPtr->weightFixed( 1, -1, CoupEW::vf2af2(1) );
    inFluxPtr->weightFixed( 2, -2, CoupEW::vf2af2(2) );
  } 

  // Multiply by colour factor 1/3.
  inFluxPtr->weightInvCol();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2ffbar2HZ::sigmaHat() { 

  // Evaluate differential cross section.
  double sigma = (M_PI / sH2) * 8. * pow2(alpEM * thetaWRat)
    * (tH * uH - s3 * s4 + 2. * sH * s4) / (pow2(sH - mZS) + mwZS);

  // Answer, leaving out flavour-dependent pieces stored above.
  return CONVERT2MB * sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2ffbar2HZ::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, 25, 23);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 10) setColAcol( 1, 0, 0, 1, 0, 0);
  else               setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//*********

// Evaluate weight for decay angles.

double Sigma2ffbar2HZ::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Hand over Higgs decay to Z0 Z0 or W+ W- to 4 fermions.
  if (process[process[iResBeg].mother1()].idAbs() == 25) 
       return weightHiggsDecay( process, iResBeg, iResEnd);

  // Handle decay of Z0 created along with Higgs.
  if (iResBeg != 5 || iResEnd != 7) return 1.;

  // Order so that fbar(1) f(2) -> H() f'(3) fbar'(4). 
  int i1       = (process[3].id() < 0) ? 3 : 4;
  int i2       = 7 - i1; 
  int i3       = process[6].daughter1();
  int i4       = process[6].daughter2();
  if (process[i3].id() < 0) swap( i3, i4); 

  // Find left- and righthanded couplings of fermion pairs.
  int    idAbs = process[i1].idAbs(); 
  double liS   = pow2( CoupEW::lf(idAbs) );
  double riS   = pow2( CoupEW::rf(idAbs) );
  idAbs        = process[i3].idAbs(); 
  double lfS   = pow2( CoupEW::lf(idAbs) );
  double rfS   = pow2( CoupEW::rf(idAbs) );

  // Evaluate relevant four-products.
  double pp13  = process[i1].p() * process[i3].p();
  double pp14  = process[i1].p() * process[i4].p();
  double pp23  = process[i2].p() * process[i3].p();
  double pp24  = process[i2].p() * process[i4].p();

  // Weight and maximum.
  double wt    = (liS * lfS + riS * rfS) * pp13 * pp24
               + (liS * rfS + riS * lfS) * pp14 * pp23;
  double wtMax = (liS + riS) * (lfS + rfS) * (pp13 + pp14) * (pp23 + pp24);

  // Done.
  return wt / wtMax;

}

//**************************************************************************

// Sigma2ffbar2HW class.
// Cross section for f fbar -> H0 W+- (H0 SM Higgs). 

//*********

// Initialize process. 
  
void Sigma2ffbar2HW::initProc() {

  // Store W+- mass and width for propagator. Common coupling factor.
  mW   = ParticleDataTable::m0(24);
  widW = ParticleDataTable::mWidth(24);
  mWS  = mW*mW;
  mwWS = pow2(mW * widW);
  thetaWRat = 1. / (4. * CoupEW::sin2thetaW());  

} 

//*********

// Initialize parton-flux object. 
  
void Sigma2ffbar2HW::initFlux() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxffbarChg();

  // Multiply by squared CKM matrix elements and colour factor 1/3.
  inFluxPtr->weightCKM2();
  inFluxPtr->weightInvCol();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2ffbar2HW::sigmaHat() { 

  //  Evaluate differential cross section.
  double sigma = (M_PI / sH2) * 2. * pow2(alpEM * thetaWRat)
    * (tH * uH - s3 * s4 + 2. * sH * s4) / (pow2(sH - mWS) + mwWS);

  // Answer, leaving out flavour-dependent pieces stored above.
  return CONVERT2MB * sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2ffbar2HW::setIdColAcol() {

  // Sign of outgoing W. 
  int sign = 1 - 2 * (abs(id1)%2);
  if (id1 < 0) sign = -sign;
  setId( id1, id2, 25, 24 * sign);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 10) setColAcol( 1, 0, 0, 1, 0, 0);
  else               setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//*********

// Evaluate weight for decay angles.

double Sigma2ffbar2HW::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // Hand over Higgs decay to Z0 Z0 or W+ W- to 4 fermions.
  if (process[process[iResBeg].mother1()].idAbs() == 25) 
       return weightHiggsDecay( process, iResBeg, iResEnd);

  // Handle decay of W+- created along with Higgs.
  if (iResBeg != 5 || iResEnd != 7) return 1.;

  // Order so that fbar(1) f(2) -> H() f'(3) fbar'(4). 
  int i1       = (process[3].id() < 0) ? 3 : 4;
  int i2       = 7 - i1; 
  int i3       = process[6].daughter1();
  int i4       = process[6].daughter2();
  if (process[i3].id() < 0) swap( i3, i4); 

  // Evaluate relevant four-products.
  double pp13  = process[i1].p() * process[i3].p();
  double pp14  = process[i1].p() * process[i4].p();
  double pp23  = process[i2].p() * process[i3].p();
  double pp24  = process[i2].p() * process[i4].p();

  // Weight and maximum.
  double wt    = pp13 * pp24;
  double wtMax = (pp13 + pp14) * (pp23 + pp24);

  // Done.
  return wt / wtMax;

}


//**************************************************************************

// Sigma2qg2Hq class.
// Cross section for q g -> H0 q (H0 SM Higgs). 

//*********

// Initialize process. 
  
void Sigma2qg2Hq::initProc() {

  m2W       = pow2( ParticleDataTable::m0(24) );
  thetaWRat = 1. / (24. * CoupEW::sin2thetaW());  
  
} 

//*********

// Initialize parton-flux object. 
  
void Sigma2qg2Hq::initFlux() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxqg();

} 

//*********

// Evaluate sigmaHat(sHat). 
// Note: owing to the fermion-mass dependence,
// the incoming partial width is put in inFlux. 

double Sigma2qg2Hq::sigmaHat() { 

  // Initial values.
  double mHat  = sqrt(sH);

  // Common part of cross section.
  double sigCom = (M_PI / sH2) * alpS * alpEM * thetaWRat;

  // Loop over incoming quark species. Fix and running mass.
  for (int idAbs = 1; idAbs < 6; ++idAbs) {
    double m2Fix = pow2( ParticleDataTable::m0(idAbs) );
    double m2Run = pow2( ParticleDataTable::mRun(idAbs, mHat) );

    // Flavour-specific part of cross section.. 
    double sigFlav = (m2Run/m2W) * ( sH / (m2Fix - uH) 
      + 2. * m2Fix * (s3 - uH) / pow2(m2Fix - uH) 
      + (m2Fix - uH) / sH - 2. * m2Fix / (m2Fix - uH) 
      + 2. * (s3 - uH)  * (s3 - m2Fix - sH) / ((m2Fix - uH) * sH) );
    
    // Store flavour-specific part.
    inFluxPtr->weightInState( idAbs, 21, sigFlav, true, true, false);
  }

  // Done.
  return CONVERT2MB * sigCom;  

}

//*********

// Select identity, colour and anticolour.

void Sigma2qg2Hq::setIdColAcol() {

  // Flavour set up for q g -> H0 q.
  int idq = (id2 == 21) ? id1 : id2;
  setId( id1, id2, 25, idq);

  // tH defined between f and f': must swap tHat <-> uHat if q g in.
  swapTU = (id2 == 21); 

  // Colour flow topologies. Swap when antiquarks.
  if (id2 == 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else           setColAcol( 2, 1, 1, 0, 0, 0, 2, 0);
  if (idq < 0) swapColAcol();

}

//*********

// Evaluate weight for Z0 Z0 or W+W- decay angles in Higgs decay.

double Sigma2qg2Hq::weightDecay( Event& process, int iResBeg,
  int iResEnd) {

  // For Higgs decay hand over to standard routine, else done.
  if (process[process[iResBeg].mother1()].idAbs() == 25) 
       return weightHiggsDecay( process, iResBeg, iResEnd);
  else return 1.; 

}

//**************************************************************************

} // end namespace Pythia8
