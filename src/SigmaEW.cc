// Function definitions (not found in the header) for the 
// electroweak simulation classes (including photon processes). 
// Copyright C 2006 Torbjorn Sjostrand

#include "SigmaEW.h"

namespace Pythia8 {

//**************************************************************************

// Sigma2qg2qgamma class.
// Cross section for q g -> q gamma.

//*********

// Initialize process, especially parton-flux object. 
  
void Sigma2qg2qgamma::initProc() {

  // Set up for q qbar initial state, e2-weighted.
  inFluxPtr = new InFluxqg();
  inFluxPtr->weightCharge2();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qg2qgamma::sigmaHat() {  

  // Calculate kinematics dependence.
  sigUS = (1./3.) * (sH2+uH2)/(-sH*uH);

  // Answer.
  return CONVERT2MB * (M_PI/sH2) * alpS * alpEM * sigUS;  

}

//*********

// Select identity, colour and anticolour.

void Sigma2qg2qgamma::setIdColAcol() {

  // Construct outgoing flavours.
  id3 = (id1 == 21) ? 22 : id1;
  id4 = (id2 == 21) ? 22 : id2;
  setId( id1, id2, id3, id4);

  // Colour flow topology. Swap if first is gluon, or when antiquark.
  setColAcol( 1, 0, 2, 1, 2, 0, 0, 0);
  if (id1 == 21) swapCol1234();
  if (id1 < 0 || id2 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2qqbar2ggamma class.
// Cross section for q qbar -> g gamma.

//*********

// Initialize process, especially parton-flux object. 
  
void Sigma2qqbar2ggamma::initProc() {

  // Set up for q qbar initial state, e2-weighted.
  inFluxPtr = new InFluxqqbarSame();
  inFluxPtr->weightCharge2();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qqbar2ggamma::sigmaHat() { 

  // Calculate kinematics dependence.
  double sigTU = (8./9.) * (tH2+uH2)/(tH*uH);

  // Answer.
  return CONVERT2MB * (M_PI/sH2) * alpS * alpEM * sigTU;

}

//*********

// Select identity, colour and anticolour.

void Sigma2qqbar2ggamma::setIdColAcol() {

  // Outgoing flavours trivial.
  setId( id1, id2, 21, 22);

  // One colour flow topology. Swap if first is antiquark.
  setColAcol( 1, 0, 0, 2, 1, 2, 0, 0);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2gg2ggamma class.
// Cross section for g g -> g gamma.
// Proceeds through a quark box, by default using 5 massless quarks.

//*********

// Initialize process, especially parton-flux object. 
  
void Sigma2gg2ggamma::initProc() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxgg();

  // Calculate charge factor from the allowed quarks in the box. 
  int nQuarkInLoop = Settings::mode("SigmaProcess:nQuarkInLoop");
  chargeSum = - 1./3. + 2./3. - 1./3.;
  if (nQuarkInLoop >= 4) chargeSum += 2./3.;
  if (nQuarkInLoop >= 5) chargeSum -= 1./3.;
  if (nQuarkInLoop >= 6) chargeSum += 2./3.;

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2gg2ggamma::sigmaHat() { 

  // Logarithms of Mandelstam variable ratios.
  double logST = log( -sH / tH );
  double logSU = log( -sH / uH );
  double logTU = log(  tH / uH );

  // Real and imaginary parts of separate amplitudes.
  double b0stuRe = 1. + (tH - uH) / sH * logTU 
    + 0.5 * (tH2 + uH2) / sH2 * (pow2(logTU) + pow2(M_PI));   
  double b0stuIm = 0.;
  double b0tsuRe = 1. + (sH - uH) / tH * logSU 
    + 0.5 * (sH2 + uH2) / tH2 * pow2(logSU);
  double b0tsuIm = -M_PI * ( (sH - uH) / tH + (sH2 + uH2) / tH2 * logSU);
  double b0utsRe = 1. + (sH - tH) / uH * logST 
    + 0.5 * (sH2 + tH2) / uH2 * pow2(logST);
  double b0utsIm = -M_PI * ( (sH - tH) / uH + (sH2 + tH2) / uH2 * logST);
  double b1stuRe = -1.;
  double b1stuIm = 0.;
  double b2stuRe = -1.;
  double b2stuIm = 0.;

  // Calculate kinematics dependence.
  double sigBox = pow2(b0stuRe) + pow2(b0stuIm) + pow2(b0tsuRe) 
    + pow2(b0tsuIm) + pow2(b0utsRe) + pow2(b0utsIm) + 4. * pow2(b1stuRe) 
    + 4. * pow2(b1stuIm) + pow2(b2stuRe) + pow2(b2stuIm);
  
  // Answer.
  return CONVERT2MB * (5. / (192. * M_PI * sH2)) * pow2(chargeSum) 
    * pow3(alpS) * alpEM * sigBox;

}

//*********

// Select identity, colour and anticolour.

void Sigma2gg2ggamma::setIdColAcol() {

  // Flavours and colours are trivial.
  setId( id1, id2, 21, 22);
  setColAcol( 1, 2, 2, 3, 1, 3, 0, 0);
  if (Rndm::flat() > 0.5) swapColAcol();

}

//**************************************************************************

// Sigma2qqbar2gammagamma class.
// Cross section for q qbar -> gamma gamma.

//*********

// Initialize process, especially parton-flux object. 
  
void Sigma2qqbar2gammagamma::initProc() {

  // Set up for q qbar initial state, e4-weighted.
  inFluxPtr = new InFluxqqbarSame();
  inFluxPtr->weightCharge2();
  inFluxPtr->weightCharge2();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qqbar2gammagamma::sigmaHat() { 

  // Calculate kinematics dependence. Colour factor for quarks in.
  sigTU = 2. * (tH2+uH2)/(tH*uH);
  double colFac = 1./3.;

  // Answer contains factor 1/2 from identical photons.
  return CONVERT2MB * (M_PI/sH2) * pow2(alpEM) * 0.5 * sigTU * colFac;

}

//*********

// Select identity, colour and anticolour.

void Sigma2qqbar2gammagamma::setIdColAcol() {

  // Outgoing flavours trivial.
  setId( id1, id2, 22, 22);

  // One colour flow topology. Swap if first is antiquark.
  setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}
//**************************************************************************

// Sigma2gg2gammagamma class.
// Cross section for g g -> gamma gamma.
// Proceeds through a quark box, by default using 5 massless quarks.

//*********

// Initialize process, especially parton-flux object. 
  
void Sigma2gg2gammagamma::initProc() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxgg();

  // Calculate charge factor from the allowed quarks in the box. 
  int nQuarkInLoop = Settings::mode("SigmaProcess:nQuarkInLoop");
  charge2Sum = 1./9. + 4./9. + 1./9.;
  if (nQuarkInLoop >= 4) charge2Sum += 4./9.;
  if (nQuarkInLoop >= 5) charge2Sum += 1./9.;
  if (nQuarkInLoop >= 6) charge2Sum += 4./9.;

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2gg2gammagamma::sigmaHat() { 

  // Logarithms of Mandelstam variable ratios.
  double logST = log( -sH / tH );
  double logSU = log( -sH / uH );
  double logTU = log(  tH / uH );

  // Real and imaginary parts of separate amplitudes.
  double b0stuRe = 1. + (tH - uH) / sH * logTU 
    + 0.5 * (tH2 + uH2) / sH2 * (pow2(logTU) + pow2(M_PI));   
  double b0stuIm = 0.;
  double b0tsuRe = 1. + (sH - uH) / tH * logSU 
    + 0.5 * (sH2 + uH2) / tH2 * pow2(logSU);
  double b0tsuIm = -M_PI * ( (sH - uH) / tH + (sH2 + uH2) / tH2 * logSU);
  double b0utsRe = 1. + (sH - tH) / uH * logST 
    + 0.5 * (sH2 + tH2) / uH2 * pow2(logST);
  double b0utsIm = -M_PI * ( (sH - tH) / uH + (sH2 + tH2) / uH2 * logST);
  double b1stuRe = -1.;
  double b1stuIm = 0.;
  double b2stuRe = -1.;
  double b2stuIm = 0.;

  // Calculate kinematics dependence.
  double sigBox = pow2(b0stuRe) + pow2(b0stuIm) + pow2(b0tsuRe) 
    + pow2(b0tsuIm) + pow2(b0utsRe) + pow2(b0utsIm) + 4. * pow2(b1stuRe) 
    + 4. * pow2(b1stuIm) + pow2(b2stuRe) + pow2(b2stuIm);

  // Answer contains factor 1/2 from identical photons.
  return CONVERT2MB * (0.5 / (16. * M_PI * sH2)) * pow2(charge2Sum) 
    * pow2(alpS) * pow2(alpEM) * sigBox;

}

//*********

// Select identity, colour and anticolour.

void Sigma2gg2gammagamma::setIdColAcol() {

  // Flavours and colours are trivial.
  setId( id1, id2, 22, 22);
  setColAcol( 1, 2, 2, 1, 0, 0, 0, 0);

}

//**************************************************************************

// Sigma2ff2ff9gmZ class.
// Cross section for f f' -> f f' via t-channel gamma*/Z0 exchange
// (f is quark or lepton). 

//*********

// Initialize process, especially parton-flux object. 
  
void Sigma2ff2ff9gmZ::initProc() {

  // Set up for f f' initial state.
  inFluxPtr = new InFluxff();

  // Multiply by spin factor 2 for neutrinos.
  inFluxPtr->weightNeutrinoSpin();

  // Store Z0 mass for propagator. Common coupling factor.
  gmZmode = Settings::mode("SigmaProcess:gmZmode");
  mZ = ParticleDataTable::m0(23);
  mZS = mZ*mZ;
  thetaWRat = 1. / (16. * CoupEW::sin2thetaW() * CoupEW::cos2thetaW());  

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 
// Currently only done for quarks??

double Sigma2ff2ff9gmZ::sigmaHat() {

  // Cross section part common for all incoming flavours.
  double sigma0 = (M_PI / sH2) * pow2(alpEM);

  // Kinematical functions for gamma-gamma, gammaZ and Z-Z parts.   
  double sigmagmgm = 2. * (sH2 + uH2) / tH2;
  double sigmagmZ = 4. * thetaWRat * sH2 / (tH * (tH - mZS));
  double sigmaZZ = 2. * pow2(thetaWRat) * sH2 / pow2(tH - mZS);
  if (gmZmode == 1) {sigmagmZ = 0.; sigmaZZ = 0.;}
  if (gmZmode == 2) {sigmagmgm = 0.; sigmagmZ = 0.;} 
  
  // Loop over incoming down- or up-type quarks and read out couplings.
  double e1, v1, a1, e2, v2, a2, sigmaFlav;
  for (int combi = 0; combi < 3; ++combi) {
    int id1 = (combi < 2) ? 1 : 2; 
    int id2 = (combi < 1) ? 1 : 2; 
    e1 = CoupEW::ef(id1);
    v1 = CoupEW::vf(id1);
    a1 = CoupEW::af(id1);
    e2 = CoupEW::ef(id2);
    v2 = CoupEW::vf(id2);
    a2 = CoupEW::af(id2);
    
    // Evaluate and store values for same-sign quarks.
    sigmaFlav = sigmagmgm * pow2(e1 * e2) + sigmagmZ * e1 * e2 
      * (v1 * v2 * (1. + uH2 / sH2) + a1 * a2 * (1. - uH2 / sH2)) 
      + sigmaZZ * ((v1*v1 + a1*a1) * (v2*v2 + a2*a2) * (1. + uH2 / sH2)
      + 4. * v1 * a1 * v2 * a2 * (1. - uH2 / sH2));
    inFluxPtr->weightInState(  id1,  id2, sigmaFlav);
    
    // Evaluate and store values for opposite-sign quarks.
    sigmaFlav = sigmagmgm * pow2(e1 * e2) + sigmagmZ * e1 * e2 
      * (v1 * v2 * (1. + uH2 / sH2) + a1 * a2 * (1. - uH2 / sH2)) 
      + sigmaZZ * ((v1*v1 + a1*a1) * (v2*v2 + a2*a2) * (1. + uH2 / sH2)
      + 4. * v1 * a1 * v2 * a2 * (1. - uH2 / sH2));
    inFluxPtr->weightInState(  id1, -id2, sigmaFlav);
  }

  // Answer, leaving out flavour-dependent pieces stored above.
  return CONVERT2MB * sigma0;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2ff2ff9gmZ::setIdColAcol() {

  // Trivial flavours: out = in.
  setId( id1, id2, id1, id2);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 10 && abs(id2) < 10 && id1*id2 > 0) 
    setColAcol( 1, 0, 2, 0, 1, 0, 2, 0);
  else if (abs(id1) < 10 && abs(id2) < 10)
    setColAcol( 1, 0, 0, 2, 1, 0, 0, 2); 
  else if (abs(id1) < 10) setColAcol( 1, 0, 0, 0, 1, 0, 0, 0); 
  else if (abs(id2) < 10) setColAcol( 0, 0, 1, 0, 0, 0, 1, 0); 
  else setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if ( (abs(id1) < 10 && id1 < 0) || (abs(id1) > 10 && id2 < 0) ) 
    swapColAcol();

}

//**************************************************************************

// Sigma2ff2ff9W class.
// Cross section for f_1 f_2 -> f_3 f_4 via t-channel W+- exchange
// (f is quark or lepton). 

//*********

// Initialize process, especially parton-flux object. 
  
void Sigma2ff2ff9W::initProc() {

  // Set up for f f' initial state.
  inFluxPtr = new InFluxff();

  // Multiply by sum of relevant squared CKM matrix elements.
  inFluxPtr->weightCKM2sum(3);

  // Multiply by spin factor 2 for neutrinos.
  inFluxPtr->weightNeutrinoSpin();

  // Store W+- mass for propagator. Common coupling factor.
  mW = ParticleDataTable::m0(24);
  mWS = mW*mW;
  thetaWRat = 1. / (4. * CoupEW::sin2thetaW());  

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 
// Currently only done for quarks??

double Sigma2ff2ff9W::sigmaHat() {

  // Cross section part common for all incoming flavours.
  double sigma0 = (M_PI / sH2) * pow2(alpEM * thetaWRat)
    * 4. * sH2 / pow2(tH - mWS);
    
  // Store values for allowed quark combinations.
  inFluxPtr->weightInState(  1,  2, 1.);
  inFluxPtr->weightInState(  1, -1, uH2 / sH2);
  inFluxPtr->weightInState(  2, -2, uH2 / sH2);

  // Answer, leaving out flavour-dependent pieces stored above.
  return CONVERT2MB * sigma0;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2ff2ff9W::setIdColAcol() {

  // Pick out-flavours by relative CKM weights.
  int id3 = VCKM::V2pick(id1);
  int id4 = VCKM::V2pick(id2);
  setId( id1, id2, id3, id4);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 10 && abs(id2) < 10 && id1*id2 > 0) 
    setColAcol( 1, 0, 2, 0, 1, 0, 2, 0);
  else if (abs(id1) < 10 && abs(id2) < 10)
    setColAcol( 1, 0, 0, 2, 1, 0, 0, 2); 
  else if (abs(id1) < 10) setColAcol( 1, 0, 0, 0, 1, 0, 0, 0); 
  else if (abs(id2) < 10) setColAcol( 0, 0, 1, 0, 0, 0, 1, 0); 
  else setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if ( (abs(id1) < 10 && id1 < 0) || (abs(id1) > 10 && id2 < 0) ) 
    swapColAcol();

}


//**************************************************************************

// Sigma2qq2Qq9W class.
// Cross section for q q' -> Q q" via t-channel W+- exchange. 
// Related to Sigma2ff2ffViaW class, but with massive matrix elements.

//*********

// Initialize process, especially parton-flux object. 
  
void Sigma2qq2Qq9W::initProc() {

  // Set up for f f' initial state.
  inFluxPtr = new InFluxff();

  // Multiply by sum of relevant squared CKM matrix elements.
  inFluxPtr->weightCKM2sum(4, idNew);

  // Store W+- mass for propagator. Common coupling factor.
  mW = ParticleDataTable::m0(24);
  mWS = mW*mW;
  thetaWRat = 1. / (4. * CoupEW::sin2thetaW());  

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qq2Qq9W::sigmaHat() {

  // Cross section part common for all incoming flavours.
  double sigma0 = (M_PI / sH2) * pow2(alpEM * thetaWRat)
    * 4. / pow2(tH - mWS);
    
  // Store values for allowed quark combinations.
  inFluxPtr->weightInState(  1,  2, sH * (sH - s3));
  if (idNew%2 == 0) {
    inFluxPtr->weightInState(  1, -1, uH * (uH - s3));
    inFluxPtr->weightInState(  2, -2, 0.);
  } else {
    inFluxPtr->weightInState(  2, -2, uH * (uH - s3));
    inFluxPtr->weightInState(  1, -1, 0.);
  }

  // Answer, leaving out flavour-dependent pieces stored above.
  return CONVERT2MB * sigma0;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2qq2Qq9W::setIdColAcol() {

  // For topologies like d dbar -> (t/c/u) (t/c/u)bar pick side.
  int id1Abs = abs(id1);
  int id2Abs = abs(id2);
  int side = 1;
  if ( (id1Abs + idNew)%2 == 1 && (id2Abs + idNew)%2 == 1 ) {
    double prob1 = VCKM::V2id(id1Abs, idNew) * VCKM::V2sum(id2Abs);
    double prob2 = VCKM::V2id(id2Abs, idNew) * VCKM::V2sum(id1Abs);
    if (prob2 > Rndm::flat() * (prob1 + prob2)) side = 2;
  } 
  else if ((id2Abs + idNew)%2 == 1) side = 2;

  // Pick out-flavours by relative CKM weights.
  int id3, id4; 
  if (side == 1) {
    // q q' -> t q" : correct order from start.
    id3 = (id1 > 0) ? idNew : -idNew;
    id4 = VCKM::V2pick(id2);
    setId( id1, id2, id3, id4);
  } else {
    // q q' -> q" t : stored as t q" so swap tHat <-> uHat.
    swapTU = true;   
    id3 = VCKM::V2pick(id1);
    id4 = (id2 > 0) ? idNew : -idNew;
    setId( id1, id2, id4, id3);
  }

  // Colour flow topologies. Swap when antiquarks.
  if (id1*id2 > 0) setColAcol( 1, 0, 2, 0, 1, 0, 2, 0);
  else             setColAcol( 1, 0, 0, 2, 1, 0, 0, 2); 
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// Sigma1ffbar2gmZ class.
// Cross section for f fbar -> gamma*/Z0 (f is quark or lepton). 

//*********

// Initialize process, especially parton-flux object. 
  
void Sigma1ffbar2gmZ::initProc() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxffbarSame();

  // Multiply by colour factor 1/3.
  inFluxPtr->weightInvCol();

} 

//*********

// Evaluate sigmaHat(sHat). 
// Note: owing to the interference structure in this process,
// it is not possible to decouple sigmaHat from the incoming parton 
// composition. Part of the cross section is therefore put in inFlux. 

double Sigma1ffbar2gmZ::sigmaHat() { 

  // Store cross section separately for down- and up-type incoming flavours.
  GmZRes.width(sqrt(sH));
  double sigmaD = GmZRes.sigma(1);
  double sigmaU = GmZRes.sigma(2);
  inFluxPtr->weightInState( 1, -1, sigmaD);
  inFluxPtr->weightInState( 2, -2, sigmaU);

  // Answer here is only the conversion factor, rest in inFlux.
  return CONVERT2MB;    

}

//*********

// Select identity, colour and anticolour.

void Sigma1ffbar2gmZ::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, 23);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 10) setColAcol( 1, 0, 0, 1, 0, 0);
  else               setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}


//************************************

// Sigma1ffbar2W class.
// Cross section for f fbar' -> W+- (f is quark or lepton). 

//*********

// Initialize process, especially parton-flux object. 
  
void Sigma1ffbar2W::initProc() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxffbarChg();

  // Multiply by squared CKM matrix elements and colour factor 1/3.
  inFluxPtr->weightCKM2();
  inFluxPtr->weightInvCol();

} 

//*********

// Evaluate sigmaHat(sHat). 

double Sigma1ffbar2W::sigmaHat() {

  // Find cross section part common for all incoming flavours.
  WRes.width(sqrt(sH));
  double sigma = WRes.sigma(1);

  // Answer.
  return CONVERT2MB * sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma1ffbar2W::setIdColAcol() {

  // Sign of outgoing W.
  int sign = 1 - 2 * (abs(id1)%2);
  if (id1 < 0) sign = -sign;
  setId( id1, id2, 24 * sign);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 10) setColAcol( 1, 0, 0, 1, 0, 0);
  else               setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2ffbar2ZW class.
// Cross section for f fbar' -> W+ W- (f is quark or lepton). 

//*********

// Initialize process, especially parton-flux object. 
  
void Sigma2ffbar2ZW::initProc() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxffbarChg();

  // Multiply by squared CKM matrix elements and colour factor 1/3.
  inFluxPtr->weightCKM2();
  inFluxPtr->weightInvCol();

  // Store W+- mass and width for propagator. 
  mW = ParticleDataTable::m0(24);
  widW = ParticleDataTable::mWidth(24);
  mWS = mW*mW;
  mwWS = pow2(mW * widW);

  // Left-handed couplings for up- and downtype quarks.
  lu = CoupEW::lf(2);
  ld = CoupEW::lf(1); 

  // Common weak coupling factor.
  thetaWRat = 1. / (4. * CoupEW::cos2thetaW());  
  thetaWpt = (9. - 8. * CoupEW::sin2thetaW()) / 4.;
  thetaWmm = (8. * CoupEW::sin2thetaW() - 6.) / 4.;

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2ffbar2ZW::sigmaHat() {

  // Evaluate cross section.
  double resBW = 1. / (pow2(sH - mWS) + mwWS);
  double sigma = (M_PI / sH2) * 0.5 * pow2(alpEM / CoupEW::sin2thetaW());
  sigma *= sH * resBW * (thetaWpt * pT2 + thetaWmm * (s3 + s4))
    + (sH - mWS) * resBW * sH * (pT2 - s3 - s4) * (lu / tH - ld / uH)
    + thetaWRat * sH * pT2 * ( lu*lu / tH2 + ld*ld / uH2 )
    + 2. * thetaWRat * sH * (s3 + s4) * lu * ld / (tH * uH);  

  // Answer, leaving out flavour-dependent pieces stored above.
  return CONVERT2MB * sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2ffbar2ZW::setIdColAcol() {

  // Sign of outgoing W. 
  int sign = 1 - 2 * (abs(id1)%2);
  if (id1 < 0) sign = -sign;
  setId( id1, id2, 23, 24 * sign);

  // tHat is defined between (f, W-) or (fbar, W+), 
  // so OK for u/ubar on side 1, but must swap tHat <-> uHat if d/dbar.   
  if (abs(id1)%2 == 1) swapTU = true;

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 10) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else               setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2ffbar2WW class.
// Cross section for f fbar -> W+ W- (f is quark or lepton). 

//*********

// Initialize process, especially parton-flux object. 
  
void Sigma2ffbar2WW::initProc() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxffbarSame();

  // Multiply by colour factor 1/3.
  inFluxPtr->weightInvCol();

  // Store Z0 mass and width for propagator. Common coupling factor.
  mZ = ParticleDataTable::m0(23);
  widZ = ParticleDataTable::mWidth(23);
  mZS = mZ*mZ;
  mwZS = pow2(mZ * widZ);
  thetaWRat = 1. / (4. * CoupEW::sin2thetaW());  

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2ffbar2WW::sigmaHat() {

  // Cross section part common for all incoming flavours.
  double sigma = (M_PI / sH2) * pow2(alpEM);

  // Z0 propagator and gamma*/Z0 interference.
  double Zprop = sH2 / (pow2(sH - mZS) + mwZS);
  double Zint = Zprop * (1. - mZS / sH);

  // Common coupling factors (g = gamma*, Z = Z0, f = t-channel fermion).
  double cgg = 0.5;
  double cgZ = thetaWRat * Zint;
  double cZZ = 0.5 * pow2(thetaWRat) * Zprop;
  double cfg = thetaWRat;
  double cfZ = pow2(thetaWRat) * Zint;
  double cff = pow2(thetaWRat);

  // Kinematical functions.   
  double rat34 = sH * (2. * (s3 + s4) + pT2) / (s3 * s4);
  double lambdaS = pow2(sH - s3 - s4) - 4.* s3 * s4;
  double intA = (sH - s3 - s4) * rat34 / sH;
  double intB = 4. * (s3 + s4 - pT2);
  double gSS = (lambdaS * rat34 + 12. * sH * pT2) / sH2;
  double gTT = rat34 + 4. * sH * pT2 / tH2;
  double gST = intA + intB / tH;
  double gUU = rat34 + 4. * sH * pT2 / uH2;
  double gSU = intA + intB / uH;   
   
  // Answer for down-type quarks.
  double ei = CoupEW::ef(1);
  double vi = CoupEW::vf(1); 
  double ai = CoupEW::af(1); 
  double sigmaD = (cgg * ei*ei + cgZ * ei * vi + cZZ * (vi*vi + ai*ai)) * gSS
    + (cfg * ei + cfZ * (vi + ai)) * gST + cff * gTT;
  
  // Answer for up-type quarks.
  ei = CoupEW::ef(2);
  vi = CoupEW::vf(2); 
  ai = CoupEW::af(2); 
  double sigmaU = (cgg * ei*ei + cgZ * ei * vi + cZZ * (vi*vi + ai*ai)) * gSS
    - (cfg * ei + cfZ * (vi + ai)) * gSU + cff * gUU;

  // Store cross section separately for down- and up-type incoming flavours.
  inFluxPtr->weightInState( 1, -1, sigmaD);
  inFluxPtr->weightInState( 2, -2, sigmaU);

  // Answer, leaving out flavour-dependent pieces stored above.
  return CONVERT2MB * sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2ffbar2WW::setIdColAcol() {

  // Sign of outgoing W's. W+ always first???
  int sign = (id1 < 0) ? 1 : -1;
  setId( id1, id2, 24 * sign, -24 * sign);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 10) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else               setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2qqbar2Wg class.
// Cross section for q qbar' -> W+- g. 

//*********

// Initialize process, especially parton-flux object. 
  
void Sigma2qqbar2Wg::initProc() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxffbarChg();

  // Multiply by squared CKM matrix elements.
  inFluxPtr->weightCKM2();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qqbar2Wg::sigmaHat() {

  // Cross section part common for all incoming flavours.
  double sigma = (M_PI / sH2) * (alpEM * alpS / CoupEW::sin2thetaW())
    * (2./9.) * (tH2 + uH2 + 2. * sH * s3) / (tH * uH);

  // Answer.
  return CONVERT2MB * sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2qqbar2Wg::setIdColAcol() {

  // Sign of outgoing W.
  int sign = 1 - 2 * (abs(id1)%2);
  if (id1 < 0) sign = -sign;
  setId( id1, id2, 24 * sign, 21);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 2, 0, 0, 1, 2);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2qg2Wq class.
// Cross section for q g -> W+- q'. 

//*********

// Initialize process, especially parton-flux object. 
  
void Sigma2qg2Wq::initProc() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxqg();

  // Multiply by sum of relevant squared CKM matrix elements.
  inFluxPtr->weightCKM2sum(1);

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qg2Wq::sigmaHat() {

  // Cross section part common for all incoming flavours (3 = W).
  double sigma = (M_PI / sH2) * (alpEM * alpS / CoupEW::sin2thetaW())
    * (1./12.) * (sH2 + uH2 + 2. * tH * s3) / (-sH * uH);

  // Answer.
  return CONVERT2MB * sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2qg2Wq::setIdColAcol() {

  // Sign of outgoing W. Flavour of outgoing quark.
  int idq = (id2 == 21) ? id1 : id2;
  int sign = 1 - 2 * (abs(idq)%2);
  if (idq < 0) sign = -sign;
  id4 = VCKM::V2pick(idq);

  // Flavour set up for q g -> W q.
  setId( id1, id2, 24 * sign, id4);

  // tH defined between f and f': must swap tHat <-> uHat if q g in.
  swapTU = (id2 == 21) ? true : false; 

  // Colour flow topologies. Swap when antiquarks.
  if (id2 == 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else           setColAcol( 2, 1, 1, 0, 0, 0, 2, 0);
  if (idq < 0) swapColAcol();

}

//**************************************************************************

// Sigma2ffbar2Wgm class.
// Cross section for f fbar' -> W+- gamma. 

//*********

// Initialize process, especially parton-flux object. 
  
void Sigma2ffbar2Wgm::initProc() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxffbarChg();

  // Multiply by squared CKM matrix elements and colour factor 1/3.
  inFluxPtr->weightCKM2();
  inFluxPtr->weightInvCol();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2ffbar2Wgm::sigmaHat() {

  // Cross section part common for all incoming flavours.
  double sigma = (M_PI / sH2) * (alpEM*alpEM / CoupEW::sin2thetaW())
    * 0.5 * (tH2 + uH2 + 2. * sH * s3) / (tH * uH);

  // Extrafactor different for q qbar' and e nu instate. ??
  sigma *= pow2( uH / (tH + uH) - 1./3.);
  // sigma *= pow2( tH / (tH + uH) );

  // Answer.
  return CONVERT2MB * sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2ffbar2Wgm::setIdColAcol() {

  // Sign of outgoing W.
  int sign = 1 - 2 * (abs(id1)%2);
  if (id1 < 0) sign = -sign;
  setId( id1, id2, 24 * sign, 21);

  // tH defined between (f,W-) or (fbar',W+).
  swapTU = (sign * id1 > 0) ? true : false;

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 10) setColAcol( 1, 0, 0, 2, 0, 0, 1, 2);
  else               setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

} // end namespace Pythia8
