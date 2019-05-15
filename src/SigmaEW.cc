// Function definitions (not found in the header) for the 
// electroweak simulation classes (including photon processes). 
// Copyright C 2007 Torbjorn Sjostrand

#include "SigmaEW.h"

namespace Pythia8 {

//**************************************************************************

// Sigma2qg2qgamma class.
// Cross section for q g -> q gamma.

//*********

// Initialize parton-flux object. 
  
void Sigma2qg2qgamma::initFlux() {

  // Set up for q qbar initial state, e2-weighted.
  inFluxPtr = new InFluxqg();
  inFluxPtr->weightCharge2();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qg2qgamma::sigmaHat() {  

  // Calculate kinematics dependence.
  sigUS = (1./3.) * (sH2+uH2)/(-sH*uH);

  // Flavours and thereby charge already fixed for Multiple Interactions.
  if (id12IsSet) { 
    int idNow    = (id2 == 21) ? id1 : id2;    
    double eNow  = CoupEW::ef(idNow);    
    sigUS       *= pow2(eNow);
  }

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

// Initialize parton-flux object. 
  
void Sigma2qqbar2ggamma::initFlux() {

  // Set up for q qbar initial state, e2-weighted.
  inFluxPtr = new InFluxqqbarSame();
  inFluxPtr->weightCharge2();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qqbar2ggamma::sigmaHat() { 

  // Calculate kinematics dependence.
  double sigTU = (8./9.) * (tH2+uH2)/(tH*uH);

  // Flavours and thereby charge already fixed for Multiple Interactions.
  if (id12IsSet) { 
    double eNow  = CoupEW::ef( abs(id1) );    
    sigTU       *= pow2(eNow);
  }

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

  // Calculate charge factor from the allowed quarks in the box. 
  int nQuarkInLoop = Settings::mode("SigmaProcess:nQuarkInLoop");
  chargeSum                         = - 1./3. + 2./3. - 1./3.;
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

// Initialize parton-flux object. 
  
void Sigma2qqbar2gammagamma::initFlux() {

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

  // Flavours and thereby charge already fixed for Multiple Interactions.
  if (id12IsSet) { 
    double eNow  = CoupEW::ef( abs(id1) );    
    sigTU       *= pow4(eNow);
  }

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

// Initialize process. 
  
void Sigma2gg2gammagamma::initProc() {

  // Calculate charge factor from the allowed quarks in the box. 
  int nQuarkInLoop = Settings::mode("SigmaProcess:nQuarkInLoop");
  charge2Sum                         = 1./9. + 4./9. + 1./9.;
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

// Sigma2ff2fftgmZ class.
// Cross section for f f' -> f f' via t-channel gamma*/Z0 exchange
// (f is quark or lepton). 

//*********

// Initialize process. 
  
void Sigma2ff2fftgmZ::initProc() {

  // Store Z0 mass for propagator. Common coupling factor.
  gmZmode   = Settings::mode("SigmaProcess:gmZmode");
  mZ        = ParticleDataTable::m0(23);
  mZS       = mZ*mZ;
  thetaWRat = 1. / (16. * CoupEW::sin2thetaW() * CoupEW::cos2thetaW());  

} 

//*********

// Initialize parton-flux object. 
  
void Sigma2ff2fftgmZ::initFlux() {

  // Set up for f f' initial state.
  inFluxPtr = new InFluxff();

    // Multiply by spin factor 2 for neutrinos.
  inFluxPtr->weightNeutrinoSpin();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 
// Currently only done for quarks??

double Sigma2ff2fftgmZ::sigmaHat() {

  // Cross section part common for all incoming flavours.
  double sigma0 = (M_PI / sH2) * pow2(alpEM);

  // Kinematical functions for gamma-gamma, gamma-Z and Z-Z parts.   
  double sigmagmgm = 2. * (sH2 + uH2) / tH2;
  double sigmagmZ  = 4. * thetaWRat * sH2 / (tH * (tH - mZS));
  double sigmaZZ   = 2. * pow2(thetaWRat) * sH2 / pow2(tH - mZS);
  if (gmZmode == 1) {sigmagmZ = 0.; sigmaZZ = 0.;}
  if (gmZmode == 2) {sigmagmgm = 0.; sigmagmZ = 0.;} 
  
  // Fix flavours for lepton beams or multiple interactions, else 
  // three possibilites in Loop over incoming down- or up-type quarks
  int id1Now, id2Now;
  double e1, v1, a1, e2, v2, a2, sigmaFlav;
  int combiMax = (hasLeptonBeams || id12IsSet) ? 1 : 3; 
  for (int combi = 0; combi < combiMax; ++combi) {
    if (hasLeptonBeams) {
      id1Now = abs(idA);
      id2Now = abs(idB);
    } else if (id12IsSet) {
      id1Now = abs(id1);
      id2Now = abs(id2);
    } else {
      id1Now = (combi < 2) ? 1 : 2; 
      id2Now = (combi < 1) ? 1 : 2; 
    }

    // Read out couplings for chosen flavour combination.
    e1 = CoupEW::ef(id1Now);
    v1 = CoupEW::vf(id1Now);
    a1 = CoupEW::af(id1Now);
    e2 = CoupEW::ef(id2Now);
    v2 = CoupEW::vf(id2Now);
    a2 = CoupEW::af(id2Now);
    
    // Evaluate and store values for same-sign quarks/leptons.
    sigmaFlav = sigmagmgm * pow2(e1 * e2) + sigmagmZ * e1 * e2 
      * (v1 * v2 * (1. + uH2 / sH2) + a1 * a2 * (1. - uH2 / sH2)) 
      + sigmaZZ * ((v1*v1 + a1*a1) * (v2*v2 + a2*a2) * (1. + uH2 / sH2)
      + 4. * v1 * a1 * v2 * a2 * (1. - uH2 / sH2));
    if (hasLeptonBeams && idA * idB > 0) sigma0 *= sigmaFlav; 
    else if (id12IsSet && id1 * id2 > 0) sigma0 *= sigmaFlav;
    else if (!hasLeptonBeams && !id12IsSet) 
      inFluxPtr->weightInState(  id1Now,  id2Now, sigmaFlav);
    
    // Evaluate and store values for opposite-sign quarks/leptons.
    sigmaFlav = sigmagmgm * pow2(e1 * e2) + sigmagmZ * e1 * e2 
      * (v1 * v2 * (1. + uH2 / sH2) + a1 * a2 * (1. - uH2 / sH2)) 
      + sigmaZZ * ((v1*v1 + a1*a1) * (v2*v2 + a2*a2) * (1. + uH2 / sH2)
      + 4. * v1 * a1 * v2 * a2 * (1. - uH2 / sH2));
    if (hasLeptonBeams && idA * idB < 0) sigma0 *= sigmaFlav; 
    else if (id12IsSet && id1 * id2 < 0) sigma0 *= sigmaFlav;
    else if (!hasLeptonBeams && !id12IsSet) 
      inFluxPtr->weightInState(  id1Now, -id2Now, sigmaFlav);
  }

  // Answer, leaving out flavour-dependent pieces stored above.
  return CONVERT2MB * sigma0;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2ff2fftgmZ::setIdColAcol() {

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

// Sigma2ff2fftW class.
// Cross section for f_1 f_2 -> f_3 f_4 via t-channel W+- exchange
// (f is quark or lepton). 

//*********

// Initialize process. 
  
void Sigma2ff2fftW::initProc() {

  // Store W+- mass for propagator. Common coupling factor.
  mW        = ParticleDataTable::m0(24);
  mWS       = mW*mW;
  thetaWRat = 1. / (4. * CoupEW::sin2thetaW());  

} 

//*********

// Initialize parton-flux object. 
  
void Sigma2ff2fftW::initFlux() {

  // Set up for f f' initial state.
  inFluxPtr = new InFluxff();

  // Multiply by sum of relevant squared CKM matrix elements.
  if (!hasLeptonBeams) inFluxPtr->weightCKM2sum(3);

  // Multiply by spin factor 2 for neutrinos.
  inFluxPtr->weightNeutrinoSpin();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 
// Currently only done for quarks??

double Sigma2ff2fftW::sigmaHat() {

  // Cross section part common for all incoming flavours.
  double sigma0 = (M_PI / sH2) * pow2(alpEM * thetaWRat)
    * 4. * sH2 / pow2(tH - mWS);

  // Simple answer for lepton beams.
  if (hasLeptonBeams) {
    if (abs(idA)%2 == abs(idB)%2) sigma0 *= uH2 / sH2;
    
  // Find value for known incoming flavour.
  } else if (id12IsSet) {
    if (abs(id1)%2 == abs(id2)%2) sigma0 *= uH2 / sH2;
    sigma0 *= VCKM::V2sum(id1) *  VCKM::V2sum(id2);

  // Else store values for allowed quark combinations.
  } else {
    inFluxPtr->weightInState(  1,  2, 1.);
    inFluxPtr->weightInState(  1, -1, uH2 / sH2);
    inFluxPtr->weightInState(  2, -2, uH2 / sH2);
  }

  // Answer, leaving out flavour-dependent pieces stored above.
  return CONVERT2MB * sigma0;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2ff2fftW::setIdColAcol() {

  // Pick out-flavours by relative CKM weights.
  int id3 = VCKM::V2pick(id1);
  int id4 = VCKM::V2pick(id2);
  setId( id1, id2, id3, id4);

  // Colour flow topologies. Swap when antiquarks.
  if      (abs(id1) < 10 && abs(id2) < 10 && id1*id2 > 0) 
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

// Sigma2qq2QqtW class.
// Cross section for q q' -> Q q" via t-channel W+- exchange. 
// Related to Sigma2ff2ffViaW class, but with massive matrix elements.

//*********

// Initialize process. 
  
void Sigma2qq2QqtW::initProc() {

  // Store W+- mass for propagator. Common coupling factor.
  mW        = ParticleDataTable::m0(24);
  mWS       = mW*mW;
  thetaWRat = 1. / (4. * CoupEW::sin2thetaW());  

} 

//*********

// Initialize parton-flux object. 
  
void Sigma2qq2QqtW::initFlux() {

  // Set up for f f' initial state.
  inFluxPtr = new InFluxff();

  // Multiply by sum of relevant squared CKM matrix elements.
  inFluxPtr->weightCKM2sum(4, idNew);


} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qq2QqtW::sigmaHat() {

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

void Sigma2qq2QqtW::setIdColAcol() {

  // For topologies like d dbar -> (t/c/u) (t/c/u)bar pick side.
  int id1Abs = abs(id1);
  int id2Abs = abs(id2);
  int side   = 1;
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

  // Colour flow topologies. Swap when antiquarks on side 1.
  if      (side == 1 && id1 * id2 > 0) setColAcol( 1, 0, 2, 0, 1, 0, 2, 0);
  else if              (id1 * id2 > 0) setColAcol( 1, 0, 2, 0, 2, 0, 1, 0);
  else if (side == 1)                  setColAcol( 1, 0, 0, 2, 1, 0, 0, 2); 
  else                                 setColAcol( 1, 0, 0, 2, 0, 2, 1, 0); 
  if (id1 < 0) swapColAcol();

}

//*********

// Evaluate weight for decay angles of W in top decay.

double Sigma2qq2QqtW::weightDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // For top decay hand over to standard routine, else done.
  if (idNew == 6 && process[process[iResBeg].mother1()].idAbs() == 6) 
       return weightTopDecay( process, iResBeg, iResEnd);
  else return 1.; 

}

//**************************************************************************

// Sigma1ffbar2gmZ class.
// Cross section for f fbar -> gamma*/Z0 (f is quark or lepton). 

//*********

// Initialize process. 
  
void Sigma1ffbar2gmZ::initProc() {

  // Allow to pick only gamma* or Z0 part of full gamma*/Z0 expression.
  gmZmode     = Settings::mode("SigmaProcess:gmZmode");

  // Store Z0 mass and width for propagator. 
  mRes        = ParticleDataTable::m0(23);
  GammaRes    = ParticleDataTable::mWidth(23);
  m2Res       = mRes*mRes;
  GamMRat     = GammaRes / mRes;
  thetaWRat   = 1. / (16. * CoupEW::sin2thetaW() * CoupEW::cos2thetaW());

  // Set pointer to particle properties and decay table.
  particlePtr = ParticleDataTable::particleDataPtr(23);

} 

//*********

// Initialize parton-flux object. 
  
void Sigma1ffbar2gmZ::initFlux() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxffbarSame();

  // Multiply by colour factor 1/3.
  inFluxPtr->weightInvCol();

} 

//*********

// Evaluate sigmaHat(sHat). 
// Note: owing to the interference structure in this process,
// it is not possible to decouple sigmaHat from the incoming parton 
// composition. Most of the cross section is therefore put in inFlux. 

double Sigma1ffbar2gmZ::sigmaHat() { 

  // Resonance mass. Common coupling factors.
  double mHat  = sqrt(sH);
  double colQ = 3. * (1. + alpS / M_PI);

  // Reset quantities to sum. Declare variables in loop.
  double gamSum = 0.;
  double intSum = 0.;
  double resSum = 0.;
  double mf, m2Rat, psvec, psaxi, betaf, ef2, efvf, vf2af2, colf;

  // Loop over all Z0 decay channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    int idAbs = abs( particlePtr->decay[i].product(0) );

    // Only contributions from three fermion generations, except top.
    if ( (idAbs > 0 && idAbs < 6) || ( idAbs > 10 && idAbs < 17)) {
      mf = ParticleDataTable::m0(idAbs);

      // Check that above threshold. Phase space.
      if (mHat > 2. * mf + MASSMARGIN) {
        m2Rat = pow2(mf / mHat);
        betaf = sqrtpos(1. - 4. * m2Rat); 
        psvec = betaf * (1. + 2. * m2Rat);
        psaxi  = pow3(betaf);

        // Combine phase space with couplings.
        ef2    = CoupEW::ef2(idAbs) * psvec;
        efvf   = CoupEW::efvf(idAbs) * psvec;
        vf2af2 = CoupEW::vf2(idAbs) * psvec + CoupEW::af2(idAbs) * psaxi; 
        colf   = (idAbs < 6) ? colQ : 1.;

        // Store sum of combinations. For outstate only open channels.
        if (particlePtr->decay[i].onMode() > 0) {
          gamSum += colf * ef2;
          intSum += colf * efvf;
          resSum += colf * vf2af2;
	}

      // End loop over fermions.
      }
    }
  }

  // Calculate prefactors for gamma/interference/Z0 cross section terms.
  gamProp = 4. * M_PI * pow2(alpEM) / (3. * sH); 
  intProp = gamProp * 2. * thetaWRat * sH * (sH - m2Res)
          / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  resProp = gamProp * pow2(thetaWRat * sH) 
          / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );

  // Optionally only keep gamma* or Z0 term.
  if (gmZmode == 1) {intProp = 0.; resProp = 0.;}
  if (gmZmode == 2) {gamProp = 0.; intProp = 0.;}

  // Loop over incoming fermion flavour species.
  for (int i = 0; i < 2; ++i) {
    int idAbs = i + 1;
    if (hasLeptonBeams && i == 0) idAbs = abs(idA);
    else if (hasLeptonBeams) continue;

    // Combine gamma, interference and Z0 parts.
    double sigma =  CoupEW::ef2(idAbs)    * gamProp * gamSum 
                  + CoupEW::efvf(idAbs)   * intProp * intSum
                  + CoupEW::vf2af2(idAbs) * resProp * resSum;

    // Send flavour-dependent parts to inFlux.
    inFluxPtr->weightInState( idAbs, -idAbs, sigma);
  }

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

//*********

// Evaluate weight for gamma*/Z0 decay angle.
  
double Sigma1ffbar2gmZ::weightDecay( Event& process, int iResBeg, int iResEnd) {

  // Z should sit in entry 5.
  if (iResBeg != 5 || iResEnd != 6) return 1.;

  // Couplings for in- and out-flavours.
  int idInAbs  = process[3].idAbs();
  double ei    = CoupEW::ef(idInAbs);
  double vi    = CoupEW::vf(idInAbs);
  double ai    = CoupEW::af(idInAbs);
  int idOutAbs = process[6].idAbs();
  double ef    = CoupEW::ef(idOutAbs);
  double vf    = CoupEW::vf(idOutAbs);
  double af    = CoupEW::af(idOutAbs);

  // Phase space factors. (One power of beta left out in formulae.)
  double mf    = process[6].m();
  double m2Rat = mf*mf / sH;
  double betaf = sqrtpos(1. - 4. * m2Rat); 

  // Coefficients of angular expression.
  double coefTran = ei*ei * gamProp * ef*ef + ei * vi * intProp * ef * vf
    + (vi*vi + ai*ai) * resProp * (vf*vf + pow2(betaf) * af*af);
  double coefLong = 4. * m2Rat * ( ei*ei * gamProp * ef*ef 
    + ei * vi * intProp * ef * vf + (vi*vi + ai*ai) * resProp * vf*vf );
  double coefAsym = betaf * ( ei * ai * intProp * ef * af 
    + 4. * vi * ai * resProp * vf * af );

  // Flip asymmetry for in-fermion + out-antifermion.
  if (process[3].id() * process[6].id() < 0) coefAsym = -coefAsym;

  // Reconstruct decay angle and weight for it.
  double cosThe = 2. * process[3].p() * (process[7].p() - process[6].p())
    / (sH * betaf);
  double wtMax = 2. * (coefTran + abs(coefAsym));
  double wt    = coefTran * (1. + pow2(cosThe)) 
     + coefLong * (1. - pow2(cosThe)) + 2. * coefAsym * cosThe;

  // Done.
  return (wt / wtMax);

}

//************************************

// Sigma1ffbar2W class.
// Cross section for f fbar' -> W+- (f is quark or lepton). 

//*********

// Initialize process. 
  
void Sigma1ffbar2W::initProc() {

  // Store W+- mass and width for propagator. 
  mRes     = ParticleDataTable::m0(24);
  GammaRes = ParticleDataTable::mWidth(24);
  m2Res    = mRes*mRes;
  GamMRat  = GammaRes / mRes;

} 

//*********

// Initialize parton-flux object. 
  
void Sigma1ffbar2W::initFlux() {

  // Set up for f fbar initial state.
  inFluxPtr = new InFluxffbarChg();

  // Multiply by squared CKM matrix elements and colour factor 1/3.
  inFluxPtr->weightCKM2();
  inFluxPtr->weightInvCol();

} 

//*********

// Evaluate sigmaHat(sHat). 

double Sigma1ffbar2W::sigmaHat() {

  // Incoming width, excluding channel-dependent initFlux() factors.
  double mHat  = sqrt(sH);
  double widthIn = WRes.widthIn(mHat);
 
  // Set up Breit-Wigner. Width out only includes open channels. 
  double sigBW    = 12. * M_PI / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );    
  double widthOut = WRes.width( mHat, true);

  // Done.
  return CONVERT2MB * widthIn * sigBW * widthOut;    

}

//*********

// Select identity, colour and anticolour.

void Sigma1ffbar2W::setIdColAcol() {

  // Sign of outgoing W.
  int sign          = 1 - 2 * (abs(id1)%2);
  if (id1 < 0) sign = -sign;
  setId( id1, id2, 24 * sign);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 10) setColAcol( 1, 0, 0, 1, 0, 0);
  else               setColAcol( 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//*********

// Evaluate weight for W decay angle.
  
double Sigma1ffbar2W::weightDecay( Event& process, int iResBeg, int iResEnd) {

  // W should sit in entry 5.
  if (iResBeg != 5 || iResEnd != 6) return 1.;

  // Phase space factors.
  double m2Rat1 = pow2(process[6].m()) / sH;
  double m2Rat2 = pow2(process[7].m()) / sH;
  double betaf = sqrtpos( pow2(1. - m2Rat1 - m2Rat2) - 4. * m2Rat1 * m2Rat2); 
   
  // Sign of asymmetry.
  double eps  = (process[3].id() * process[6].id() > 0) ? 1. : -1.;

  // Reconstruct decay angle and weight for it.
  double cosThe = 2. * process[3].p() * (process[7].p() - process[6].p())
    / (sH * betaf);
  double wtMax = 4.;
  double wt  = pow2(1. + betaf * eps * cosThe) - pow2(m2Rat1 - m2Rat2); 
 
  // Done.
  return (wt / wtMax);

}

//**************************************************************************

// Sigma2ffbar2ffbarsgm class.
// Cross section f fbar -> gamma* -> f' fbar', for multiple interactions.

//*********

// Initialize parton-flux object. 
  
void Sigma2ffbar2ffbarsgm::initFlux() {

  // Set up for f fbar initial state.
  inFluxPtr = new InFluxffbarSame();

  // Multiply by colour factor 1/3 and charge factor.
  inFluxPtr->weightInvCol();
  inFluxPtr->weightCharge2();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2ffbar2ffbarsgm::sigmaHat() { 

  // Pick new flavour. Allow three leptons and five quarks.
  double colQ     = 1. + (alpS / M_PI);
  double flavWt   = 3. + colQ * 11. / 3.;
  double flavRndm = Rndm::flat() * flavWt;
  if (flavRndm < 3.) {
    if      (flavRndm < 1.) idNew = 11;
    else if (flavRndm < 2.) idNew = 13;
    else                    idNew = 15; 
  } else { 
    flavRndm = 3. * (flavWt - 3.) / colQ;
    if      (flavRndm <  4.) idNew = 2;
    else if (flavRndm <  8.) idNew = 4;
    else if (flavRndm <  9.) idNew = 1;
    else if (flavRndm < 10.) idNew = 3;
    else                     idNew = 5; 
  }
  double mNew  = ParticleDataTable::m0(idNew);
  double m2New = mNew*mNew;

  // Calculate kinematics dependence. Give correct mass factors for 
  // tHat, uHat defined as if massless kinematics, d(sigma)/d(Omega)
  // = beta (1 + cos^2(theta) + (1 - beta^2) sin^2(theta)).
  // Special case related to phase space form in multiple interactions.
  double sigS = 0.;
  if (sH > 4. * m2New) {
    double beta = sqrt(1. - 4. * m2New / sH);
    sigS = flavWt * beta * (2.* (tH2 + uH2) 
      + 4. * (1. - beta * beta) * tH * uH) / sH2; 
  }

  // Flavours and thereby charge/colour already fixed for MI.
  if (id12IsSet) { 
    double eNow = CoupEW::ef( abs(id1) );    
    sigS *= pow2(eNow);
    if (abs(id1) < 9) sigS /= 3.;
  }

  // Answer is proportional to number of outgoing flavours.
  return CONVERT2MB * (M_PI/sH2) * pow2(alpEM) * sigS;  

}

//*********

// Select identity, colour and anticolour.

void Sigma2ffbar2ffbarsgm::setIdColAcol() {

  // Set outgoing flavours.
  id3 = (id1 > 0) ? idNew : -idNew;
  setId( id1, id2, id3, -id3);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9 && idNew < 9) setColAcol( 1, 0, 0, 1, 2, 0, 0, 2);
  else if (abs(id1) < 9)         setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else if (idNew < 9)            setColAcol( 0, 0, 0, 0, 1, 0, 0, 1);
  else                           setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2ffbar2FFbarsgmZ class.
// Cross section f fbar -> gamma*/Z0 -> F Fbar.

//*********

// Initialize process. 
  
void Sigma2ffbar2FFbarsgmZ::initProc() {

  // Allow to pick only gamma* or Z0 part of full gamma*/Z0 expression.
  gmZmode     = Settings::mode("SigmaProcess:gmZmode");

  // Store Z0 mass and width for propagator. 
  mRes        = ParticleDataTable::m0(23);
  GammaRes    = ParticleDataTable::mWidth(23);
  m2Res       = mRes*mRes;
  GamMRat     = GammaRes / mRes;
  thetaWRat   = 1. / (16. * CoupEW::sin2thetaW() * CoupEW::cos2thetaW());

  // Store couplings of F.
  ef          = CoupEW::ef(idNew);
  vf          = CoupEW::vf(idNew);
  af          = CoupEW::af(idNew);

} 

//*********

// Initialize parton-flux object. 
  
void Sigma2ffbar2FFbarsgmZ::initFlux() {

  // Set up for f fbar initial state.
  inFluxPtr = new InFluxffbarSame();

  // Multiply by colour factor 1/3.
  inFluxPtr->weightInvCol();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2ffbar2FFbarsgmZ::sigmaHat() { 

  // Check that above threshold.
  double mHat    = sqrt(sH);
  if (mHat < m3 + m4 + MASSMARGIN) return 0.;

  // Define average F, Fbar mass so same beta. Phase space.
  double s34Avg  = 0.5 * (s3 + s4) - 0.25 * pow2(s3 - s4) / sH; 
  double m2Rat   = s34Avg / sH;
  double betaf   = sqrtpos(1. - 4. * m2Rat);
 
  // Reconstruct decay angle so can reuse 2 -> 1 cross section.
  double cosThe  = (tH - uH) / (betaf * sH); 

  // Calculate prefactors for gamma/interference/Z0 cross section terms.
  double gamProp = M_PI * pow2(alpEM) / sH2; 
  double intProp = gamProp * 2. * thetaWRat * sH * (sH - m2Res)
                 / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  double resProp = gamProp * pow2(thetaWRat * sH) 
                 / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );

  // Optionally only keep gamma* or Z0 term.
  if (gmZmode == 1) {intProp = 0.; resProp = 0.;}
  if (gmZmode == 2) {gamProp = 0.; intProp = 0.;}

  // Loop over incoming fermion flavour species.
  for (int i = 0; i < 2; ++i) {
    int idAbs = i + 1;
    if (hasLeptonBeams && i == 0) idAbs = abs(idA);
    else if (hasLeptonBeams) continue;

    // Couplings for in-flavours.
    double ei    = CoupEW::ef(idAbs);
    double vi    = CoupEW::vf(idAbs);
    double ai    = CoupEW::af(idAbs);

    // Coefficients of angular expression.
    double coefTran = ei*ei * gamProp * ef*ef + ei * vi * intProp * ef * vf
      + (vi*vi + ai*ai) * resProp * (vf*vf + pow2(betaf) * af*af);
    double coefLong = 4. * m2Rat * ( ei*ei * gamProp * ef*ef 
      + ei * vi * intProp * ef * vf + (vi*vi + ai*ai) * resProp * vf*vf );
    double coefAsym = betaf * ( ei * ai * intProp * ef * af 
      + 4. * vi * ai * resProp * vf * af );

    // Combine gamma, interference and Z0 parts.
    double sigma    = coefTran * (1. + pow2(cosThe)) 
     + coefLong * (1. - pow2(cosThe)) + 2. * coefAsym * cosThe; 

    // Send flavour-dependent parts to inFlux.
    inFluxPtr->weightInState( idAbs, -idAbs, sigma);
  }

  // Answer here is only the colour factor, rest in inFlux.
  double colF   = (idNew < 9) ? 3. * (1. + alpS / M_PI) : 1.;
  return CONVERT2MB * colF;    
   
}

//*********

// Select identity, colour and anticolour.

void Sigma2ffbar2FFbarsgmZ::setIdColAcol() {

  // Set outgoing flavours.
  id3 = (id1 > 0) ? idNew : -idNew;
  setId( id1, id2, id3, -id3);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9 && idNew < 9) setColAcol( 1, 0, 0, 1, 2, 0, 0, 2);
  else if (abs(id1) < 9)         setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else if (idNew < 9)            setColAcol( 0, 0, 0, 0, 1, 0, 0, 1);
  else                           setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//*********

// Evaluate weight for decay angles of W in top decay.

double Sigma2ffbar2FFbarsgmZ::weightDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // For top decay hand over to standard routine, else done.
  if (idNew == 6 && process[process[iResBeg].mother1()].idAbs() == 6) 
       return weightTopDecay( process, iResBeg, iResEnd);
  else return 1.; 

}

//**************************************************************************

// Sigma2ffbar2FfbarsW class.
// Cross section f fbar' -> W+- -> F fbar".

//*********

// Initialize process. 
  
void Sigma2ffbar2FfbarsW::initProc() {

  // Store W+- mass and width for propagator. 
  mRes      = ParticleDataTable::m0(24);
  GammaRes  = ParticleDataTable::mWidth(24);
  m2Res     = mRes*mRes;
  GamMRat   = GammaRes / mRes;
  thetaWRat = 1. / (12. * CoupEW::sin2thetaW());

  // Normally assume partner massless, but for t want to use b mass.
  idPartner = 0;
  if (idNew == 6) idPartner = 5;

  // Sum of CKM weights for quarks.
  V2New     = (idNew < 9) ? VCKM::V2sum(idNew) : 1.;
} 

//*********

// Initialize parton-flux object. 
  
void Sigma2ffbar2FfbarsW::initFlux() {

  // Set up for f fbar initial state.
  inFluxPtr = new InFluxffbarChg();

  // Multiply by squared CKM matrix elements and colour factor 1/3.
  inFluxPtr->weightCKM2();
  inFluxPtr->weightInvCol();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2ffbar2FfbarsW::sigmaHat() { 

  // Check that above threshold.
  double mHat   = sqrt(sH);
  if (mHat < m3 + m4 + MASSMARGIN) return 0.;

  // Phase space factors.
  double m2Rat1 = s3 / sH;
  double m2Rat2 = s4 / sH;
  double betaf  = sqrtpos( pow2(1. - m2Rat1 - m2Rat2) - 4. * m2Rat1 * m2Rat2); 
 
  // Reconstruct decay angle so can reuse 2 -> 1 cross section.
  double cosThe = (tH - uH) / (betaf * sH); 
 
  // Set up Breit-Wigner and in- and out-widths.
  double sigBW  = 9. * M_PI * pow2(alpEM * thetaWRat) 
                / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );    
  double colF   = (idNew < 9) ? 3. * (1. + alpS / M_PI) * V2New : 1.;

  // Angular dependence.
  double wt     = pow2(1. + betaf * cosThe) - pow2(m2Rat1 - m2Rat2); 

  // Done.
  return CONVERT2MB * sigBW * colF * wt;    
   
}

//*********

// Select identity, colour and anticolour.

void Sigma2ffbar2FfbarsW::setIdColAcol() {

  // Set outgoing flavours.
  id3 = idNew;
  id4 = VCKM::V2pick(idNew);
  if (idNew%2 == 0) {
    int idInUp = (abs(id1)%2 == 0) ? id1 : id2;
    if (idInUp > 0) id4 = -id4;
    else            id3 = -id3;
  } else {
    int idInDn = (abs(id1)%2 == 1) ? id1 : id2;
    if (idInDn > 0) id4 = -id4;
    else            id3 = -id3;
  }
  setId( id1, id2, id3, id4);

  // Swap tHat and uHat for fbar' f -> F f".
  if (id1 * id3 < 0) swapTU = true;   

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9 && idNew < 9) setColAcol( 1, 0, 0, 1, 2, 0, 0, 2);
  else if (abs(id1) < 9)         setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else if (idNew < 9)            setColAcol( 0, 0, 0, 0, 1, 0, 0, 1);
  else                           setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapCol12();
  if (id3 < 0) swapCol34();

}

//*********

// Evaluate weight for decay angles of W in top decay.

double Sigma2ffbar2FfbarsW::weightDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // For top decay hand over to standard routine, else done.
  if (idNew == 6 && process[process[iResBeg].mother1()].idAbs() == 6) 
       return weightTopDecay( process, iResBeg, iResEnd);
  else return 1.; 

}

//**************************************************************************

// Sigma2ffbargmZWgmZW class.
// Collects common methods for f fbar ->  gamma*/Z0/W+- gamma*/Z0/W-+.

//*********

// Calculate and store internal products.

void Sigma2ffbargmZWgmZW::setupProd( Event& process, int i1, int i2, 
  int i3, int i4, int i5, int i6) {

  // Store incoming and outgoing momenta,
  pRot[1] = process[i1].p();
  pRot[2] = process[i2].p();
  pRot[3] = process[i3].p();
  pRot[4] = process[i4].p();
  pRot[5] = process[i5].p();
  pRot[6] = process[i6].p();

  // Do random rotation to avoid accidental zeroes in HA expressions. 
  bool smallPT = false;
  do {
    smallPT = false;
    double theta = acos(2. * Rndm::flat() - 1.);
    double phi   = 2. * M_PI * Rndm::flat();
    for (int i = 1; i <= 6; ++i) { 
      pRot[i].rot( theta, phi);
      if (pRot[i].pT2() < 1e-4 * pRot[i].pAbs2()) smallPT = true;
    }
  } while (smallPT); 

  // Calculate internal products.
  for (int i = 1; i < 6; ++i) {
    for (int j = i + 1; j <= 6; ++j) { 
      hA[i][j] = 
          sqrt( (pRot[i].e() - pRot[i].pz()) * (pRot[j].e() + pRot[j].pz()) 
        / pRot[i].pT2() ) * complex( pRot[i].px(), pRot[i].py() ) 
        - sqrt( (pRot[i].e() + pRot[i].pz()) * (pRot[j].e() - pRot[j].pz()) 
        / pRot[j].pT2() ) * complex( pRot[j].px(), pRot[j].py() ); 
      hC[i][j] = conj( hA[i][j] );
      if (i <= 2) {
        hA[i][j] *= complex( 0., 1.);
        hC[i][j] *= complex( 0., 1.);
      }
      hA[j][i] = - hA[i][j]; 
      hC[j][i] = - hC[i][j]; 
    }
  }

}

//*********

// Evaluate the F function of Gunion and Kunszt.

complex Sigma2ffbargmZWgmZW::fGK(int j1, int j2, int j3, int j4, int j5, 
  int j6) {
 
  return 4. * hA[j1][j3] * hC[j2][j6] 
         * ( hA[j1][j5] * hC[j1][j4] + hA[j3][j5] * hC[j3][j4] ); 

}

//*********

// Evaluate the Xi function of Gunion and Kunszt.

double Sigma2ffbargmZWgmZW::xiGK( double tHnow, double uHnow) {

  return - 4. * s3 * s4 + tHnow * (3. * tHnow + 4. * uHnow) 
         + tHnow * tHnow * ( tHnow * uHnow / (s3 * s4)
           - 2. * (1. / s3 + 1./s4) * (tHnow + uHnow)  
           + 2. * (s3 / s4 + s4 / s3) );

}

//*********

// Evaluate the Xj function of Gunion and Kunszt.

double Sigma2ffbargmZWgmZW::xjGK( double tHnow, double uHnow) {

  return 8. * pow2(s3 + s4) - 8. * (s3 + s4) * (tHnow + uHnow)
         - 6. * tHnow * uHnow - 2. * tHnow * uHnow * ( tHnow * uHnow 
           / (s3 * s4) - 2. * (1. / s3 + 1. / s4) * (tHnow + uHnow) 
           + 2. * (s3 / s4 + s4 / s3) );

}

//**************************************************************************

// Sigma2ffbar2gmZgmZ class.
// Cross section for f fbar -> gamma*/Z0 gamma*/Z0 (f is quark or lepton). 

//*********

// Initialize process. 
  
void Sigma2ffbar2gmZgmZ::initProc() {

  // Allow to pick only gamma* or Z0 part of full gamma*/Z0 expression.
  gmZmode     = Settings::mode("SigmaProcess:gmZmode");

  // Store Z0 mass and width for propagator. 
  mRes        = ParticleDataTable::m0(23);
  GammaRes    = ParticleDataTable::mWidth(23);
  m2Res       = mRes*mRes;
  GamMRat     = GammaRes / mRes;
  thetaWRat   = 1. / (16. * CoupEW::sin2thetaW() * CoupEW::cos2thetaW());

  // Set pointer to particle properties and decay table.
  particlePtr = ParticleDataTable::particleDataPtr(23);

} 

//*********

// Initialize parton-flux object. 
  
void Sigma2ffbar2gmZgmZ::initFlux() {

  // Set up for f fbar initial state.
  inFluxPtr = new InFluxffbarSame();

  // Multiply by colour factor 1/3.
  inFluxPtr->weightInvCol();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2ffbar2gmZgmZ::sigmaHat() {

  // Cross section part common for all incoming flavours.
  double sigma = (M_PI / sH2) * pow2(alpEM) * 0.5   
               * ( (tH2 + uH2 + 2. * (s3 + s4) * sH) / (tH * uH)
                 - s3 * s4 * (1./tH2 + 1./uH2) );

  // Common coupling factors at the resonance masses
  double alpEM3 = alphaEM.alphaEM(s3);
  double alpS3  = alphaS.alphaS(s3);
  double colQ3  = 3. * (1. + alpS3 / M_PI);
  double alpEM4 = alphaEM.alphaEM(s4);
  double alpS4  = alphaS.alphaS(s4);
  double colQ4  = 3. * (1. + alpS4 / M_PI);

  // Reset quantities to sum. Declare variables in loop.
  gamSum3 = 0.;
  intSum3 = 0.;
  resSum3 = 0.;
  gamSum4 = 0.;
  intSum4 = 0.;
  resSum4 = 0.;
  double mf, m2Rat, psvec, psaxi, betaf, ef2, efvf, vf2af2, colf;

  // Loop over all Z0 decay channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    int idAbs = abs( particlePtr->decay[i].product(0) );

    // Only contributions from three fermion generations, except top.
    if ( (idAbs > 0 && idAbs < 6) || ( idAbs > 10 && idAbs < 17)) {
      mf = ParticleDataTable::m0(idAbs);

      // First Z0: check that above threshold. Phase space.
      if (m3 > 2. * mf + MASSMARGIN) {
        m2Rat = pow2(mf / m3);
        betaf = sqrtpos(1. - 4. * m2Rat); 
        psvec = betaf * (1. + 2. * m2Rat);
        psaxi  = pow3(betaf);

        // First Z0: combine phase space with couplings.
        ef2    = CoupEW::ef2(idAbs) * psvec;
        efvf   = CoupEW::efvf(idAbs) * psvec;
        vf2af2 = CoupEW::vf2(idAbs) * psvec + CoupEW::af2(idAbs) * psaxi; 
        colf   = (idAbs < 6) ? colQ3 : 1.;

        // First Z0: store sum of combinations for open outstate channels.
        if (particlePtr->decay[i].onMode() > 0) {
          gamSum3 += colf * ef2;
          intSum3 += colf * efvf;
          resSum3 += colf * vf2af2;
	}
      }

      // Second Z0: check that above threshold. Phase space.
      if (m4 > 2. * mf + MASSMARGIN) {
        m2Rat = pow2(mf / m4);
        betaf = sqrtpos(1. - 4. * m2Rat); 
        psvec = betaf * (1. + 2. * m2Rat);
        psaxi  = pow3(betaf);

        // Second Z0: combine phase space with couplings.
        ef2    = CoupEW::ef2(idAbs) * psvec;
        efvf   = CoupEW::efvf(idAbs) * psvec;
        vf2af2 = CoupEW::vf2(idAbs) * psvec + CoupEW::af2(idAbs) * psaxi; 
        colf   = (idAbs < 6) ? colQ4 : 1.;

        // Second Z0: store sum of combinations for open outstate channels.
        if (particlePtr->decay[i].onMode() > 0) {
          gamSum4 += colf * ef2;
          intSum4 += colf * efvf;
          resSum4 += colf * vf2af2;
	}
      }

    // End loop over fermions.
    }
  }

  // First Z0: calculate prefactors for gamma/interference/Z0 terms.
  gamProp3 = 4. * alpEM3 / (3. * M_PI * s3); 
  intProp3 = gamProp3 * 2. * thetaWRat * s3 * (s3 - m2Res)
           / ( pow2(s3 - m2Res) + pow2(s3 * GamMRat) );
  resProp3 = gamProp3 * pow2(thetaWRat * s3) 
           / ( pow2(s3 - m2Res) + pow2(s3 * GamMRat) );

  // First Z0: optionally only keep gamma* or Z0 term.
  if (gmZmode == 1) {intProp3 = 0.; resProp3 = 0.;}
  if (gmZmode == 2) {gamProp3 = 0.; intProp3 = 0.;}

  // Second Z0: calculate prefactors for gamma/interference/Z0 terms.
  gamProp4 = 4. * alpEM4 / (3. * M_PI * s4); 
  intProp4 = gamProp4 * 2. * thetaWRat * s4 * (s4 - m2Res)
           / ( pow2(s4 - m2Res) + pow2(s4 * GamMRat) );
  resProp4 = gamProp4 * pow2(thetaWRat * s4) 
           / ( pow2(s4 - m2Res) + pow2(s4 * GamMRat) );

  // Second Z0: optionally only keep gamma* or Z0 term.
  if (gmZmode == 1) {intProp4 = 0.; resProp4 = 0.;}
  if (gmZmode == 2) {gamProp4 = 0.; intProp4 = 0.;}

  // Loop over incoming fermion flavour species.
  for (int i = 0; i < 2; ++i) {
    int idAbs = i + 1;
    if (hasLeptonBeams && i == 0) idAbs = abs(idA);
    else if (hasLeptonBeams) continue;

    // Charge/2, left- and righthanded couplings for in-fermion.
    double ei = 0.5 * CoupEW::ef(idAbs);
    double li =       CoupEW::lf(idAbs);
    double ri =       CoupEW::rf(idAbs);

    // Combine left/right gamma, interference and Z0 parts for each Z0.
    double left3  = ei * ei * gamProp3 * gamSum3 
                  + ei * li * intProp3 * intSum3 
                  + li * li * resProp3 * resSum3;
    double right3 = ei * ei * gamProp3 * gamSum3 
                  + ei * ri * intProp3 * intSum3 
                  + ri * ri * resProp3 * resSum3;
    double left4  = ei * ei * gamProp4 * gamSum4 
                  + ei * li * intProp4 * intSum4 
                  + li * li * resProp4 * resSum4;
    double right4 = ei * ei * gamProp4 * gamSum4 
                  + ei * ri * intProp4 * intSum4 
                  + ri * ri * resProp4 * resSum4; 

    // Combine left- and right-handed couplings for the two Z0's.
    double wtFlav = left3 * left4 + right3 * right4;    

    // Correct for the running-width Z0 propagators weight in PhaseSpace. 
    wtFlav /= (runBW3 * runBW4);

    // Send flavour-dependent parts to inFlux.
    inFluxPtr->weightInState( idAbs, -idAbs, wtFlav);
  }

  // Answer for hadrons, leaving out flavour-dependent pieces stored above.
  return CONVERT2MB * sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2ffbar2gmZgmZ::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, 23, 23);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 10) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else               setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//*********

// Evaluate correlated decay flavours of the two gamma*/Z0.
// Unique complication, caused by gamma*/Z0 mix different left/right.

double Sigma2ffbar2gmZgmZ::weightDecayFlav( Event& process) {

  // Order so that fbar(1) f(2) -> f'(3) fbar'(4) f"(5) fbar"(6).
  i1 = (process[3].id() < 0) ? 3 : 4;
  i2 = 7 - i1; 
  i3 = (process[7].id() > 0) ? 7 : 8;
  i4 = 15 - i3;
  i5 = (process[9].id() > 0) ? 9 : 10;
  i6 = 19 - i5;

  // Charge/2, left- and righthanded couplings for in- and out-fermions.
  int idAbs = process[i1].idAbs();
  double ei = 0.5 * CoupEW::ef(idAbs);
  double li =       CoupEW::lf(idAbs);
  double ri =       CoupEW::rf(idAbs);
  idAbs     = process[i3].idAbs();
  double e3  = 0.5 * CoupEW::ef(idAbs);
  double l3  =       CoupEW::lf(idAbs);
  double r3  =       CoupEW::rf(idAbs);
  idAbs      = process[i5].idAbs();
  double e4  = 0.5 * CoupEW::ef(idAbs);
  double l4  =       CoupEW::lf(idAbs);
  double r4  =       CoupEW::rf(idAbs);

  // Left- and righthanded couplings combined with propagators.
  c3LL = ei * ei * gamProp3 * e3 * e3
       + ei * li * intProp3 * e3 * l3
       + li * li * resProp3 * l3 * l3; 
  c3LR = ei * ei * gamProp3 * e3 * e3
       + ei * li * intProp3 * e3 * r3
       + li * li * resProp3 * r3 * r3; 
  c3RL = ei * ei * gamProp3 * e3 * e3
       + ei * ri * intProp3 * e3 * l3
       + ri * ri * resProp3 * l3 * l3; 
  c3RR = ei * ei * gamProp3 * e3 * e3
       + ei * ri * intProp3 * e3 * r3
       + ri * ri * resProp3 * r3 * r3; 
  c4LL = ei * ei * gamProp4 * e4 * e4
       + ei * li * intProp4 * e4 * l4
       + li * li * resProp4 * l4 * l4; 
  c4LR = ei * ei * gamProp4 * e4 * e4
       + ei * li * intProp4 * e4 * r4
       + li * li * resProp4 * r4 * r4; 
  c4RL = ei * ei * gamProp4 * e4 * e4
       + ei * ri * intProp4 * e4 * l4
       + ri * ri * resProp4 * l4 * l4; 
  c4RR = ei * ei * gamProp4 * e4 * e4
       + ei * ri * intProp4 * e4 * r4
       + ri * ri * resProp4 * r4 * r4; 

  // Flavour weight and maximum.
  flavWt = (c3LL + c3LR) * (c4LL + c4LR) + (c3RL + c3RR) * (c4RL + c4RR);
  double flavWtMax = (c3LL + c3LR + c3RL + c3RR) * (c4LL + c4LR + c4RL + c4RR); 

  // Done.
  return flavWt / flavWtMax;

}

//*********

// Evaluate weight for decay angles of the two gamma*/Z0.

double Sigma2ffbar2gmZgmZ::weightDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // Two resonance decays, but with common weight.
  if (iResBeg != 5 || iResEnd != 7) return 1.;

  // Set up four-products and internal products.
  setupProd( process, i1, i2, i3, i4, i5, i6); 

  // Flip tHat and uHat if first incoming is fermion. 
  double tHres   = tH;
  double uHres   = uH;
  if (process[3].id() > 0) swap( tHres, uHres);

  // Kinematics factors (norm(x) = |x|^2). 
  double fGK135 = norm( fGK( 1, 2, 3, 4, 5, 6) / tHres
                      + fGK( 1, 2, 5, 6, 3, 4) / uHres );
  double fGK145 = norm( fGK( 1, 2, 4, 3, 5, 6) / tHres
		      + fGK( 1, 2, 5, 6, 4, 3) / uHres );
  double fGK136 = norm( fGK( 1, 2, 3, 4, 6, 5) / tHres
		      + fGK( 1, 2, 6, 5, 3, 4) / uHres );
  double fGK146 = norm( fGK( 1, 2, 4, 3, 6, 5) / tHres
		      + fGK( 1, 2, 6, 5, 4, 3) / uHres );
  double fGK253 = norm( fGK( 2, 1, 5, 6, 3, 4) / tHres
      	              + fGK( 2, 1, 3, 4, 5, 6) / uHres );
  double fGK263 = norm( fGK( 2, 1, 6, 5, 3, 4) / tHres
      	              + fGK( 2, 1, 3, 4, 6, 5) / uHres );
  double fGK254 = norm( fGK( 2, 1, 5, 6, 4, 3) / tHres
      	  	      + fGK( 2, 1, 4, 3, 5, 6) / uHres );
  double fGK264 = norm( fGK( 2, 1, 6, 5, 4, 3) / tHres
      	    	      + fGK( 2, 1, 4, 3, 6, 5) / uHres );

  // Weight and maximum.
  double wt     = c3LL * c4LL * fGK135 + c3LR * c4LL * fGK145
                + c3LL * c4LR * fGK136 + c3LR * c4LR * fGK146
                + c3RL * c4RL * fGK253 + c3RR * c4RL * fGK263
                + c3RL * c4RR * fGK254 + c3RR * c4RR * fGK264;
  double wtMax  = 16. * s3 * s4 * flavWt 
    * ( (tHres*tHres + uHres*uHres + 2. * sH * (s3 + s4)) / (tHres * uHres)
      - s3 * s4 * (1. / (tHres*tHres) + 1. / (uHres*uHres)) ); 

  // Done.
  return wt / wtMax;

}

//**************************************************************************

// Sigma2ffbar2ZW class.
// Cross section for f fbar' -> W+ W- (f is quark or lepton). 

//*********

// Initialize process. 
  
void Sigma2ffbar2ZW::initProc() {

  // Store W+- mass and width for propagator. 
  mW   = ParticleDataTable::m0(24);
  widW = ParticleDataTable::mWidth(24);
  mWS  = mW*mW;
  mwWS = pow2(mW * widW);

  // Left-handed couplings for up/nu- and down/e-type quarks.
  lun   = (hasLeptonBeams) ? CoupEW::lf(12) : CoupEW::lf(2);
  lde   = (hasLeptonBeams) ? CoupEW::lf(11) : CoupEW::lf(1); 

  // Common weak coupling factor.
  sin2thetaW = CoupEW::sin2thetaW();
  cos2thetaW = CoupEW::cos2thetaW();
  thetaWRat  = 1. / (4. * cos2thetaW);  
  thetaWpt   = (9. - 8. * sin2thetaW) / 4.;
  thetaWmm   = (8. * sin2thetaW - 6.) / 4.;

} 

//*********

// Initialize parton-flux object. 
  
void Sigma2ffbar2ZW::initFlux() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxffbarChg();

  // Multiply by squared CKM matrix elements and colour factor 1/3.
  inFluxPtr->weightCKM2();
  inFluxPtr->weightInvCol();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2ffbar2ZW::sigmaHat() {

  // Evaluate cross section.
  double resBW = 1. / (pow2(sH - mWS) + mwWS);
  double sigma = (M_PI / sH2) * 0.5 * pow2(alpEM / sin2thetaW);
  sigma *= sH * resBW * (thetaWpt * pT2 + thetaWmm * (s3 + s4))
    + (sH - mWS) * resBW * sH * (pT2 - s3 - s4) * (lun / tH - lde / uH)
    + thetaWRat * sH * pT2 * ( lun*lun / tH2 + lde*lde / uH2 )
    + 2. * thetaWRat * sH * (s3 + s4) * lun * lde / (tH * uH);  

  // Answer, leaving out flavour-dependent pieces stored above.
  // (Protect against slightly negative cross sections, probably 
  // caused by addition of width to the W propagator.)
  return CONVERT2MB * max(0., sigma);    

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

//*********

// Evaluate weight for Z0 and W+- decay angles.

double Sigma2ffbar2ZW::weightDecay( Event& process, int iResBeg, int iResEnd) {

  // Two resonance decays, but with common weight.
  if (iResBeg != 5 || iResEnd != 7) return 1.;

  // Order so that fbar(1) f(2) -> f'(3) fbar'(4) f"(5) fbar"(6)
  // with f' fbar' from W+- and f" fbar" from Z0 (note flip Z0 <-> W+-).
  int i1 = (process[3].id() < 0) ? 3 : 4;
  int i2 = 7 - i1; 
  int i3 = (process[9].id() > 0) ? 9 : 10;
  int i4 = 19 - i3;
  int i5 = (process[7].id() > 0) ? 7 : 8;
  int i6 = 15 - i5;

  // Set up four-products and internal products.
  setupProd( process, i1, i2, i3, i4, i5, i6); 

  // Swap tHat and uHat if incoming fermion is downtype. 
  double tHres   = tH;
  double uHres   = uH;
  if (process[i2].id()%2 == 1) swap( tHres, uHres);

  //  Couplings of incoming (anti)fermions and outgoing from Z0.
  int idAbs     = process[i1].idAbs();
  double ai     = CoupEW::af(idAbs); 
  double li1    = CoupEW::lf(idAbs); 
  idAbs         = process[i2].idAbs();
  double li2    = CoupEW::lf(idAbs); 
  idAbs         = process[i5].idAbs();
  double l4     = CoupEW::lf(idAbs); 
  double r4     = CoupEW::rf(idAbs); 

  // W propagator/interference factor.
  double Wint   = cos2thetaW * (sH - mWS) / (pow2(sH - mWS) + mwWS);

  // Combinations of couplings and kinematics (norm(x) = |x|^2).
  double aWZ    = li2 / tHres - 2. * Wint * ai;
  double bWZ    = li1 / uHres + 2. * Wint * ai;
  double fGK135 = norm( aWZ * fGK( 1, 2, 3, 4, 5, 6) 
		      + bWZ * fGK( 1, 2, 5, 6, 3, 4) );
  double fGK136 = norm( aWZ * fGK( 1, 2, 3, 4, 6, 5) 
		      + bWZ * fGK( 1, 2, 6, 5, 3, 4) );
  double xiT    = xiGK( tHres, uHres);
  double xiU    = xiGK( uHres, tHres);
  double xjTU   = xjGK( tHres, uHres);

  // Weight and maximum weight.
  double wt     = l4*l4 * fGK135 + r4*r4 * fGK136;
  double wtMax  = 4. * s3 * s4 * (l4*l4 + r4*r4) 
                * (aWZ * aWZ * xiT + bWZ * bWZ * xiU + aWZ * bWZ * xjTU);

  // Done.
  return wt / wtMax;

}

//**************************************************************************

// Sigma2ffbar2WW class.
// Cross section for f fbar -> W- W+ (f is quark or lepton). 

//*********

// Initialize process. 
  
void Sigma2ffbar2WW::initProc() {

  // Store Z0 mass and width for propagator. Common coupling factor.
  mZ        = ParticleDataTable::m0(23);
  widZ      = ParticleDataTable::mWidth(23);
  mZS       = mZ*mZ;
  mwZS      = pow2(mZ * widZ);
  thetaWRat = 1. / (4. * CoupEW::sin2thetaW());  

} 

//*********

// Initialize parton-flux object. 
  
void Sigma2ffbar2WW::initFlux() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxffbarSame();

  // Multiply by colour factor 1/3.
  inFluxPtr->weightInvCol();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2ffbar2WW::sigmaHat() {

  // Cross section part common for all incoming flavours.
  double sigma   = (M_PI / sH2) * pow2(alpEM);

  // Z0 propagator and gamma*/Z0 interference.
  double Zprop   = sH2 / (pow2(sH - mZS) + mwZS);
  double Zint    = Zprop * (1. - mZS / sH);

  // Common coupling factors (g = gamma*, Z = Z0, f = t-channel fermion).
  double cgg     = 0.5;
  double cgZ     = thetaWRat * Zint;
  double cZZ     = 0.5 * pow2(thetaWRat) * Zprop;
  double cfg     = thetaWRat;
  double cfZ     = pow2(thetaWRat) * Zint;
  double cff     = pow2(thetaWRat);

  // Kinematical functions.   
  double rat34   = sH * (2. * (s3 + s4) + pT2) / (s3 * s4);
  double lambdaS = pow2(sH - s3 - s4) - 4. * s3 * s4;
  double intA    = (sH - s3 - s4) * rat34 / sH;
  double intB    = 4. * (s3 + s4 - pT2);
  double gSS     = (lambdaS * rat34 + 12. * sH * pT2) / sH2;
  double gTT     = rat34 + 4. * sH * pT2 / tH2;
  double gST     = intA + intB / tH;
  double gUU     = rat34 + 4. * sH * pT2 / uH2;
  double gSU     = intA + intB / uH;   

  // Loop over incoming fermion flavour species.
  for (int i = 0; i < 2; ++i) {
    int idAbs = i + 1;
    if (hasLeptonBeams && i == 0) idAbs = abs(idA);
    else if (hasLeptonBeams) continue;
 
    // Flavour-specific couplings.
    double ei = CoupEW::ef(idAbs);
    double vi = CoupEW::vf(idAbs); 
    double ai = CoupEW::af(idAbs); 

    // Combine, with different cases for up- and down-type in-flavours.
    double wtFlav = (idAbs%2 == 1)
      ? (cgg * ei*ei + cgZ * ei * vi + cZZ * (vi*vi + ai*ai)) * gSS
        + (cfg * ei + cfZ * (vi + ai)) * gST + cff * gTT
      : (cgg * ei*ei + cgZ * ei * vi + cZZ * (vi*vi + ai*ai)) * gSS
        - (cfg * ei + cfZ * (vi + ai)) * gSU + cff * gUU;

    // Send flavour-dependent parts to inFlux.
    inFluxPtr->weightInState( idAbs, -idAbs, wtFlav );
  }

  // Answer, leaving out flavour-dependent pieces stored above.
  return CONVERT2MB * sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2ffbar2WW::setIdColAcol() {

  // Always order W- W+, i.e. W- first.
  setId( id1, id2, -24, 24);

  // tHat is defined between (f, W-) or (fbar, W+), 
  if (id1 < 0) swapTU = true;

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 10) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else               setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//*********

// Evaluate weight for W+ and W- decay angles.
  
double Sigma2ffbar2WW::weightDecay( Event& process, int iResBeg, int iResEnd) {

  // Two resonance decays, but with common weight.
  if (iResBeg != 5 || iResEnd != 7) return 1.;

  // Order so that fbar(1) f(2) -> f'(3) fbar'(4) f"(5) fbar"(6).
  // with f' fbar' from W- and f" fbar" from W+.
  int i1 = (process[3].id() < 0) ? 3 : 4;
  int i2 = 7 - i1; 
  int i3 = (process[7].id() > 0) ? 7 : 8;
  int i4 = 15 - i3;
  int i5 = (process[9].id() > 0) ? 9 : 10;
  int i6 = 19 - i5;

  // Set up four-products and internal products.
  setupProd( process, i1, i2, i3, i4, i5, i6); 

  // tHat and uHat of fbar f -> W- W+ opposite to previous convention. 
  double tHres   = uH;
  double uHres   = tH;

  //  Couplings of incoming (anti)fermion.
  int idAbs     = process[i1].idAbs();
  double ai     = CoupEW::af(idAbs); 
  double li     = CoupEW::lf(idAbs); 
  double ri     = CoupEW::rf(idAbs); 

  // gamma*/Z0 propagator/interference factor.
  double Zint   = mZS * (sH - mZS) / (pow2(sH - mZS) + mwZS);

  // Combinations of couplings and kinematics (norm(x) = |x|^2).
  double dWW    = (li * Zint + ai) / sH;
  double aWW    = dWW + 0.5 * (ai + 1.) / tHres; 
  double bWW    = dWW + 0.5 * (ai - 1.) / uHres;
  double cWW    = ri * Zint / sH; 
  double fGK135 = norm( aWW * fGK( 1, 2, 3, 4, 5, 6) 
		      - bWW * fGK( 1, 2, 5, 6, 3, 4) );
  double fGK253 = norm( cWW * ( fGK( 2, 1, 5, 6, 3, 4) 
                              - fGK( 2, 1, 3, 4, 5, 6) ) );
  double xiT    = xiGK( tHres, uHres);
  double xiU    = xiGK( uHres, tHres);
  double xjTU   = xjGK( tHres, uHres);

  // Weight and maximum weight.
  double wt     = fGK135 + fGK253;
  double wtMax  = 4. * s3 * s4
                * ( aWW * aWW * xiT + bWW * bWW * xiU - aWW * bWW * xjTU
	          + cWW * cWW * (xiT + xiU - xjTU) );

  // Done.
  return wt / wtMax;
}

//**************************************************************************

// Sigma2ffbargmZggm class.
// Collects common methods for f fbar -> gamma*/Z0 g/gamma and permutations.

//*********

// Initialize process.
  
void Sigma2ffbargmZggm::initProc() {

  // Allow to pick only gamma* or Z0 part of full gamma*/Z0 expression.
  gmZmode     = Settings::mode("SigmaProcess:gmZmode");

  // Store Z0 mass and width for propagator. 
  mRes        = ParticleDataTable::m0(23);
  GammaRes    = ParticleDataTable::mWidth(23);
  m2Res       = mRes*mRes;
  GamMRat     = GammaRes / mRes;
  thetaWRat   = 1. / (16. * CoupEW::sin2thetaW() * CoupEW::cos2thetaW());

  // Set pointer to particle properties and decay table.
  particlePtr = ParticleDataTable::particleDataPtr(23);

}

//*********

// Evaluate sum of flavour couplings times phase space.

void Sigma2ffbargmZggm::flavSum() {

  // Coupling factors for Z0 subsystem. 
  double alpSZ = alphaS.alphaS(s3);
  double colQZ = 3. * (1. + alpSZ / M_PI);

  // Reset quantities to sum. Declare variables in loop.
  gamSum = 0.;
  intSum = 0.;
  resSum = 0.;
  double mf, m2Rat, psvec, psaxi, betaf, ef2, efvf, vf2af2, colf;

  // Loop over all Z0 decay channels. 
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    int idAbs = abs( particlePtr->decay[i].product(0) );

    // Only contributions from three fermion generations, except top.
    if ( (idAbs > 0 && idAbs < 6) || ( idAbs > 10 && idAbs < 17)) {
      mf = ParticleDataTable::m0(idAbs);

      // Check that above threshold. Phase space.
      if (m3 > 2. * mf + MASSMARGIN) {
        m2Rat = pow2(mf / m3);
        betaf = sqrtpos(1. - 4. * m2Rat); 
        psvec = betaf * (1. + 2. * m2Rat);
        psaxi  = pow3(betaf);

        // Combine phase space with couplings.
        ef2    = CoupEW::ef2(idAbs) * psvec;
        efvf   = CoupEW::efvf(idAbs) * psvec;
        vf2af2 = CoupEW::vf2(idAbs) * psvec + CoupEW::af2(idAbs) * psaxi; 
        colf   = (idAbs < 6) ? colQZ : 1.;

        // Store sum of combinations. For outstate only open channels.
        if (particlePtr->decay[i].onMode() > 0) {
          gamSum += colf * ef2;
          intSum += colf * efvf;
          resSum += colf * vf2af2;
	}

      // End loop over fermions.
      }
    }
  }

  // Done. Return values in gamSum, intSum and resSum.

}

//*********

// Calculate common parts of gamma/interference/Z0 propagator terms.
  
void Sigma2ffbargmZggm::propTerm() {

  // Calculate prefactors for gamma/interference/Z0 cross section terms.
  gamProp = 4. * alpEM / (3. * M_PI * s3); 
  intProp = gamProp * 2. * thetaWRat * s3 * (s3 - m2Res)
          / ( pow2(s3 - m2Res) + pow2(s3 * GamMRat) );
  resProp = gamProp * pow2(thetaWRat * s3) 
          / ( pow2(s3 - m2Res) + pow2(s3 * GamMRat) );

  // Optionally only keep gamma* or Z0 term.
  if (gmZmode == 1) {intProp = 0.; resProp = 0.;}
  if (gmZmode == 2) {gamProp = 0.; intProp = 0.;}

}

//*********

// Evaluate weight for gamma*/Z0 decay angle.
  
double Sigma2ffbargmZggm::weightDecay( Event& process, int iResBeg, int iResEnd) {

  // Z should sit in entry 5 and one more parton in entry 6.
  if (iResBeg != 5 || iResEnd != 7) return 1.;

  // In an outgoing sense fermions are labelled f(1) fbar(2) f'(3) fbar'(4)
  // where f' fbar' come from gamma*/Z0 decay. 
  int i1, i2;
  int i3 = (process[7].id() > 0) ? 7 : 8;
  int i4 = 15 - i3;

  // Order so that fbar(1) f(2) -> gamma*/Z0 g/gamma.
  if (process[3].idAbs() < 20 && process[4].idAbs() < 20) {
    i1 = (process[3].id() < 0) ? 3 : 4;
    i2 = 7 - i1; 

  // Order so that f(2)/fbar(1)  g/gamma -> f(1)/fbar(2) f'(3) gamma*/Z0.
  } else if (process[3].idAbs() < 20) {
    i1 = (process[3].id() < 0) ? 3 : 6;
    i2 = 9 - i1; 
  } else {
    i1 = (process[4].id() < 0) ? 4 : 6;
    i2 = 10 - i1; 
  }

  // Charge/2, left- and righthanded couplings for in- and out-fermion.
  int id1Abs   = process[i1].idAbs();
  double ei    = 0.5 * CoupEW::ef(id1Abs);
  double li    =       CoupEW::lf(id1Abs);
  double ri    =       CoupEW::rf(id1Abs);
  int id3Abs   = process[i3].idAbs();
  double ef    = 0.5 * CoupEW::ef(id3Abs);
  double lf    =       CoupEW::lf(id3Abs);
  double rf    =       CoupEW::rf(id3Abs);

  // Combinations of left/right for in/out, gamma*/interference/Z0.
  double clilf = ei*ei * gamProp * ef*ef + ei*li * intProp * ef*lf 
               + li*li * resProp * lf*lf;
  double clirf = ei*ei * gamProp * ef*ef + ei*li * intProp * ef*rf 
               + li*li * resProp * rf*rf;
  double crilf = ei*ei * gamProp * ef*ef + ei*ri * intProp * ef*lf 
               + ri*ri * resProp * lf*lf;
  double crirf = ei*ei * gamProp * ef*ef + ei*ri * intProp * ef*rf 
               + ri*ri * resProp * rf*rf;

  // Evaluate four-vector products.
  double p13   = process[i1].p() * process[i3].p(); 
  double p14   = process[i1].p() * process[i4].p(); 
  double p23   = process[i2].p() * process[i3].p(); 
  double p24   = process[i2].p() * process[i4].p(); 

  // Calculate weight and its maximum.
  double wt    = (clilf + crirf) * (p13*p13 + p24*p24)
               + (clirf + crilf) * (p14*p14 + p23*p23) ;
  double wtMax = (clilf + clirf + crilf + crirf)  
               * (pow2(p13 + p14) + pow2(p23 + p24));
 
  // Done.
  return (wt / wtMax);

}

//**************************************************************************

// Sigma2qqbar2gmZg class.
// Cross section for q qbar -> gamma*/Z0 g. 

//*********

// Initialize parton-flux object. 
  
void Sigma2qqbar2gmZg::initFlux() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxqqbarSame();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qqbar2gmZg::sigmaHat() {

  // Cross section part common for all incoming flavours.
  double sigma = (M_PI / sH2) * (alpEM * alpS) 
    * (2./9.) * (tH2 + uH2 + 2. * sH * s3) / (tH * uH);

  // Calculate flavour sums for final state.
  flavSum();

  // Calculate prefactors for gamma/interference/Z0 cross section terms.
  propTerm();  

  // Loop over incoming fermion flavour species.
  for (int i = 0; i < 2; ++i) {
    int idAbs = i + 1;

    // Combine gamma, interference and Z0 parts.
    double wtFlav = CoupEW::ef2(idAbs)    * gamProp * gamSum 
                  + CoupEW::efvf(idAbs)   * intProp * intSum
                  + CoupEW::vf2af2(idAbs) * resProp * resSum;

    // Correct for the running-width Z0 propagater weight in PhaseSpace. 
    wtFlav /= runBW3;

    // Send flavour-dependent parts to inFlux.
    inFluxPtr->weightInState( idAbs, -idAbs, wtFlav);
  }

  // Flavour -independent part of answer.
  return CONVERT2MB * sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2qqbar2gmZg::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, 23, 21);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 2, 0, 0, 1, 2);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2qg2gmZq class.
// Cross section for q g -> gamma*/Z0 q. 

//*********

// Initialize parton-flux object. 
  
void Sigma2qg2gmZq::initFlux() {

  // Set up for q g initial state.
  inFluxPtr = new InFluxqg();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qg2gmZq::sigmaHat() {

  // Cross section part common for all incoming flavours.
  double sigma = (M_PI / sH2) * (alpEM * alpS) 
    * (1./12.) * (sH2 + uH2 + 2. * tH * s3) / (-sH * uH);

  // Calculate flavour sums for final state.
  flavSum();

  // Calculate prefactors for gamma/interference/Z0 cross section terms.
  propTerm();  

  // Loop over incoming fermion flavour species.
  for (int i = 0; i < 2; ++i) {
    int idAbs = i + 1;

    // Combine gamma, interference and Z0 parts.
    double wtFlav = CoupEW::ef2(idAbs)    * gamProp * gamSum 
                  + CoupEW::efvf(idAbs)   * intProp * intSum
                  + CoupEW::vf2af2(idAbs) * resProp * resSum;

    // Correct for the running-width Z0 propagater weight in PhaseSpace. 
    wtFlav /= runBW3;

    // Send flavour-dependent parts to inFlux.
    inFluxPtr->weightInState( idAbs, 21, wtFlav);
  }

  // Answer.
  return CONVERT2MB * sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2qg2gmZq::setIdColAcol() {

  // Flavour set up for q g -> gamma*/Z0 q.
  int idq = (id2 == 21) ? id1 : id2;
  setId( id1, id2, 23, idq);

  // tH defined between f and f': must swap tHat <-> uHat if q g in.
  swapTU = (id2 == 21); 

  // Colour flow topologies. Swap when antiquarks.
  if (id2 == 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else           setColAcol( 2, 1, 1, 0, 0, 0, 2, 0);
  if (idq < 0) swapColAcol();

}

//**************************************************************************

// Sigma2ffbar2gmZgm class.
// Cross section for f fbar -> gamma*/Z0 gamma. 

//*********

// Initialize parton-flux object. 
  
void Sigma2ffbar2gmZgm::initFlux() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxffbarSame();

  // Multiply by colour factor 1/3.
  inFluxPtr->weightInvCol();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2ffbar2gmZgm::sigmaHat() {

  // Cross section part common for all incoming flavours.
  double sigma = (M_PI / sH2) * (alpEM*alpEM) 
    * 0.5 * (tH2 + uH2 + 2. * sH * s3) / (tH * uH);

  // Calculate flavour sums for final state.
  flavSum();

  // Calculate prefactors for gamma/interference/Z0 cross section terms.
  propTerm();  

  // Loop over incoming fermion flavour species.
  for (int i = 0; i < 2; ++i) {
    int idAbs = i + 1;
    if (hasLeptonBeams && i == 0) idAbs = abs(idA);
    else if (hasLeptonBeams) continue;

    // Combine gamma, interference and Z0 parts.
    double wtFlav = CoupEW::ef2(idAbs)
                * ( CoupEW::ef2(idAbs)    * gamProp * gamSum 
                  + CoupEW::efvf(idAbs)   * intProp * intSum
                  + CoupEW::vf2af2(idAbs) * resProp * resSum );

    // Correct for the running-width Z0 propagater weight in PhaseSpace. 
    wtFlav /= runBW3;

    // Send flavour-dependent parts to inFlux.
    inFluxPtr->weightInState( idAbs, -idAbs, wtFlav);
  }

  // Answer.
  return CONVERT2MB * sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2ffbar2gmZgm::setIdColAcol() {

  // Flavours trivial.
  setId( id1, id2, 23, 22);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 10) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else               setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2fgm2gmZf class.
// Cross section for f gamma -> gamma*/Z0 f'. 

//*********

// Initialize parton-flux object. 
  
void Sigma2fgm2gmZf::initFlux() {

  // Set up for f gamma initial state.
  inFluxPtr = new InFluxfgm();

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2fgm2gmZf::sigmaHat() {

  // Cross section part common for all incoming flavours.
  double sigma = (M_PI / sH2) * (alpEM*alpEM) 
    * 0.5 * (sH2 + uH2 + 2. * tH * s3) / (- sH * uH);

  // Calculate flavour sums for final state.
  flavSum();

  // Calculate prefactors for gamma/interference/Z0 cross section terms.
  propTerm();  

  // Loop over incoming fermion flavour species.
  for (int i = 0; i < 2; ++i) {
    int idAbs = i + 1;

    // Combine gamma, interference and Z0 parts.
    double wtFlav = CoupEW::ef2(idAbs)
                * ( CoupEW::ef2(idAbs)    * gamProp * gamSum 
                  + CoupEW::efvf(idAbs)   * intProp * intSum
                  + CoupEW::vf2af2(idAbs) * resProp * resSum );

    // Correct for the running-width Z0 propagater weight in PhaseSpace. 
    wtFlav /= runBW3;

    // Send flavour-dependent parts to inFlux.
    inFluxPtr->weightInState( idAbs, 22, wtFlav);
  }

  // Answer.
  return CONVERT2MB * sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2fgm2gmZf::setIdColAcol() {

  // Flavour set up for q gamma -> gamma*/Z0 q.
  int idq = (id2 == 22) ? id1 : id2;
  setId( id1, id2, 23, idq);

  // tH defined between f and f': must swap tHat <-> uHat if q gamma in.
  swapTU = (id2 == 22); 

  // Colour flow topologies. Swap when antiquarks.
  if      (abs(id1) < 10) setColAcol( 1, 0, 0, 0, 0, 0, 1, 0);
  else if (abs(id2) < 10) setColAcol( 0, 0, 1, 0, 0, 0, 1, 0);
  else                    setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (idq < 0) swapColAcol();

}

//**************************************************************************

// Sigma2ffbarWggm class.
// Collects common methods for f fbar -> W+- g/gamma and permutations.

//*********

// Evaluate weight for W+- decay angle.
  
double Sigma2ffbarWggm::weightDecay( Event& process, int iResBeg, int iResEnd) {

  // W should sit in entry 5 and one more parton in entry 6.
  if (iResBeg != 5 || iResEnd != 7) return 1.;

  // In an outgoing sense fermions are labelled f(1) fbar(2) f'(3) fbar'(4)
  // where f' fbar' come from W+- decay. 
  int i1, i2;
  int i3 = (process[7].id() > 0) ? 7 : 8;
  int i4 = 15 - i3;

  // Order so that fbar(1) f(2) -> W+- g/gamma.
  if (process[3].idAbs() < 20 && process[4].idAbs() < 20) {
    i1 = (process[3].id() < 0) ? 3 : 4;
    i2 = 7 - i1; 

  // Order so that f(2)/fbar(1)  g/gamma -> f(1)/fbar(2) f'(3) W+-.
  } else if (process[3].idAbs() < 20) {
    i1 = (process[3].id() < 0) ? 3 : 6;
    i2 = 9 - i1; 
  } else {
    i1 = (process[4].id() < 0) ? 4 : 6;
    i2 = 10 - i1; 
  }

  // Evaluate four-vector products.
  double p13 = process[i1].p() * process[i3].p(); 
  double p14 = process[i1].p() * process[i4].p(); 
  double p23 = process[i2].p() * process[i3].p(); 
  double p24 = process[i2].p() * process[i4].p(); 

  // Calculate weight and its maximum.
  double wt    = pow2(p13) + pow2(p24);
  double wtMax = pow2(p13 + p14) + pow2(p23 + p24);
 
  // Done.
  return (wt / wtMax);

}

//**************************************************************************

// Sigma2qqbar2Wg class.
// Cross section for q qbar' -> W+- g. 

//*********

// Initialize parton-flux object. 
  
void Sigma2qqbar2Wg::initFlux() {

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

// Initialize parton-flux object. 
  
void Sigma2qg2Wq::initFlux() {

  // Set up for q g initial state.
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
  int idq           = (id2 == 21) ? id1 : id2;
  int sign          = 1 - 2 * (abs(idq)%2);
  if (idq < 0) sign = -sign;
  id4 = VCKM::V2pick(idq);

  // Flavour set up for q g -> W q.
  setId( id1, id2, 24 * sign, id4);

  // tH defined between f and f': must swap tHat <-> uHat if q g in.
  swapTU = (id2 == 21); 

  // Colour flow topologies. Swap when antiquarks.
  if (id2 == 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else           setColAcol( 2, 1, 1, 0, 0, 0, 2, 0);
  if (idq < 0) swapColAcol();

}

//**************************************************************************

// Sigma2ffbar2Wgm class.
// Cross section for f fbar' -> W+- gamma. 

//*********

// Initialize parton-flux object. 
  
void Sigma2ffbar2Wgm::initFlux() {

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

  // Extrafactor different for e nu and q qbar' instate.
  double chgUp = (hasLeptonBeams) ? 0. : 2./3.;
  sigma *= pow2( chgUp - tH / (tH + uH) );

  // Answer.
  return CONVERT2MB * sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2ffbar2Wgm::setIdColAcol() {

  // Sign of outgoing W.
  int sign          = 1 - 2 * (abs(id1)%2);
  if (id1 < 0) sign = -sign;
  setId( id1, id2, 24 * sign, 22);

  // tH defined between (f,W-) or (fbar',W+).
  swapTU = (sign * id1 > 0);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 10) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else               setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2fgm2Wf class.
// Cross section for f gamma -> W+- f'. 

//*********

// Initialize parton-flux object. 
  
void Sigma2fgm2Wf::initFlux() {

  // Set up for f gamma initial state.
  inFluxPtr = new InFluxfgm();

  // Multiply by sum of relevant squared CKM matrix elements.
  inFluxPtr->weightCKM2sum(1);

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2fgm2Wf::sigmaHat() {

  // Cross section part common for all incoming flavours.
  double sigma = (M_PI / sH2) * (alpEM*alpEM / CoupEW::sin2thetaW())
    * 0.5 * (sH2 + uH2 + 2. * tH * s3) / (pT2 * s3 - sH * uH);

  // Extrafactor dependent on charge of incoming fermion.
  if (hasLeptonBeams) sigma *= pow2( 1. - sH / (sH + uH) );
  else {
    inFluxPtr->weightInState( 1, 22, pow2( 1./3. - sH / (sH + uH) ) );
    inFluxPtr->weightInState( 2, 22, pow2( 2./3. - sH / (sH + uH) ) );
  }

  // Answer.
  return CONVERT2MB * sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2fgm2Wf::setIdColAcol() {

  // Sign of outgoing W. Flavour of outgoing fermion.
  int idq           = (id2 == 22) ? id1 : id2;
  int sign          = 1 - 2 * (abs(idq)%2);
  if (idq < 0) sign = -sign;
  id4 = VCKM::V2pick(idq);

  // Flavour set up for q gamma -> W q.
  setId( id1, id2, 24 * sign, id4);

  // tH defined between f and f': must swap tHat <-> uHat if q gamma in.
  swapTU = (id2 == 22); 

  // Colour flow topologies. Swap when antiquarks.
  if      (abs(id1) < 10) setColAcol( 1, 0, 0, 0, 0, 0, 1, 0);
  else if (abs(id2) < 10) setColAcol( 0, 0, 1, 0, 0, 0, 1, 0);
  else                    setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (idq < 0) swapColAcol();

}

//**************************************************************************

} // end namespace Pythia8
