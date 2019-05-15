// Function definitions (not found in the header) for the
// SigmaProcess class.
// Copyright C 2006 Torbjorn Sjostrand

#include "SigmaProcess.h"

namespace Pythia8 {

//**************************************************************************

// The SigmaProcess class.

//*********
 
// Definitions of static variables.
// (Values will be overwritten in initStatic call, so are purely dummy.)

int SigmaProcess::alphaSorder = 1;
int SigmaProcess::nQuark = 5;
double SigmaProcess::alphaSvalue = 0.1265;
double SigmaProcess::alphaEM =  0.00729735;

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// The sum of outgoing masses must not be too close to the cm energy.
const double SigmaProcess::MASSMARGIN = 0.1;

// Conversion of GeV^{-2} to mb for cross section.
const double SigmaProcess::CONVERT2MB = 0.389380; 

//*********

// Initialize static data members.

void SigmaProcess::initStatic() {

  // Parameters of alphaStrong generation .
  alphaSvalue = Settings::parameter("SigmaProcess:alphaSvalue");
  alphaSorder = Settings::mode("SigmaProcess:alphaSorder");

  // Maximum new quark flavour, alphaEM and conversion Gev^{-2} -> mb.
  nQuark = Settings::mode("SigmaProcess:nQuark");
  alphaEM = Settings::parameter("StandardModel:alphaEMfix"); 

}

//*********

// Initialize alphaStrong. 

void SigmaProcess::init() {

  // Initialize alphaStrong generation.
  alphaScalc.init( alphaSvalue, alphaSorder); 

}

//*********

// Derive uHat from sHat and tHat input; alpha_strong calculated.

void SigmaProcess::setupShTh(int id1In, int id2In, double sHatIn, 
  double tHatIn) {

  // Default values: same flavours out as in, no colours.
  id1 = id1In;
  id2 = id2In;
  id3 = id1;
  id4 = id2;
  flow = 0;

  // Standard Mandelstam variables and their squares.
  sH = sHatIn;
  tH = tHatIn;
  uH = - sH - tH;
  sH2 = sH*sH;
  tH2 = tH*tH;
  uH2 = uH*uH;

  // Calculate running alpha_strong, using scale pT2 = tH*uH/sH.
  double pT2 = tH * uH / sH;
  alphaS = alphaScalc.alphaS(pT2);  
  alphaS2 = alphaS * alphaS;  

}

//*********

// Derive uHat from sHat and tHat input; alpha_strong input.

void SigmaProcess::setupShThAs(int id1In, int id2In, double sHatIn, 
  double tHatIn, double alphaSIn) {

  // Default values: same flavours out as in, no colours.
  id1 = id1In;
  id2 = id2In;
  id3 = id1;
  id4 = id2;
  flow = 0;

  // Standard Mandelstam variables and their squares.
  sH = sHatIn;
  tH = tHatIn;
  uH = - sH - tH;
  sH2 = sH*sH;
  tH2 = tH*tH;
  uH2 = uH*uH;

  // alpha_strong given as input.
  alphaS = alphaSIn;  
  alphaS2 = alphaS * alphaS;  

}

//*********

// g + g -> g + g.

double SigmaProcess::gg2gg() {

  // Calculate.  
  double sigTS = (9./4.) * (tH2/sH2 + 2.*tH/sH + 3. + 2.*sH/tH + sH2/tH2);
  double sigUS = (9./4.) * (uH2/sH2 + 2.*uH/sH + 3. + 2.*sH/uH + sH2/uH2);
  double sigTU = (9./4.) * (tH2/uH2 + 2.*tH/uH + 3. + 2.*uH/tH + uH2/tH2);
  double sigSum = sigTS + sigUS + sigTU;

  // Three colour flow topologies, each with two orientations.
  double sigRand = sigSum * Rndm::flat();
  if (sigRand < sigTS) flow = (Rndm::flat() < 0.5) ? 12231443 : 12314234;
  else if (sigRand < sigTS + sigUS) flow = (Rndm::flat() < 0.5) 
    ? 12313442 : 12234314;
  else flow = (Rndm::flat() < 0.5) ?  12341432 : 12343214;

  // Answer contains factor 1/2 from identical gluons.
  return CONVERT2MB * (M_PI/sH2) * alphaS2 * 0.5 * sigSum;  

}

//*********

// q + g -> q + g.

double SigmaProcess::qg2qg() {
  
  // Calculate.
  double sigTS = uH2/tH2 - (4./9.) * uH/sH;
  double sigTU = sH2/tH2 - (4./9.) * sH/uH;
  double sigSum = sigTS + sigTU;
  flow = (sigSum * Rndm::flat() < sigTS) ? 10213023 : 10232013;

  return CONVERT2MB * (M_PI/sH2) * alphaS2 * sigSum;  

}

//*********

// q1 + q2 -> q1 + q2, where q1 and q2 are different.

double SigmaProcess::qq2qqDiff() {
  
  // Calculate.
  double sigT = (4./9.) * (sH2+uH2)/tH2;
  flow = 10202010;

  return CONVERT2MB * (M_PI/sH2) * alphaS2 * sigT;  

}

//*********

// q + q -> q + q, involving identical quarks.

double SigmaProcess::qq2qqSame() {
  
  // Calculate.
  double sigT = (4./9.) * (sH2+uH2)/tH2;
  double sigU = (4./9.) * (sH2+tH2)/uH2;
  double sigTU = - (8./27.) * sH2/(tH*uH);
  flow = ((sigT + sigU) * Rndm::flat() < sigT) ? 10202010 : 10201020;

  // Answer contains factor 1/2 from identical quarks.
  return CONVERT2MB * (M_PI/sH2) * alphaS2 * 0.5 * (sigT + sigU + sigTU);  

}

//*********

// q1 + qbar2 -> q1 + qbar2, where q1 and q2 are different.

double SigmaProcess::qqbar2qqbarDiff() {
  
  // Calculate.
  double sigT = (4./9.) * (sH2+uH2)/tH2;
  flow = 10012002;
  return CONVERT2MB * (M_PI/sH2) * alphaS2 * sigT;  

}

//*********

// q + qbar -> q + qbar, where q and qbar are each others antiparticles.

double SigmaProcess::qqbar2qqbarAnti() {
  
  // Calculate.
  double sigT = (4./9.) * (sH2+uH2)/tH2 - (8./27.) * uH2/(sH*tH);
  flow = 10012002;

  return CONVERT2MB * (M_PI/sH2) * alphaS2 * sigT;  

}

//*********

// q + qbar -> q' + q'bar, i.e. creation of new flavour (possibly same 
// as old). Quark masses still not properly taken into account.


double SigmaProcess::qqbar2qqbarNew() {

  // Pick new flavour.
  idNew = 1 + int( nQuark * Rndm::flat() ); 
  mNew = ParticleDataTable::m0(idNew);
  double  m2New = mNew*mNew;
  id3 = (id1 > 0) ? idNew : -idNew;
  id4 = -id3;
  
  // Calculate.
  double sigS = 0;
  if (sH > 4. * m2New) sigS = (4./9.) * (tH2+uH2)/sH2 
    * sqrt(1. - 4. * m2New / sH) * (1. +  2. * m2New / sH);
  flow = 10021002;

  return CONVERT2MB * (M_PI/sH2) * alphaS2 * nQuark * sigS;  

}

//*********

// q + qbar -> g + g.

double SigmaProcess::qqbar2gg() {

  // New flavours.
  id3 = 21;
  id4 = 21;
  
  // Calculate.
  double sigTS = (32./27.) * uH/tH - (8./3.) * uH2/sH2;
  double sigUS = (32./27.) * tH/uH - (8./3.) * tH2/sH2;
  double sigSum = sigTS + sigUS;
  flow = (sigSum * Rndm::flat() < sigTS) ? 10021332 : 10023213;

  return CONVERT2MB * (M_PI/sH2) * alphaS2 * sigSum;  

}

//*********

// g + g -> q + qbar, i.e. creation of new flavour.
// Quark masses still not properly taken into account.

double SigmaProcess::gg2qqbar() {

  // Pick new flavour.
  idNew = 1 + int( nQuark * Rndm::flat() ); 
  mNew = ParticleDataTable::m0(idNew);
  double  m2New = mNew*mNew;
  id3 = idNew;
  id4 = -idNew;
  
  // Calculate.
  double sigTS = 0.;
  double sigUS = 0.;
  if (sH > 4. * m2New) {
    // beta = sqrt(1. - 4. * m2New/sHat);
    sigTS = (1./6.) * uH/tH - (3./8.) * uH2/sH2;
    sigUS = (1./6.) * tH/uH - (3./8.) * tH2/sH2; 
  }
  double sigSum = sigTS + sigUS;
  flow = (sigSum * Rndm::flat() < sigTS) ? 12231003 : 12313002;

  return CONVERT2MB * (M_PI/sH2) * alphaS2 * (sigTS + sigUS);  

}

//*********

// q + g -> q + gamma.

double SigmaProcess::qg2qgam() {

  // New flavour.
  id4 = 22;
  
  // Calculate.
  double charge = ParticleDataTable::charge(id1);
  double sigUS = (1./3.) * (sH2+uH2)/(-sH*uH);
  flow = 10212000;

  return CONVERT2MB * (M_PI/sH2) * alphaS * alphaEM * pow2(charge) * sigUS;

}

//*********

// q + qbar -> g + gamma.

double SigmaProcess::qqbar2ggam() {

  // New flavours.
  id3 = 21;
  id4 = 22;
  
  // Calculate.
  double charge = ParticleDataTable::charge(id1);
  double sigTU = (8./9.) * (tH2+uH2)/(tH*uH);
  flow = 10021200;

  return CONVERT2MB * (M_PI/sH2) * alphaS * alphaEM * pow2(charge) * sigTU;

}

//*********

// q + qbar -> gamma + gamma.

double SigmaProcess::qqbar2gamgam() {

  // New flavours.
  id3 = 22;
  id4 = 22;
  
  // Calculate.
  double charge = ParticleDataTable::charge(id1);
  double sigTU = 2. * (tH2+uH2)/(tH*uH);
  flow = 10010000;

  return CONVERT2MB * (M_PI/sH2) * pow2(alphaEM) * pow4(charge) * sigTU;  

}

//*********

// Set up colour flow and kinematics of 2 -> 2 subprocess.
// Two versions: one feeding in all the required values, the
// other using current values. The latter can be used directly 
// if nothing happened since ME evaluation.

bool SigmaProcess::doKinematics(int id1In, int id2In, int id3In, int id4In, 
int colFlowIn, double sHatIn, double tHatIn) {
   
  // Set current values from input.
  id1 = id1In;
  id2 = id2In;
  id3 = id3In;
  id4 = id4In;
  flow = colFlowIn;
  sH = sHatIn;
  tH = tHatIn;
  uH = - sH - tH;
  
  // Call routine that does the real work.
  return doKinematics();

}

bool SigmaProcess::doKinematics() {

  // Check that masses of outgoing particles not too big.
  double m3 = ParticleDataTable::m0(id3);
  double m4 = ParticleDataTable::m0(id4);
  double eCM = sqrt(sH);
  if (m3 + m4 + MASSMARGIN > eCM) return false;

  // Unpack colours. Sometimes need to swap sides or colours<->anticolours.
  int cols[] = {flow/10000000, (flow/1000000)%10, (flow/100000)%10, 
    (flow/10000)%10, (flow/1000)%10, (flow/100)%10, (flow/10)%10, flow%10 };
  if (id1 > 20 && abs(id2) < 20) {
    swap( cols[0], cols[2]); swap( cols[1], cols[3]);  
    swap( cols[4], cols[6]); swap( cols[5], cols[7]); 
  }
  if (id1 < 0 || (id1 > 20 && id2 < 0) ) {
    swap( cols[0], cols[1]); swap( cols[2], cols[3]);  
    swap( cols[4], cols[5]); swap( cols[6], cols[7]); 
  }

  // Do kinematics of the decay.
  double eIn = 0.5 * eCM;
  double e3 = 0.5 * (sH + m3*m3 - m4*m4) / eCM;
  double e4 = 0.5 * (sH + m4*m4 - m3*m3) / eCM;
  double pAbs = sqrtpos( e3*e3 - m3*m3 );
  double cosTheta = 1. + 2. * tH / sH;
  double sinTheta = 2. * sqrtpos( tH * uH ) / sH;
  double phi = 2. * M_PI * Rndm::flat();
  double pZ = pAbs * cosTheta;
  double pX = pAbs * sinTheta * sin(phi);
  double pY = pAbs * sinTheta * cos(phi);
  double scale = eIn * sinTheta;

  // Fill particle info.
  parton[0] = Particle( id1, -31, -1, -1, 2, 3, cols[0], cols[1],
    0., 0., eIn, eIn, 0., scale);
  parton[1] = Particle( id2, -31, -1, -1, 2, 3, cols[2], cols[3],
    0., 0., -eIn, eIn, 0., scale);
  parton[2] = Particle( id3, 33, 0, 1, -1, -1, cols[4], cols[5],
    pX, pY, pZ, e3, m3, scale);
  parton[3] = Particle( id4, 33, 0, 1, -1, -1,cols[6], cols[7],
    -pX, -pY, -pZ, e4, m4, scale);

  // Done.
  return true;

}  
 
//**************************************************************************

} // end namespace Pythia8
