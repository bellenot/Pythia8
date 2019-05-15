// Function definitions (not found in the header) for the 
// InFlux and SigmaProcess classes, and classes derived from them.
// Copyright C 2006 Torbjorn Sjostrand

#include "SigmaProcess.h"

namespace Pythia8 {

//**************************************************************************

// InFlux class.
// Base class for the combined incoming parton flux.

//*********
 
// Definitions of static variables and functions.
// (Values will be overwritten in initStatic call, so are purely dummy.)

int InFlux::nQuark = 5;

//*********

// Initialize static data members.

void InFlux::initStatic() {

  // Maximum incoming quark flavour.
  nQuark = Settings::mode("inFlux:nQuark");

}

//*********

// Weight channels by (two or four) powers of e.

void InFlux::weightCharge(int ePow) {

  for (int i = 0; i < int(weightAB.size()); ++i) {
    int id = (abs(idPartonPairA[i]) < 9) ? idPartonPairA[i] 
      : idPartonPairB[i];
    double e2 = (id%2 == 0) ? 4./9. : 1./9.;
    if (ePow == 2) weightAB[i] *= e2;
    if (ePow == 4) weightAB[i] *= e2 * e2;
  }

}

//*********

// Calculate products of parton densities for allowed combinations.

double InFlux::flux(double x1, double x2, double Q2) {

  // Evaluate and store the required parton densities.
  for (int j = 0; j < int(idPartonA.size()); ++j) 
    pdfA[j] = pdfAPtr->xf( idPartonA[j], x1, Q2); 
  for (int j = 0; j < int(idPartonB.size()); ++j) 
    pdfB[j] = pdfBPtr->xf( idPartonB[j], x2, Q2); 

  // Multiply on these densities for each of the allowed channels. Sum.
  fluxwtSum = 0.;
  for (int i = 0; i < int(weightAB.size()); ++i) {
    for (int j = 0; j < int(idPartonA.size()); ++j) 
    if (idPartonPairA[i] == idPartonA[j]) {
      pdfPairA[i] = pdfA[j];
      break;
    }
    for (int j = 0; j < int(idPartonB.size()); ++j) 
    if (idPartonPairB[i] == idPartonB[j]) {
      pdfPairB[i] = pdfB[j];
      break;
    }
    fluxweightAB[i] = pdfPairA[i] * pdfPairB[i] * weightAB[i];
    fluxwtSum += fluxweightAB[i];
  }
 
  // Done.
  return fluxwtSum;

}

//*********

// Pick one of the possible channels according to their weight.

void InFlux::pick() {

  // Pick channel. Extract channel flavours.
  double fluxwtRand =  fluxwtSum * Rndm::flat();
  for (int i = 0; i < int(fluxweightAB.size()); ++i) {
    fluxwtRand -= fluxweightAB[i];
    if (fluxwtRand <= 0.) {
      idNow1 = idPartonPairA[i];
      idNow2 = idPartonPairB[i];
      pdfNow1 = pdfPairA[i];
      pdfNow2 = pdfPairB[i];
      break;
    }
  }

}

//**************************************************************************

// InFluxgg class.
// Derived class for a g g incoming state.

//*********

// Fill arrays for a g g incoming state at initialization.

void InFluxgg::initChannels() {

  // The beams individually.
  idPartonA.push_back(21);
  idPartonB.push_back(21);
  pdfA.push_back(0.);
  pdfB.push_back(0.);
 
  // The combined channels.
  idPartonPairA.push_back(21);
  idPartonPairB.push_back(21);
  pdfPairA.push_back(0.);
  pdfPairB.push_back(0.);
  weightAB.push_back(1.);
  fluxweightAB.push_back(1.);

}

//**************************************************************************

// InFluxqg class.
// Derived class for q g incoming states.

//*********

// Fill arrays for a q g incoming state at initialization.

void InFluxqg::initChannels() {

  // The beams individually.
  for (int i = -nQuark; i <= nQuark; ++i) {
    int id = (i == 0) ? 21 : i;
    idPartonA.push_back(id);
    idPartonB.push_back(id);
    pdfA.push_back(0.);
    pdfB.push_back(0.);
  }

  // The combined channels.
  for (int id = -nQuark; id <= nQuark; ++id) 
  if (id != 0) {
    idPartonPairA.push_back(id);
    idPartonPairB.push_back(21);
    pdfPairA.push_back(0.);
    pdfPairB.push_back(0.);
    weightAB.push_back(1.);
    fluxweightAB.push_back(0.);
    idPartonPairA.push_back(21);
    idPartonPairB.push_back(id);
    pdfPairA.push_back(0.);
    pdfPairB.push_back(0.);
    weightAB.push_back(1.);
    fluxweightAB.push_back(1.);
  }

}

//**************************************************************************

// InFluxqqbarqqDiff class.
// Derived class for q qbbar' or q q' (qbar qbar') incoming states, 
// with q' != q.

//*********

// Fill arrays for a q q' incoming states at initialization.

void InFluxqqbarqqDiff::initChannels() {

  // The beams individually.
  for (int id = -nQuark; id <= nQuark; ++id) 
  if (id != 0) {
    idPartonA.push_back(id);
    idPartonB.push_back(id);
    pdfA.push_back(0.);
    pdfB.push_back(0.);
  }

  // The combined channels.
  for (int id1 = -nQuark; id1 <= nQuark; ++id1) 
  if (id1 != 0) 
  for (int id2 = -nQuark; id2 <= nQuark; ++id2) 
  if (id2 != 0 && id2 != id1) {
    idPartonPairA.push_back(id1);
    idPartonPairB.push_back(id2);
    pdfPairA.push_back(0.);
    pdfPairB.push_back(0.);
    weightAB.push_back(1.);
    fluxweightAB.push_back(1.);
  }

}

//**************************************************************************

// InFluxqqDiff class.
// Derived class for q q' (qbar qbar') incoming states, with q' != q.

//*********

// Fill arrays for a q q' incoming states at initialization.

void InFluxqqDiff::initChannels() {

  // The beams individually.
  for (int id = -nQuark; id <= nQuark; ++id) 
  if (id != 0) {
    idPartonA.push_back(id);
    idPartonB.push_back(id);
    pdfA.push_back(0.);
    pdfB.push_back(0.);
  }

  // The combined channels.
  for (int id1 = -nQuark; id1 <= nQuark; ++id1) 
  if (id1 != 0) 
  for (int id2 = -nQuark; id2 <= nQuark; ++id2) 
  if (id2 != 0 && id2 != id1 && id1 * id2 > 0) {
    idPartonPairA.push_back(id1);
    idPartonPairB.push_back(id2);
    pdfPairA.push_back(0.);
    pdfPairB.push_back(0.);
    weightAB.push_back(1.);
    fluxweightAB.push_back(1.);
  }

}
//**************************************************************************

// InFluxqqSame class.
// Derived class for q q (qbar qbar) incoming states.

//*********

// Fill arrays for a q q incoming states at initialization.

void InFluxqqSame::initChannels() {

  // The beams individually.
  for (int id = -nQuark; id <= nQuark; ++id) 
  if (id != 0) {
    idPartonA.push_back(id);
    idPartonB.push_back(id);
    pdfA.push_back(0.);
    pdfB.push_back(0.);
  }

  // The combined channels.
  for (int id = -nQuark; id <= nQuark; ++id) 
  if (id != 0) {
    idPartonPairA.push_back(id);
    idPartonPairB.push_back(id);
    pdfPairA.push_back(0.);
    pdfPairB.push_back(0.);
    weightAB.push_back(1.);
    fluxweightAB.push_back(1.);
  }

}

//**************************************************************************

// InFluxqqbarDiff class.
// Derived class for q qbar' incoming states.

//*********

// Fill arrays for a q qbar incoming state at initialization.

void InFluxqqbarDiff::initChannels() {

  // The beams individually.
  for (int id = -nQuark; id <= nQuark; ++id) 
  if (id != 0) {
    idPartonA.push_back(id);
    idPartonB.push_back(id);
    pdfA.push_back(0.);
    pdfB.push_back(0.);
  }

  // The combined channels.
  for (int id1 = -nQuark; id1 <= nQuark; ++id1) 
  if (id1 != 0) 
  for (int id2 = -nQuark; id2 <= nQuark; ++id2) 
  if (id2 != 0 && id2 != -id1 && id1 * id2 < 0) {
    idPartonPairA.push_back(id1);
    idPartonPairB.push_back(id2);
    pdfPairA.push_back(0.);
    pdfPairB.push_back(0.);
    weightAB.push_back(1.);
    fluxweightAB.push_back(1.);
  }

}

//**************************************************************************

// InFluxqqbarSame class.
// Derived class for q qbar antiparticle incoming states.

//*********

// Fill arrays for a q qbar incoming state at initialization.

void InFluxqqbarSame::initChannels() {

  // The beams individually.
  for (int id = -nQuark; id <= nQuark; ++id) 
  if (id != 0) {
    idPartonA.push_back(id);
    idPartonB.push_back(id);
    pdfA.push_back(0.);
    pdfB.push_back(0.);
  }

  // The combined channels.
  for (int id = -nQuark; id <= nQuark; ++id) 
  if (id != 0) {
    idPartonPairA.push_back(id);
    idPartonPairB.push_back(-id);
    pdfPairA.push_back(0.);
    pdfPairB.push_back(0.);
    weightAB.push_back(1.);
    fluxweightAB.push_back(1.);
  }

}

//**************************************************************************

// The SigmaProcess class.
// Base class for cross sections.

//*********
 
// Definitions of static variables and functions.
// (Values will be overwritten in initStatic call, so are purely dummy.)

int SigmaProcess::alphaSorder = 1;
int SigmaProcess::nQuark = 3;
double SigmaProcess::alphaSvalue = 0.1265;
double SigmaProcess::alphaEMfix =  0.00729735;
AlphaStrong SigmaProcess::alphaScalc;

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Conversion of GeV^{-2} to mb for cross section.
const double SigmaProcess::CONVERT2MB = 0.389380; 

// The sum of outgoing masses must not be too close to the cm energy.
const double SigmaProcess::MASSMARGIN = 0.1;

//*********

// Initialize static data members.

void SigmaProcess::initStatic() {

  // Parameters of alphaStrong generation .
  alphaSvalue = Settings::parameter("SigmaProcess:alphaSvalue");
  alphaSorder = Settings::mode("SigmaProcess:alphaSorder");

  // Initialize alphaStrong generation.
  alphaScalc.init( alphaSvalue, alphaSorder); 

  // Maximum new quark flavour, alphaEM.
  nQuark = Settings::mode("SigmaProcess:nQuark");
  alphaEMfix = Settings::parameter("StandardModel:alphaEMfix"); 

}

//**************************************************************************

// The Sigma2Process class.
// Base class for resolved 2 -> 2 cross sections; derived from SigmaProcess.

//*********

// Input and complement kinematics for resolved 2 -> 2 process. 

bool Sigma2Process::set2Kin( double x1in, double x2in, double sHin, double tHin) {

  // Incoming parton momentum fractions.
  x1 = x1in;
  x2 = x2in;

  // Standard Mandelstam variables and their squares.
  sH = sHin;
  tH = tHin;
  uH = (masslessKin) ? -(sH + tH) : m3S + m4S -(sH + tH); 
  sH2 = sH * sH;
  tH2 = tH * tH;
  uH2 = uH * uH;

  // Calculate squared transverse momentum.
  pT2 = (masslessKin) ?  tH * uH / sH : (tH * uH - m3S * m4S) / sH;

  // Use pT2 as renormalization scale, generalized to m_T3 * m_T4
  // for massive case. Evaluate alpha_strong.
  Q2RenH = (masslessKin) ? pT2 : sqrt((pT2 + m3S) * (pT2 + m4S));
  alpS = alphaScalc.alphaS(Q2RenH);  

  // Use pT2 as factorization scale, also in massive case.
  Q2FacH = pT2;

  // Done.
  return true;

}

//*********

// As above, special kinematics for multiple interactions. 

bool Sigma2Process::set2KinMI( int id1in, int id2in, double x1in, double x2in,
  double sHin, double tHin, double uHin, double alpSin) {
 
  // Incoming flavours and x values.
  id1 = id1in;
  id2 = id2in;
  x1 = x1in;
  x2 = x2in;

  // Standard Mandelstam variables and their squares.
  sH = sHin;
  tH = tHin;
  uH = uHin; 
  sH2 = sH * sH;
  tH2 = tH * tH;
  uH2 = uH * uH;

  // Strong coupling.
  alpS = alpSin;

  //  Done.
  return true;

}

//*********

// Perform kinematics for a Multiple Interaction.

bool Sigma2Process::final2KinMI() {

  // Have to set flavours and colours.
  setIdColAcol();

  // Check that masses of outgoing particles not too big.
  m3 = ParticleDataTable::m0(idH[3]);
  m4 = ParticleDataTable::m0(idH[4]);
  double eCM = sqrt(sH);
  if (m3 + m4 + MASSMARGIN > eCM) return false;

  // Do kinematics of the decay.
  double eIn = 0.5 * eCM;
  double e3 = 0.5 * (sH + m3*m3 - m4*m4) / eCM;
  double e4 = 0.5 * (sH + m4*m4 - m3*m3) / eCM;
  double pAbs = sqrtpos( e3*e3 - m3*m3 );
  double cosTheta = 1. + 2. * tH / sH;
  double sinTheta = 2. * sqrtpos( tH * uH ) / sH;
  thetaH = atan2( sinTheta, cosTheta);
  phiH = 2. * M_PI * Rndm::flat();
  double pZ = pAbs * cosTheta;
  double pX = pAbs * sinTheta * sin(phiH);
  double pY = pAbs * sinTheta * cos(phiH);
  double scale = eIn * sinTheta;

  // Fill particle info.
  parton[1] = Particle( idH[1], -31, 0, 0, 3, 4, colH[1], acolH[1],
    0., 0., eIn, eIn, 0., scale);
  parton[2] = Particle( idH[2], -31, 0, 0, 3, 4, colH[2], acolH[2],
    0., 0., -eIn, eIn, 0., scale);
  parton[3] = Particle( idH[3],  33, 1, 2, 0, 0, colH[3], acolH[3],
    pX, pY, pZ, e3, m3, scale);
  parton[4] = Particle( idH[4],  33, 1, 2, 0, 0, colH[4], acolH[4],
    -pX, -pY, -pZ, e4, m4, scale);

  // Boost particles from subprocess rest frame to event rest frame.
  double betaZ = (x1 - x2) / (x1 + x2);
  for (int i = 1; i <= 4; ++i) parton[i].bst(0., 0., betaZ);

  // Done.
  return true;

}  

//**************************************************************************

// Sigma0AB2AB class.
// Cross section for elastic scattering A B -> A B.

//*********

// Select identity, colour and anticolour.

void Sigma0AB2AB::setIdColAcol() {

  // Flavours and colours are trivial. 
  setId( idA, idB, idA, idB);
  setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
}

//**************************************************************************

// Sigma0AB2XB class.
// Cross section for single diffractive scattering A B -> X B.

//*********

// Select identity, colour and anticolour.

void Sigma0AB2XB::setIdColAcol() {

  // Flavours and colours are trivial. 
  int idX = 10* (abs(idA) / 10) + 9900000; 
  if (idA < 0) idX = -idX;
  setId( idA, idB, idX, idB);
  setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);

}

//**************************************************************************

// Sigma0AB2AX class.
// Cross section for single diffractive scattering A B -> A X.

//*********

// Select identity, colour and anticolour.

void Sigma0AB2AX::setIdColAcol() {

  // Flavours and colours are trivial. 
  int idX = 10* (abs(idB) / 10) + 9900000; 
  if (idB < 0) idX = -idX;
  setId( idA, idB, idA, idX);
  setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);

}

//**************************************************************************

// Sigma0AB2XX class.
// Cross section for double diffractive scattering A B -> X X.

//*********

// Select identity, colour and anticolour.

void Sigma0AB2XX::setIdColAcol() {

  // Flavours and colours are trivial. 
  int idX1 = 10* (abs(idA) / 10) + 9900000; 
  if (idA < 0) idX1 = -idX1;
  int idX2 = 10* (abs(idB) / 10) + 9900000; 
  if (idB < 0) idX2 = -idX2;
  setId( idA, idB, idX1, idX2);
  setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);

}

//**************************************************************************

// Sigma2gg2gg class.
// Cross section for g g -> g g.

//*********

// Evaluate d(sigmaHat)/d(tHat).

double Sigma2gg2gg::sigmaHat() {

  // Calculate kinematics dependence.
  sigTS = (9./4.) * (tH2/sH2 + 2.*tH/sH + 3. + 2.*sH/tH + sH2/tH2);
  sigUS = (9./4.) * (uH2/sH2 + 2.*uH/sH + 3. + 2.*sH/uH + sH2/uH2);
  sigTU = (9./4.) * (tH2/uH2 + 2.*tH/uH + 3. + 2.*uH/tH + uH2/tH2);
  sigSum = sigTS + sigUS + sigTU;

  // Answer contains factor 1/2 from identical gluons.
  return CONVERT2MB * (M_PI/sH2) * pow2(alpS) * 0.5 * sigSum;  

}

//*********

// Initialize parton flux object. 
  
void Sigma2gg2gg::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for g g initial state and parton densities.
  inFluxPtr = new InFluxgg();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Select identity, colour and anticolour.

void Sigma2gg2gg::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, 21, 21);

  // Three colour flow topologies, each with two orientations.
  double sigRand = sigSum * Rndm::flat();
  if (sigRand < sigTS) setColAcol( 1, 2, 2, 3, 1, 4, 4, 3);
  else if (sigRand < sigTS + sigUS) 
                       setColAcol( 1, 2, 3, 1, 3, 4, 4, 2);
  else                 setColAcol( 1, 2, 3, 4, 1, 4, 3, 2); 
  if (Rndm::flat() > 0.5) swapColAcol();

}

//**************************************************************************

// SigmaHgg2qqbar class.
// Cross section for g g -> q qbar (q = u, d, s, i.e. almost massless).

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2gg2qqbar::sigmaHat() { 

  // Pick new flavour.
  idNew = 1 + int( nQuark * Rndm::flat() ); 
  mNew = ParticleDataTable::m0(idNew);
  m2New = mNew*mNew;
  
  // Calculate kinematics dependence.
  sigTS = 0.;
  sigUS = 0.;
  if (sH > 4. * m2New) {
    sigTS = (1./6.) * uH/tH - (3./8.) * uH2/sH2;
    sigUS = (1./6.) * tH/uH - (3./8.) * tH2/sH2; 
  }
  sigSum = sigTS + sigUS;

  // Answer is proportional to number of outgoing flavours.
  return CONVERT2MB * (M_PI/sH2) * pow2(alpS) * nQuark * sigSum;  

}

//*********

// Initialize parton flux object. 
  
void Sigma2gg2qqbar::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for g g initial state and parton densities.
  inFluxPtr = new InFluxgg();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Select identity, colour and anticolour.

void Sigma2gg2qqbar::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idNew, -idNew);

  // Two colour flow topologies.
  double sigRand = sigSum * Rndm::flat();
  if (sigRand < sigTS) setColAcol( 1, 2, 2, 3, 1, 0, 0, 3);
  else                 setColAcol( 1, 2, 3, 1, 3, 0, 0, 2); 

}

//**************************************************************************

// Sigma2qg2qg class.
// Cross section for q g -> q g.

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qg2qg::sigmaHat() { 

  // Calculate kinematics dependence.
  sigTS = uH2/tH2 - (4./9.) * uH/sH;
  sigTU = sH2/tH2 - (4./9.) * sH/uH;
  sigSum = sigTS + sigTU;

  // Answer.
  return CONVERT2MB * (M_PI/sH2) * pow2(alpS) * sigSum;  

}

//*********

// Initialize parton flux object. 
  
void Sigma2qg2qg::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxqg();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Select identity, colour and anticolour.

void Sigma2qg2qg::setIdColAcol() {

  // Outgoing = incoming flavours.
  setId( id1, id2, id1, id2);

  // Two colour flow topologies. Swap if first is gluon, or when antiquark.
  double sigRand = sigSum * Rndm::flat();
  if (sigRand < sigTS) setColAcol( 1, 0, 2, 1, 3, 0, 2, 3);
  else                 setColAcol( 1, 0, 2, 3, 2, 0, 1, 3); 
  if (id1 == 21) swap1234();
  if (id1 < 0 || id2 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2qq2qqDiff class.
// Cross section for q qbar' -> q qbar' or q q' -> q q' 
// (qbar qbar' -> qbar qbar'), q' != q.

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qq2qqDiff::sigmaHat() { 

  // Calculate kinematics dependence.
  sigT = (4./9.) * (sH2+uH2)/tH2;

  // Answer.
  return CONVERT2MB * (M_PI/sH2) * pow2(alpS) * sigT;  

}

//*********

// Initialize parton flux object. 
  
void Sigma2qq2qqDiff::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxqqbarqqDiff();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Select identity, colour and anticolour.

void Sigma2qq2qqDiff::setIdColAcol() {

  // Outgoing = incoming flavours.
  setId( id1, id2, id1, id2);

  // Colour flow topologies. Swap when antiquarks.
  if (id1 * id2 > 0) setColAcol( 1, 0, 2, 0, 2, 0, 1, 0);
  else               setColAcol( 1, 0, 0, 1, 2, 0, 0, 2);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2qq2qqSame class.
// Cross section for q q -> q q (qbar qbar -> qbar qbar).

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qq2qqSame::sigmaHat() { 

  // Calculate kinematics dependence.
  sigT = (4./9.) * (sH2+uH2)/tH2;
  sigU = (4./9.) * (sH2+tH2)/uH2;
  sigTU = - (8./27.) * sH2/(tH*uH);
  sigSum = sigT + sigU + sigTU;

  // Answer contains factor 1/2 from identical quarks.
  return CONVERT2MB * (M_PI/sH2) * pow2(alpS) * 0.5 * sigSum;  

}

//*********

// Initialize parton flux object. 
  
void Sigma2qq2qqSame::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxqqSame();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Select identity, colour and anticolour.

void Sigma2qq2qqSame::setIdColAcol() {

  // Outgoing = incoming flavours.
  setId( id1, id2, id1, id2);

  // Two colour flow topologies. Swap if first is antiquark.
  double sigRand = (sigT + sigU) * Rndm::flat();
  if (sigRand < sigT) setColAcol( 1, 0, 2, 0, 2, 0, 1, 0);
  else                setColAcol( 1, 0, 2, 0, 1, 0, 2, 0); 
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2qqbar2qqbarSame class.
// Cross section q qbar -> q qbar.

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qqbar2qqbarSame::sigmaHat() {

  // Calculate kinematics dependence.
  sigT = (4./9.) * (sH2+uH2)/tH2 - (8./27.) * uH2/(sH*tH);

  // Answer.
  return CONVERT2MB * (M_PI/sH2) * pow2(alpS) * sigT;  

}

//*********

// Initialize parton flux object. 
  
void Sigma2qqbar2qqbarSame::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxqqbarSame();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Select identity, colour and anticolour.

void Sigma2qqbar2qqbarSame::setIdColAcol() {

  // Outgoing = incoming flavours.
  setId( id1, id2, id1, id2);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 1, 2, 0, 0, 2);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2qqbar2qqbarNew class.
// Cross section q qbar -> q' qbar'.

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qqbar2qqbarNew::sigmaHat() { 

  // Pick new flavour.
  idNew = 1 + int( nQuark * Rndm::flat() ); 
  mNew = ParticleDataTable::m0(idNew);
  m2New = mNew*mNew;

  // Calculate kinematics dependence.
  sigS = 0.;
  if (sH > 4. * m2New) sigS = (4./9.) * (tH2+uH2)/sH2; 

  // Answer is proportional to number of outgoing flavours.
  return CONVERT2MB * (M_PI/sH2) * pow2(alpS) * nQuark * sigS;  

}

//*********

// Initialize parton flux object. 
  
void Sigma2qqbar2qqbarNew::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxqqbarSame();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Select identity, colour and anticolour.

void Sigma2qqbar2qqbarNew::setIdColAcol() {

  // Set outgoing flavours ones.
  id3 = (id1 > 0) ? idNew : -idNew;
  setId( id1, id2, id3, -id3);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 2, 1, 0, 0, 2);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2qqbar2gg class.
// Cross section for q qbar -> g g.

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qqbar2gg::sigmaHat() { 

  // Calculate kinematics dependence.
  sigTS = (32./27.) * uH/tH - (8./3.) * uH2/sH2;
  sigUS = (32./27.) * tH/uH - (8./3.) * tH2/sH2;
  sigSum = sigTS + sigUS;

  // Answer contains factor 1/2 from identical gluons.
  return CONVERT2MB * (M_PI/sH2) * pow2(alpS) * 0.5 * sigSum;  

}

//*********

// Initialize parton flux object. 
  
void Sigma2qqbar2gg::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxqqbarSame();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Select identity, colour and anticolour.

void Sigma2qqbar2gg::setIdColAcol() {

  // Outgoing flavours trivial.
  setId( id1, id2, 21, 21);

  // Two colour flow topologies. Swap if first is antiquark.
  double sigRand = sigSum * Rndm::flat();
  if (sigRand < sigTS) setColAcol( 1, 0, 0, 2, 1, 3, 3, 2);
  else                 setColAcol( 1, 0, 0, 2, 3, 2, 1, 3); 
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2gg2QQbar class.
// Cross section g g -> Q Qbar (Q = c or b).

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2gg2QQbar::sigmaHat() { 

  // Modified Mandelstam variables for massive kinematics.
  double tHQ = tH - m2New;
  double uHQ = uH - m2New; 
  double tHQ2 = tHQ * tHQ;
  double uHQ2 = uHQ * uHQ;

  // Calculate kinematics dependence.
  sigTS = 0.;
  sigUS = 0.;
  if (sH > 4. * m2New) {
    double tumHQ = tHQ * uHQ - m2New * sH;
    sigTS = ( uHQ / tHQ - 2.25 * uHQ2 / sH2 + 4.5 * m2New * tumHQ 
      / ( sH * tHQ2) + 0.5 * m2New * (tHQ + m2New) / tHQ2 
      - m2New*m2New / (sH * tHQ) ) / 6.;
    sigUS = ( tHQ / uHQ - 2.25 * tHQ2 / sH2 + 4.5 * m2New * tumHQ 
      / ( sH * uHQ2) + 0.5 * m2New * (uHQ + m2New) / uHQ2 
      - m2New*m2New / (sH * uHQ) ) / 6.;
  }
  sigSum = sigTS + sigUS;

  // Answer.
  return CONVERT2MB * (M_PI/sH2) * pow2(alpS) * sigSum;  

}

//*********

// Initialize parton flux object. 
  
void Sigma2gg2QQbar::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxgg();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

  // Set up masses.
  mNew = ParticleDataTable::m0(idNew);
  m2New = mNew*mNew;
  mH[3] = mNew;
  mH[4] = mNew;
  m3 = mNew;
  m3S = m2New;
  m4 = mNew;
  m4S = m2New;

} 

//*********

// Select identity, colour and anticolour.

void Sigma2gg2QQbar::setIdColAcol() {

  // Flavours are trivial.
  setId( id1, id2, idNew, -idNew);

  // Two colour flow topologies.
  double sigRand = sigSum * Rndm::flat();
  if (sigRand < sigTS) setColAcol( 1, 2, 2, 3, 1, 0, 0, 3);
  else                 setColAcol( 1, 2, 3, 1, 3, 0, 0, 2); 

}

//**************************************************************************

// Sigma2qqbar2QQbar class.
// Cross section q qbar -> Q Qbar (Q = c or b).

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qqbar2QQbar::sigmaHat() { 

  // Modified Mandelstam variables and their squares.
  double tHQ = tH - m2New;
  double uHQ = uH - m2New; 
  double tHQ2 = tHQ * tHQ;
  double uHQ2 = uHQ * uHQ;

  // Calculate kinematics dependence.
  sigS = 0.;
  if (sH > 4. * m2New) sigS = (4./9.) * ((tHQ2 + uHQ2) / sH2 
    + 2. * m2New / sH); 

  // Answer.
  return CONVERT2MB * (M_PI/sH2) * pow2(alpS) * sigS;  

}

//*********

// Initialize parton flux object. 
  
void Sigma2qqbar2QQbar::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxqqbarSame();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

  // Set up masses.
  mNew = ParticleDataTable::m0(idNew);
  m2New = mNew*mNew;
  mH[3] = mNew;
  mH[4] = mNew;
  m3 = mNew;
  m3S = m2New;
  m4 = mNew;
  m4S = m2New;

} 

//*********

// Select identity, colour and anticolour.

void Sigma2qqbar2QQbar::setIdColAcol() {

  // Set outgoing flavours.
  id3 = (id1 > 0) ? idNew : -idNew;
  setId( id1, id2, id3, -id3);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 2, 1, 0, 0, 2);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2qg2qgamma class.
// Cross section for q g -> q gamma.

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qg2qgamma::sigmaHat() {  

  // Calculate kinematics dependence.
  sigUS = (1./3.) * (sH2+uH2)/(-sH*uH);

  // Answer.
  return CONVERT2MB * (M_PI/sH2) * alpS * alpEM * sigUS;  

}

//*********

// Initialize parton flux object. 
  
void Sigma2qg2qgamma::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state, e2-weighted, and parton densities.
  inFluxPtr = new InFluxqg();
  inFluxPtr->init( pdfAPtr, pdfBPtr);
  inFluxPtr->weightCharge(2);

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
  if (id1 == 21) swap1234();
  if (id1 < 0 || id2 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2qqbar2ggamma class.
// Cross section for q qbar -> g gamma.

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qqbar2ggamma::sigmaHat() { 

  // Calculate kinematics dependence.
  double sigTU = (8./9.) * (tH2+uH2)/(tH*uH);

  // Answer.
  return CONVERT2MB * (M_PI/sH2) * alpS * alpEM * sigTU;

}

//*********

// Initialize parton flux object. 
  
void Sigma2qqbar2ggamma::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state, e2-weighted, and parton densities.
  inFluxPtr = new InFluxqqbarSame();
  inFluxPtr->init( pdfAPtr, pdfBPtr);
  inFluxPtr->weightCharge(2);

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

// Initialize parton flux object. 
  
void Sigma2gg2ggamma::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxgg();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

  // Calculate charge factor from the allowed quarks in the box. 
  int nQuarkInLoop = Settings::mode("SigmaProcess:nQuarkInLoop");
  chargeSum = - 1./3. + 2./3. - 1./3.;
  if (nQuarkInLoop >= 4) chargeSum += 2./3.;
  if (nQuarkInLoop >= 5) chargeSum -= 1./3.;
  if (nQuarkInLoop >= 6) chargeSum += 2./3.;

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

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qqbar2gammagamma::sigmaHat() { 

  // Calculate kinematics dependence. Colour factor for quarks in.
  sigTU = 2. * (tH2+uH2)/(tH*uH);
  double colFac = 1./3.;

  // Answer contains factor 1/2 from identical photons.
  return CONVERT2MB * (M_PI/sH2) * pow2(alpEM) * 0.5 * sigTU * colFac;

}

//*********

// Initialize parton flux object. 
  
void Sigma2qqbar2gammagamma::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state, e4-weighted, and parton densities.
  inFluxPtr = new InFluxqqbarSame();
  inFluxPtr->init( pdfAPtr, pdfBPtr);
  inFluxPtr->weightCharge(4);

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

// Initialize parton flux object. 
  
void Sigma2gg2gammagamma::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxgg();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

  // Calculate charge factor from the allowed quarks in the box. 
  int nQuarkInLoop = Settings::mode("SigmaProcess:nQuarkInLoop");
  charge2Sum = 1./9. + 4./9. + 1./9.;
  if (nQuarkInLoop >= 4) charge2Sum += 4./9.;
  if (nQuarkInLoop >= 5) charge2Sum += 1./9.;
  if (nQuarkInLoop >= 6) charge2Sum += 4./9.;

} 

//*********

// Select identity, colour and anticolour.

void Sigma2gg2gammagamma::setIdColAcol() {

  // Flavours and colours are trivial.
  setId( id1, id2, 22, 22);
  setColAcol( 1, 2, 2, 1, 0, 0, 0, 0);

}


//**************************************************************************

} // end namespace Pythia8
