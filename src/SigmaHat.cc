// Function definitions (not found in the header) for the 
// InFlux, SigmaHat and SigmaHgg2gg classes.
// Copyright C 2006 Torbjorn Sjostrand

#include "SigmaHat.h"

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

// InFluxNone class.
// Derived class for dummy cases, e.g. elastic/diffractive scattering.

//*********

// Fill arrays for an empty incoming state at initialization.

void InFluxNone::initChannels() {

  // The beams individually.
  idPartonA.push_back(0);
  idPartonB.push_back(0);
  pdfA.push_back(0.);
  pdfB.push_back(0.);
 
  // The combined channels.
  idPartonPairA.push_back(0);
  idPartonPairB.push_back(0);
  pdfPairA.push_back(0.);
  pdfPairB.push_back(0.);
  weightAB.push_back(1.);
  fluxweightAB.push_back(1.);

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

// The SigmaHat class.
// Base class for cross sections.

//*********
 
// Definitions of static variables and functions.
// (Values will be overwritten in initStatic call, so are purely dummy.)

int SigmaHat::alphaSorder = 1;
int SigmaHat::nQuark = 3;
double SigmaHat::alphaSvalue = 0.1265;
double SigmaHat::alphaEM =  0.00729735;
AlphaStrong SigmaHat::alphaScalc;

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Conversion of GeV^{-2} to mb for cross section.
const double SigmaHat::CONVERT2MB = 0.389380; 

//*********

// Initialize static data members.

void SigmaHat::initStatic() {

  // Parameters of alphaStrong generation .
  alphaSvalue = Settings::parameter("SigmaHat:alphaSvalue");
  alphaSorder = Settings::mode("SigmaHat:alphaSorder");

  // Initialize alphaStrong generation.
  alphaScalc.init( alphaSvalue, alphaSorder); 

  // Maximum new quark flavour, alphaEM.
  nQuark = Settings::mode("SigmaHat:nQuark");
  alphaEM = Settings::parameter("StandardModel:alphaEMfix"); 

}

//*********

// Store info of use during the generation.

void SigmaHat::initInfo(Info& info) {

  // Beam identities and masses.
  idA = info.idA();
  idB = info.idB();
  mA = info.mA();
  mB = info.mB();

}

//**************************************************************************

// SigmaHAB2AB class.
// Cross section for elastic scattering A B -> A B.

//*********

// Initialize parton flux object - dummy in this case. 
  
void SigmaHAB2AB::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  inFluxPtr = new InFluxNone();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Select identity, colour and anticolour.

void SigmaHAB2AB::setIdColAcol() {

  // Flavours and colours are trivial. 
  setId( idA, idB, idA, idB);
  setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);

}

//**************************************************************************

// SigmaHAB2XB class.
// Cross section for single diffractive scattering A B -> X B.

//*********

// Initialize parton flux object - dummy in this case. 
  
void SigmaHAB2XB::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  inFluxPtr = new InFluxNone();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Select identity, colour and anticolour.

void SigmaHAB2XB::setIdColAcol() {

  // Flavours and colours are trivial. 
  int idX = 10* (abs(idA) / 10) + 9900000; 
  if (idA < 0) idX = -idX;
  setId( idA, idB, idX, idB);
  setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);

}

//**************************************************************************

// SigmaHAB2AX class.
// Cross section for single diffractive scattering A B -> A X.

//*********

// Initialize parton flux object - dummy in this case. 
  
void SigmaHAB2AX::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  inFluxPtr = new InFluxNone();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Select identity, colour and anticolour.

void SigmaHAB2AX::setIdColAcol() {

  // Flavours and colours are trivial. 
  int idX = 10* (abs(idB) / 10) + 9900000; 
  if (idB < 0) idX = -idX;
  setId( idA, idB, idA, idX);
  setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);

}

//**************************************************************************

// SigmaHAB2XX class.
// Cross section for double diffractive scattering A B -> X X.

//*********

// Initialize parton flux object - dummy in this case. 
  
void SigmaHAB2XX::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  inFluxPtr = new InFluxNone();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Select identity, colour and anticolour.

void SigmaHAB2XX::setIdColAcol() {

  // Flavours and colours are trivial. 
  int idX1 = 10* (abs(idA) / 10) + 9900000; 
  if (idA < 0) idX1 = -idX1;
  int idX2 = 10* (abs(idB) / 10) + 9900000; 
  if (idB < 0) idX2 = -idX2;
  setId( idA, idB, idX1, idX2);
  setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);

}

//**************************************************************************

// SigmaHgg2gg class.
// Cross section for g g -> g g.

//*********

// Initialize parton flux object. 
  
void SigmaHgg2gg::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for g g initial state and parton densities.
  inFluxPtr = new InFluxgg();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double SigmaHgg2gg::sigma2( double x1, double x2, double sH, double tH) {

  // Standard Mandelstam variables and their squares.
  double uH = -(sH + tH); 
  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

  // Calculate kinematics dependence.
  sigTS = (9./4.) * (tH2/sH2 + 2.*tH/sH + 3. + 2.*sH/tH + sH2/tH2);
  sigUS = (9./4.) * (uH2/sH2 + 2.*uH/sH + 3. + 2.*sH/uH + sH2/uH2);
  sigTU = (9./4.) * (tH2/uH2 + 2.*tH/uH + 3. + 2.*uH/tH + uH2/tH2);
  sigSum = sigTS + sigUS + sigTU;

  // Calculate running alpha_strong, using scale pT2 = tH*uH/sH.
  double pT2 = tH * uH / sH;
  double alphaS = alphaScalc.alphaS(pT2);  
  double alphaS2 = alphaS * alphaS;  

  // Also use pT2 as PDF scale. Evaluate incoming parton densities.
  q2pdf = pT2;
  double pdf = inFluxPtr->flux( x1, x2, q2pdf);

  // Answer contains factor 1/2 from identical gluons.
  return CONVERT2MB * (M_PI/sH2) * alphaS2 * 0.5 * sigSum * pdf;  

}

//*********

// Select identity, colour and anticolour.

void SigmaHgg2gg::setIdColAcol() {

  // Flavours are trivial.
  setId( 21, 21, 21, 21);

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

// Initialize parton flux object. 
  
void SigmaHgg2qqbar::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for g g initial state and parton densities.
  inFluxPtr = new InFluxgg();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double SigmaHgg2qqbar::sigma2( double x1, double x2, double sH, double tH) {

  // Standard Mandelstam variables and their squares.
  double uH = -(sH + tH); 
  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

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
  double sigSum = sigTS + sigUS;

  // Calculate running alpha_strong, using scale pT2 = tH*uH/sH.
  double pT2 = tH * uH / sH;
  double alphaS = alphaScalc.alphaS(pT2);  
  double alphaS2 = alphaS * alphaS;  

  // Also use pT2 as PDF scale. Evaluate incoming parton densities.
  q2pdf = pT2;
  double pdf = inFluxPtr->flux( x1, x2, q2pdf);

  // Answer is proportional to number of outgoing flavours.
  return CONVERT2MB * (M_PI/sH2) * alphaS2 * nQuark * sigSum * pdf;  

}

//*********

// Select identity, colour and anticolour.

void SigmaHgg2qqbar::setIdColAcol() {

  // Flavours are trivial.
  setId( 21, 21, idNew, -idNew);

  // Two colour flow topologies.
  double sigRand = sigSum * Rndm::flat();
  if (sigRand < sigTS) setColAcol( 1, 2, 2, 3, 1, 0, 0, 3);
  else                 setColAcol( 1, 2, 3, 1, 3, 0, 0, 2); 

}

//**************************************************************************

// SigmaHqg2qg class.
// Cross section for q g -> q g.

//*********

// Initialize parton flux object. 
  
void SigmaHqg2qg::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxqg();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double SigmaHqg2qg::sigma2( double x1, double x2, double sH, double tH) {

  // Standard Mandelstam variables and their squares.
  double uH = -(sH + tH); 
  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

  // Calculate kinematics dependence.
  sigTS = uH2/tH2 - (4./9.) * uH/sH;
  sigTU = sH2/tH2 - (4./9.) * sH/uH;
  sigSum = sigTS + sigTU;

  // Calculate running alpha_strong, using scale pT2 = tH*uH/sH.
  double pT2 = tH * uH / sH;
  double alphaS = alphaScalc.alphaS(pT2);  
  double alphaS2 = alphaS * alphaS;  

  // Also use pT2 as PDF scale. Evaluate incoming parton densities.
  q2pdf = pT2;
  double pdf = inFluxPtr->flux( x1, x2, q2pdf);

  // Answer.
  return CONVERT2MB * (M_PI/sH2) * alphaS2 * sigSum * pdf;  

}

//*********

// Select identity, colour and anticolour.

void SigmaHqg2qg::setIdColAcol() {

  // Extract incoming = outgoing flavours.
  int id1 = inFluxPtr->id1();
  int id2 = inFluxPtr->id2();
  setId( id1, id2, id1, id2);

  // Two colour flow topologies. Swap if first is gluon, or when antiquark.
  double sigRand = sigSum * Rndm::flat();
  if (sigRand < sigTS) setColAcol( 1, 0, 2, 1, 3, 0, 2, 3);
  else                 setColAcol( 1, 0, 2, 3, 2, 0, 1, 3); 
  if (id1 == 21) swap1234();
  if (id1 < 0 || id2 < 0) swapColAcol();

}

//**************************************************************************

// SigmaHqq2qqDiff class.
// Cross section for q qbar' -> q qbar' or q q' -> q q' 
// (qbar qbar' -> qbar qbar'), q' != q.

//*********

// Initialize parton flux object. 
  
void SigmaHqq2qqDiff::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxqqbarqqDiff();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double SigmaHqq2qqDiff::sigma2( double x1, double x2, double sH, 
  double tH) {

  // Standard Mandelstam variables and their squares.
  double uH = -(sH + tH); 
  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

  // Calculate kinematics dependence.
  sigT = (4./9.) * (sH2+uH2)/tH2;

  // Calculate running alpha_strong, using scale pT2 = tH*uH/sH.
  double pT2 = tH * uH / sH;
  double alphaS = alphaScalc.alphaS(pT2);  
  double alphaS2 = alphaS * alphaS;  

  // Also use pT2 as PDF scale. Evaluate incoming parton densities.
  q2pdf = pT2;
  double pdf = inFluxPtr->flux( x1, x2, q2pdf);

  // Answer.
  return CONVERT2MB * (M_PI/sH2) * alphaS2 * sigT * pdf;  

}

//*********

// Select identity, colour and anticolour.

void SigmaHqq2qqDiff::setIdColAcol() {

  // Extract incoming = outgoing flavours.
  int id1 = inFluxPtr->id1();
  int id2 = inFluxPtr->id2();
  setId( id1, id2, id1, id2);

  // Colour flow topologies. Swap when antiquarks.
  if (id1 * id2 > 0) setColAcol( 1, 0, 2, 0, 2, 0, 1, 0);
  else               setColAcol( 1, 0, 0, 1, 2, 0, 0, 2);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// SigmaHqq2qqSame class.
// Cross section for q q -> q q (qbar qbar -> qbar qbar).

//*********

// Initialize parton flux object. 
  
void SigmaHqq2qqSame::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxqqSame();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double SigmaHqq2qqSame::sigma2( double x1, double x2, double sH, 
  double tH) {

  // Standard Mandelstam variables and their squares.
  double uH = -(sH + tH); 
  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

  // Calculate kinematics dependence.
  sigT = (4./9.) * (sH2+uH2)/tH2;
  sigU = (4./9.) * (sH2+tH2)/uH2;
  sigTU = - (8./27.) * sH2/(tH*uH);
  sigSum = sigT + sigU + sigTU;

  // Calculate running alpha_strong, using scale pT2 = tH*uH/sH.
  double pT2 = tH * uH / sH;
  double alphaS = alphaScalc.alphaS(pT2);  
  double alphaS2 = alphaS * alphaS;  

  // Also use pT2 as PDF scale. Evaluate incoming parton densities.
  q2pdf = pT2;
  double pdf = inFluxPtr->flux( x1, x2, q2pdf);

  // Answer contains factor 1/2 from identical quarks.
  return CONVERT2MB * (M_PI/sH2) * alphaS2 * 0.5 * sigSum * pdf;  

}

//*********

// Select identity, colour and anticolour.

void SigmaHqq2qqSame::setIdColAcol() {

  // Extract incoming = outgoing flavours.
  int id = inFluxPtr->id1();
  setId( id, id, id, id);

  // Two colour flow topologies. Swap if first is antiquark.
  double sigRand = (sigT + sigU) * Rndm::flat();
  if (sigRand < sigT) setColAcol( 1, 0, 2, 0, 2, 0, 1, 0);
  else                setColAcol( 1, 0, 2, 0, 1, 0, 2, 0); 
  if (id < 0) swapColAcol();

}

//**************************************************************************

// SigmaHqqbar2qqbarSame class.
// Cross section q qbar -> q qbar.

//*********

// Initialize parton flux object. 
  
void SigmaHqqbar2qqbarSame::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxqqbarSame();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double SigmaHqqbar2qqbarSame::sigma2( double x1, double x2, double sH, 
  double tH) {

  // Standard Mandelstam variables and their squares.
  double uH = -(sH + tH); 
  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

  // Calculate kinematics dependence.
  sigT = (4./9.) * (sH2+uH2)/tH2 - (8./27.) * uH2/(sH*tH);

  // Calculate running alpha_strong, using scale pT2 = tH*uH/sH.
  double pT2 = tH * uH / sH;
  double alphaS = alphaScalc.alphaS(pT2);  
  double alphaS2 = alphaS * alphaS;  

  // Also use pT2 as PDF scale. Evaluate incoming parton densities.
  q2pdf = pT2;
  double pdf = inFluxPtr->flux( x1, x2, q2pdf);

  // Answer.
  return CONVERT2MB * (M_PI/sH2) * alphaS2 * sigT * pdf;  

}

//*********

// Select identity, colour and anticolour.

void SigmaHqqbar2qqbarSame::setIdColAcol() {

  // Extract incoming = outgoing flavours.
  int id1 = inFluxPtr->id1();
  int id2 = -id1;
  setId( id1, id2, id1, id2);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 1, 2, 0, 0, 2);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// SigmaHqqbar2qqbarNew class.
// Cross section q qbar -> q' qbar'.

//*********

// Initialize parton flux object. 
  
void SigmaHqqbar2qqbarNew::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxqqbarSame();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double SigmaHqqbar2qqbarNew::sigma2( double x1, double x2, double sH, 
  double tH) {

  // Standard Mandelstam variables and their squares.
  double uH = -(sH + tH); 
  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

  // Pick new flavour.
  idNew = 1 + int( nQuark * Rndm::flat() ); 
  mNew = ParticleDataTable::m0(idNew);
  m2New = mNew*mNew;

  // Calculate kinematics dependence.
  sigS = 0.;
  if (sH > 4. * m2New) sigS = (4./9.) * (tH2+uH2)/sH2; 

  // Calculate running alpha_strong, using scale pT2 = tH*uH/sH.
  double pT2 = tH * uH / sH;
  double alphaS = alphaScalc.alphaS(pT2);  
  double alphaS2 = alphaS * alphaS;  

  // Also use pT2 as PDF scale. Evaluate incoming parton densities.
  q2pdf = pT2;
  double pdf = inFluxPtr->flux( x1, x2, q2pdf);

  // Answer is proportional to number of outgoing flavours.
  return CONVERT2MB * (M_PI/sH2) * alphaS2 * nQuark * sigS * pdf;  

}

//*********

// Select identity, colour and anticolour.

void SigmaHqqbar2qqbarNew::setIdColAcol() {

  // Extract incoming flavours, set outgoing ones.
  int id1 = inFluxPtr->id1();
  int id3 = (id1 > 0) ? idNew : -idNew;
  setId( id1, -id1, id3, -id3);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 2, 1, 0, 0, 2);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// SigmaHqqbar2gg class.
// Cross section for q qbar -> g g.

//*********

// Initialize parton flux object. 
  
void SigmaHqqbar2gg::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxqqbarSame();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double SigmaHqqbar2gg::sigma2( double x1, double x2, double sH, double tH) {

  // Standard Mandelstam variables and their squares.
  double uH = -(sH + tH); 
  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

  // Calculate kinematics dependence.
  sigTS = (32./27.) * uH/tH - (8./3.) * uH2/sH2;
  sigUS = (32./27.) * tH/uH - (8./3.) * tH2/sH2;
  sigSum = sigTS + sigUS;

  // Calculate running alpha_strong, using scale pT2 = tH*uH/sH.
  double pT2 = tH * uH / sH;
  double alphaS = alphaScalc.alphaS(pT2);  
  double alphaS2 = alphaS * alphaS;  

  // Also use pT2 as PDF scale. Evaluate incoming parton densities.
  q2pdf = pT2;
  double pdf = inFluxPtr->flux( x1, x2, q2pdf);

  // Answer contains factor 1/2 from identical gluons.
  return CONVERT2MB * (M_PI/sH2) * alphaS2 * 0.5 * sigSum * pdf;  

}

//*********

// Select identity, colour and anticolour.

void SigmaHqqbar2gg::setIdColAcol() {

  // Extract incoming flavours; outgoing ones trivial.
  int id1 = inFluxPtr->id1();
  int id2 = inFluxPtr->id2();
  setId( id1, id2, 21, 21);

  // Two colour flow topologies. Swap if first is antiquark.
  double sigRand = sigSum * Rndm::flat();
  if (sigRand < sigTS) setColAcol( 1, 0, 0, 2, 1, 3, 3, 2);
  else                 setColAcol( 1, 0, 0, 2, 3, 2, 1, 3); 
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// SigmaHqqbar2QQbar class.
// Cross section q qbar -> Q Qbar (Q = c or b).

//*********

// Initialize parton flux object. 
  
void SigmaHqqbar2QQbar::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxqqbarSame();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double SigmaHqqbar2QQbar::sigma2( double x1, double x2, double sH, 
  double tH) {

  // Set up masses.
  mNew = ParticleDataTable::m0(idNew);
  m2New = mNew*mNew;

  // Standard Mandelstam variables and their squares.
  double uH = 2. * m2New - (sH + tH); 
  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

  // Calculate kinematics dependence.
  sigS = 0.;
  if (sH > 4. * m2New) sigS = (4./9.) * ((tH2 + uH2) / sH2 
    + 2. * m2New / sH); 

  // Calculate running alpha_strong, using scale mT2 = m2 + pT2.
  double pT2 = (tH * uH - pow2(m2New)) / sH;
  double mT2 = m2New + pT2;
  double alphaS = alphaScalc.alphaS(mT2);  
  double alphaS2 = alphaS * alphaS;  

  // Also use mT2 as PDF scale. Evaluate incoming parton densities.
  q2pdf = mT2;
  double pdf = inFluxPtr->flux( x1, x2, q2pdf);
  // Should put = 0 for inflavour = outflavour ??

  // Answer is proportional to number of outgoing flavours.
  return CONVERT2MB * (M_PI/sH2) * alphaS2 * sigS * pdf;  

}

//*********

// Select identity, colour and anticolour.

void SigmaHqqbar2QQbar::setIdColAcol() {

  // Extract incoming flavours, set outgoing ones.
  int id1 = inFluxPtr->id1();
  int id3 = (id1 > 0) ? idNew : -idNew;
  setId( id1, -id1, id3, -id3);

  // Colour flow topologies. Swap when antiquarks.
  setColAcol( 1, 0, 0, 2, 1, 0, 0, 2);
  if (id1 < 0) swapColAcol();

}


//**************************************************************************

// SigmaHqg2qgamma class.
// Cross section for q g -> q gamma.

//*********

// Initialize parton flux object. 
  
void SigmaHqg2qgamma::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state, e2-weighted, and parton densities.
  inFluxPtr = new InFluxqg();
  inFluxPtr->init( pdfAPtr, pdfBPtr);
  inFluxPtr->weightCharge(2);

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double SigmaHqg2qgamma::sigma2( double x1, double x2, double sH, 
  double tH) {

  // Standard Mandelstam variables and their squares.
  double uH = -(sH + tH); 
  double sH2 = sH * sH;
  double uH2 = uH * uH;

  // Calculate kinematics dependence.
  sigUS = (1./3.) * (sH2+uH2)/(-sH*uH);

  // Calculate running alpha_strong, using scale pT2 = tH*uH/sH.
  double pT2 = tH * uH / sH;
  double alphaS = alphaScalc.alphaS(pT2);  

  // Also use pT2 as PDF scale. Evaluate incoming parton densities.
  q2pdf = pT2;
  double pdf = inFluxPtr->flux( x1, x2, q2pdf);

  // Answer.
  return CONVERT2MB * (M_PI/sH2) * alphaS * alphaEM * sigUS * pdf;  

}

//*********

// Select identity, colour and anticolour.

void SigmaHqg2qgamma::setIdColAcol() {

  // Extract incoming flavours 
  int id1 = inFluxPtr->id1();
  int id2 = inFluxPtr->id2();

  // Construct outgoing flavours.
  int id3 = (id1 == 21) ? 22 : id1;
  int id4 = (id2 == 21) ? 22 : id2;
  setId( id1, id2, id3, id4);

  // Colour flow topology. Swap if first is gluon, or when antiquark.
  setColAcol( 1, 0, 2, 1, 2, 0, 0, 0);
  if (id1 == 21) swap1234();
  if (id1 < 0 || id2 < 0) swapColAcol();

}

//**************************************************************************

// SigmaHqqbar2ggamma class.
// Cross section for q qbar -> g gamma.

//*********

// Initialize parton flux object. 
  
void SigmaHqqbar2ggamma::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state, e2-weighted, and parton densities.
  inFluxPtr = new InFluxqqbarSame();
  inFluxPtr->init( pdfAPtr, pdfBPtr);
  inFluxPtr->weightCharge(2);

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double SigmaHqqbar2ggamma::sigma2( double x1, double x2, double sH, 
  double tH) {

  // Standard Mandelstam variables and their squares.
  double uH = -(sH + tH); 
  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

  // Calculate kinematics dependence.
  double sigTU = (8./9.) * (tH2+uH2)/(tH*uH);

  // Calculate running alpha_strong, using scale pT2 = tH*uH/sH.
  double pT2 = tH * uH / sH;
  double alphaS = alphaScalc.alphaS(pT2);  

  // Also use pT2 as PDF scale. Evaluate incoming parton densities.
  q2pdf = pT2;
  double pdf = inFluxPtr->flux( x1, x2, q2pdf);

  // Answer.
  return CONVERT2MB * (M_PI/sH2) * alphaS * alphaEM * sigTU * pdf;

}

//*********

// Select identity, colour and anticolour.

void SigmaHqqbar2ggamma::setIdColAcol() {

  // Extract incoming flavours; outgoing ones trivial.
  int id1 = inFluxPtr->id1();
  int id2 = inFluxPtr->id2();
  setId( id1, id2, 21, 22);

  // One colour flow topology. Swap if first is antiquark.
  setColAcol( 1, 0, 0, 2, 1, 2, 0, 0);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// SigmaHgg2ggamma class.
// Cross section for g g -> g gamma.
// Proceeds through a quark box, by default using 5 massless quarks.


//*********

// Initialize parton flux object. 
  
void SigmaHgg2ggamma::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxgg();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

  // Calculate charge factor from the allowed quarks in the box. 
  int nQuarkInLoop = Settings::mode("SigmaHat:nQuarkInLoop");
  chargeSum = - 1./3. + 2./3. - 1./3.;
  if (nQuarkInLoop >= 4) chargeSum += 2./3.;
  if (nQuarkInLoop >= 5) chargeSum -= 1./3.;
  if (nQuarkInLoop >= 6) chargeSum += 2./3.;

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double SigmaHgg2ggamma::sigma2( double x1, double x2, double sH, 
  double tH) {

  // Standard Mandelstam variables and their squares.
  double uH = -(sH + tH); 
  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

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
  

  // Calculate running alpha_strong, using scale pT2 = tH*uH/sH.
  double pT2 = tH * uH / sH;
  double alphaS = alphaScalc.alphaS(pT2);  

  // Also use pT2 as PDF scale. Evaluate incoming parton densities.
  q2pdf = pT2;
  double pdf = inFluxPtr->flux( x1, x2, q2pdf);

  // Answer.
  return CONVERT2MB * (5. / (192. * M_PI * sH2)) * pow2(chargeSum) 
    * pow3(alphaS) * alphaEM * sigBox * pdf;

}

//*********

// Select identity, colour and anticolour.

void SigmaHgg2ggamma::setIdColAcol() {

  // Flavours and colours are trivial.
  setId( 21, 21, 21, 22);
  setColAcol( 1, 2, 2, 3, 1, 3, 0, 0);
  if (Rndm::flat() > 0.5) swapColAcol();

}

//**************************************************************************

// SigmaHqqbar2gammagamma class.
// Cross section for q qbar -> gamma gamma.

//*********

// Initialize parton flux object. 
  
void SigmaHqqbar2gammagamma::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state, e4-weighted, and parton densities.
  inFluxPtr = new InFluxqqbarSame();
  inFluxPtr->init( pdfAPtr, pdfBPtr);
  inFluxPtr->weightCharge(4);

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double SigmaHqqbar2gammagamma::sigma2( double x1, double x2, double sH, 
  double tH) {

  // Standard Mandelstam variables and their squares.
  double uH = -(sH + tH); 
  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

  // Calculate kinematics dependence. Colour factor for quarks in.
  sigTU = 2. * (tH2+uH2)/(tH*uH);
  double colFac = 1./3.;

  // Calculate running alpha_strong, using scale pT2 = tH*uH/sH.
  double pT2 = tH * uH / sH;
  double alphaEM2 = alphaEM * alphaEM;

  // Also use pT2 as PDF scale. Evaluate incoming parton densities.
  q2pdf = pT2;
  double pdf = inFluxPtr->flux( x1, x2, q2pdf);

  // Answer contains factor 1/2 from identical photons.
  return CONVERT2MB * (M_PI/sH2) * alphaEM2 * 0.5 * sigTU * colFac * pdf;

}

//*********

// Select identity, colour and anticolour.

void SigmaHqqbar2gammagamma::setIdColAcol() {

  // Extract incoming flavours; outgoing ones trivial.
  int id1 = inFluxPtr->id1();
  int id2 = inFluxPtr->id2();
  setId( id1, id2, 22, 22);

  // One colour flow topology. Swap if first is antiquark.
  setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}
//**************************************************************************

// SigmaHgg2gammagamma class.
// Cross section for g g -> gamma gamma.
// Proceeds through a quark box, by default using 5 massless quarks.


//*********

// Initialize parton flux object. 
  
void SigmaHgg2gammagamma::initFlux( PDF* pdfAPtr, PDF* pdfBPtr) {

  // Set up for q qbar initial state and parton densities.
  inFluxPtr = new InFluxgg();
  inFluxPtr->init( pdfAPtr, pdfBPtr);

  // Calculate charge factor from the allowed quarks in the box. 
  int nQuarkInLoop = Settings::mode("SigmaHat:nQuarkInLoop");
  charge2Sum = 1./9. + 4./9. + 1./9.;
  if (nQuarkInLoop >= 4) charge2Sum += 4./9.;
  if (nQuarkInLoop >= 5) charge2Sum += 1./9.;
  if (nQuarkInLoop >= 6) charge2Sum += 4./9.;

} 

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double SigmaHgg2gammagamma::sigma2( double x1, double x2, double sH, 
  double tH) {

  // Standard Mandelstam variables and their squares.
  double uH = -(sH + tH); 
  double sH2 = sH * sH;
  double tH2 = tH * tH;
  double uH2 = uH * uH;

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
  

  // Calculate running alpha_strong, using scale pT2 = tH*uH/sH.
  double pT2 = tH * uH / sH;
  double alphaS = alphaScalc.alphaS(pT2);  

  // Also use pT2 as PDF scale. Evaluate incoming parton densities.
  q2pdf = pT2;
  double pdf = inFluxPtr->flux( x1, x2, q2pdf);

  // Answer contains factor 1/2 from identical photons.
  return CONVERT2MB * (0.5 / (16. * M_PI * sH2)) * pow2(charge2Sum) 
    * pow2(alphaS) * pow2(alphaEM) * sigBox * pdf;

}

//*********

// Select identity, colour and anticolour.

void SigmaHgg2gammagamma::setIdColAcol() {

  // Flavours and colours are trivial.
  setId( 21, 21, 22, 22);
  setColAcol( 1, 2, 2, 1, 0, 0, 0, 0);

}

//**************************************************************************

} // end namespace Pythia8
