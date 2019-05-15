// Function definitions (not found in the header) for the 
// InFlux class, and classes derived from it.
// Copyright C 2006 Torbjorn Sjostrand

#include "InFlux.h"

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

// Weight channels by (two or four) powers of CKM matrix elements.
// Currently assume already known instate is q qbar'.

void InFlux::weightCKM(int ckmPow) {

  for (int i = 0; i < int(weightAB.size()); ++i) {

    // Remove impossible combinations.
    int idA = idPartonPairA[i];
    int idB = idPartonPairB[i]; 
    if (idA * idB > 0 || (idA + idB)%2 == 0) weightAB[i] *= 0.;

    // Find up-type and down-type quarks, and read off CKM value.
    else { 
      int idAabs = abs(idA);
      int idBabs = abs(idB);
      if (idAabs%2 == 1) swap(idAabs, idBabs);
      double ckm2 = VCKM::V2( (idAabs/2), (idBabs+1)/2 );
      if (ckmPow == 2) weightAB[i] *= ckm2;
      if (ckmPow == 4) weightAB[i] *= ckm2 * ckm2;
    }
  }

}

//*********

// Weight by inverse colour factor in annihilation graphs: 
// 1/3 for q qbar, 1/8 for g g.

void InFlux::weightInvCol() {

  for (int i = 0; i < int(weightAB.size()); ++i) {

    // Incoming q qbar pair or g g pair.
    int idA = idPartonPairA[i];
    int idB = idPartonPairB[i]; 
    if (idA * idB < 0 && abs(idA) < 10 && abs(idB) << 10) weightAB[i] /= 3.;
    if (idA == 21 && idB == 21)  weightAB[i] /= 8.;
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

// InFluxffbarChg class.
// Derived class for f fbar' incoming states with net charge +-1.

//*********

// Fill arrays for a f fbar' incoming state at initialization.
// Still needs to do for lepton beams??

void InFluxffbarChg::initChannels() {

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
  if (id2 != 0 && id1 * id2 < 0 && (abs(id1) + abs(id2))%2 == 1) {
    idPartonPairA.push_back(id1);
    idPartonPairB.push_back(id2);
    pdfPairA.push_back(0.);
    pdfPairB.push_back(0.);
    weightAB.push_back(1.);
    fluxweightAB.push_back(1.);
  }

}

//**************************************************************************

} // end namespace Pythia8
