// Function definitions (not found in the header) for the 
// InFlux class, and classes derived from it.
// Copyright C 2007 Torbjorn Sjostrand

#include "InFlux.h"

namespace Pythia8 {

//**************************************************************************

// InFlux class.
// Base class for the combined incoming parton flux.

//*********
 
// Definitions of static variables and functions.
// (Values will be overwritten in initStatic call, so are purely dummy.)

int  InFlux::nQuark = 5;
bool InFlux::showChannels = false;

// Pointers to the parton densities of the incoming beams.
BeamParticle* InFlux::beamAPtr; 
BeamParticle* InFlux::beamBPtr; 

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Map fermions onto their first-generation equivalents.
const int InFlux::MAPTOFIRST[40] = {0, 0, -12, -11, -12, -11, -12, -11, 
  -12, -11, 0, 0, -2, -1, -2, -1, -2, -1, -2, -1, 0, 1, 2, 1, 2, 1, 2, 
  1, 2, 0, 0, 11, 12, 11, 12, 11, 12, 11, 12, 0}; 

//*********

// Initialize static data members.

void InFlux::initStatic() {

  // Maximum incoming quark flavour.
  nQuark       = Settings::mode("inFlux:nQuark");

  // Print or not initialization data.
  showChannels = Settings::flag("inFlux:showChannels");

}

//*********

// Weight channels by two powers of e.

void InFlux::weightCharge2() {

  for (int i = 0; i < sizePair(); ++i) {
    int idAabs = abs(inPair[i].idA);
    int idBabs = abs(inPair[i].idB); 
    int id = (idAabs < 19) ? idAabs : idBabs;
    if (id < 9) inPair[i].fixWeight *= (id%2 == 0) ? 4./9. : 1./9.;
    else if (id%2 == 0) inPair[i].fixWeight = 0.; 
  }

}

//*********

// Weight channels by two powers of CKM matrix elements.
// Currently assume already known instate is q qbar'??

void InFlux::weightCKM2() {

  for (int i = 0; i < sizePair(); ++i) {

    // Remove impossible combinations.
    int idA = inPair[i].idA;
    int idB = inPair[i].idB; 
    if (idA * idB > 0 || (idA + idB)%2 == 0) inPair[i].fixWeight = 0.;

    // Read off CKM value.
    else inPair[i].fixWeight *= VCKM::V2id(idA, idB);
  }

}

//*********

// Weight channels by sum of squared CKM matrix element 
// for initial flavour(s). Three possibilities implemented:
// mode = 1 : CKM weight on one incoming side.
//      = 2 : decoupled CKM weights on both incoming sides.
//      = 3 : coupled CKM weights on both incoming sides, as 
//          consistent with t-channel W exhange (e.g. W+W- -> h0). 
//      = 4 : as 3, but with specified flavour produced on one side. 


void InFlux::weightCKM2sum(int mode, int idQ) {

  // Loop over existing inchannels and read instate.
  for (int i = 0; i < sizePair(); ++i) {
    int idA = inPair[i].idA;
    int idB = inPair[i].idB;
    int idAabs = abs(idA);
    int idBabs = abs(idB);

    // mode 1: only f g and f gamma allowed.
    if (mode == 1) {
      int idMin = min( idAabs, idBabs);
      int idMax = max( idAabs, idBabs);
      if (idMax !=21 && idMax != 22) inPair[i].fixWeight = 0.;
      else inPair[i].fixWeight *= VCKM::V2sum(idMin);
    
    // mode 2: decoupled f f' incoming state.
    } else if (mode == 2) {
      if (idAabs > 20 || idBabs > 20) inPair[i].fixWeight = 0.;   
      else inPair[i].fixWeight *= VCKM::V2sum(idA) *  VCKM::V2sum(idB);

    // mode 3: coupled f f' incoming state (like t-channel W).  
    } else if (mode == 3) {
      if (idAabs > 20 || idBabs > 20) inPair[i].fixWeight = 0.; 
      else if ( (idAabs%2 == idBabs%2 && idA * idB < 0)
        || (idAabs%2 != idBabs%2 && idA * idB > 0) ) 
        inPair[i].fixWeight *= VCKM::V2sum(idA) *  VCKM::V2sum(idB);
      else inPair[i].fixWeight = 0.;   

    // mode 4: coupled f f' incoming state, on one side to specified flavour.  
    } else if (mode == 4) {
      if (idAabs > 20 || idBabs > 20) inPair[i].fixWeight = 0.; 
      else if ( (idA * idB < 0 && (idAabs + idBabs)%2 == 0)
        || (idA * idB > 0 && (idAabs + idBabs)%2 == 1) ) {
        if ( (idQ + idAabs)%2 == 1 && (idQ + idBabs)%2 == 1 )
          inPair[i].fixWeight *= ( VCKM::V2id(idA, idQ) *  VCKM::V2sum(idB)
            + VCKM::V2sum(idA) * VCKM::V2id(idB, idQ) );
        else if ( (idQ + idAabs)%2 == 1) 
          inPair[i].fixWeight *= VCKM::V2id(idA, idQ) * VCKM::V2sum(idB);
        else if ( (idQ + idBabs)%2 == 1) 
          inPair[i].fixWeight *= VCKM::V2sum(idA) * VCKM::V2id(idB, idQ);
        else inPair[i].fixWeight = 0.; 
      }
      else inPair[i].fixWeight = 0.;   
    }
  }

}

//*********

// Weight by inverse colour factor in annihilation graphs: 
// 1/3 for q qbar, 1/8 for g g.

void InFlux::weightInvCol() {

  for (int i = 0; i < sizePair(); ++i) {

    // Incoming q qbar pair or g g pair.
    int idA = inPair[i].idA;
    int idB = inPair[i].idB; 
    if (idA * idB < 0 && abs(idA) < 10 && abs(idB) < 10) 
      inPair[i].fixWeight /= 3.;
    if (idA == 21 && idB == 21)  inPair[i].fixWeight /= 8.;
  }

}

//*********

// Weight by a factor of 2 for neutrinos: since only one spin state
// no spin averaging factor to be applied in some processes.

void InFlux::weightNeutrinoSpin() {

  for (int i = 0; i < sizePair(); ++i) {

    // Look for incoming neutrinos.
    int idAabs = abs(inPair[i].idA);
    if (idAabs == 12 || idAabs == 14 || idAabs == 16 || idAabs == 18)
      inPair[i].fixWeight *= 2.;
    int idBabs = abs(inPair[i].idB);
    if (idBabs == 12 || idBabs == 14 || idBabs == 16 || idBabs == 18)
      inPair[i].fixWeight *= 2.;
  }

}

//*********

// Weight by channel-specific factors.

void InFlux::weightFixed(int id1, int id2, double nowWeight, 
  bool flipSide, bool conjugate, bool allGen) {

  // Define charge-conjugate instate. Check whether to be used. 
  int id1Anti = (abs(id1) < 20) ? -id1 : id1;
  int id2Anti = (abs(id2) < 20) ? -id2 : id2;
  bool checkAnti = conjugate && ( (id1Anti != id1) || (id2Anti != id2) ); 

  // Loop over channels and read flavours.
  for (int i = 0; i < sizePair(); ++i) {
    int idA = inPair[i].idA;
    int idB = inPair[i].idB; 

    // Default option: treat all generations same way as first one.
    if (allGen) {
      if (abs(idA) < 20) idA = MAPTOFIRST[idA + 20];    
      if (abs(idB) < 20) idB = MAPTOFIRST[idB + 20];    
    }

    // If channels match then multiply weights by input number.
    if ( idA == id1 && idB == id2 ) 
      inPair[i].fixWeight *= nowWeight;
    else if ( flipSide && idA == id2 && idB == id1 ) 
      inPair[i].fixWeight *= nowWeight;
    else if ( checkAnti &&  idA == id1Anti && idB == id2Anti ) 
      inPair[i].fixWeight *= nowWeight;
    else if ( flipSide && checkAnti && idA == id2Anti && idB == id1Anti ) 
      inPair[i].fixWeight *= nowWeight;
  }

}

//*********

// Remove empty channels and optionally list remaining ones.

void InFlux::checkChannels(string processName, ostream& os) {

  // Remove empty channels, i.e. fix weight = 0, from list.
  vector<InPair>::iterator iterPair;
  for (iterPair = inPair.begin(); iterPair != inPair.end(); ++iterPair) 
  if (iterPair->fixWeight <= 0.) {
    inPair.erase(iterPair); 
    --iterPair;
  }   

  // Done if not to list flavours and channels.
  if (!showChannels) return;

  // List process name and allowed incoming flavours in beams.
  os << "\n InFlux initialization of process " << processName << "\n";
  os << " Allowed flavours in beam A:";
  for (int i = 0; i < int(inBeamA.size()); ++i) 
    os << "  " << inBeamA[i].id;
  os << "\n"; 
  os << " Allowed flavours in beam B:";
  for (int i = 0; i < int(inBeamB.size()); ++i) 
    os << "  " << inBeamB[i].id;
  os << "\n"; 

  // List allowed channels and fixed initialization weights.
  os << " Allowed channels of two incoming flavours with associated "
     << "fixed initial weight:\n";
  os << "    idA idB     weight    idA idB     weight    idA idB     "
     << "weight    idA idB     weight \n" << scientific << setprecision(3);
  int nList = inPair.size();
  int nLine = (nList + 3) / 4;
  for (int iList = 0; iList < 4 * nLine; ++iList) {
    if (iList > 0 && iList%4 == 0) os << "\n";
    int i = nLine * iList - (4 * nLine - 1) * (iList/4);
    if (i < nList) os << setw(7) << inPair[i].idA << setw(4) 
      << inPair[i].idB << setw(11) << inPair[i].fixWeight;
  }
  os << endl;

}

//*********

// Weight by channel-specific factors.

void InFlux::weightInState(int id1, int id2, double nowWeight, 
  bool flipSide, bool conjugate, bool allGen) {

  // Define charge-conjugate instate. Check whether to be used. 
  int id1Anti = (abs(id1) < 20) ? -id1 : id1;
  int id2Anti = (abs(id2) < 20) ? -id2 : id2;
  bool checkAnti = conjugate && ( (id1Anti != id1) || (id2Anti != id2) ); 

  // Loop over channels and read flavours.
  for (int i = 0; i < sizePair(); ++i) {
    int idA = inPair[i].idA;
    int idB = inPair[i].idB; 

    // Default option: treat all generations same way as first one.
    if (allGen) {
      if (abs(idA) < 20) idA = MAPTOFIRST[idA + 20];    
      if (abs(idB) < 20) idB = MAPTOFIRST[idB + 20];    
    }

    // If channels match then multiply weights by input number.
    if ( idA == id1 && idB == id2 ) 
      inPair[i].varWeight = nowWeight;
    else if ( flipSide && idA == id2 && idB == id1 ) 
      inPair[i].varWeight = nowWeight;
    else if ( checkAnti &&  idA == id1Anti && idB == id2Anti ) 
      inPair[i].varWeight = nowWeight;
    else if ( flipSide && checkAnti && idA == id2Anti && idB == id1Anti ) 
      inPair[i].varWeight = nowWeight;
  }

}

//*********

// Calculate products of parton densities for allowed combinations.

double InFlux::flux(double x1, double x2, double Q2) {

  // Evaluate and store the required parton densities.
  for (int j = 0; j < sizeBeamA(); ++j) 
    inBeamA[j].pdf = beamAPtr->xfHard( inBeamA[j].id, x1, Q2); 
  for (int j = 0; j < sizeBeamB(); ++j) 
    inBeamB[j].pdf = beamBPtr->xfHard( inBeamB[j].id, x2, Q2); 

  // Multiply on these densities for each of the allowed channels. Sum.
  fluxwtSum = 0.;
  for (int i = 0; i < sizePair(); ++i) {
    for (int j = 0; j < sizeBeamA(); ++j) 
    if (inPair[i].idA == inBeamA[j].id) {
      inPair[i].pdfA = inBeamA[j].pdf;
      break;
    }
    for (int j = 0; j < sizeBeamB(); ++j) 
    if (inPair[i].idB == inBeamB[j].id) {
      inPair[i].pdfB = inBeamB[j].pdf;
      break;
    }
    inPair[i].fluxWeight = inPair[i].fixWeight * inPair[i].varWeight 
      * inPair[i].pdfA * inPair[i].pdfB;
    fluxwtSum += inPair[i].fluxWeight;
  }
 
  // Done.
  return fluxwtSum;

}

//*********

// Pick one of the possible channels according to their weight.

void InFlux::pick() {

  // Pick channel. Extract channel flavours.
  double fluxwtRand =  fluxwtSum * Rndm::flat();
  for (int i = 0; i < sizePair(); ++i) {
    fluxwtRand -= inPair[i].fluxWeight;
    if (fluxwtRand <= 0.) {
      idNow1 = inPair[i].idA;
      idNow2 = inPair[i].idB;
      pdfNow1 = inPair[i].pdfA; 
      pdfNow2 = inPair[i].pdfB; 
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
  addBeamA(21);
  addBeamB(21);
 
  // The combined channels.
  addPair(21, 21);

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
    addBeamA(id);
    addBeamB(id);
  }

  // The combined channels.
  for (int id = -nQuark; id <= nQuark; ++id) 
  if (id != 0) {
    addPair(id, 21);
    addPair(21, id);
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
    addBeamA(id);
    addBeamB(id);
  }

  // The combined channels.
  for (int id1 = -nQuark; id1 <= nQuark; ++id1) 
  if (id1 != 0) 
  for (int id2 = -nQuark; id2 <= nQuark; ++id2) 
  if (id2 != 0 && id2 != id1) 
    addPair(id1, id2);

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
    addBeamA(id);
    addBeamB(id);
  }

  // The combined channels.
  for (int id1 = -nQuark; id1 <= nQuark; ++id1) 
  if (id1 != 0) 
  for (int id2 = -nQuark; id2 <= nQuark; ++id2) 
  if (id2 != 0 && id2 != id1 && id1 * id2 > 0) 
    addPair(id1, id2);

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
    addBeamA(id);
    addBeamB(id);
  }

  // The combined channels.
  for (int id = -nQuark; id <= nQuark; ++id) 
  if (id != 0) 
    addPair(id, id);

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
    addBeamA(id);
    addBeamB(id);
  }

  // The combined channels.
  for (int id1 = -nQuark; id1 <= nQuark; ++id1) 
  if (id1 != 0) 
  for (int id2 = -nQuark; id2 <= nQuark; ++id2) 
  if (id2 != 0 && id2 != -id1 && id1 * id2 < 0) 
    addPair(id1, id2);

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
    addBeamA(id);
    addBeamB(id);
  }

  // The combined channels.
  for (int id = -nQuark; id <= nQuark; ++id) 
  if (id != 0) 
    addPair(id, -id);

}

//**************************************************************************

// InFluxff class.
// Derived class for f f', f fbar', fbar fbar' incoming state.

//*********

// Fill arrays for a f f' incoming state at initialization.

void InFluxff::initChannels() {

  // If beams are leptons then they are also the colliding partons.
  if ( beamAPtr->isLepton() && beamBPtr->isLepton() ) {
    addBeamA( beamAPtr->id() );
    addBeamB( beamBPtr->id() );
    addPair( beamAPtr->id(), beamBPtr->id() );
    return;
  }

  // The beams individually.
  for (int id = -nQuark; id <= nQuark; ++id) 
  if (id != 0) {
    addBeamA(id);
    addBeamB(id);
  }

  // The combined channels.
  for (int id1 = -nQuark; id1 <= nQuark; ++id1) 
  if (id1 != 0) 
  for (int id2 = -nQuark; id2 <= nQuark; ++id2) 
  if (id2 != 0) 
    addPair(id1, id2);

}

//**************************************************************************

// InFluxffbarSame class.
// Derived class for f fbar antiparticle incoming state.

//*********

// Fill arrays for a f fbar incoming state at initialization.

void InFluxffbarSame::initChannels() {

  // If beams are antiparticle pair and leptons then also colliding partons.
  if ( beamAPtr->id() + beamBPtr->id() == 0 && beamAPtr->isLepton() ) {
    addBeamA( beamAPtr->id() );
    addBeamB( beamBPtr->id() );
    addPair( beamAPtr->id(), beamBPtr->id() );
    return;
  }

  // Else assume both to be hadrons, for better or worse.
  // The beams individually.
  for (int id = -nQuark; id <= nQuark; ++id) 
  if (id != 0) {
    addBeamA(id);
    addBeamB(id);
  }

  // The combined channels.
  for (int id = -nQuark; id <= nQuark; ++id) 
  if (id != 0) 
    addPair(id, -id);

}

//**************************************************************************

// InFluxffbarChg class.
// Derived class for f fbar' incoming states with net charge +-1.

//*********

// Fill arrays for a f fbar' incoming state at initialization.

void InFluxffbarChg::initChannels() {

  // If beams are leptons then also colliding partons.
  if ( beamAPtr->isLepton() && beamBPtr->isLepton() 
    && abs ( ParticleDataTable::chargeType(beamAPtr->id()) 
           + ParticleDataTable::chargeType(beamBPtr->id()) ) == 3 ) {
    addBeamA( beamAPtr->id() );
    addBeamB( beamBPtr->id() );
    addPair( beamAPtr->id(), beamBPtr->id() );
    return;
  }

  // The beams individually.
  for (int id = -nQuark; id <= nQuark; ++id) 
  if (id != 0) {
    addBeamA(id);
    addBeamB(id);
  }

  // The combined channels.
  for (int id1 = -nQuark; id1 <= nQuark; ++id1) 
  if (id1 != 0) 
  for (int id2 = -nQuark; id2 <= nQuark; ++id2) 
  if (id2 != 0 && id1 * id2 < 0 && (abs(id1) + abs(id2))%2 == 1) 
    addPair(id1, id2);

}

//**************************************************************************

// InFluxfgm class.
// Derived class for a f gamma incoming state.

//*********

// Fill arrays for a gamma gamma incoming state at initialization.

void InFluxfgm::initChannels() {

  // Fermion from incoming side A.
  if ( beamAPtr->isLepton() ) {
    addBeamA( beamAPtr->id() );
    addPair( beamAPtr->id(), 22);
  } else {  
    for (int id = -nQuark; id <= nQuark; ++id) 
    if (id != 0) {
      addBeamA(id);
      addPair(id, 22);
    }
  }

  // Fermion from incoming side B.
  if ( beamBPtr->isLepton() ) {
    addBeamB( beamBPtr->id() );
    addPair( 22, beamBPtr->id());
  } else {  
    for (int id = -nQuark; id <= nQuark; ++id) 
    if (id != 0) {
      addBeamB(id);
      addPair(22, id);
    }
  }

  // Photons in the beams.
  addBeamA(22);
  addBeamB(22);

}

//**************************************************************************

// InFluxgmgm class.
// Derived class for a gamma gamma incoming state.

//*********

// Fill arrays for a gamma gamma incoming state at initialization.

void InFluxgmgm::initChannels() {

  // The beams individually.
  addBeamA(22);
  addBeamB(22);
 
  // The combined channels.
  addPair(22, 22);

}

//**************************************************************************

} // end namespace Pythia8
