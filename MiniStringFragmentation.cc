// Function definitions (not found in the header) for the .
// MiniStringFragmentation class
// Copyright © 2005 Torbjörn Sjöstrand

#include "MiniStringFragmentation.h"

namespace Pythia8 {

//**************************************************************************

// The MiniStringFragmentation class.

//*********
 
// Definitions of static variables.
// (Values will be overwritten in initStatic call, so are purely dummy.)

int MiniStringFragmentation::nTryMass = 2;
double MiniStringFragmentation::sigma = 0.35;
double MiniStringFragmentation::sigma2Had = 0.245;
double MiniStringFragmentation::bLund = 0.58;

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// After one-body fragmentation failed, try two-body once more. 
const int MiniStringFragmentation::NTRYLASTRESORT = 100;

// To avoid division by zero one must have sigma > 0.
const double MiniStringFragmentation::SIGMAMIN = 0.01;

// Loop try to combine available endquarks to valid hadron. 
const int MiniStringFragmentation::NTRYFLAV = 10;

//*********

// Initialize parameters of ministring fragmentation.

void MiniStringFragmentation::initStatic() {

  // Initialize the MiniStringFragmentation class proper.
  nTryMass = Settings::mode("MiniStringFragmentation:nTry");
  sigma = Settings::parameter("StringPT:sigma");
  sigma2Had = 2. * pow2( max( SIGMAMIN, sigma) );

  // Initialize the b parameter of the z spectrum, used when joining jets.
  bLund = Settings::parameter("StringZ:bLund");

}

//*********

// Do the fragmentation: driver routine.
  
bool MiniStringFragmentation::fragment(int iSub, ColConfig& colConfig, 
  Event& event) {

  // Read in info on system to be treated.
  id1 = event[ colConfig.front(iSub) ].id();
  id2 = event[ colConfig.back(iSub) ].id(); 
  pSum = colConfig[iSub].pSum;
  mSum = colConfig[iSub].mass;
  m2Sum = mSum*mSum;
  isClosed = colConfig[iSub].isClosed;
  iParton.resize(0);
  for (int i = 0; i < colConfig.iPartonSize(iSub); ++i) 
    iParton.push_back(colConfig.iParton(iSub, i));  

  // First try to produce two particles from the system.
  if (ministring2two( nTryMass, event)) return true;  

  // If this fails, then form one hadron and shuffle momentum.
  if (ministring2one( iSub, colConfig, event)) return true;  

  // If also this fails, then try harder to produce two particles.
  if (ministring2two( NTRYLASTRESORT, event)) return true;  

  // Else complete failure.
  return false;

}

//*********

  // Attempt to produce two particles from the ministring.
  
bool MiniStringFragmentation::ministring2two( int nTry, Event& event) {

  // Properties of the produced hadrons.
  int idHad1 = 0;
  int idHad2 = 0;
  double mHad1 = 0.;
  double mHad2 = 0.;
  double mHadSum = 0.;

  // Allow a few attempts to find a particle pair with low enough masses.
  for (int iTry = 0; iTry < nTry; ++iTry) {

    // For closed gluon loop need to pick an initial flavour.
    if (isClosed) do {
      int idTry = (Rndm::flat() < 0.5) ? 1 : 2;
      idTry = flavSel.pick( idTry);
      id1 = flavSel.pick( idTry);
      id2 = -id1;
    } while (id1 == 0);
   
    // Create a new q qbar flavour to form two hadrons. 
    // Start from a diquark, if any.
    do {
      int id3 = (abs(id1) > 8 || (abs(id2) < 9 && Rndm::flat() < 0.5) )
        ? flavSel.pick( id1) : -flavSel.pick( id2);
      idHad1 = flavSel.combine( id1, id3);
      idHad2 = flavSel.combine( id2, -id3);
    } while (idHad1 == 0 || idHad2 == 0);

    // Check whether the mass sum fits inside the available phase space.  
    mHad1 = ParticleDataTable::mass(idHad1);
    mHad2 = ParticleDataTable::mass(idHad2);
    mHadSum = mHad1 + mHad2;
    if (mHadSum < mSum) break;
  } 
  if (mHadSum >= mSum) return false;

  // Define an effective two-parton string, by splitting intermediate
  // gluon momenta in proportion to their closeness to either endpoint.
  Vec4 pSum1 = event[ iParton.front() ].p();
  Vec4 pSum2 = event[ iParton.back() ].p();
  if (iParton.size() > 2) {
    Vec4 pEnd1 = pSum1;
    Vec4 pEnd2 = pSum2;
    Vec4 pEndSum = pEnd1 + pEnd2; 
    for (int i = 1; i < int(iParton.size()) - 1 ; ++i) {
      Vec4 pNow = event[ iParton[i] ].p();
      double ratio = (pEnd2 * pNow) / (pEndSum * pNow);
      pSum1 += ratio * pNow;
      pSum2 += (1. - ratio) * pNow;
    }
  }

  // Set up a string region based on the two effective endpoints.
  StringRegion region;
  region.setUp( pSum1, pSum2);

  // Generate an isotropic decay in the ministring rest frame, 
  // suppressed at large pT by a fragmentation pT Gaussian.
  double pAbs2 = 0.25 * ( pow2(m2Sum - mHad1*mHad1 - mHad2*mHad2)
    - pow2(2. * mHad1 * mHad2) ) / m2Sum; 
  double pT2 = 0.;
  do {
    double cosTheta = Rndm::flat();
    if (sigma < SIGMAMIN) cosTheta = 1.;
    pT2 = (1. - pow2(cosTheta)) * pAbs2;
  } while ( exp( -pT2 / sigma2Had) < Rndm::flat() ); 

  // Construct the forward-backward asymmetry of the two particles.
  double mT21 = mHad1*mHad1 + pT2;
  double mT22 = mHad2*mHad2 + pT2;
  double lambda = sqrtpos( pow2(m2Sum  - mT21 - mT22) - 4. * mT21 * mT22 );
  double probReverse = 1. / (1. + exp( min( 50., bLund * lambda) ) ); 

  // Construct kinematics, as viewed in the transverse rest frame. 
  double xpz1 = 0.5 * lambda/ m2Sum;
  if (probReverse > Rndm::flat()) xpz1 = -xpz1; 
  double xmDiff = (mT21 - mT22) / m2Sum;
  double xe1 = 0.5 * (1. + xmDiff);
  double xe2 = 0.5 * (1. - xmDiff ); 

  // Distribute pT isotropically in angle.
  double phi = 2. * M_PI * Rndm::flat();
  double pT = sqrt(pT2);
  double px = pT * cos(phi);
  double py = pT * sin(phi);

  // Translate this into kinematics in the string frame.
  Vec4 pHad1 = region.pHad( xe1 + xpz1, xe1 - xpz1,  px,  py);
  Vec4 pHad2 = region.pHad( xe2 - xpz1, xe2 + xpz1, -px, -py);

  // Add produced particles to the event record.
  int iFirst = event.append( idHad1, 82, iParton.front(), iParton.back(), 
    0, 0, 0, 0, pHad1, mHad1);
  int iLast = event.append( idHad2, 82, iParton.front(), iParton.back(), 
    0, 0, 0, 0, pHad2, mHad2);

  // Set decay vertex when this is displaced.
  if (event[iParton.front()].hasVertex()) {
    Vec4 vDec = event[iParton.front()].vDec();
    event[iFirst].vProd( vDec );
    event[iLast].vProd( vDec );
  }

  // Set lifetime of hadrons.
  event[iFirst].tau( event[iFirst].tau0() * Rndm::exp() );
  event[iLast].tau( event[iLast].tau0() * Rndm::exp() );

  // Mark original partons as hadronized and set their daughter range.
  for (int i = 0; i < int(iParton.size()); ++i) {
    event[ iParton[i] ].statusNeg();
    event[ iParton[i] ].daughters(iFirst, iLast);
  }    

  // Successfully done.
  return true;

}

//*********

// Attempt to produce one particles from a ministring.
// Current algorithm: find the system with largest invariant mass
// relative to the existing one, and boost that system appropriately.
// Try more sophisticated alternatives later?? (Z0 mass shifted??)
// Also, if problems, attempt several times to obtain closer mass match??
  
bool MiniStringFragmentation::ministring2one( int iSub, 
  ColConfig& colConfig, Event& event) {

  // Cannot handle qq + qbarqbar system. 
  if (abs(id1) > 100 && abs(id2) > 100) return false;

  // For closed gluon loop need to pick an initial flavour.
  if (isClosed) do {
    int idTry = (Rndm::flat() < 0.5) ? 1 : 2;
    idTry = flavSel.pick( idTry);
    id1 = flavSel.pick( idTry);
    id2 = -id1;
  } while (abs(id1) > 100);

  // Select hadron flavour from available quark flavours.
  int idHad = 0;
  for (int iTryFlav = 0; iTryFlav < NTRYFLAV; ++iTryFlav) {
    idHad = flavSel.combine( id1, id2);
    if (idHad != 0) break;
  } 
  if (idHad == 0) return false;

  // Find mass.  
  double mHad = ParticleDataTable::mass(idHad);
  
  // Find the untreated parton system which combines to the largest 
  // squared mass above mimimum required. 
  int iMax = -1;
  double deltaM2 = mHad*mHad - mSum*mSum; 
  double delta2Max = 0.;
  for (int iRec = iSub + 1; iRec < colConfig.size(); ++iRec) {
    double delta2Rec = 2. * (pSum * colConfig[iRec].pSum) - deltaM2 
      - 2. * mHad * colConfig[iRec].mass; 
    if (delta2Rec > delta2Max) { iMax = iRec; delta2Max = delta2Rec;}
  }
  if (iMax == -1) return false;  

  // Construct kinematics of the hadron and recoiling system. 
  Vec4& pRec = colConfig[iMax].pSum;
  double mRec = colConfig[iMax].mass;
  double vecProd = pSum * pRec; 
  double coefOld = mSum*mSum + vecProd;
  double coefNew = mHad*mHad + vecProd;
  double coefRec = mRec*mRec + vecProd;
  double coefSum = coefOld + coefNew;
  double sHat = coefOld + coefRec;
  double root = sqrtpos( (pow2(coefSum) - 4. * sHat * mHad*mHad)
    / (pow2(vecProd) - pow2(mSum * mRec)) );
  double k2 = 0.5 * (coefOld * root - coefSum) / sHat;
  double k1 = (coefRec * k2 + 0.5 * deltaM2) / coefOld;
  Vec4 pHad = (1. + k1) * pSum - k2 * pRec;
  Vec4 pRecNew = (1. + k2) * pRec - k1 * pSum;
  
  // Add the produced particle to the event record.
  int iHad = event.append( idHad, 81, iParton.front(), iParton.back(), 
    0, 0, 0, 0, pHad, mHad);

  // Set decay vertex when this is displaced.
  if (event[iParton.front()].hasVertex()) {
    Vec4 vDec = event[iParton.front()].vDec();
    event[iHad].vProd( vDec );
  }

  // Set lifetime of hadron.
  event[iHad].tau( event[iHad].tau0() * Rndm::exp() );

  // Mark original partons as hadronized and set their daughter range.
  for (int i = 0; i < int(iParton.size()); ++i) {
    event[ iParton[i] ].statusNeg();
    event[ iParton[i] ].daughters(iHad, iHad);
  }    
   
  // Copy down recoiling system, with boosted momentum. Update current partons.
  RotBstMatrix M;
  M.bst(pRec, pRecNew); 
  for (int i = 0; i < colConfig[iMax].size(); ++i) {
    int iOld = colConfig.iParton(iMax, i);
    // Do not touch negative iOld = beginning of new junction leg.
    if (iOld >= 0) {
      int iNew = event.copy(iOld, 72);
      event[iNew].rotbst(M);
      colConfig.iParton(iMax, i) = iNew;
    }
  }
  colConfig[iMax].pSum = pRecNew;
  colConfig[iMax].isCollected = true;

  // Successfully done.
  return true;

}

//**************************************************************************

} // end namespace Pythia8
