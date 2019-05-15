// Function definitions (not found in the header) for the 
// supersymmetry simulation classes. 
// Copyright C 2007 Torbjorn Sjostrand

#include "SigmaSUSY.h"

namespace Pythia8 {

//**************************************************************************

// Sigma2qqbar2chi0chi0 class.
// Cross section for q qbar -> chi0_i chi0_j, i.e. neutralino pairs . 

//*********

// Initialize process. 
  
void Sigma2qqbar2chi0chi0::initProc() {

  // Construct neutralino id codes from ordering indices.
  id3                  = 1000022; 
  if (id3chi == 2) id3 = 1000023; 
  if (id3chi == 3) id3 = 1000025; 
  if (id3chi == 4) id3 = 1000035; 
  id4                  = 1000022; 
  if (id4chi == 2) id4 = 1000023; 
  if (id4chi == 3) id4 = 1000025; 
  if (id4chi == 4) id4 = 1000035; 

  // Construct name of process. 
  nameSave = "q qbar -> " + ParticleDataTable::name(id3) + " " 
    + ParticleDataTable::name(id4);

  // Set up couplings
  sin2W = CoupEW::sin2thetaW();
  mZ    = ParticleDataTable::m0(23);
  wZ    = ParticleDataTable::mWidth(23);

  // For future use, when full mixing is implemented:
  // Shorthand for SUSY couplings
  // By default, use the tan(beta) given in MINPAR(3)
  // If not found, use the running one in HMIX (or EXTPAR?)
  //double tanb = slha->minpar.exists(3) ? slha->minpar(3) : slha->hmix(2);
  //double sinW=sqrt(sin2W);
  //double cosW=sqrt(1.0-sin2W);
  //double sinb = sqrt(max(0.0,1.0-cosb*cosb));
  //double cosb = sqrt( 1.0 / (1.0 + tanb*tanb) );
  //SusyLesHouches::matrixblock<6> Su(slha->usqmix);
  //SusyLesHouches::matrixblock<6> Sd(slha->dsqmix);
  //SusyLesHouches::matrixblock<6> imSu(slha->imusqmix);
  //SusyLesHouches::matrixblock<6> imSd(slha->imdsqmix);  

  // Local complex copies of neutralino mixing matrix entries. 
  complex ni1( slha->nmix(id3chi,1), slha->imnmix(id3chi,1) );
  complex nj1( slha->nmix(id4chi,1), slha->imnmix(id4chi,1) );
  complex ni2( slha->nmix(id3chi,2), slha->imnmix(id3chi,2) );
  complex nj2( slha->nmix(id4chi,2), slha->imnmix(id4chi,2) );
  complex ni3( slha->nmix(id3chi,3), slha->imnmix(id3chi,3) );
  complex nj3( slha->nmix(id4chi,3), slha->imnmix(id4chi,3) );
  complex ni4( slha->nmix(id3chi,4), slha->imnmix(id3chi,4) );
  complex nj4( slha->nmix(id4chi,4), slha->imnmix(id4chi,4) );

  // Change to positive mass convention.
  complex iRot( 0., 1.);
  if (slha->mass(id3) < 0.) {
    ni1 *= iRot;
    ni2 *= iRot;
    ni3 *= iRot;
    ni4 *= iRot;
  };
  if (slha->mass(id4) < 0.) {
    nj1 *= iRot;
    nj2 *= iRot;
    nj3 *= iRot;
    nj4 *= iRot;
  };

  // Z chi_i chi_j 
  OL = -0.5 * ni3 * conj(nj3) + 0.5 * ni4 * conj(nj4);
  OR =  0.5 * conj(ni3) * nj3 - 0.5 * conj(ni4) * nj4;

  // Loop over possible (diagonal) incoming flavours  
  // Z q_{iq} q_{iq} (def with extra factor 2 compared to [Okun])
  for (int iq = 1; iq <= 5; ++iq) {
    LqqZ[iq] = CoupEW::lf(iq);
    RqqZ[iq] = CoupEW::rf(iq);
  }

}
//*********

// Initialize parton-flux object. 
  
void Sigma2qqbar2chi0chi0::initFlux() {

  // Set up for q qbar initial state.
  inFluxPtr = new InFluxqqbarSame();

  // Multiply by colour factor 1/3.
  inFluxPtr->weightInvCol();

  // Identical-fermion factor
  if (id3chi == id4chi) inFluxPtr->weightFixed(0.5);

}

//*********

// Evaluate d(sigmaHat)/d(tHat). 

double Sigma2qqbar2chi0chi0::sigmaHat() {

  // Common flavour-independent factor
  double sigma = (M_PI / sH2) * pow2(alpEM) / pow2(sin2W * (1 - sin2W)); 

  // Auxiliary factors for use below
  double ui = uH - s3;
  double uj = uH - s4;
  double ti = tH - s3;
  double tj = tH - s4;
  double sz = sH - pow2(mZ);
  double d  = pow2(sz) + pow2(mZ * wZ);
  complex propZ( sz / d, mZ * wZ / d);
 
  // Flavour-dependent factors. Sum over diagonal in-flavours.
  // (NB: will eventually need also class for non-diagonal)
  for (int iq = 1; iq <= 2; ++iq) {

    // Flavour-dependent kinematics-dependent couplings
    complex QuLL = LqqZ[iq] * OL/2.0 * propZ;
    complex QtLL = LqqZ[iq] * OR/2.0 * propZ;
    complex QuRR = RqqZ[iq] * OR/2.0 * propZ;
    complex QtRR = RqqZ[iq] * OL/2.0 * propZ;
    complex QuLR = 0.0;
    complex QtLR = 0.0;
    complex QuRL = 0.0;
    complex QtRL = 0.0;
  
    // Compute matrix element weight
    double weight = 0;
    // Average over separate helicity contributions
    // LL (ha = -1, hb = +1) (divided by 4 for average)            
    weight += norm(QuLL) * ui * uj + norm(QtLL) * ti * tj
      + 2 * real(conj(QuLL) * QtLL) * m3 * m4 * sH;
    // RR (ha =  1, hb = -1) (divided by 4 for average)        
    weight += norm(QtRR) * ti * tj + norm(QuRR) * ui * uj  
      + 2 * real(conj(QuRR) * QtRR) * m3 * m4 * sH;
    // RL (ha =  1, hb =  1) (divided by 4 for average)        
    weight += norm(QuRL) * ui * uj + norm(QtRL) * ti * tj
      - real(conj(QuRL) * QtRL) * (uH * tH - s3 * s4);
    // LR (ha = -1, hb = -1) (divided by 4 for average)        
    weight += norm(QuLR) * ui * uj + norm(QtLR) * ti * tj
      - real(conj(QuLR) * QtLR) * (uH * tH - s3 * s4);
    
    // Apply weight to inFlux. First-generation results automatically
    // extend also to second and third (so far, with no squark contribution).  
    inFluxPtr->weightInState( iq, -iq, weight);
  
  }
  // Answer, leaving out flavour-dependent pieces stored above.
  return CONVERT2MB * sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2qqbar2chi0chi0::setIdColAcol() {

  // Set flavours.
  setId( id1, id2, id3, id4);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 10) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else               setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

} // end namespace Pythia8
