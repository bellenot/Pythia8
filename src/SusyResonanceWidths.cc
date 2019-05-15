// ResonanceWidths.cc is a part of the PYTHIA event generator.
// Copyright (C) 2010 
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for 
// the ResonanceWidths class and classes derived from it.

#include "SusyResonanceWidths.h"
#include "SusyCouplings.h"
#include "PythiaComplex.h"

namespace Pythia8 {

//==========================================================================

// The SUSYResonanceWidths Class
// Derived class for SUSY resonances

const bool SUSYResonanceWidths::DEBUG = false;

//--------------------------------------------------------------------------

bool SUSYResonanceWidths::init(Info* infoPtrIn, Settings* settingsPtrIn,
   ParticleData* particleDataPtrIn, Couplings* couplingsPtrIn) {

  // Save pointers.
  infoPtr         = infoPtrIn;
  settingsPtr     = settingsPtrIn;
  particleDataPtr = particleDataPtrIn;
  coupSUSYPtr     = (couplingsPtrIn->isSUSY ? (CoupSUSY *) couplingsPtrIn: 0 );

  // No initialization necessary for SM-only
  if(!couplingsPtrIn->isSUSY) return true;

  // Minimal decaying-resonance width. Minimal phase space for meMode = 103.
  minWidth     = settingsPtr->parm("ResonanceWidths:minWidth");
  minThreshold = settingsPtr->parm("ResonanceWidths:minThreshold");

  // Pointer to particle species.
  particlePtr  = particleDataPtr->particleDataEntryPtr(idRes);
  if (particlePtr == 0) {
    infoPtr->errorMsg("Error in ResonanceWidths::init:"
      " unknown resonance identity code");   
    return false;
  }  

  // Generic particles should not have meMode < 100, but allow 
  // some exceptions: not used Higgses and not used Technicolor.
  if (idRes == 35 || idRes == 36 || idRes == 37 
    || idRes/1000000 == 3) isGeneric = false;

  // Resonance properties: antiparticle, mass, width
  hasAntiRes   = particlePtr->hasAnti();
  mRes         = particlePtr->m0();
  GammaRes     = particlePtr->mWidth();
  m2Res        = mRes*mRes;

  // For very narrow resonances assign fictitious small width.
  if (GammaRes < minWidth) GammaRes = 0.1 * minWidth;  
  GamMRat      = GammaRes / mRes;

  // Secondary decay chains by default all on.
  openPos      = 1.;
  openNeg      = 1.;

  // Allow option where on-shell width is forced to current value.
  doForceWidth = particlePtr->doForceWidth();
  forceFactor  = 1.;

  // Check that resonance OK.
  if (particlePtr == 0) infoPtr->errorMsg("Error in ResonanceWidths::init:"
      " unknown resonance identity code");   

  // Calculate various common prefactors for the current mass.
  mHat          = mRes;

  // Initialize constants used for a resonance.

  initConstants();
  calcPreFac(true);

  // Reset quantities to sum. Declare variables inside loop.
  double widTot = 0.; 
  double widPos = 0.;
  double widNeg = 0.;
  int    idNow, idAnti;
  double openSecPos, openSecNeg;

  // Loop over all decay channels. Basic properties of channel.
  for (int i = 0; i < particlePtr->sizeChannels(); ++i) {
    iChannel    = i;
    onMode      = particlePtr->channel(i).onMode();
    meMode      = particlePtr->channel(i).meMode();
    mult        = particlePtr->channel(i).multiplicity();
    widNow      = 0.;

    // Warn if not relevant meMode.
    if ( meMode < 0 || meMode > 103 || (isGeneric && meMode < 100) ) { 
      infoPtr->errorMsg("Error in ResonanceWidths::init:"
        " resonance meMode not acceptable"); 
    }

    // Check if decay table was read in via SLHA
    bool hasDecayTable = false;
    for(unsigned int iDec = 1; iDec < (coupSUSYPtr->slhaPtr)->decays.size(); iDec++)
      hasDecayTable = ((coupSUSYPtr->slhaPtr)->decays[iDec].getId() == abs(idRes));

    // Calculation of SUSY particle widths
    if (meMode == 103 && GammaRes > 0. &&
	(!settingsPtr->flag("SLHA:useDecayTable") || !hasDecayTable)) {
      // Read out information on channel: primarily use first two. 
      id1       = particlePtr->channel(i).product(0);
      id2       = particlePtr->channel(i).product(1);
      id1Abs    = abs(id1);
      id2Abs    = abs(id2);
       
      // Order first two in descending order of absolute values.
      if (id2Abs > id1Abs) {swap( id1, id2); swap( id1Abs, id2Abs);}

      // Allow for third product to be treated in derived classes.
      if (mult > 2) { 
        id3     = particlePtr->channel(i).product(2);
        id3Abs  = abs(id3);
        
        // Also order third into descending order of absolute values.
        if (id3Abs > id2Abs) {swap( id2, id3); swap( id2Abs, id3Abs);}
        if (id2Abs > id1Abs) {swap( id1, id2); swap( id1Abs, id2Abs);}
      }

      // Read out masses. Calculate two-body phase space.
      mf1       = particleDataPtr->m0(id1Abs);
      mf2       = particleDataPtr->m0(id2Abs);
      mr1       = pow2(mf1 / mHat);
      mr2       = pow2(mf2 / mHat);
      ps        = (mHat < mf1 + mf2 + MASSMARGIN) ? 0. 
                : sqrtpos( pow2(1. - mr1 - mr2) - 4. * mr1 * mr2 );
      if (mult > 2) {      
        mf3     = particleDataPtr->m0(id3Abs);
        mr3     = pow2(mf3 / mHat);
      }

      // Let derived class calculate width for channel provided.
      calcWidth(true);
    }

    // Channels with meMode >= 100 are calculated based on stored values.
    else widNow = GammaRes * particlePtr->channel(i).bRatio();
   
    // Find secondary open fractions of partial width.
    openSecPos  = 1.;
    openSecNeg  = 1.;
    if (widNow > 0.) for (int j = 0; j < mult; ++j) {
      idNow     = particlePtr->channel(i).product(j);
      idAnti    = (particleDataPtr->hasAnti(idNow)) ? -idNow : idNow;
      // Secondary widths not yet initialized for heavier states,
      // so have to assume unit open fraction there.
      if (idNow == 23 || abs(idNow) == 24 
	|| particleDataPtr->m0(abs(idNow)) < mRes) {
        openSecPos *= particleDataPtr->resOpenFrac(idNow); 
        openSecNeg *= particleDataPtr->resOpenFrac(idAnti);
      } 
    }

    // Store partial widths and secondary open fractions.
    particlePtr->channel(i).onShellWidth(widNow); 
    particlePtr->channel(i).openSec( idRes, openSecPos);  
    particlePtr->channel(i).openSec(-idRes, openSecNeg);  

    // Update sum over all channnels and over open channels only.
    widTot     += widNow;    
    if (onMode == 1 || onMode == 2) widPos += widNow * openSecPos;
    if (onMode == 1 || onMode == 3) widNeg += widNow * openSecNeg;
  }

  // If no decay channels are open then set particle stable and done.
  if (widTot < minWidth) { 
    particlePtr->setMayDecay(false, false);
    particlePtr->setMWidth(0., false);
    for (int i = 0; i < particlePtr->sizeChannels(); ++i) 
      particlePtr->channel(i).bRatio( 0., false);
    return true;
  }

  // Normalize branching ratios to unity.
  double bRatio;
  for (int i = 0; i < particlePtr->sizeChannels(); ++i) {
    bRatio      = particlePtr->channel(i).onShellWidth() / widTot;
    particlePtr->channel(i).bRatio( bRatio, false);
  }

  // Optionally force total width by rescaling of all partial ones.
  if (doForceWidth) {
    forceFactor = GammaRes / widTot;
    for (int i = 0; i < particlePtr->sizeChannels(); ++i)
      particlePtr->channel(i).onShellWidthFactor( forceFactor);
  } 

  // Else update newly calculated partial width.
  else {
    particlePtr->setMWidth(widTot, false);
    GammaRes    = widTot;
  }

  // Updated width-to-mass ratio. Secondary widths for open.
  GamMRat       = GammaRes / mRes;  
  openPos       = widPos / widTot;
  openNeg       = widNeg / widTot;

  // Done.
  return true;

}  



//==========================================================================

// The ResonanceSquark class
// Derived class for Squark resonances

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceSquark::initConstants() {

  // Locally stored properties and couplings.
  alpS  = coupSUSYPtr->alphaS(mHat * mHat );
  alpEM = coupSUSYPtr->alphaEM(mHat * mHat);
  s2W   = coupSUSYPtr->sin2W;
}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceSquark::calcPreFac(bool) {

  // Common coupling factors.
  preFac = 1.0 / (s2W * pow(mHat,3));

}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceSquark::calcWidth(bool) {

  // Squark type -- in u_i/d_i and generation
  int ksusy = 1000000;
  bool idown = (abs(idRes)%2 == 0 ? false : true);
  int isq = (abs(idRes)/ksusy == 2) ? 
    (abs(idRes)%10+1)/2 + 3: (abs(idRes)%10+1)/2;
  int isqgen = (abs(idRes)%10 + 1)/2;

  // Check that above threshold.
  if (ps == 0.) return;
  kinFac = (mHat * mHat - mf1 * mf1 -mf2 * mf2);
  lambda = lam(mHat*mHat, mf1*mf1, mf2*mf2);

  double fac = 0.0 , wid = 0.0;

  //Case 1: RPV decay
  if(id1Abs < 7 && id2Abs < 7){

    //Temporary till colour assignments sorted out
    //widNow = 0.0;
    //return;

    // Quark generations
    int iq1 = (id1Abs + 1)/2;
    int iq2 = (id2Abs + 1)/2;

    // Check for RPV UDD couplings
    if(!coupSUSYPtr->isUDD) {widNow = 0; return;}

    // ~q -> q_i + q_j

    fac = mHat;
    if(idown) {
      if ((id1Abs+id2Abs)%2 == 1){
	if(id1Abs%2==1)
	  for(int isq2 = 4; isq2 < 7; isq2++)
	    wid = coupSUSYPtr->rvUDD[iq2][iq1][isqgen] *
	      norm(coupSUSYPtr->Rdsq[isq][isq2]);
	else
	  for(int isq2 = 4; isq2 < 7; isq2++)
	    wid = coupSUSYPtr->rvUDD[iq1][iq2][isqgen] *
	      norm(coupSUSYPtr->Rdsq[isq][isq2]);
      }
    }
    else {
      if ((id1Abs+id2Abs)%2 != 0) widNow = 0.0;
      else
	for(int isq2 = 4; isq2 < 7; isq2++)
	  wid = coupSUSYPtr->rvUDD[isq][iq1][iq2] *
	    norm(coupSUSYPtr->Rusq[isq][isq2]);
    }
  }

  //Case 2: quark + gaugino (higgsino)
  else if (id1Abs > ksusy && id2Abs < 7) {
    
    int iq = (id2Abs + 1)/2;

    // ~q -> ~g + q
    if(id1Abs == 1000021 && idRes%10 == id2Abs) {
      fac = 2.0 / 3.0 * alpS *  preFac * sqrt(lambda);
      if(idown)
	wid = kinFac * (norm(coupSUSYPtr->LsddG[isq][iq])+norm(coupSUSYPtr->RsddG[isq][iq]))
	    - 4.0 * mf1 * mf2 * real(coupSUSYPtr->LsddG[isq][iq]*coupSUSYPtr->RsddG[isq][iq]);
      else
	wid = kinFac * (norm(coupSUSYPtr->LsuuG[isq][iq])+norm(coupSUSYPtr->RsuuG[isq][iq]))
	    - 4.0 * mf1 * mf2 * real(coupSUSYPtr->LsuuG[isq][iq]*coupSUSYPtr->RsuuG[isq][iq]);

    } 
    else 
      for(int i=1; i<6 ; i++){
	// ~q -> ~chi0 + q
	if(coupSUSYPtr->idNeut(i)==id1Abs && idRes%2 == id2Abs%2){
	  fac = alpEM *  preFac * sqrt(lambda)/ (4.0 * (1 - s2W));
	  if(idown)
	    wid = kinFac * (norm(coupSUSYPtr->LsddX[isq][iq][i]) + norm(coupSUSYPtr->RsddX[isq][iq][i]))
		- 4.0 * mf1 * mf2 * real(coupSUSYPtr->LsddX[isq][iq][i]*coupSUSYPtr->RsddX[isq][iq][i]);
	  else
	    wid = kinFac * (norm(coupSUSYPtr->LsuuX[isq][iq][i]) + norm(coupSUSYPtr->RsuuX[isq][iq][i]))
		- 4.0 * mf1 * mf2 * real(coupSUSYPtr->LsuuX[isq][iq][i]*coupSUSYPtr->RsuuX[isq][iq][i]);
	}

	// ~q -> chi- + q
	else if (i < 3 && coupSUSYPtr->idChar(i)==id1Abs && idRes%2 != id2Abs%2){

	  fac = alpEM *  preFac * sqrt(lambda)/ (4.0 * (1 - s2W));
	  if(idown)
	    wid = kinFac * (norm(coupSUSYPtr->LsduX[isq][iq][i]) + norm(coupSUSYPtr->RsduX[isq][iq][i]))
		- 4.0 * mf1 * mf2 * real(coupSUSYPtr->LsduX[isq][iq][i]*coupSUSYPtr->RsduX[isq][iq][i]);
	  else
	    wid = kinFac * (norm(coupSUSYPtr->LsudX[isq][iq][i]) + norm(coupSUSYPtr->RsudX[isq][iq][i]))
		- 4.0 * mf1 * mf2 * real(coupSUSYPtr->LsudX[isq][iq][i]*coupSUSYPtr->RsudX[isq][iq][i]);
	}
      }
  }
  
  //Case 3: ~q_i -> ~q_j + Z/W
  else if (id1Abs > ksusy && id1Abs%100 < 7 && (id2Abs == 23 || id2Abs == 24)){

    fac = alpEM * preFac/(16.0 * pow2(particleDataPtr->m0(id2Abs)) * (1.0 - s2W))
      * pow(lambda, 1.5);
    
    int isq2 = (id1Abs/ksusy == 2) ? (id1Abs%10+1)/2 + 3: (id1Abs%10+1)/2;

    if(id2Abs == 23 && id1Abs%2 == idRes%2){
      if(idown)
	wid = norm(coupSUSYPtr->LsdsdZ[isq][isq2] + coupSUSYPtr->RsdsdZ[isq][isq2]);
      else
	wid = norm(coupSUSYPtr->LsusuZ[isq][isq2] + coupSUSYPtr->RsusuZ[isq][isq2]);
    }
    else if (id2Abs == 24 && id1Abs%2 != idRes%2){
      if(idown)
	wid = norm(coupSUSYPtr->LsusdW[isq2][isq]);
      else
	wid = norm(coupSUSYPtr->LsusdW[isq][isq2]);
    }
  }

  widNow = fac * wid * ps;
  if(DEBUG) cout<<idRes<<":: id1:"<<id1Abs<<" id2:"<<id2Abs<<" Width: "<<widNow<<endl;
  return;
	
}

//==========================================================================

// The ResonanceGluino class
// Derived class for Gluino resonances

//--------------------------------------------------------------------------

// Initialize constants.

void ResonanceGluino::initConstants() {

  // Locally stored properties and couplings.
  alpS  = coupSUSYPtr->alphaS(mHat * mHat );
  return;
}

//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.

void ResonanceGluino::calcPreFac(bool) {
  // Common coupling factors.
  preFac = alpS /( 8.0 * pow(mHat,3));
  return;
}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.

void ResonanceGluino::calcWidth(bool) {


  widNow = 0.0;
  if(ps == 0.) return;
  kinFac = (mHat * mHat - mf1 * mf1 + mf2 * mf2);
  lambda = lam(mHat*mHat, mf1*mf1 , mf2*mf2);

  if(id1Abs > 1000000 && (id1Abs % 100) < 7 && id2Abs < 7) {

    int isq = (abs(id1Abs)/1000000 == 2) ? 
      (abs(id1Abs)%10+1)/2 + 3: (abs(id1Abs)%10+1)/2;
    bool idown = (id2Abs%2 == 1);
    int iq = (id2Abs + 1)/2;

    // ~g -> ~q + q
    if(idown){
      widNow = kinFac * (norm(coupSUSYPtr->LsddG[isq][iq])+norm(coupSUSYPtr->RsddG[isq][iq]))
	+ 4.0 * mHat * mf2 * real( coupSUSYPtr->LsddG[isq][iq] 
						  * conj(coupSUSYPtr->RsddG[isq][iq]));
    }
    else{
      widNow = kinFac * (norm(coupSUSYPtr->LsuuG[isq][iq])+norm(coupSUSYPtr->RsuuG[isq][iq]))
	+ 4.0 * mHat * mf2 * real( coupSUSYPtr->LsuuG[isq][iq] 
						  * conj(coupSUSYPtr->RsuuG[isq][iq]));

    }
    widNow = widNow * preFac * ps * sqrt(lambda);
    if(DEBUG) {
      cout<<"Gluino:: id1:"<<id1Abs<<" id2:"<<id2Abs<<" Width: ";
      cout<<scientific<<widNow<<endl;
    }
    return;
  }
}

//==========================================================================

//  Class ResonanceNeut
//  Derived class for Neutralino Resonances

//--------------------------------------------------------------------------


void ResonanceNeut::initConstants(){
  
  alpEM = coupSUSYPtr->alphaEM(mHat * mHat);
  s2W   = coupSUSYPtr->sin2W;
}
 
//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.
void  ResonanceNeut::calcPreFac(bool){

  // Common coupling factors.
  preFac = alpEM / (8.0 * s2W * pow(mHat,3));
  return;
}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.
void  ResonanceNeut::calcWidth(bool){

  widNow = 0.0;

  if(ps ==0.) return;
  kinFac = mHat * mHat - mf1 * mf1 + mf2 * mf2;
  kinFac2 = pow(mHat,4) + pow(mf1,4) - 2.0 * pow(mf2,4) 
    + pow2(mHat) * pow2(mf2) + pow2(mf1) * pow2(mf2) - pow2(mHat) * pow2(mf1);
  lambda = lam(mHat*mHat, mf1*mf1, mf2*mf2);

  // Stable lightest neutralino
  if(idRes == 1000022) return;

  double fac = 0.0;
  int iNeut1 = typeNeut(idRes);
  int iNeut2 = typeNeut(id1Abs);
  int iChar1 = typeChar(id1Abs);

  if(iNeut2>0 && id2Abs == 23){
    // ~chi0_i -> chi0_j + Z
    fac = kinFac2 * (norm(coupSUSYPtr->OLpp[iNeut1][iNeut2]) + norm(coupSUSYPtr->ORpp[iNeut1][iNeut2]));
    fac -= 12.0 * mHat * mf1 * pow2(mf2) *
      real(coupSUSYPtr->OLpp[iNeut1][iNeut2] * conj(coupSUSYPtr->ORpp[iNeut1][iNeut2]));
    fac /= pow2(mf2) * (1.0 - s2W);
  }
  else if(iChar1>0 && id2Abs==24){
    // ~chi0_i -> chi+_j + W- (or c.c.)

    fac = kinFac2 * (norm(coupSUSYPtr->OL[iNeut1][iChar1]) + norm(coupSUSYPtr->OR[iNeut1][iChar1]));
    fac -= 12.0 * mHat * mf1 * pow2(mf2) * 
      real(coupSUSYPtr->OL[iNeut1][iChar1] * conj(coupSUSYPtr->OR[iNeut1][iChar1]));
    fac /= pow2(mf2);
  }
  else if(id1Abs > 1000000 && id1Abs%100 < 7 && id2Abs < 7){
    // ~chi0_k -> ~q + q
    bool idown = (id1Abs%2 == 1);
    int iq = (id2Abs + 1 )/ 2;
    int isq = (abs(idRes)/1000000 == 2) ? 
      (abs(idRes)%10+1)/2 + 3: (abs(idRes)%10+1)/2;

    if(idown){
      fac = kinFac * (norm(coupSUSYPtr->LsddX[isq][iq][iNeut1]) + norm(coupSUSYPtr->RsddX[isq][iq][iNeut1]));
      fac += 4.0 * mHat * mf2 * 
	real(coupSUSYPtr->LsddX[isq][iq][iNeut1] * conj(coupSUSYPtr->RsddX[isq][iq][iNeut1]));
    }
    else{
      fac = kinFac * (norm(coupSUSYPtr->LsuuX[isq][iq][iNeut1]) + norm(coupSUSYPtr->RsuuX[isq][iq][iNeut1]));
      fac += 4.0 * mHat * mf2 * sqrt(lambda) *
	real(coupSUSYPtr->LsuuX[isq][iq][iNeut1] * conj(coupSUSYPtr->RsuuX[isq][iq][iNeut1]));
    }
    fac *= 2.0/(1 - s2W);
  }

  widNow = fac * preFac * ps * sqrt(lambda);
  if(DEBUG) {
    cout<<idRes<<":: id1:"<<id1Abs<<" id2:"<<id2Abs<<" Width: ";
    cout<<scientific<<widNow<<endl;
  }
  return;
}

//==========================================================================

//  Class ResonanceChar
//  Derived class for Neutralino Resonances

//--------------------------------------------------------------------------


void ResonanceChar::initConstants(){

  alpEM = coupSUSYPtr->alphaEM(mHat * mHat);
  s2W   = coupSUSYPtr->sin2W;
  return;
}
 
//--------------------------------------------------------------------------

// Calculate various common prefactors for the current mass.
void  ResonanceChar::calcPreFac(bool){

  preFac = alpEM / (8.0 * s2W * pow(mHat,3));
  return;
}

//--------------------------------------------------------------------------

// Calculate width for currently considered channel.
void  ResonanceChar::calcWidth(bool){

  widNow = 0.0;
  if(ps == 0.) return;

  double fac = 0.0;
  kinFac = mHat * mHat - mf1 * mf1 + mf2 * mf2;
  kinFac2 = pow(mHat,4) + pow(mf1,4) - 2.0 * pow(mf2,4) 
    + pow2(mHat) * pow2(mf2) + pow2(mf1) * pow2(mf2) - pow2(mHat) * pow2(mf1);
  lambda = lam(mHat*mHat , mf1*mf1 , mf2*mf2);

  int idChar1 = typeChar(idRes);
  int idChar2 = typeChar(id1Abs);
  int idNeut1 = typeNeut(id1Abs);

  if(idChar2>0 && id2Abs == 23){
    // ~chi_i -> chi_j + Z
    fac = kinFac2 * (norm(coupSUSYPtr->OLp[idChar1][idChar2]) 
		     + norm(coupSUSYPtr->ORp[idChar1][idChar2]));
    fac -= 12.0 * mHat * mf1 * pow2(mf2) *
      real(coupSUSYPtr->OLp[idChar1][idChar2] 
	   * conj(coupSUSYPtr->ORp[idChar1][idChar2]));
    fac /= pow2(mf2) * (1.0 - s2W);
  }
  else if(idNeut1>0 && id2Abs==24){
    // ~chi_i -> chi0_j + W- (or c.c.)

    fac = kinFac2 * (norm(coupSUSYPtr->OL[idNeut1][idChar1]) + norm(coupSUSYPtr->OR[idNeut1][idChar1]));
    fac -= 12.0 * mHat * mf1 * pow2(mf2) *
      real(coupSUSYPtr->OL[idNeut1][idChar1] * conj(coupSUSYPtr->OR[idNeut1][idChar1]));
    fac /= pow2(mf2);
  }
  else if(id1Abs > 1000000 && id1Abs%100 < 7 && id2Abs < 7){
    // ~chi0_k -> ~q + q
    bool idown = (id1Abs%2 == 1);
    int iq = (id2Abs + 1 )/ 2;
    int isq = (abs(idRes)/1000000 == 2) ? 
      (abs(idRes)%10+1)/2 + 3: (abs(idRes)%10+1)/2;

    if(idown){
      fac = kinFac * (norm(coupSUSYPtr->LsduX[isq][iq][idChar1]) + norm(coupSUSYPtr->RsduX[isq][iq][idChar1]));
      fac += 4.0 * mHat * mf2 * 
	real(coupSUSYPtr->LsduX[isq][iq][idChar1] * conj(coupSUSYPtr->RsduX[isq][iq][idChar1]));
    }
    else{
      fac = kinFac * (norm(coupSUSYPtr->LsudX[isq][iq][idChar1]) + norm(coupSUSYPtr->RsudX[isq][iq][idChar1]));
      fac += 4.0 * mHat * mf2 * 
	real(coupSUSYPtr->LsudX[isq][iq][idChar1] * conj(coupSUSYPtr->RsudX[isq][iq][idChar1]));
    }
    fac *= 2.0/(1 - s2W);
  }

  widNow = fac * preFac * ps * sqrt(lambda) ;
  if(DEBUG) {
    cout<<idRes<<":: id1:"<<id1Abs<<" id2:"<<id2Abs<<" Width: ";
    cout<<scientific<<widNow<<endl;
  }
  return;
}


//==========================================================================

//Return neutralino code; zero if not a neutralino

int SUSYResonanceWidths::typeNeut(int idPDG) {
  int type = 0;
  int idAbs = abs(idPDG);
  if(idAbs == 1000022) type = 1;
  else if(idAbs == 1000023) type = 2;
  else if(idAbs == 1000025) type = 3;
  else if(idAbs == 1000035) type = 4;
  else if(coupSUSYPtr->isNMSSM && idAbs == 1000045) type = 5;
  return type;

}  


//--------------------------------------------------------------------------

//Check whether particle is a Chargino

int SUSYResonanceWidths::typeChar(int idPDG) {
  int type = 0;
  if(abs(idPDG) == 1000024) type = 1;
  else if (abs(idPDG) == 1000037)type = 2;
  return type;
}  

//--------------------------------------------------------------------------

// Function for Kallen function

double SUSYResonanceWidths::lam(double x, double y, double z){
  
  double val = x*x + y*y + z*z - 2.0* (x*y + y*z + z*x);
  return val;

}


//--------------------------------------------------------------------------



} //end namespace Pythia8


