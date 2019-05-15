// Function definitions (not found in the header) for the UserHooks class.
// Copyright C 2007 Torbjorn Sjostrand

#include "UserHooks.h"

namespace Pythia8 {
 
//**************************************************************************

// The UserHooks class.

//*********

// multiplySigmaBy allows the user to introduce a multiplicative factor 
// that modifies the cross section of a hard process. Since it is called
// from before the event record is generated in full, the normal analysis
// does not work. The code here provides a rather extensive summary of
// which methods actually do work. It is a convenient starting point for 
// writing your own derived routine.

double UserHooks::multiplySigmaBy( const SigmaProcess* sigmaProcessPtr, 
  const PhaseSpace* phaseSpacePtr, bool inEvent) {

  // Process code, necessary when some to be treated differently.
  //int code       = sigmaProcessPtr->code();

  // Final multiplicity, i.e. whether 2 -> 1 or 2 -> 2.
  //int nFinal     = sigmaProcessPtr->nFinal();

  // Incoming x1 and x2 to the hard collision, and factorization scale.
  //double x1      = phaseSpacePtr->x1();
  //double x2      = phaseSpacePtr->x2();
  //double Q2Fac   = sigmaProcessPtr->Q2Fac();

  // Renormalization scale and assumed alpha_strong and alpha_EM.
  //double Q2Ren   = sigmaProcessPtr->Q2Ren();
  //double alphaS  = sigmaProcessPtr->alphaSH();
  //double alphaEM = sigmaProcessPtr->alphaEMH();
  
  // Subprocess mass-square.
  //double sHat = phaseSpacePtr->sHat();

  // Now methods only relevant for 2 -> 2.
  //if (nFinal == 2) {
    
    // Mandelstam variables and hard-process pT.
    //double tHat  = phaseSpacePtr->tHat();
    //double uHat  = phaseSpacePtr->uHat();
    //double pTHat = phaseSpacePtr->pTHat();
  
    // Masses of the final-state particles. (Here 0 for light quarks.)
    //double m3    = sigmaProcessPtr->m(3);
    //double m4    = sigmaProcessPtr->m(4);
  //}

  // Dummy statement to avoid compiler warnings.
  return ((inEvent && sigmaProcessPtr->code() == 0 
    && phaseSpacePtr->sHat() < 0.) ? 0. : 1.);

}

//*********

// subEvent extracts currently resolved partons in the hard process.

void UserHooks::subEvent(const Event& event, bool isHardest) {

  // Reset work event to be empty. 
  workEvent.clear();  

  // Find which subsystem to study.
  int iSys = 0;
  if (!isHardest) iSys = event.sizeSystems() - 1;

  // Loop through all the final partons of the given subsystem.
  for (int i = 2; i < event.sizeSystem(iSys); ++i) {
    int iOld = event.getInSystem( iSys, i);

    // Copy partons to work event.
    int iNew = workEvent.append( event[iOld]); 

    // No mothers. Position in full event as daughters.  
    workEvent[iNew].mothers( 0, 0);
    workEvent[iNew].daughters( iOld, iOld);
  }
 
}
 
//**************************************************************************

// The SuppressSmallPT class, derived from UserHooks.

double SuppressSmallPT::multiplySigmaBy( const SigmaProcess* sigmaProcessPtr, 
  const PhaseSpace* phaseSpacePtr, bool inEvent) {

  // Dummy to avoid compiler warnings.
  inEventNow = inEvent;
        
  // Only modify 2 -> 2 processes.
  int nFinal = sigmaProcessPtr->nFinal();
  if (nFinal != 2) return 1.;

  // Extract pT and modify by pT^4 / (pT^2 + pT0^2)^2 
  double pT2 = pow2( phaseSpacePtr->pTHat() );
  return pow2( pT2 / (pT02 + pT2) ); 

}

 
//**************************************************************************

} // end namespace Pythia8
