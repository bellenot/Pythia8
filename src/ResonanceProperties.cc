// Function definitions (not found in the header) for 
// the ResonanceProperties class and classes derived from it.
// Copyright C 2006 Torbjorn Sjostrand

#include "ResonanceProperties.h"

namespace Pythia8 {

//**************************************************************************

// The Resonance Properties class.
// Base class for the various resonances.

//*********
 
// Definitions of static variables and functions.
// (Values will be overwritten in initStatic call, so are purely dummy.)

int ResonanceProperties::alphaSorder = 1;
double ResonanceProperties::alphaSvalue = 0.1265;
AlphaStrong ResonanceProperties::alphaScalc;
AlphaEM ResonanceProperties::alphaEMcalc;

// The sum of product masses must not be too close to the resonance mass.
const double ResonanceProperties::MASSMARGIN = 0.01;

//*********

// Initialize static data members.

void ResonanceProperties::initStatic() {

  // Parameters of alphaStrong generation .
  alphaSvalue = Settings::parm("SigmaProcess:alphaSvalue");
  alphaSorder = Settings::mode("SigmaProcess:alphaSorder");

  // Initialize alphaStrong generation.
  alphaScalc.init( alphaSvalue, alphaSorder); 

}
 
//**************************************************************************

// The ResonanceGmZ class.
// Derived class for gamma*/Z0 properties.

//*********
 
// Definitions of static variables and functions.
int ResonanceGmZ::idRes = 23;
ParticleDataEntry* ResonanceGmZ::particlePtr;
int ResonanceGmZ::gmZmode = 0; 
double ResonanceGmZ::mRes, ResonanceGmZ::GammaRes, ResonanceGmZ::m2Res, 
  ResonanceGmZ::GamMRat, ResonanceGmZ::thetaWRat; 

//*********

// Initialize static data members.

void ResonanceGmZ::initStatic() {

  // Set pointer to particle properties and decay table.
  particlePtr = ParticleDataTable::particleDataPtr(idRes);

  // Z0 properties.
  mRes = ParticleDataTable::m0(idRes);
  GammaRes = ParticleDataTable::mWidth(idRes);
  m2Res = mRes*mRes;
  GamMRat = GammaRes / mRes;
  thetaWRat = 1. / (16. * CoupEW::sin2thetaW() * CoupEW::cos2thetaW());

  // Allow to pick only gamma* or Z0 part of full gamma*/Z0 expression.
  gmZmode = Settings::mode("SigmaProcess:gmZmode");

}

//*********

// Calculate the total width and store phase-space-weighted coupling sums.

double ResonanceGmZ::width(double mH) {

  // Common coupling factors.
  double sH = mH*mH;
  double colQ = 3. * (1. + alphaScalc.alphaS(sH) / M_PI);
  double alpEM = alphaEMcalc.alphaEM(sH); 

  // Reset quantities to sum. Declare others.
  for (int i = 0; i < 4; ++i) widSum[i] = 0.; 
  double mf, m2Rat, psvec, psax, betaf, colf, ef2, efvf, vf2af2;

  // Loop over three fermion generations, except top.
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    int idAbs = abs( particlePtr->decay[i].product(0) );
    if ( (idAbs > 0 && idAbs < 6) || ( idAbs > 10 && idAbs < 17)) {
      mf = ParticleDataTable::m0(idAbs);

      // Check that above threshold. Phase space.
      if (mH > 2. * mf + MASSMARGIN) {
        m2Rat = pow2(mf / mH);
        betaf = sqrtpos(1. - 4. * m2Rat); 
        psvec = betaf * (1. + 2. * m2Rat);
        psax = pow3(betaf);

        // Combine phase space with couplings.
        colf = (idAbs < 6) ? colQ : 1.;
        ef2 = CoupEW::ef2(idAbs) * psvec;
        efvf = CoupEW::ef(idAbs) * CoupEW::vf(idAbs) * psvec;
        vf2af2 = CoupEW::vf2(idAbs) * psvec + CoupEW::af2(idAbs) * psax; 

        // Store sum of combinations. For outstate only open channels.
        widSum[0] += colf * vf2af2;
        if (particlePtr->decay[i].onMode() > 0) {
          widSum[1] += colf * ef2;
          widSum[2] += colf * efvf;
          widSum[3] += colf * vf2af2;
	}

      // End loop over fermions.
      }
    }
  }

  // Calculate prefactors for width and for gamma/interference/Z0 parts.
  widNorm = alpEM * mH * thetaWRat / 3.;
  sigNorm = 4. * M_PI * pow2(alpEM) / (3. * sH); 
  gamNorm = sigNorm;
  intNorm = sigNorm * 2. * thetaWRat * sH * (sH - m2Res)
    / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  resNorm = sigNorm * pow2(thetaWRat * sH) 
    / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );

  // Optionally only keep gamma* or Z0 term.
  if (gmZmode == 1) {intNorm = 0.; resNorm = 0.;}
  if (gmZmode == 2) {gamNorm = 0.; intNorm = 0.;}

  // Done.
  return widNorm * widSum[0];
  
}

//*********

// Use information stored by width(...) above to calculate cross section
// for a given incoming flavour to any (allowed) final flavour.
// Does not include colour factor 1/3 for incoming quark pair.

double ResonanceGmZ::sigma(int id) {

  // Index in coupling tables. Only accept standard fermions.
  int idAbs = abs(id);
  if (idAbs == 0 || (idAbs > 5 && idAbs < 11) || idAbs > 16) return 0.;

  // Combine gamma, interference and Z0 parts for the given inflavour.
  double sig =  CoupEW::ef2(idAbs) * gamNorm * widSum[1] 
    + CoupEW::ef(idAbs) * CoupEW::vf(idAbs) * intNorm * widSum[2]
    + (CoupEW::vf2(idAbs) + CoupEW::af2(idAbs)) * resNorm * widSum[3];
  return sig;

}

//*********

// Select decay products for given mass and incoming flavour.
  
DecayChannel& ResonanceGmZ::dynamicDecay( double mH, int idIn) {

  // Couplings for the incoming fermion.
  int idInAbs = abs(idIn);
  double ei2 = CoupEW::ef2(idInAbs);
  double eivi = CoupEW::ef(idInAbs) * CoupEW::vf(idInAbs);
  double vi2ai2 = CoupEW::vf2(idInAbs) + CoupEW::af2(idInAbs); 

  // Calculate prefactors for interference and resonance part.
  double sH = mH*mH;
  gamNorm = 1.;
  intNorm = 2. * thetaWRat * sH * (sH - m2Res)
    / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  resNorm = pow2(thetaWRat * sH) 
    / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );

  // Optionally only keep gamma* or Z0 term.
  if (gmZmode == 1) {intNorm = 0.; resNorm = 0.;}
  if (gmZmode == 2) {gamNorm = 0.; intNorm = 0.;}

  // Common coupling factor for quarks.
  double colQ = 3. * (1. + alphaScalc.alphaS(sH) / M_PI);

  // Declare quantities to use.
  double mf, m2Rat, psvec, psax, betaf, colf, ef2, efvf, vf2af2,
    relativeSigma;

  // Loop over three fermion generations, except top and not open ones.
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    relativeSigma = 0.;
    int idAbs = abs( particlePtr->decay[i].product(0) );
    if ( ( (idAbs > 0 && idAbs < 6) || ( idAbs > 10 && idAbs < 17) ) 
      && particlePtr->decay[i].onMode() > 0) {
      mf = ParticleDataTable::m0(idAbs);
 
      // Check that above threshold. Phase space.
      if (mH > 2. * mf + MASSMARGIN) {
        m2Rat = pow2(mf / mH);
        betaf = sqrtpos(1. - 4. * m2Rat); 
        psvec = betaf * (1. + 2. * m2Rat);
        psax = pow3(betaf);

        // Combine phase space with couplings.
        colf = (idAbs < 6) ? colQ : 1.;
        ef2 = CoupEW::ef2(idAbs) * psvec;
        efvf = CoupEW::ef(idAbs) * CoupEW::vf(idAbs) * psvec;
        vf2af2 = CoupEW::vf2(idAbs) * psvec + CoupEW::af2(idAbs) * psax; 

        // Combine instate, propagator and outstate to relative sigma.
        relativeSigma = colf * (ei2 * gamNorm * ef2 + eivi * intNorm * efvf
          + vi2ai2 * resNorm * vf2af2);

      // End loop over fermions. Store result.
      }
    }
    particlePtr->decay[i].dynamicBR(relativeSigma); 
  }
  
  // Pick one channel according to relative sigmas calculated above.
  return particlePtr->decay.dynamicPick();
}

//*********

// Select decay angle for given mass and incoming/outgoing flavours.
  
double ResonanceGmZ::cosTheta( double mH, int idIn, int idOut) {

  // Couplings for in- and out-flavours.
  int idInAbs = abs(idIn);
  double ei = CoupEW::ef(idInAbs);
  double vi = CoupEW::vf(idInAbs);
  double ai = CoupEW::af(idInAbs);
  int idOutAbs = abs(idOut);
  double ef = CoupEW::ef(idOutAbs);
  double vf = CoupEW::vf(idOutAbs);
  double af = CoupEW::af(idOutAbs);

  // Phase space factors. (One power of beta left out in formulae.)
  double mf = ParticleDataTable::m0(idOutAbs);
  double m2Rat = pow2(mf / mH);
  double betaf = sqrtpos(1. - 4. * m2Rat); 

  // Calculate prefactors for interference and resonance part.
  double sH = mH*mH;
  gamNorm = 1.;
  intNorm = 2. * thetaWRat * sH * (sH - m2Res)
    / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  resNorm = pow2(thetaWRat * sH) 
    / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );

  // Optionally only keep gamma* or Z0 term.
  if (gmZmode == 1) {intNorm = 0.; resNorm = 0.;}
  if (gmZmode == 2) {gamNorm = 0.; intNorm = 0.;}

  // Coefficients of angular expression.
  double coefTran = ei*ei * gamNorm * ef*ef + ei * vi * intNorm * ef * vf
    + (vi*vi + ai*ai) * resNorm * (vf*vf + pow2(betaf) * af*af);
  double coefLong = m2Rat * ( ei*ei * gamNorm * ef*ef 
    + ei * vi * intNorm * ef * vf + (vi*vi + ai*ai) * resNorm * vf*vf );
  double coefAsym = betaf * ( ei * ai * intNorm * ef * af 
    + 4. * vi * ai * resNorm * vf * af );

  // Flip asymmetry for in-fermion + out-antifermion.
  if (idIn * idOut < 0) coefAsym = -coefAsym;

  // Do random choice of decay angle and reweight it.
  double wtMax = 2. * (coefTran + abs(coefAsym));
  double wtNow, cosThe;
  do {  
    cosThe = 2. * Rndm::flat() - 1.;
    wtNow = coefTran * (1. + pow2(cosThe)) + coefLong * (1. - pow2(cosThe))  
      + 2. * coefAsym * cosThe;
  } while (wtNow < wtMax * Rndm::flat());

  // Done.
  return cosThe;  
}
 
//**************************************************************************

// The ResonanceW class.
// Derived class for W+- properties.

//*********
 
// Definitions of static variables and functions.
int ResonanceW::idRes = 24;
ParticleDataEntry* ResonanceW::particlePtr;
double ResonanceW::mRes, ResonanceW::GammaRes, ResonanceW::m2Res, 
  ResonanceW::GamMRat, ResonanceW::thetaWRat;

//*********

// Initialize static data members.

void ResonanceW::initStatic() {

  // Set pointer to particle properties and decay table.
  particlePtr = ParticleDataTable::particleDataPtr(idRes);

  // W properties.
  mRes = ParticleDataTable::m0(idRes);
  GammaRes = ParticleDataTable::mWidth(idRes);
  m2Res = mRes*mRes;
  GamMRat = GammaRes / mRes;
  thetaWRat = 1. / (12. * CoupEW::sin2thetaW());

}

//*********

// Calculate the total width and store phases-space-weighted coupling sums.

double ResonanceW::width(double mH) {

  // Common coupling factors.
  double sH = mH*mH;
  double colQ = 3. * (1. + alphaScalc.alphaS(sH) / M_PI);
  double alpEM = alphaEMcalc.alphaEM(sH); 

  // Reset quantities to sum. Declare others.
  widSum = 0.; 
  double widOpen = 0.;
  int id1Abs, id2Abs;
  double mf1, mf2, m2Rat1, m2Rat2, ps, widNow;

  // Loop over three fermion generations, except top.
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    id1Abs = abs( particlePtr->decay[i].product(0) );
    id2Abs = abs( particlePtr->decay[i].product(1) );
    if ( (id1Abs > 0 && id1Abs < 6 && id2Abs > 0 && id2Abs < 6) 
      || (id1Abs > 10 && id1Abs < 17 && id2Abs > 10 && id2Abs < 17)) {
      mf1 = ParticleDataTable::m0(id1Abs);
      mf2 = ParticleDataTable::m0(id2Abs);

      // Check that above threshold. Phase space.
      if (mH > mf1 + mf2 + MASSMARGIN) {
        m2Rat1 = pow2(mf1 / mH);
        m2Rat2 = pow2(mf2 / mH);
        ps = (1. - 0.5 * (m2Rat1 + m2Rat2) - 0.5 * pow2(m2Rat1 - m2Rat2))
          * sqrtpos( pow2(1. - m2Rat1 - m2Rat2) - 4. * m2Rat1 * m2Rat2 ); 

        // Combine phase space with colour factor annd CKM couplings.
        widNow = ps;
        if (id1Abs < 6) widNow *= colQ * VCKM::V2id(id1Abs, id2Abs);
        widSum += widNow;

        // For cross section need width of open channels only.
        if (particlePtr->decay[i].onMode() > 0) widOpen += widNow;


      // End loop over fermions.
      }
    }
  }

  // Include prefactors for width. 
  widSum *= alpEM * mH * thetaWRat;
  widOpen *= alpEM * mH * thetaWRat;

  // Prepare cross section info for width() routine.
  double sigIn = alpEM * sH * thetaWRat;
  double sigProp = 12. / ( pow2(sH - m2Res) + pow2(sH * GamMRat) );
  double sigOut = sH * widOpen / mRes; 
  sigOpen = (M_PI/sH) * sigIn * sigProp * sigOut;

  // Done.
  return widSum;
  
}

//*********

// Use information stored by width(...) above to calculate cross section
// for a given incoming flavour to any (allowed) final flavour.
// Does not include colour factor 1/3 for incoming quark pair.

double ResonanceW::sigma(int id) {

  // No dependence on incoming flavour, except check it is OK.
  int idAbs = abs(id);
  if (idAbs == 0 || (idAbs > 5 && idAbs < 11) || idAbs > 16) return 0.;  

  // Answer.
  return sigOpen;

}

//*********

// Select decay products for given mass (incoming flavour dummy).
  
DecayChannel& ResonanceW::dynamicDecay( double mH, int idIn) {

  // Dummy to avoid compiler warnings.
  idSave = idIn;

  // Common coupling factor for quarks.
  double sH = mH*mH;
  double colQ = 3. * (1. + alphaScalc.alphaS(sH) / M_PI);

  // Declare quantities to use.
  int id1Abs, id2Abs;
  double mf1, mf2, m2Rat1, m2Rat2, ps, relativeSigma;

  // Loop over three fermion generations, except top and not open ones.
  for (int i = 0; i < particlePtr->decay.size(); ++i) {
    relativeSigma = 0.;
    id1Abs = abs( particlePtr->decay[i].product(0) );
    id2Abs = abs( particlePtr->decay[i].product(1) );
    if ( ( (id1Abs > 0 && id1Abs < 6 && id2Abs > 0 && id2Abs < 6) 
      || (id1Abs > 10 && id1Abs < 17 && id2Abs > 10 && id2Abs < 17) ) 
      && particlePtr->decay[i].onMode() > 0)  {
      mf1 = ParticleDataTable::m0(id1Abs);
      mf2 = ParticleDataTable::m0(id2Abs);

      // Check that above threshold. Phase space.
      if (mH > mf1 + mf2 + MASSMARGIN) {
        m2Rat1 = pow2(mf1 / mH);
        m2Rat2 = pow2(mf2 / mH);
        ps = (1. - 0.5 * (m2Rat1 + m2Rat2) - 0.5 * pow2(m2Rat1 - m2Rat2))
          * sqrtpos( pow2(1. - m2Rat1 - m2Rat2) - 4. * m2Rat1 * m2Rat2 ); 

        // Combine phase space with colour factor and CKM couplings.
        if (id1Abs > 10) relativeSigma = ps;
        else relativeSigma = colQ * VCKM::V2id(id1Abs, id2Abs) * ps;

      // End loop over fermions. Store result.
      }
    }
    particlePtr->decay[i].dynamicBR(relativeSigma); 
  }
  
  // Pick one channel according to relative sigmas calculated above.
  return particlePtr->decay.dynamicPick();
}

//*********

// Select decay angle for given mass and incoming/outgoing flavours.
  
double ResonanceW::cosTheta( double mH, int idIn, int idOut1, int idOut2) {

  // Sign of asymmetry.
  double asymSign = (idIn * idOut1 > 0) ? 1. : -1.;

  // Phase space factors. 
  double mf1 = ParticleDataTable::m0(abs(idOut1));
  double mf2 = ParticleDataTable::m0(abs(idOut2));
  double m2Rat1 = pow2(mf1 / mH);
  double m2Rat2 = pow2(mf2 / mH);
  double betaf = sqrtpos( pow2(1. - m2Rat1 - m2Rat2) - 4. * m2Rat1 * m2Rat2 ); 
  double m2RatDiff = pow2(m2Rat1 - m2Rat2);

  // Do random choice of decay angle and reweight it.
  double wtNow, cosThe;
  do {  
    cosThe = 2. * Rndm::flat() - 1.;
    wtNow = pow2(1. + asymSign * betaf * cosThe) - m2RatDiff;
  } while (wtNow < 4. * Rndm::flat());

  // Done.
  return cosThe;  
}

//**************************************************************************

} // end namespace Pythia8
