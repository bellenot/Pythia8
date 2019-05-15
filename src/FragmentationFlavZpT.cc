// Function definitions (not found in the header) for the 
// StringFlav, StringZ and StringPT classes.
// Copyright C 2006 Torbjorn Sjostrand

#include "FragmentationFlavZpT.h"

namespace Pythia8 {

//**************************************************************************

// The StringFlav class.

//*********
 
// Definitions of static variables.
// (Values will be overwritten in initStatic call, so are purely dummy.)

double StringFlav::probQQtoQ = 0.1;
double StringFlav::probStoU = 0.3;
double StringFlav::probSQtoQQ = 0.4;
double StringFlav::probQQ1toQQ0 = 0.05;
double StringFlav::probQandQQ = 1.1;
double StringFlav::probQandS = 2.3;
double StringFlav::probQandSinQQ = 2.12;
double StringFlav::probQQ1corr = 0.15;
double StringFlav::probQQ1corrInv = 6.667;
double StringFlav::probQQ1norm = 0.13;
double StringFlav::mesonUspin1 = 0.5;
double StringFlav::mesonSspin1 = 0.6;
double StringFlav::mesonCspin1 = 0.7;
double StringFlav::mesonBspin1 = 0.75;
double StringFlav::mesonMix1[2][4] 
  = {{ 0., 0.5, 0.5, 0.},  {0., 0.5, 0.5, 0.} };
double StringFlav::mesonMix2[2][4] 
  = { {0., 0.75, 0.75, 0.5}, {0., 1., 1., 0.} };
double StringFlav::suppressEta = 1.0;
double StringFlav::suppressEtaPrime = 0.4;

// Clebsch-Gordan coefficients for baryon octet and decuplet are
// fixed once and for all, so only weighted sum needs to be edited.
double StringFlav::baryonClebsch12[6] 
  = { 0.75, 0.5, 0., 0.1667, 0.0833, 0.1667};
double StringFlav::baryonClebsch32[6] 
  = { 0.,  0.,  1., 0.3333, 0.6667, 0.3333};
double StringFlav::baryonClebschSum[6] 
  = { 0.75, 0.5, 1., 0.5, 0.75, 0.5};

//*********

// Initialize static data members of the flavour generation.

void StringFlav::initStatic() {

  // Basic parameters for generation of new flavour.
  probQQtoQ = Settings::parm("StringFlav:probQQtoQ");
  probStoU = Settings::parm("StringFlav:probStoU");
  probSQtoQQ = Settings::parm("StringFlav:probSQtoQQ");
  probQQ1toQQ0 = Settings::parm("StringFlav:probQQ1toQQ0");

  // Parameters derived from above.
  probQandQQ = 1. + probQQtoQ;
  probQandS = 2. + probStoU;
  probQandSinQQ = 2. + probSQtoQQ * probStoU;
  probQQ1corr = 3. * probQQ1toQQ0;
  probQQ1corrInv = 1. / probQQ1corr;
  probQQ1norm = probQQ1corr / (1. + probQQ1corr);

  // Parameters for meson production.
  mesonUspin1 = Settings::parm("StringFlav:mesonUspin1");
  mesonSspin1 = Settings::parm("StringFlav:mesonSspin1");
  mesonCspin1 = Settings::parm("StringFlav:mesonCspin1");
  mesonBspin1 = Settings::parm("StringFlav:mesonBspin1");
  suppressEta = Settings::parm("StringFlav:suppressEta");
  suppressEtaPrime = Settings::parm("StringFlav:suppressEtaPrime");

  // Set parameters for uubar - ddbar - ssbar meson mixing.
  for (int spin = 0; spin < 2; ++spin) { 
    // Mixing angle for pseudoscalar and vector multiplets.
    double alpha;
    if (spin == 0) {
      double thetaPS = Settings::parm("StringFlav:thetaPS");
      alpha = 90. - (thetaPS + 54.7);
    } else { 
      double thetaV = Settings::parm("StringFlav:thetaV");
      alpha = thetaV + 54.7;
    } 
    alpha *= M_PI / 180.;
    // Fill in (spin, flavour)-dependent probability of producing
    // the lightest or the lightest two mesons of the nonet. 
    mesonMix1[spin][1] = 0.5;
    mesonMix2[spin][1] = 0.5 * (1. + pow2(sin(alpha)));
    mesonMix1[spin][2] = mesonMix1[spin][1];
    mesonMix2[spin][2] = mesonMix2[spin][1];
    mesonMix1[spin][3] = 0.;
    mesonMix2[spin][3] = pow2(cos(alpha));
  }

  // Sum of baryon octet and decuplet weights.
  double suppressDecuplet = Settings::parm("StringFlav:suppressDecuplet");
  for (int i = 0; i < 6; ++i) baryonClebschSum[i]
    = baryonClebsch12[i] + suppressDecuplet * baryonClebsch32[i];

}

//*********

// Pick a new flavour (including diquarks) given an incoming one.

int StringFlav::pick(int idOld) {

  // Choose whether to generate a meson or a baryon.
  int idOldAbs = abs(idOld);
  bool doBaryon = (idOldAbs <9 && probQandQQ * Rndm::flat() < probQQtoQ) 
    ? true : false; 

  // Single quark for new meson or baryon.
  if (!doBaryon) {
    double rndmFlav = probQandS * Rndm::flat();
    int idNewAbs = 1;
    if (rndmFlav > 1.) idNewAbs = 2;
    if (rndmFlav > 2.) idNewAbs = 3;
    return ( (idOld > 0 && idOld < 9) || idOld < -8 ) ? -idNewAbs : idNewAbs; 
  }

  // Constituent flavours in diquark.
  int idNew1, idNew2, spin;
  double weight;
  do {
    double rndmFlav = probQandSinQQ * Rndm::flat();
    idNew1 = 1;
    if (rndmFlav > 1.) idNew1 = 2;
    if (rndmFlav > 2.) idNew1 = 3;
    rndmFlav = probQandSinQQ * Rndm::flat();
    idNew2 = 1;
    if (rndmFlav > 1.) idNew2 = 2;
    if (rndmFlav > 2.) idNew2 = 3;
    // Diquark spin 1 or 0; reweight to get correct mix.
    spin = (idNew1 >= idNew2) ? 1 : 0;
    weight = (spin == 0) ? probQQ1corrInv : probQQ1corr; 
  } while (weight < Rndm::flat());    
  // Combined diquark code.
  int idNewAbs = 1000 * max(idNew1, idNew2) + 100 * min(idNew1, idNew2) 
    + 2 * spin + 1; 
  return (idOld > 0) ? idNewAbs : -idNewAbs; 

}

//*********

// Combine two flavours (including diquarks) to produce a hadron.
// The weighting of the combination may fail, giving output 0.

int StringFlav::combine(int id1, int id2) {

  // Recognize largest and smallest flavour.
  int id1Abs = abs(id1);
  int id2Abs = abs(id2);
  int idMax = max(id1Abs, id2Abs);
  int idMin = min(id1Abs, id2Abs);
 
  // Form meson: pick spin state and preliminary code.
  if (idMax < 9) {
    int spin;
    if (idMax < 3) spin = (Rndm::flat() < mesonUspin1) ? 1 : 0;
    else if (idMax == 3) spin = (Rndm::flat() < mesonSspin1) ? 1 : 0;
    else if (idMax == 4) spin = (Rndm::flat() < mesonCspin1) ? 1 : 0;
    else spin = (Rndm::flat() < mesonBspin1) ? 1 : 0;
    int idMeson = 100 * idMax + 10 * idMin + 2 * spin + 1;

    // For nondiagonal mesons distinguish particle/antiparticle.
    if (id1Abs != id2Abs) {
      int sign = (idMax%2 == 0) ? 1 : -1;
      if ( (idMax == id1Abs && id1 < 0) 
        || (idMax == id2Abs && id2 < 0) ) sign = -sign;
      idMeson *= sign;  

    // For light diagonal mesons include uubar - ddbar - ssbar mixing.
    } else if (idMax < 4) {
      double rMix = Rndm::flat();
      if (rMix < mesonMix1[spin][idMax]) idMeson = 110 + 2 * spin + 1;
      else if (rMix < mesonMix2[spin][idMax]) idMeson = 220 + 2 * spin + 1;
      else idMeson = 330 + 2 * spin + 1;

      // Additional suppression of eta and eta' may give failure.
      if (idMeson == 221 && suppressEta < Rndm::flat()) return 0;
      if (idMeson == 331 && suppressEtaPrime < Rndm::flat()) return 0;
    }

    // Finished for mesons.
    return idMeson;
  }

  // SU(6) factors for baryon production may give failure.
  int idQQ1 = idMax / 1000;
  int idQQ2 = (idMax / 100) % 10;
  int spinQQmod = idMax % 10;
  int spinFlav = spinQQmod - 1;
  if (spinFlav == 2 && idQQ1 != idQQ2) spinFlav = 4;
  if (idMin != idQQ1 && idMin != idQQ2) spinFlav++;   
  if (baryonClebschSum[spinFlav] < Rndm::flat()) return 0;

  // Order quarks to form baryon. Pick spin.
  int idOrd1 = max( idMin, max( idQQ1, idQQ2) ); 
  int idOrd3 = min( idMin, min( idQQ1, idQQ2) ); 
  int idOrd2 = idMin + idQQ1 + idQQ2 - idOrd1 - idOrd3;
  int spinMod = (baryonClebschSum[spinFlav] *Rndm::flat() 
    < baryonClebsch12[spinFlav]) ? 2 : 4;
  
  // Distinguish Lambda- and Sigma-like. 
  bool LambdaLike = false;
  if (spinMod == 2 && idOrd1 > idOrd2 && idOrd2 > idOrd3) {
    if (spinQQmod == 1 && idOrd1 == idMin) LambdaLike = true;
    else if (spinQQmod == 1) LambdaLike = (Rndm::flat() < 0.25) 
      ? false : true;
    else if (idOrd1 != idMin) LambdaLike = (Rndm::flat() < 0.75)  
      ? false : true;    
  }

  // Form baryon code and return with sign.  
  int idBaryon = (LambdaLike) 
    ? 1000 * idOrd1 + 100 * idOrd3 + 10 * idOrd2 + spinMod
    : 1000 * idOrd1 + 100 * idOrd2 + 10 * idOrd3 + spinMod;
   return (id1 > 0) ? idBaryon : -idBaryon;

}

//*********

// Combine two quarks to produce a diquark. 
// Normally according to production compostion, but nonvanishing idHad
// means diquark from known hadron content, so use SU(6) wave fucntion.

int StringFlav::makeDiquark(int id1, int id2, int idHad) {

  // Initial values.
  int idMin = min( abs(id1), abs(id2));
  int idMax = max( abs(id1), abs(id2));
  int spin = 1;

  // Select spin of diquark formed from two valence quarks in proton.
  // (More hadron cases??)
  if (abs(idHad) == 2212) {
    if (idMin == 1 && idMax == 2 && Rndm::flat() < 0.75) spin = 0;  

  // Else select spin of diquark according to production composition.
  } else {
    if (idMin != idMax && Rndm::flat() > probQQ1norm) spin = 0; 
  }   

  // Combined diquark code.
  int idNewAbs = 1000 * idMax + 100 * idMin + 2 * spin + 1; 
  return (id1 > 0) ? idNewAbs : -idNewAbs; 

}
 
//**************************************************************************

// The StringZ class.

//*********
 
// Definitions of static variables.
// (Values will be overwritten in initStatic call, so are purely dummy.)
 
bool StringZ::usePetersonC = false;
bool StringZ::usePetersonB = false;
bool StringZ::usePetersonH = false;
double StringZ::mc2 = 2.25;
double StringZ::mb2 = 25.0;
double StringZ::aLund = 0.3;
double StringZ::bLund = 0.58;
double StringZ::aExtraDiquark = 0.5;
double StringZ::rFactC = 1.0;
double StringZ::rFactB = 1.0;
double StringZ::rFactH = 1.0;
double StringZ::epsilonC = 0.05;
double StringZ::epsilonB = 0.005;
double StringZ::epsilonH = 0.005;

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// When a or c are close to special cases, default to these.
const double StringZ::CFROMUNITY = 0.01;
const double StringZ::AFROMZERO = 0.02;
const double StringZ::AFROMC = 0.01;

// Do not take exponent of too large or small number.
const double StringZ::EXPMAX = 50.; 

//*********

// Initialize static data members of the string z selection.

void StringZ::initStatic() {

  // c and b quark masses.
  mc2 = pow2( ParticleDataTable::m0(4)); 
  mb2 = pow2( ParticleDataTable::m0(5)); 

  // Paramaters of Lund/Bowler symmetric fragmentation function.
  aLund = Settings::parm("StringZ:aLund");
  bLund = Settings::parm("StringZ:bLund");
  aExtraDiquark = Settings::parm("StringZ:aExtraDiquark");
  rFactC = Settings::parm("StringZ:bExtraC");
  rFactB = Settings::parm("StringZ:bExtraB");
  rFactH = Settings::parm("StringZ:bExtraH");

  // Flags and parameters of Peterson/SLAC fragmentation function.
  usePetersonC = Settings::flag("StringZ:usePetersonC");
  usePetersonB = Settings::flag("StringZ:usePetersonB");
  usePetersonH = Settings::flag("StringZ:usePetersonH");
  epsilonC = Settings::parm("StringZ:epsilonC");
  epsilonB = Settings::parm("StringZ:epsilonB");
  epsilonH = Settings::parm("StringZ:epsilonH");

}

//*********

// Generate the fraction z that the next hadron will take, 
// using either Lund/Bowler or, for heavy, Peterson/SLAC functions.
// Note: for a heavy new coloured particle we assume pT negligible.

double StringZ::zFrag( int idOld, int idNew, double mT2) {

  // Find if old or new flavours correspond to diquarks.
  int idOldAbs = abs(idOld);
  int idNewAbs = abs(idNew);
  bool isOldDiquark = (idOldAbs > 1000 && idOldAbs < 10000) ? true : false;
  bool isNewDiquark = (idNewAbs > 1000 && idNewAbs < 10000) ? true : false;

  // Find heaviest quark in fragmenting parton/diquark.
  int idFrag = idOldAbs;
  if (isOldDiquark) idFrag = max( idOldAbs / 1000, (idOldAbs / 100) % 10);
  
  // Use Peterson where explicitly requested for heavy flavours.
  if (idFrag == 4 && usePetersonC) return zPeterson( epsilonC);
  if (idFrag == 5 && usePetersonB) return zPeterson( epsilonB);
  if (idFrag >  5 && usePetersonH) {
    double epsilon = epsilonH * mb2 / mT2; 
    return zPeterson( epsilon);
  }

  // Shape parameters of Lund symmetric fragmentation function.
  double aShape = aLund;
  if (isOldDiquark) aShape += aExtraDiquark;
  double bShape = bLund * mT2;
  double cShape = 1.;
  if (isOldDiquark) cShape -= aExtraDiquark;
  if (isNewDiquark) cShape += aExtraDiquark;
  if (idFrag == 4) cShape += rFactC * bLund * mc2;
  if (idFrag == 5) cShape += rFactB * bLund * mb2;
  if (idFrag >  5) cShape += rFactH * bLund * mT2;
  return zLund( aShape, bShape, cShape);

}

//*********

// Generate a random z according to the Lund/Bowler symmetric
// fragmentation function f(z) = (1 -z)^a * exp(-b/z) / z^c.
// Normalized so that f(z_max) = 1  it can also be written as
// f(z) = exp( a * ln( (1 - z) / (1 - z_max) ) + b * (1/z_max - 1/z) 
//           + c * ln(z_max/z) ).  

double StringZ::zLund( double a, double b, double c) {

  // Special cases for c = 1, a = 0 and a = c. 
  bool cIsUnity = (abs( c - 1.) < CFROMUNITY) ? true : false;
  bool aIsZero = (a < AFROMZERO) ? true : false;
  bool aIsC = (abs(a - c) < AFROMC) ? true : false;

  // Determine position of maximum.
  double zMax;
  if (aIsZero) zMax = (c > b) ? b / c : 1.; 
  else if (aIsC) zMax = b / (b + c);
  else { zMax = 0.5 * (b + c - sqrt( pow2(b - c) + 4. * a * b)) / (c - a);
         if (zMax > 0.9999 && b > 100.) zMax = min(zMax, 1. - a / b); }   
        
  // Subdivide z range if distribution very peaked near either endpoint.
  bool peakedNearZero = (zMax < 0.1) ? true : false;
  bool peakedNearUnity = (zMax > 0.85 && b > 1.) ? true : false;

  // Find integral of trial function everywhere bigger than f. 
  // (Dummy start values.)
  double fIntLow = 1.; 
  double fIntHigh = 1.; 
  double fInt = 2.; 
  double zDiv = 0.5; 
  double zDivC = 0.5;
  // When z_max is small use that f(z)
  //   < 1     for z < z_div = 2.75 * z_max,
  //   < (z_div/z)^c for z > z_div (=> logarithm for c = 1, else power).   
  if (peakedNearZero) {
    zDiv = 2.75 * zMax;
    fIntLow = zDiv; 
    if (cIsUnity) fIntHigh = -zDiv * log(zDiv);
    else { zDivC = pow( zDiv, 1. - c);
           fIntHigh = zDiv * (1. - 1./zDivC) / (c - 1.);} 
    fInt = fIntLow + fIntHigh;
  // When z_max large use that f(z)
  //   < exp( b * (z - z_div) ) for z < z_div with z_div messy expression,
  //   < 1   for z > z_div.
  // To simplify expressions the integral is extended to z =  -infinity.
  } else if (peakedNearUnity) {
    double rcb = sqrt(4. + pow2(c / b));
    zDiv = rcb - 1./zMax - (c / b) * log( zMax * 0.5 * (rcb + c / b) );  
    if (!aIsZero) zDiv += (a/b) * log(1. - zMax);
    zDiv = min( zMax, max(0., zDiv));
    fIntLow = 1. / b;
    fIntHigh = 1. - zDiv; 
    fInt = fIntLow + fIntHigh;
  }

  // Choice of z, preweighted for peaks at low or high z. (Dummy start values.)
  double z = 0.5;
  double fPrel = 1.; 
  double fVal = 1.; 
  do { 
    // Choice of z flat good enough for distribution peaked in the middle;
    // if not this z can be reused as a random number in general.
    z = Rndm::flat();
    fPrel = 1.;
    // When z_max small use flat below z_div and 1/z^c above z_div.
    if (peakedNearZero) {
      if (fInt * Rndm::flat() < fIntLow) z = zDiv * z;
      else if (cIsUnity) {z = pow( zDiv, z); fPrel = zDiv / z;}
      else { z = pow( zDivC + (1. - zDivC) * z, 1. / (1. - c) );
             fPrel = pow( zDiv / z, c); }
    // When z_max large use exp( b * (z -z_div) ) below z_div and flat above it.
    } else if (peakedNearUnity) {
      if (fInt * Rndm::flat() < fIntLow) { 
        z = zDiv + log(z) / b;
        fPrel = exp( b * (z - zDiv) ); 
      } else z = zDiv + (1. - zDiv) * z; 
    }  

    // Evaluate actual f(z) (if in physical range) and correct.
    if (z > 0 && z < 1) {
      double fExp = b * (1. / zMax - 1. / z)+ c * log(zMax / z);
      if (!aIsZero) fExp += a * log( (1. - z) / (1. - zMax) );
      fVal = exp( max( -EXPMAX, min( EXPMAX, fExp) ) ) ;
    } else fVal = 0.;
  } while (fVal < Rndm::flat() * fPrel);

  // Done.
  return z;

}

//*********

// Generate a random z according to the Peterson/SLAC formula
// f(z) = 1 / ( z * (1 - 1/z - epsilon/(1-z))^2 )
//      = z * (1-z)^2 / ((1-z)^2 + epsilon * z)^2.

double StringZ::zPeterson( double epsilon) {

  double z, fVal;

  // For large epsilon pick z flat and reject, 
  // knowing that 4 * epsilon * f(z) < 1 everywhere.
  if (epsilon > 0.01) { 
    do { 
      z = Rndm::flat();
      fVal = 4. * epsilon * z * pow2(1. - z) 
        / pow2( pow2(1. - z) + epsilon * z);
    } while (fVal < Rndm::flat());
    return z; 
  } 
  
  // Else split range, using that 4 * epsilon * f(z) 
  //   < 4 * epsilon / (1 - z)^2 for 0 < z < 1 - 2 * sqrt(epsilon)
  //   < 1                       for 1 - 2 * sqrt(epsilon) < z < 1
  double epsRoot = sqrt(epsilon);
  double epsComb = 0.5 / epsRoot - 1.;
  double fIntLow = 4. * epsilon * epsComb;
  double fInt = fIntLow + 2. * epsRoot;
  do { 
    if (Rndm::flat() * fInt < fIntLow) {
      z = 1. - 1. / (1. + Rndm::flat() * epsComb);
      fVal = z * pow2( pow2(1. - z) / (pow2(1. - z) + epsilon * z) );
    } else {
      z = 1. - 2. * epsRoot * Rndm::flat();
      fVal = 4. * epsilon * z * pow2(1. - z) 
        / pow2( pow2(1. - z) + epsilon * z);
    }
  } while (fVal < Rndm::flat());
  return z; 

} 
 
//**************************************************************************

// The StringPT class.

//*********
 
// Definitions of static variables.
// (Values will be overwritten in init call, so are purely dummy.)

double StringPT::sigmaQ = 0.25;
double StringPT::enhancedFraction = 0.01;
double StringPT::enhancedWidth = 2.;

//*********

// Initialize static data members of the string pT selection.

void StringPT::initStatic() {

  // Parameters of the pT width and enhancement.
  sigmaQ = Settings::parm("StringPT:sigma") / sqrt(2.);
  enhancedFraction = Settings::parm("StringPT:enhancedFraction");
  enhancedWidth = Settings::parm("StringPT:enhancedWidth");

}

//*********

// Generate Gaussian pT such that <p_x^2> = <p_x^2> = sigma^2 = width^2/2,
// but with small fraction multiplied up to a broader spectrum.

double StringPT::pxy() {

  double pxy = sigmaQ * Rndm::gauss();
  if (Rndm::flat() < enhancedFraction) pxy *= enhancedWidth;
  return pxy;

}
  
//**************************************************************************

} // end namespace Pythia8
