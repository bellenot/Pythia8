// Function definitions (not found in the header) for the
// MultipleInteractions class.
// Copyright C 2006 Torbjorn Sjostrand
  
#include "MultipleInteractions.h"

namespace Pythia8 {

//**************************************************************************

// The MultipleInteractions class.

//*********
 
// Definitions of static variables.
// (Values will be overwritten in initStatic call, so are purely dummy.)

double MultipleInteractions::alphaSvalue = 0.127;
int MultipleInteractions::alphaSorder = 1; 
double MultipleInteractions::Kfactor = 1.0; 
double MultipleInteractions::pT0Ref= 3.0; 
double MultipleInteractions::ecmRef = 1800.; 
double MultipleInteractions::ecmPow = 0.16; 
double MultipleInteractions::pTmin = 0.2; 
int MultipleInteractions::bProfile = 2; 
double MultipleInteractions::coreRadius = 0.4; 
double MultipleInteractions::coreFraction = 0.5; 
double MultipleInteractions::expPow = 1.; 
int MultipleInteractions::nQuark = 5; 
int MultipleInteractions::nSample = 1000; 

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Naive upper estimate of cross section too pessimistic, so reduce by this.
const double MultipleInteractions::SIGMAFUDGE = 0.7; 

// Number of bins that the dpT2 / (pT2 + r * pT20)^2 range is split into.
const int MultipleInteractions::NBINS = 100;

// The r value above, picked to allow a flatter correct/trial cross section.
const double MultipleInteractions::RPT20 = 0.25;

// Reduce pT0 by factor pT0STEP if sigmaInt < SIGMASTEP * sigmaND.
const double MultipleInteractions::PT0STEP = 0.9;
const double MultipleInteractions::SIGMASTEP = 1.1;

// Refuse too low expPow in impact parameter profile.
const double MultipleInteractions::EXPPOWMIN = 0.4; 

// Basic step size for b integration; sometimes modified.
const double MultipleInteractions::BSTEP = 0.02;

// Stop b integration when integrand dropped enough.
const double MultipleInteractions::BMAX = 1e-6;

// Do not allow too large argument to exp function.
const double MultipleInteractions::EXPMAX = 50.;

// Convergence criterion for k iteration.
const double MultipleInteractions::KCONVERGE = 1e-5;

// Conversion of GeV^{-2} to mb for cross section.
const double MultipleInteractions::CONVERT2MB = 0.389380; 

//*********

// Initialize static data members.

void MultipleInteractions::initStatic() {

  //  Parameters of alphaStrong and cross section generation.
  alphaSvalue = Settings::parameter("MultipleInteractions:alphaSvalue");
  alphaSorder = Settings::mode("MultipleInteractions:alphaSorder");
  Kfactor = Settings::parameter("MultipleInteractions:Kfactor");

  // Regularization of QCD evolution for pT -> 0. 
  pT0Ref = Settings::parameter("MultipleInteractions:pT0Ref");
  ecmRef = Settings::parameter("MultipleInteractions:ecmRef");
  ecmPow = Settings::parameter("MultipleInteractions:ecmPow");
  pTmin = Settings::parameter("MultipleInteractions:pTmin");

  // Impact parameter profile.
  bProfile = Settings::mode("MultipleInteractions:bProfile");
  coreRadius = Settings::parameter("MultipleInteractions:coreRadius");
  coreFraction = Settings::parameter("MultipleInteractions:coreFraction");
  expPow = max( EXPPOWMIN, 
    Settings::parameter("MultipleInteractions:expPow"));

  // Various other parameters. 
  nQuark = Settings::mode("MultipleInteractions:nQuark");
  nSample = Settings::mode("MultipleInteractions:nSample");

}

//*********

// Initialize the generation process for given beams.

bool MultipleInteractions::init( BeamParticle& beamA, BeamParticle& beamB) {

  // Initialize alpha_strong generation.
  alphaS.init( alphaSvalue, alphaSorder); 

  // Initialize matrix-element calculation.
  dSigmaDt1.init();
  dSigmaDt2.init();

  // Calculate invariant mass of system. Set current pT0 scale.
  sCM = m2( beamA.p(), beamB.p());
  eCM = sqrt(sCM);
  pT0 = pT0Ref * pow(eCM / ecmRef, ecmPow);

  // Get the total inelastic and nondiffractive cross section. Output.
  bool canDoMI = sigmaTot.init( beamA.id(), beamB.id(), eCM);
  if (!canDoMI) return false;
  sigmaND = sigmaTot.sigmaND();
  cout << "\n *-------  PYTHIA Multiple Interactions Initialization  --"
       << "-----* \n"
       << " |                                                        "
       << "     | \n"
       << " |                 sigmaNonDiffractive = " << fixed 
       << setprecision(2) << setw(7) << sigmaND << " mb            | \n"
       << " |                                                        "
       << "     | \n";

  // The pT0 value may need to be decreased, if sigmaInt < sigmaND.
  for ( ; ; ) { 

    // Derived pT kinematics combinations.
    pT20 = pT0*pT0;
    pT2min = pTmin*pTmin;
    pTmax = 0.5*eCM;
    pT2max = pTmax*pTmax;
    pT20R = RPT20 * pT20;
    pT20minR = pT2min + pT20R;
    pT20maxR = pT2max + pT20R;
    pT20min0maxR = pT20minR * pT20maxR;
    pT2maxmin = pT2max - pT2min;   

    // Provide upper estimate of interaction rate d(Prob)/d(pT2).
    upperEnvelope( beamA, beamB);

    // Integrate the parton-parton interaction cross section.
    jetCrossSection( beamA, beamB);

    // Sufficiently big SigmaInt or reduce pT0. Output.
    if (sigmaInt > SIGMASTEP * sigmaND) break; 
    cout << " |  pT0 = "  << setw(5) << pT0 << " gives sigmaInteraction = " 
         << setw(7) << sigmaInt << " mb: rejected  | \n";
    pT0 *= PT0STEP;
  }
  cout << " |  pT0 = " << setw(5) << pT0 << " gives sigmaInteraction = " 
       << setw(7) << sigmaInt << " mb: accepted  | \n"
       << " |                                                        "
       << "     | \n"
       << " *-------  End PYTHIA Multiple Interactions Initialization"
       << "  ---* " << endl;

  // Calculate factor relating matter overlap and interaction rate.
  overlapFactor();

  // Done.
  return true;
}

//*********

// Determine constant in d(Prob)/d(pT2) < const / (pT2 + r * pT20)^2.  

void MultipleInteractions::upperEnvelope( BeamParticle& beamA, 
  BeamParticle& beamB) {  

  // Initially determine constant in jet cross section upper estimate 
  // d(sigma_approx)/d(pT2) < const / (pT2 + r * pT20)^2. 
  pT4dSigmaMax = 0.;
  
  // Loop thorough allowed pT range logarithmically evenly.
  for (int iPT = 0; iPT < NBINS; ++iPT) {
    double pT = pTmin * pow( pTmax/pTmin, (iPT + 0.5) / NBINS);
    pT2 = pT*pT;
    pT2shift = pT2 + pT20;
    xT = 2. * pT / eCM;

    // Evaluate parton density sums at x1 = x2 = xT.
    double xPDF1sumMax = (9./4.) * beamA.xf(21, xT, pT2shift);
    for (int id = 1; id <= nQuark; ++id) 
      xPDF1sumMax += beamA.xf(id, xT, pT2shift) 
        + beamA.xf(-id, xT, pT2shift);
    double xPDF2sumMax = (9./4.) * beamB.xf(21, xT, pT2shift);
    for (int id = 1; id <= nQuark; ++id)
      xPDF2sumMax += beamB.xf(id, xT, pT2shift) 
        + beamB.xf(-id, xT, pT2shift);

    // Evaluate alpha_strong, matrix element and phase space volume.
    double alpS = alphaS.alphaS(pT2shift);
    double dSigmaPartonApprox = CONVERT2MB * Kfactor * 0.5 * M_PI 
      * pow2(alpS / pT2shift);
    double yMax = log(1./xT + sqrt(1./(xT*xT) - 1.));
    double volumePhSp = pow2(2. * yMax);
  
    // Final comparison to determine upper estimate.
    double dSigmaApproxNow = SIGMAFUDGE * xPDF1sumMax * xPDF2sumMax 
      * dSigmaPartonApprox * volumePhSp;
    double pT4dSigmaNow = pow2(pT2 + pT20R) * dSigmaApproxNow;
    if ( pT4dSigmaNow > pT4dSigmaMax) pT4dSigmaMax = pT4dSigmaNow;
  } 
  
  // Get wanted constant by dividing by the nondiffractive cross section.   
  pT4dProbMax = pT4dSigmaMax / sigmaND;

}

//*********

// Integrate the parton-parton interaction cross section,
// using stratified Monte Carlo sampling.
// Store result in pT bins for use as Sudakov form factors.

void MultipleInteractions::jetCrossSection( BeamParticle& beamA, 
  BeamParticle& beamB) {

  // Common factor from bin size in dpT2 / (pT2 + r * pT20)^2 and statistics.   
  double sigmaFactor = (1. / pT20minR - 1. / pT20maxR) / (NBINS * nSample);
  
  // Loop through allowed pT range evenly in dpT2 / (pT2 + r * pT20)^2.
  sigmaInt = 0.; 
  double dSigmaMax = 0.;
  sudExpPT[NBINS] = 0.;
  for (int iPT = NBINS - 1; iPT >= 0; --iPT) {
    double sigmaSum = 0.;

    // In each pT bin sample a number of random pT values.
    for (int iSample = 0; iSample < nSample; ++iSample) {
      double mappedPT2 = 1. - (iPT + Rndm::flat()) / NBINS;
      pT2 = pT20min0maxR / (pT20minR + mappedPT2 * pT2maxmin) - pT20R;

      // Evaluate cross section dSigma/dpT2 in phase space point.
      double dSigma = sigmaPT2( beamA, beamB);
 
      // Multiply by (pT2 + r * pT20)^2 to compensate for pT sampling. Sum.
      dSigma *=  pow2(pT2 + pT20R);
      sigmaSum += dSigma; 
      if (dSigma > dSigmaMax) dSigmaMax = dSigma;      
    }

    // Store total cross section and exponent of Sudakov.
    sigmaSum *= sigmaFactor;
    sigmaInt += sigmaSum;
    sudExpPT[iPT] =  sudExpPT[iPT + 1] + sigmaSum / sigmaND;

  // End of loop over pT values.
  } 

  // Update upper estimate of differential cross section. Done.
  if (dSigmaMax  > pT4dSigmaMax) {
    pT4dSigmaMax = dSigmaMax;
    pT4dProbMax = dSigmaMax / sigmaND;
  }

}  

//*********

// Calculate the actual cross section for initialization decision.
// Lightly modified copy of acceptPT2; e.g. no previous interactions.

double MultipleInteractions::sigmaPT2( BeamParticle& beamA, 
  BeamParticle& beamB) {

  // Derive shifted pT2 and rapidity limits from chosen pT2.
  pT2shift = pT2 + pT20;
  xT = 2. * sqrt(pT2) / eCM;
  xT2 = xT*xT;   
  double yMax = log(1./xT + sqrt(1./xT2 - 1.));

  // Select rapidities y3 and y4 of the two produced partons.
  double y3 = yMax * (2. * Rndm::flat() - 1.);
  double y4 = yMax * (2. * Rndm::flat() - 1.);
  y = 0.5 * (y3 + y4);

  // Reject some events at large rapidities to improve efficiency.
  // (Don't have to evaluate PDF's and ME's.)
  double WTy = (1. - pow2(y3/yMax)) * (1. - pow2(y4/yMax));
  if (WTy < Rndm::flat()) return 0.; 

  // Failure if x1 or x2 exceed unity.
  x1 = 0.5 * xT * (exp(y3) + exp(y4));
  x2 = 0.5 * xT * (exp(-y3) + exp(-y4));
  if (x1 > 1. || x2 > 1.) return 0.; 
  tau = x1 * x2;

  // Evaluate parton densities at actual x1 and x2.
  double xPDF1[21];
  double xPDF1sum = 0.;
  for (int id = -nQuark; id <= nQuark; ++id) {
    if (id == 0) xPDF1[10] = (9./4.) * beamA.xf(21, x1, pT2shift);
    else xPDF1[id+10] = beamA.xf(id, x1, pT2shift);
    xPDF1sum += xPDF1[id+10];
  }
  double xPDF2[21];
  double xPDF2sum = 0.;
  for (int id = -nQuark; id <= nQuark; ++id) {
    if (id == 0) xPDF2[10] = (9./4.) * beamB.xf(21, x2, pT2shift);
    else xPDF2[id+10] = beamB.xf(id, x2, pT2shift);
    xPDF2sum += xPDF2[id+10];
  }

  // Select incoming flavours according to actual PDF's.
  id1 = -nQuark-1;
  double temp = xPDF1sum * Rndm::flat();
  do { temp -= xPDF1[(++id1) + 10]; } 
  while (temp > 0. && id1 < nQuark);
  if (id1 == 0) id1 = 21; 
  id2 = -nQuark-1;
  temp = xPDF2sum * Rndm::flat();
  do { temp -= xPDF2[(++id2) + 10]; } 
  while (temp > 0. && id2 < nQuark);  
  if (id2 == 0) id2 = 21; 

  // Prepare to generate differential cross sections.
  sHat = tau * sCM;
  double root = sqrtpos(1. - xT2 / tau);
  double alpS = alphaS.alphaS(pT2shift);
  double dSigma1, dSigma2;

  // Loop over two symmetrical configurations (tHat <-> uHat).
  // (Not necessary, but makes for better MC efficiency.)
  for (int i = 1; i < 3; ++i) {
    SigmaProcess& dSigmaDtNow = (i == 1) ? dSigmaDt1 : dSigmaDt2;
    double tHatNow = (i == 1) ? -0.5 * sHat * (1. - root) 
      : -0.5 * sHat * (1. + root);
    dSigmaDtNow.setupShThAs(id1, id2, sHat, tHatNow, alpS);  
    double& dSigma = (i ==1) ? dSigma1 : dSigma2; 

    // Evaluate cross section based on flavours and kinematics.
    // Factor 4./9. per incoming gluon to compensate for preweighting.
    // g + g -> g + g.  
    if (id1 == 21 && id2 == 21) dSigma = (16./81.) * dSigmaDtNow.gg2gg();
    // q + g -> q + g.   
    else if (id1 == 21 || id2 == 21) dSigma = (4./9.) * dSigmaDtNow.qg2qg();
    // q + q -> q + q, involving identical quarks.
    else if (id1 == id2) dSigma = dSigmaDtNow.qq2qqSame();
    // q1 + q2 -> q1 + q2, where q1 and q2 are different.  
    else if (id1 * id2 > 0) dSigma = dSigmaDtNow.qq2qqDiff();
    // q + qbar -> q + qbar, where q and qbar are each others antiparticles.
    else if (id1 == -id2) dSigma = dSigmaDtNow.qqbar2qqbarAnti();
    // q1 + q2bar -> q1 + q2bar, where q1 and q2 are different.  
    else dSigma = dSigmaDtNow.qqbar2qqbarDiff();

    // Include possibility of K factor, common to all cross sections.
    dSigma *= Kfactor; 
  }

  // Form average dSigmaHat/dt, representing dSigmaHat/dpT2 (up to factors).
  double dSigmaPartonCorr = (dSigma1 + dSigma2) / 2.;

  // Combine cross section, pdf's and phase space integral.
  double volumePhSp = pow2(2. * yMax) / WTy;
  double dSigmaCorr = dSigmaPartonCorr * xPDF1sum * xPDF2sum * volumePhSp;

  // Dampen cross section at small pT values; part of formalism.
  dSigmaCorr *= pow2(pT2 / pT2shift);

  // Done
  return dSigmaCorr;
}

//*********

// Calculate factor relating matter overlap and interaction rate,
// i.e. k in <n_interaction(b)> = k * overlap(b) (neglecting
// n_int = 0 corrections and energy-momentum conservation effects).

void MultipleInteractions::overlapFactor() {

  // Initial values for iteration. Step size of b integration.
  double nAvg = sigmaInt / sigmaND;
  double kNow = 0.5;
  int stepDir = 1;
  double deltaB = BSTEP;
  if (bProfile == 2) deltaB *= min( 0.5, 2.5 * coreRadius); 
  if (bProfile == 3) deltaB *= max(1., pow(2. / expPow, 1. / expPow)); 
  
  // Further variables, with dummy initial values.
  double nNow = 0.;
  double kLow = 0.;
  double nLow = 0.;
  double kHigh = 0.;
  double nHigh = 0.;
  double overlapNow = 0.;
  double probNow = 0.; 
  double overlapInt = 0.5;
  double probInt = 0.; 
  double probOverlapInt = 0.;
  double normPi = 1. / (2. * M_PI);

  // Some common combinations for double Gaussian, as shorthand.
  fractionA = pow2(1. - coreFraction);
  fractionB = 2. * coreFraction * (1. - coreFraction);
  fractionC = pow2(coreFraction); 
  radius2B = 0.5 * (1. + pow2(coreRadius));
  radius2C = pow2(coreRadius);

  // First close k into an interval by binary steps,
  // then find k by successive interpolation.  
  do {
    if (stepDir == 1) kNow *= 2.;
    else if (stepDir == -1) kNow *= 0.5;
    else kNow = kLow + (nAvg - nLow) * (kHigh - kLow) / (nHigh - nLow);

    // Overlap trivial if no impact parameter dependence.
    if (bProfile <= 0 || bProfile > 3) {
      probInt = 0.5 * M_PI * (1. - exp(-kNow));
      probOverlapInt = probInt / M_PI;

    // Else integrate overlap over impact parameter.
    } else { 

      // Reset integrals.
      overlapInt = (bProfile == 3) ? 0. : 0.5;
      probInt = 0.; 
      probOverlapInt = 0.;

      // Step in b space.
      double b = -0.5 * deltaB;
      double bArea = 0.;
      do {
        b += deltaB;
        bArea = 2. * M_PI * b * deltaB;

        // Evaluate overlap at current b value.
        if (bProfile == 1) { 
          overlapNow = normPi * exp( -b*b);
	} else if (bProfile == 2) {          
          overlapNow = normPi * ( fractionA * exp( -min(EXPMAX, b*b))
            + fractionB * exp( -min(EXPMAX, b*b / radius2B)) / radius2B
            + fractionC * exp( -min(EXPMAX, b*b / radius2C)) / radius2C );
	} else {
          overlapNow = normPi * exp( -pow( b, expPow));  
          overlapInt += bArea * overlapNow;
	}

        // Calculate interaction probability and integrate.
        probNow = 1. - exp( -min(EXPMAX, M_PI * kNow * overlapNow));
        probInt += bArea * probNow;
        probOverlapInt += bArea * overlapNow * probNow;

      // Continue out in b until overlap too small.
      } while (b < 1. || b * probNow > BMAX); 
    }      
 
    // Ratio of b-integrated k * overlap / (1 - exp( - k * overlap)).
    nNow = M_PI * kNow * overlapInt / probInt;

    // Replace lower or upper limit of k.
    if (nNow < nAvg) {
      kLow = kNow;
      nLow = nNow;
      if (stepDir == -1) stepDir = 0;
    } else {
      kHigh = kNow;
      nHigh = nNow;
      if (stepDir == 1) stepDir = -1;
    } 

  // Continue iteration until convergence.
  } while (abs(nNow - nAvg) > KCONVERGE * nAvg);

  // Save relevant final numbers. Done.
  double avgOverlap  = probOverlapInt / probInt; 
  zeroIntCorr = probOverlapInt / overlapInt; 
  normOverlap = normPi * zeroIntCorr / avgOverlap;

}

//*********

// Prepare system for evolution.

void MultipleInteractions::prepare(double pTscale) {  

  // Pick interaction rate enhancement related to impact parameter.
  selectEnhance(pTscale);

}

//*********

// Pick interaction rate enhancement related to impact parameter.

void MultipleInteractions::selectEnhance(double pTscale) {

  // Default, valid for bProfile = 0. Also initial Sudakov.
  enhanceB = zeroIntCorr;
  if (bProfile <= 0 || bProfile > 3) return; 
  double pT2scale = pTscale*pTscale;
  double sudakov = 1.;

  // Begin loop over pT-dependent rejection of b value.
  do {

    // Flat enhancement distribution for simple Gaussian.
    if (bProfile == 1) {
      enhanceB = normOverlap * Rndm::flat();  

    // For double Gaussian go the way via b, according to each Gaussian.
    } else if (bProfile == 2) {
      double bType = Rndm::flat();  
      double b2 = -log( Rndm::flat() );
      if (bType < fractionA) ;
      else if (bType < fractionA + fractionB) b2 *= radius2B;
      else b2 *= radius2C; 
      enhanceB = normOverlap * ( fractionA * exp( -min(EXPMAX, b2))
        + fractionB * exp( -min(EXPMAX, b2 / radius2B)) / radius2B
        + fractionC * exp( -min(EXPMAX, b2 / radius2C)) / radius2C ); 

    // For exp( - b^expPow) transform to variable c = b^expPow so that
    // f(b) = b * exp( - b^expPow) -> f(c) = c^p * exp(-c) with p = expP. 
    // For  exp( - b^expPow) with expPow >= 2 <=> - 1 < p < 0: 
    // f(c) < c^p for c < 1,  f(c) < exp(-c) for c > 1.  
    } else if (bProfile == 3 && expPow > 1.999) {
      double expNow = max(2., expPow);
      double expP = 2. / expNow - 1.;
      double prob1 = expNow / (2. * exp(-1.) + expNow);
      double b2mod = 0.;
      double accept = 0.;
      do { 
        if (Rndm::flat() < prob1) {
          b2mod = pow( Rndm::flat(), 0.5 * expNow);
          accept = exp(-b2mod);
        } else {
          b2mod = 1. - log( Rndm::flat() );
          accept = pow( b2mod, expP);    
        } 
      } while (accept < Rndm::flat());
      enhanceB = normOverlap * exp(-b2mod);  

    // For exp( - b^expPow) with expPow < 2 <=> expP > 0: 
    // f(c) < p^p exp(-p) for c < 2p,  (2p)^p exp(-p-c/2) for c > 2p.
    } else if (bProfile == 3) {
      double expP = 2. / expPow - 1.;
      double prob1 = expP / (expP + pow(2., expP) * exp( - expP));
      double b2mod = 0.;
      double accept = 0.;  
      do {
        if (Rndm::flat() < prob1) {
          b2mod = 2. * expP * Rndm::flat();
          accept = pow( b2mod / expP, expP) * exp(expP - b2mod);
        } else {
          b2mod = 2. * (expP - log( Rndm::flat() )); 
          accept = pow(0.5 * b2mod / expP, expP) * exp(expP - 0.5 * b2mod);
        } 
      } while (accept < Rndm::flat()); 
      enhanceB = normOverlap *exp(-b2mod);  
    }

    // Evaluate "Sudakov form factor" for not having a harder interaction
    // at the selected b value, given the pT scale of the event.
    double xBin = (pT2scale - pT2min) * pT20maxR 
      / (pT2maxmin * (pT2scale + pT20R)); 
    xBin = max(1e-6, min(NBINS - 1e-6, NBINS * xBin) );
    int iBin = int(xBin);
    sudakov = exp( - (sudExpPT[iBin] 
      + (xBin - iBin) * (sudExpPT[iBin + 1] - sudExpPT[iBin]) ));  
  } while (sudakov < Rndm::flat());

  // Done.
}

//*********

// Select next pT in downwards evolution.

double MultipleInteractions::pTnext( BeamParticle& beamA, 
  BeamParticle& beamB, double pTbegAll, double pTendAll) {

  // Pick a pT using a quick-and-dirty cross section estimate.
  double pT2end = pow2( max(pTmin, pTendAll) );
  pT2 = pow2(pTbegAll);
  do pT2 = fastPT2(pT2);

  // Pick complete kinematics, evaluate correction factors to tell 
  // whether to keep the event.
  while (pT2 > pT2end && !acceptPT2( beamA, beamB)); 

  // Done.
  return (pT2 > pT2end) ? sqrt(pT2) : 0.;

}

//*********

// Pick a trial next pT, based on a simple upper estimate of the
// d(sigma)/d(pT2) spectrum.

double MultipleInteractions::fastPT2( double pT2beg) {

  // Use d(Prob)/d(pT2) < pT4dProbMax / (pT2 + r * pT20)^2. 
  double pT20begR = pT2beg + pT20R;
  double pT4dProbMaxNow = pT4dProbMax * enhanceB; 
  double pT2try = pT4dProbMaxNow * pT20begR 
    / (pT4dProbMaxNow - pT20begR * log(Rndm::flat())) - pT20R;

  // Save cross section associated with ansatz above. Done.
  dSigmaApprox = pT4dSigmaMax / pow2(pT2try + pT20R);
  return pT2try;

}

//*********

// Calculate the actual cross section to decide whether fast choice is OK.
// Select flavours and kinematics for interaction at given pT.
// Weighting may fail (return false) or succeed (return true).

bool MultipleInteractions::acceptPT2( BeamParticle& beamA, 
  BeamParticle& beamB) {
 
  // Derive shifted pT2 and rapidity limits from chosen pT2.
  pT2shift = pT2 + pT20;
  xT = 2. * sqrt(pT2) / eCM;
  xT2 = xT*xT;   
  double yMax = log(1./xT + sqrt(1./xT2 - 1.));

  // Select rapidities y3 and y4 of the two produced partons.
  double y3 = yMax * (2. * Rndm::flat() - 1.);
  double y4 = yMax * (2. * Rndm::flat() - 1.);
  y = 0.5 * (y3 + y4);

  // Reject some events at large rapidities to improve efficiency.
  // (Don't have to evaluate PDF's and ME's.)
  double WTy = (1. - pow2(y3/yMax)) * (1. - pow2(y4/yMax));
  if (WTy < Rndm::flat()) return false; 

  // Failure if x1 or x2 exceed what is left in respective beam.
  x1 = 0.5 * xT * (exp(y3) + exp(y4));
  x2 = 0.5 * xT * (exp(-y3) + exp(-y4));
  if (x1 > beamA.xMax() || x2 > beamB.xMax()) return false; 
  tau = x1 * x2;

  // Evaluate parton densities at actual x1 and x2.
  // Note that these densities take into account previous interactions.
  double xPDF1[21];
  double xPDF1sum = 0.;
  for (int id = -nQuark; id <= nQuark; ++id) {
    if (id == 0) xPDF1[10] = (9./4.) * beamA.xfMI(21, x1, pT2shift);
    else xPDF1[id+10] = beamA.xfMI(id, x1, pT2shift);
    xPDF1sum += xPDF1[id+10];
  }
  double xPDF2[21];
  double xPDF2sum = 0.;
  for (int id = -nQuark; id <= nQuark; ++id) {
    if (id == 0) xPDF2[10] = (9./4.) * beamB.xfMI(21, x2, pT2shift);
    else xPDF2[id+10] = beamB.xfMI(id, x2, pT2shift);
    xPDF2sum += xPDF2[id+10];
  }

  // Select incoming flavours according to actual PDF's.
  id1 = -nQuark-1;
  double temp = xPDF1sum * Rndm::flat();
  do { temp -= xPDF1[(++id1) + 10]; } 
  while (temp > 0. && id1 < nQuark);
  if (id1 == 0) id1 = 21; 
  id2 = -nQuark-1;
  temp = xPDF2sum * Rndm::flat();
  do { temp -= xPDF2[(++id2) + 10]; } 
  while (temp > 0. && id2 < nQuark);  
  if (id2 == 0) id2 = 21; 

  // Prepare to generate differential cross sections.
  sHat = tau * sCM;
  double root = sqrtpos(1. - xT2 / tau);
  double alpS = alphaS.alphaS(pT2shift);
  double dSigma1, dSigma2;

  // Loop over two symmetrical configurations (tHat <-> uHat).
  // (Not necessary, but makes for better MC efficiency.)
  for (int i = 1; i < 3; ++i) {
    SigmaProcess& dSigmaDtNow = (i == 1) ? dSigmaDt1 : dSigmaDt2;
    double tHatNow = (i == 1) ? -0.5 * sHat * (1. - root) 
      : -0.5 * sHat * (1. + root);
    dSigmaDtNow.setupShThAs(id1, id2, sHat, tHatNow, alpS);  
    double& dSigma = (i ==1) ? dSigma1 : dSigma2; 

    // Evaluate cross section based on flavours and kinematics.
    // Factor 4./9. per incoming gluon to compensate for preweighting.
    // g + g -> g + g.  
    if (id1 == 21 && id2 == 21) dSigma = (16./81.) * dSigmaDtNow.gg2gg();
    // q + g -> q + g.   
    else if (id1 == 21 || id2 == 21) dSigma = (4./9.) * dSigmaDtNow.qg2qg();
    // q + q -> q + q, involving identical quarks.
    else if (id1 == id2) dSigma = dSigmaDtNow.qq2qqSame();
    // q1 + q2 -> q1 + q2, where q1 and q2 are different.  
    else if (id1 * id2 > 0) dSigma = dSigmaDtNow.qq2qqDiff();
    // q + qbar -> q + qbar, where q and qbar are each others antiparticles.
    else if (id1 == -id2) dSigma = dSigmaDtNow.qqbar2qqbarAnti();
    // q1 + q2bar -> q1 + q2bar, where q1 and q2 are different.  
    else dSigma = dSigmaDtNow.qqbar2qqbarDiff();

    // Include possibility of K factor, common to all cross sections.
    dSigma *= Kfactor; 
  }

  // Form average dSigmaHat/dt, representing dSigmaHat/dpT2 (up to factors).
  double dSigmaPartonCorr = (dSigma1 + dSigma2) / 2.;

  // Combine cross section, pdf's and phase space integral.
  double volumePhSp = pow2(2. * yMax) / WTy;
  double dSigmaCorr = dSigmaPartonCorr * xPDF1sum * xPDF2sum * volumePhSp;

  // Dampen cross section at small pT values; part of formalism.
  dSigmaCorr *= pow2(pT2 / pT2shift);

  // Compare with approximate cross section used to pick pT value.
  double WTacc = dSigmaCorr / dSigmaApprox;
  if (WTacc < Rndm::flat()) return false; 

  // Pick one of two possible tHat values for acceptable pT.
  dSigmaDtSel = ((dSigma1+ dSigma2) * Rndm::flat() < dSigma1) 
    ? &dSigmaDt1 : &dSigmaDt2;

  // Done.
  return true;
}

//*********

// Set up the kinematics of the 2 -> 2 scattering process,
// and store the scattering in the event record.

bool MultipleInteractions::scatter( BeamParticle& beamA, 
  BeamParticle& beamB, Event& event) {

  // Set up kinematics in scattering subsystem.
  // Abort if consistent kinematics could not be found.
  if (!dSigmaDtSel->doKinematics()) return false;

  // Info for boost to event CM frame, history and colours.
  double betaZ = (x1 - x2) / (x1 + x2);
  int motherOffset = event.size();
  int colOffset = event.lastColTag();

  // Loop over four partons and update info as above.
  for (int i = 0; i < 4; ++i) {
    Particle parton = dSigmaDtSel->getParton(i);
    parton.bst(0., 0., betaZ);
    if (i < 2 ) parton.mothers( i + 1, 0);  
    else parton.mothers( motherOffset, motherOffset + 1);
    if (i < 2 ) parton.daughters( motherOffset + 2, motherOffset + 3);
    else parton.daughters( 0, 0);
    int col = parton.col();
    if (col > 0) parton.col( col + colOffset);
    int acol = parton.acol();
    if (acol > 0) parton.acol( acol + colOffset);

    // Put the partons into the event record.
    event.append(parton);
  }

  // Add scattered partons to list in beam remnants.
  int iA = beamA.append( motherOffset, id1, x1);
  int iB = beamB.append( motherOffset + 1, id2, x2);

  // Find whether incoming partons are valence or sea, so prepared for ISR.
  beamA.xfISR( iA, id1, x1, pT2);
  beamA.pickValSeaComp(); 
  beamB.xfISR( iB, id2, x2, pT2);
  beamB.pickValSeaComp(); 

  // Done.
  return true;
} 

//**************************************************************************

} // end namespace Pythia8
