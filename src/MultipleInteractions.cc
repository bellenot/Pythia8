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
double MultipleInteractions::alphaEMfix =  0.00729735;
double MultipleInteractions::Kfactor = 1.0; 
double MultipleInteractions::pT0Ref= 2.2; 
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

// Factorization scale pT2 by default, but could be shifted to pT2 + pT02.
// (A priori possible, but problems for flavour threshold interpretation.)
const bool MultipleInteractions::SHIFTFACSCALE = false;

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

// Define low-b region by interaction rate above given number.
const double MultipleInteractions::PROBATLOWB = 0.6;

// Basic step size for b integration; sometimes modified.
const double MultipleInteractions::BSTEP = 0.01;

// Stop b integration when integrand dropped enough.
const double MultipleInteractions::BMAX = 1e-8;

// Do not allow too large argument to exp function.
const double MultipleInteractions::EXPMAX = 50.;

// Convergence criterion for k iteration.
const double MultipleInteractions::KCONVERGE = 1e-7;

// Conversion of GeV^{-2} to mb for cross section.
const double MultipleInteractions::CONVERT2MB = 0.389380; 

//*********

// Initialize static data members.

void MultipleInteractions::initStatic() {

  //  Parameters of alphaStrong and cross section generation.
  alphaSvalue = Settings::parameter("MultipleInteractions:alphaSvalue");
  alphaSorder = Settings::mode("MultipleInteractions:alphaSorder");
  alphaEMfix = Settings::parameter("StandardModel:alphaEMfix"); 
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

bool MultipleInteractions::init( BeamParticle* beamAPtrIn, 
  BeamParticle* beamBPtrIn, bool reInit) {

  // Do not initialize if already done, and reinitialization not asked for.
  if (isInit && !reInit) return true;

  // Store input pointers for future use. 
  beamAPtr = beamAPtrIn;
  beamBPtr = beamBPtrIn;

  // Some common combinations for double Gaussian, as shorthand.
  if (bProfile == 2) {
    fracA = pow2(1. - coreFraction);
    fracB = 2. * coreFraction * (1. - coreFraction);
    fracC = pow2(coreFraction); 
    radius2B = 0.5 * (1. + pow2(coreRadius));
    radius2C = pow2(coreRadius);

  // Some common combinations for exp(b^pow), as shorthand.
  } else if (bProfile == 3) {
    lowPow = (expPow < 2.) ? true : false;
    expRev = 2. / expPow - 1.;
  } 

  // Initialize alpha_strong generation.
  alphaS.init( alphaSvalue, alphaSorder); 

  // Attach matrix-element calculation objects.
  if (!reInit) {
    sigma2gg2ggT = new Sigma2gg2gg();
    sigma2gg2ggU = new Sigma2gg2gg();
    sigma2qg2qgT = new Sigma2qg2qg();
    sigma2qg2qgU = new Sigma2qg2qg();
    sigma2qq2qqSameT = new Sigma2qq2qqSame();
    sigma2qq2qqSameU = new Sigma2qq2qqSame();
    sigma2qqbar2qqbarSameT = new Sigma2qqbar2qqbarSame();
    sigma2qqbar2qqbarSameU = new Sigma2qqbar2qqbarSame();
    sigma2qq2qqDiffT = new Sigma2qq2qqDiff();
    sigma2qq2qqDiffU = new Sigma2qq2qqDiff();
  }

  // Calculate invariant mass of system. Set current pT0 scale.
  sCM = m2( beamAPtr->p(), beamBPtr->p());
  eCM = sqrt(sCM);
  pT0 = pT0Ref * pow(eCM / ecmRef, ecmPow);

  // Get the total inelastic and nondiffractive cross section. Output.
  bool canDoMI = sigmaTot.init( beamAPtr->id(), beamBPtr->id(), eCM);
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
    upperEnvelope();

    // Integrate the parton-parton interaction cross section.
    jetCrossSection();

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
  overlapInit();

  // Done.
  isInit = true;
  return true;
}

//*********

// Select first = hardest pT in minbias process.
// Requires separate treatment at low and high b values

void MultipleInteractions::pTfirst() {  

  // Pick impact parameter and thereby interaction rate enhancement.
  overlapFirst();
  bSetInFirst = true;
  double WTacc;

  // At low b values evolve downwards with Sudakov. 
  if (atLowB) {
    pT2 = pT2max;
    do {

      // Pick a pT using a quick-and-dirty cross section estimate.
      pT2 = fastPT2(pT2);

      // If fallen below lower cutoff then need to restart at top.
      if (pT2 < pT2min) {
        pT2 = pT2max;
        WTacc = 0.;

      // Else pick complete kinematics and evaluate cross-section correction.
      } else WTacc = sigmaPT2(true) / dSigmaApprox;
    
    // Loop until acceptable pT and acceptable kinematics.
    } while (WTacc < Rndm::flat() || !dSigmaDtSel->final2KinMI()); 

  // At high b values make preliminary pT choice without Sudakov factor.
  } else {
    do {
      pT2 = pT20min0maxR / (pT20minR + Rndm::flat() * pT2maxmin) - pT20R; 

      // Evaluate upper estimate of cross section for this pT2 choice.  
      dSigmaApprox = pT4dSigmaMax / pow2(pT2 + pT20R);

      // Pick complete kinematics and evaluate cross-section correction.
      WTacc = sigmaPT2(true) / dSigmaApprox;

      // Evaluate and include Sudakov factor.
      WTacc *= sudakov( pT2, enhanceB);
    
    // Loop until acceptable pT and acceptable kinematics.
    } while (WTacc < Rndm::flat() || !dSigmaDtSel->final2KinMI()); 
  }
  
}

//*********

// Set up kinematics for first = hardest pT in minbias process.

void MultipleInteractions::setupFirstSys( Info* infoPtr, Event& process) { 

  // Remove any partons of previous failed interactions.
  if (process.size() > 3) {
    process.popBack( process.size() - 3);
    process.initColTag();
  }

  // Loop over four partons and offset info relative to subprocess itself.
  int colOffset = process.lastColTag();
  for (int i = 1; i <= 4; ++i) {
    Particle parton = dSigmaDtSel->getParton(i);
    if (i <= 2 ) parton.mothers( i, 0);  
    else parton.mothers( 3, 4);
    if (i <= 2 ) parton.daughters( 5, 6);
    else parton.daughters( 0, 0);
    int col = parton.col();
    if (col > 0) parton.col( col + colOffset);
    int acol = parton.acol();
    if (acol > 0) parton.acol( acol + colOffset);

    // Put the partons into the event record.
    process.append(parton);
  }

  // Set scale from which to begin evolution.
  process.scale(  sqrt(pT2Fac) );

  // Info on subprocess - specific to mimimum-bias events.
  string nameSub = dSigmaDtSel->name();
  int codeSub = dSigmaDtSel->code();
  int nFinalSub = dSigmaDtSel->nFinal();
  infoPtr->setSubType( nameSub, codeSub, nFinalSub);

  // Further standard info on process.
  infoPtr->setPDFalpha( id1, id2, xPDF1now, xPDF2now, pT2Fac, alpEM, alpS, 
    pT2Ren);
  double m3 = dSigmaDtSel->m(3);
  double m4 = dSigmaDtSel->m(4); 
  double theta = dSigmaDtSel->theta(); 
  double phi = dSigmaDtSel->phi(); 
  infoPtr->setKin( x1, x2, sHat, tHat, uHat, sqrt(pT2), m3, m4, theta, phi);

}

//*********

// Select next pT in downwards evolution.

double MultipleInteractions::pTnext( double pTbegAll, double pTendAll) {

  // Pick a pT using a quick-and-dirty cross section estimate.
  double WTacc;
  double pT2end = pow2( max(pTmin, pTendAll) );
  pT2 = pow2(pTbegAll);
  do {
    pT2 = fastPT2(pT2);

    // Pick complete kinematics and evaluate cross-section correction.
    WTacc = sigmaPT2(false) / dSigmaApprox;
 
    // Decide whether to keep the event.
    if (pT2 < pT2end) return 0.;
  } while (WTacc < Rndm::flat() || !dSigmaDtSel->final2KinMI()); 

  // Done.
  return sqrt(pT2);

}

//*********

// Set up the kinematics of the 2 -> 2 scattering process,
// and store the scattering in the event record.

void MultipleInteractions::scatter( Event& event) {

  // Loop over four partons and offset info relative to subprocess itself.
  int motherOffset = event.size();
  int colOffset = event.lastColTag();
  for (int i = 1; i <= 4; ++i) {
    Particle parton = dSigmaDtSel->getParton(i);
    if (i <= 2 ) parton.mothers( i, 0);  
    else parton.mothers( motherOffset, motherOffset + 1);
    if (i <= 2 ) parton.daughters( motherOffset + 2, motherOffset + 3);
    else parton.daughters( 0, 0);
    int col = parton.col();
    if (col > 0) parton.col( col + colOffset);
    int acol = parton.acol();
    if (acol > 0) parton.acol( acol + colOffset);

    // Put the partons into the event record.
    event.append(parton);
  }

  // Add scattered partons to list in beam remnants.
  int iA = beamAPtr->append( motherOffset, id1, x1);
  int iB = beamBPtr->append( motherOffset + 1, id2, x2);

  // Find whether incoming partons are valence or sea, so prepared for ISR.
  beamAPtr->xfISR( iA, id1, x1, pT2);
  beamAPtr->pickValSeaComp(); 
  beamBPtr->xfISR( iB, id2, x2, pT2);
  beamBPtr->pickValSeaComp(); 

  // Done.
} 

//*********

// Determine constant in d(Prob)/d(pT2) < const / (pT2 + r * pT20)^2.  

void MultipleInteractions::upperEnvelope() {  

  // Initially determine constant in jet cross section upper estimate 
  // d(sigma_approx)/d(pT2) < const / (pT2 + r * pT20)^2. 
  pT4dSigmaMax = 0.;
  
  // Loop thorough allowed pT range logarithmically evenly.
  for (int iPT = 0; iPT < NBINS; ++iPT) {
    double pT = pTmin * pow( pTmax/pTmin, (iPT + 0.5) / NBINS);
    pT2 = pT*pT;
    pT2shift = pT2 + pT20;
    pT2Ren = pT2shift;
    pT2Fac = (SHIFTFACSCALE) ? pT2shift : pT2;
    xT = 2. * pT / eCM;

    // Evaluate parton density sums at x1 = x2 = xT.
    double xPDF1sumMax = (9./4.) * beamAPtr->xf(21, xT, pT2Fac);
    for (int id = 1; id <= nQuark; ++id) 
      xPDF1sumMax += beamAPtr->xf(id, xT, pT2Fac) 
        + beamAPtr->xf(-id, xT, pT2Fac);
    double xPDF2sumMax = (9./4.) * beamBPtr->xf(21, xT, pT2Fac);
    for (int id = 1; id <= nQuark; ++id)
      xPDF2sumMax += beamBPtr->xf(id, xT, pT2Fac) 
        + beamBPtr->xf(-id, xT, pT2Fac);

    // Evaluate alpha_strong, matrix element and phase space volume.
    alpS = alphaS.alphaS(pT2Ren);
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

void MultipleInteractions::jetCrossSection() {

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
      double dSigma = sigmaPT2(true);

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

// Evaluate "Sudakov form factor" for not having a harder interaction
// at the selected b value, given the pT scale of the event.

double MultipleInteractions::sudakov(double pT2sud, double enhance) {

  // Find bin the pT2 scale falls in.
  double xBin = (pT2sud - pT2min) * pT20maxR 
    / (pT2maxmin * (pT2sud + pT20R)); 
  xBin = max(1e-6, min(NBINS - 1e-6, NBINS * xBin) );
  int iBin = int(xBin);

  // Interpolate inside bin. Optionally include enhancement factor.
  double sudExp = sudExpPT[iBin] + (xBin - iBin) 
    * (sudExpPT[iBin + 1] - sudExpPT[iBin]);
  return exp( -enhance * sudExp);
  
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
// Slightly different treatment for first interaction and subsequent ones.

double MultipleInteractions::sigmaPT2(bool isFirst) {
 
  // Derive shifted pT2 and rapidity limits from chosen pT2.
  pT2shift = pT2 + pT20;
  pT2Ren = pT2shift;
  pT2Fac = (SHIFTFACSCALE) ? pT2shift : pT2;
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

  // Failure if x1 or x2 exceed what is left in respective beam.
  x1 = 0.5 * xT * (exp(y3) + exp(y4));
  x2 = 0.5 * xT * (exp(-y3) + exp(-y4));
  if (isFirst) {
    if (x1 > 1. || x2 > 1.) return 0.; 
  } else {
    if (x1 > beamAPtr->xMax() || x2 > beamBPtr->xMax()) return 0.; 
  }
  tau = x1 * x2;

  // Begin evaluate parton densities at actual x1 and x2.
  double xPDF1[21];
  double xPDF1sum = 0.;
  double xPDF2[21];
  double xPDF2sum = 0.;

  // For first interaction use normal densities.
  if (isFirst) {
    for (int id = -nQuark; id <= nQuark; ++id) {
      if (id == 0) xPDF1[10] = (9./4.) * beamAPtr->xf(21, x1, pT2Fac);
      else xPDF1[id+10] = beamAPtr->xf(id, x1, pT2Fac);
      xPDF1sum += xPDF1[id+10];
    }
    for (int id = -nQuark; id <= nQuark; ++id) {
      if (id == 0) xPDF2[10] = (9./4.) * beamBPtr->xf(21, x2, pT2Fac);
      else xPDF2[id+10] = beamBPtr->xf(id, x2, pT2Fac);
      xPDF2sum += xPDF2[id+10];
    }
  
  // For subsequent interactions use rescaled densities.
  } else {
    for (int id = -nQuark; id <= nQuark; ++id) {
      if (id == 0) xPDF1[10] = (9./4.) * beamAPtr->xfMI(21, x1, pT2Fac);
      else xPDF1[id+10] = beamAPtr->xfMI(id, x1, pT2Fac);
      xPDF1sum += xPDF1[id+10];
    }
    for (int id = -nQuark; id <= nQuark; ++id) {
      if (id == 0) xPDF2[10] = (9./4.) * beamBPtr->xfMI(21, x2, pT2Fac);
      else xPDF2[id+10] = beamBPtr->xfMI(id, x2, pT2Fac);
      xPDF2sum += xPDF2[id+10];
    }
  }

  // Select incoming flavours according to actual PDF's.
  id1 = -nQuark-1;
  double temp = xPDF1sum * Rndm::flat();
  do { xPDF1now = xPDF1[(++id1) + 10]; temp -= xPDF1now; } 
  while (temp > 0. && id1 < nQuark);
  if (id1 == 0) id1 = 21; 
  id2 = -nQuark-1;
  temp = xPDF2sum * Rndm::flat();
  do { xPDF2now = xPDF2[(++id2) + 10]; temp -= xPDF2now;} 
  while (temp > 0. && id2 < nQuark);  
  if (id2 == 0) id2 = 21; 

  // Assign pointers to processes relevant for incoming flavour choice.
  // Factor 4./9. per incoming gluon to compensate for preweighting.  
  SigmaProcess *dSigmaT, *dSigmaU;
  double gluFac = 1.;

  // g + g -> g + g.  
  if (id1 == 21 && id2 == 21) {
    dSigmaT = sigma2gg2ggT;
    dSigmaU = sigma2gg2ggU;
    gluFac = 16. / 81.;

  // q + g -> q + g.   
  } else if (id1 == 21 || id2 == 21) {
    dSigmaT = sigma2qg2qgT;
    dSigmaU = sigma2qg2qgU;
    gluFac = 4. / 9.;

  // q + q -> q + q, involving identical quarks.
  } else if (id1 == id2) {
    dSigmaT = sigma2qq2qqSameT;
    dSigmaU = sigma2qq2qqSameU;

  // q + qbar -> q + qbar, pair of same flavour.
  } else if (id1 == -id2) {
    dSigmaT = sigma2qqbar2qqbarSameT;
    dSigmaU = sigma2qqbar2qqbarSameU;

  // q1 + q2 -> q1 + q2 or q1 + q2bar -> q1 + q2bar, different flavours.  
  } else {
    dSigmaT = sigma2qq2qqDiffT;
    dSigmaU = sigma2qq2qqDiffU;
  }

  // Prepare to generate differential cross sections.
  sHat = tau * sCM;
  double root = sqrtpos(1. - xT2 / tau);
  alpS = alphaS.alphaS(pT2Ren);
  alpEM = alphaEMfix;

  // Set kinematics for two symmetrical configurations (tHat <-> uHat).
  // (Not necessary, but makes for better MC efficiency.)
  tHat = -0.5 * sHat * (1. - root);
  uHat = -0.5 * sHat * (1. + root);
  dSigmaT->set2KinMI( id1, id2, x1, x2, sHat, tHat, uHat, alpS);
  dSigmaU->set2KinMI( id1, id2, x1, x2, sHat, uHat, tHat, alpS);

  // Evaluate cross sections, include possibility of K factor.
  double dSigmaNowT = Kfactor * gluFac * dSigmaT->sigmaHat();
  double dSigmaNowU = Kfactor * gluFac * dSigmaU->sigmaHat();

  // Form average dSigmaHat/dt, representing dSigmaHat/dpT2 (up to factors).
  double dSigmaPartonCorr = (dSigmaNowT + dSigmaNowU) / 2.;

  // Pick one of two possible tHat values for acceptable pT.
  dSigmaDtSel = ((dSigmaNowT + dSigmaNowU) * Rndm::flat() < dSigmaNowT) 
    ? dSigmaT : dSigmaU;
  if (dSigmaDtSel == dSigmaU) swap( tHat, uHat);

  // Combine cross section, pdf's and phase space integral.
  double volumePhSp = pow2(2. * yMax) / WTy;
  double dSigmaCorr = dSigmaPartonCorr * xPDF1sum * xPDF2sum * volumePhSp;

  // Dampen cross section at small pT values; part of formalism.
  dSigmaCorr *= pow2(pT2 / pT2shift);

  // Done.
  return dSigmaCorr;
}


//*********

// Calculate factor relating matter overlap and interaction rate,
// i.e. k in <n_interaction(b)> = k * overlap(b) (neglecting
// n_int = 0 corrections and energy-momentum conservation effects).

void MultipleInteractions::overlapInit() {

  // Initial values for iteration. Step size of b integration.
  nAvg = sigmaInt / sigmaND;
  kNow = 0.5;
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
  double bProbInt = 0.;
  normPi = 1. / (2. * M_PI);

  // Subdivision into low-b and high-b region by interaction rate.
  bool pastBDiv = false;  
  double overlapHighB = 0.;

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
      bProbInt = probInt;

    // Else integrate overlap over impact parameter.
    } else { 

      // Reset integrals.
      overlapInt = (bProfile == 3) ? 0. : 0.5;
      probInt = 0.; 
      probOverlapInt = 0.;
      bProbInt = 0.;
      pastBDiv = false;
      overlapHighB = 0.;

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
          overlapNow = normPi * ( fracA * exp( -min(EXPMAX, b*b))
            + fracB * exp( -min(EXPMAX, b*b / radius2B)) / radius2B
            + fracC * exp( -min(EXPMAX, b*b / radius2C)) / radius2C );
	} else {
          overlapNow = normPi * exp( -pow( b, expPow));  
          overlapInt += bArea * overlapNow;
	}
        if (pastBDiv) overlapHighB += bArea * overlapNow;

        // Calculate interaction probability and integrate.
        probNow = 1. - exp( -min(EXPMAX, M_PI * kNow * overlapNow));
        probInt += bArea * probNow;
        probOverlapInt += bArea * overlapNow * probNow;
        bProbInt += b * bArea * probNow;

        // Check when interaction probability has dropped sufficiently.
        if (!pastBDiv && probNow < PROBATLOWB) {
          bDiv = b + 0.5 * deltaB;
          pastBDiv = true;
        }

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

  // Save relevant final numbers for overlap values.
  double avgOverlap  = probOverlapInt / probInt; 
  zeroIntCorr = probOverlapInt / overlapInt; 
  normOverlap = normPi * zeroIntCorr / avgOverlap;
  bAvg = bProbInt / probInt;

  // Relative rates for preselection of low-b and high-b region.
  // Other useful combinations for subsequent selection.
  if (bProfile > 0 && bProfile <= 3) {
    probLowB = M_PI * bDiv*bDiv;
    double probHighB = M_PI * kNow * overlapHighB;
    if (bProfile == 1) probHighB = M_PI * kNow * 0.5 * exp( -bDiv*bDiv);
    else if (bProfile == 2) {
      fracAhigh = fracA * exp( -bDiv*bDiv);
      fracBhigh = fracB * exp( -bDiv*bDiv / radius2B);
      fracChigh = fracC * exp( -bDiv*bDiv / radius2C);
      fracABChigh = fracAhigh + fracBhigh + fracChigh;
      probHighB = M_PI * kNow * 0.5 * fracABChigh;
    } else { 
      cDiv = pow( bDiv, expPow);
      cMax = max(2. * expRev, cDiv); 
    } 
    probLowB /= (probLowB + probHighB);
  }

}

//*********

// Pick impact parameter and interaction rate enhancement beforehand,
// i.e. before even the hardest interaction for minimum-bias events. 

void MultipleInteractions::overlapFirst() {

  // Trivial values if no impact parameter dependence.
  if (bProfile <= 0 || bProfile > 3) {
    bNow = bAvg;
    enhanceB = zeroIntCorr;
    atLowB = true;
    return;
  }

  // Preliminary choice between and inside low-b and high-b regions.
  double overlapNow = 0.;
  double probAccept = 0.;
  do {

    // Treatment in low-b region: pick b flat in area.
    if (Rndm::flat() < probLowB) {
      atLowB = true;
      bNow = bDiv * sqrt(Rndm::flat());
      if (bProfile == 1) overlapNow = normPi * exp( -bNow*bNow);
      else if (bProfile == 2) overlapNow = normPi * 
        ( fracA * exp( -bNow*bNow)
        + fracB * exp( -bNow*bNow / radius2B) / radius2B
        + fracC * exp( -bNow*bNow / radius2C) / radius2C );
      else overlapNow = normPi * exp( -pow( bNow, expPow));
      probAccept = 1. - exp( -min(EXPMAX, M_PI * kNow * overlapNow));

    // Treatment in high-b region: pick b according to overlap.
    } else {

      // For simple and double Gaussian pick b according to exp(-b^2 / r^2).
      if (bProfile == 1) {
        bNow = sqrt(bDiv*bDiv - log(Rndm::flat()));
        overlapNow = normPi * exp( -min(EXPMAX, bNow*bNow));
      } else if (bProfile == 2) {
        double pickFrac = Rndm::flat() * fracABChigh; 
        if (pickFrac < fracAhigh) bNow = sqrt(bDiv*bDiv - log(Rndm::flat()));
        else if (pickFrac < fracAhigh + fracBhigh) 
          bNow = sqrt(bDiv*bDiv - radius2B * log(Rndm::flat()));
        else bNow = sqrt(bDiv*bDiv - radius2C * log(Rndm::flat()));
        overlapNow = normPi * ( fracA * exp( -min(EXPMAX, bNow*bNow))
          + fracB * exp( -min(EXPMAX, bNow*bNow / radius2B)) / radius2B
          + fracC * exp( -min(EXPMAX, bNow*bNow / radius2C)) / radius2C );

      // For exp( - b^expPow) transform to variable c = b^expPow so that
      // f(b) = b * exp( - b^expPow) -> f(c) = c^r * exp(-c) with r = expRev. 
      // case lowPow: expPow < 2 <=> r > 0: preselect according to
      // f(c) < N exp(-c/2) and then accept with N' * c^r * exp(-c/2). 
      } else if (lowPow) {
        double cNow, acceptC;
        do {      
          cNow = cDiv - 2. * log(Rndm::flat());
          acceptC = pow(cNow / cMax, expRev) * exp( -0.5 * (cNow - cMax));
        } while (acceptC < Rndm::flat());
        bNow = pow( cNow, 1. / expPow);
        overlapNow = normPi * exp( -cNow);

      // case !lowPow: expPow >= 2 <=> - 1 < r < 0: preselect according to
      // f(c) < N exp(-c) and then accept with N' * c^r. 
      } else {
        double cNow, acceptC;
        do {      
          cNow = cDiv - log(Rndm::flat());
          acceptC = pow(cNow / cDiv, expRev);
        } while (acceptC < Rndm::flat());
        bNow = pow( cNow, 1. / expPow);
        overlapNow = normPi * exp( -cNow);    
      }
      double temp = M_PI * kNow *overlapNow;
      probAccept = (1. - exp( -min(EXPMAX, temp))) / temp;    
    }

  // Confirm choice of b value. Derive enhancement factor.
  } while (probAccept < Rndm::flat());
  enhanceB = (normOverlap / normPi) * overlapNow ;  

}

//*********

// Pick impact parameter and interaction rate enhancement afterwards,
// i.e. after a hard interaction is known but before rest of MI treatment.

void MultipleInteractions::overlapNext(double pTscale) {

  // Default, valid for bProfile = 0. Also initial Sudakov.
  enhanceB = zeroIntCorr;
  if (bProfile <= 0 || bProfile > 3) return; 
  double pT2scale = pTscale*pTscale;

  // Begin loop over pT-dependent rejection of b value.
  do {

    // Flat enhancement distribution for simple Gaussian.
    if (bProfile == 1) {
      double expb2 = Rndm::flat();
      enhanceB = normOverlap * expb2;  
      bNow = sqrt( -log(expb2));

    // For double Gaussian go the way via b, according to each Gaussian.
    } else if (bProfile == 2) {
      double bType = Rndm::flat();  
      double b2 = -log( Rndm::flat() );
      if (bType < fracA) ;
      else if (bType < fracA + fracB) b2 *= radius2B;
      else b2 *= radius2C; 
      enhanceB = normOverlap * ( fracA * exp( -min(EXPMAX, b2))
        + fracB * exp( -min(EXPMAX, b2 / radius2B)) / radius2B
        + fracC * exp( -min(EXPMAX, b2 / radius2C)) / radius2C ); 
      bNow = sqrt(b2);

    // For exp( - b^expPow) transform to variable c = b^expPow so that
    // f(b) = b * exp( - b^expPow) -> f(c) = c^r * exp(-c) with r = expRev. 
    // case lowPow: expPow < 2 <=> r > 0: 
    // f(c) < r^r exp(-r) for c < 2r, < (2r)^r exp(-r-c/2) for c > 2r.
    } else if (lowPow) {
      double cNow, acceptC;
      double probLowC = expRev / (expRev + pow(2., expRev) * exp( - expRev));
      do {
        if (Rndm::flat() < probLowC) {
          cNow = 2. * expRev * Rndm::flat();
          acceptC = pow( cNow / expRev, expRev) * exp(expRev - cNow);
        } else {
          cNow = 2. * (expRev - log( Rndm::flat() )); 
          acceptC = pow(0.5 * cNow / expRev, expRev) * exp(expRev - 0.5 * cNow);
        }
      } while (acceptC < Rndm::flat()); 
      enhanceB = normOverlap *exp(-cNow);  
      bNow = pow( cNow, 1. / expPow);

    // case !lowPow: expPow >= 2 <=> - 1 < r < 0: 
    // f(c) < c^r for c < 1,  f(c) < exp(-c) for c > 1.  
    } else {
      double cNow, acceptC;
      double probLowC = expPow / (2. * exp(-1.) + expPow);
      do { 
        if (Rndm::flat() < probLowC) {
          cNow = pow( Rndm::flat(), 0.5 * expPow);
          acceptC = exp(-cNow);
        } else {
          cNow = 1. - log( Rndm::flat() );
          acceptC = pow( cNow, expRev);    
        } 
      } while (acceptC < Rndm::flat());
      enhanceB = normOverlap * exp(-cNow);  
      bNow = pow( cNow, 1. / expPow);
    }

    // Evaluate "Sudakov form factor" for not having a harder interaction.
  } while (sudakov(pT2scale, enhanceB) < Rndm::flat());

  // Done.
}

//**************************************************************************

} // end namespace Pythia8
