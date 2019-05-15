// Function definitions (not found in the header) for the 
// ProcessContainer and SetupContainers classes.
// Copyright C 2007 Torbjorn Sjostrand

#include "ProcessContainer.h"

// Internal headers for special processes.
#include "SigmaQCD.h"
#include "SigmaEW.h"
#include "SigmaOnia.h"
#include "SigmaHiggs.h"
#include "SigmaSUSY.h"

namespace Pythia8 {

//**************************************************************************

// ProcessContainer class.
// Information allowing the generation of a specific process.

//*********
 
// Definitions of static variables and functions.

// Pointer to the information object.
Info* ProcessContainer::infoPtr;

// Pointer to the resonance decay object.
ResonanceDecays* ProcessContainer::resonanceDecaysPtr;

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of event tries to check maximization finding reliability.
const int ProcessContainer::NSAMPLE = 100;

//*********

// Initialize 

bool ProcessContainer::init() {

  // Extract info about current process from SigmaProcess object.
  isMinBias  = sigmaProcessPtr->isMinBias();
  isResolved = sigmaProcessPtr->isResolved();
  isDiffA    = sigmaProcessPtr->isDiffA();
  isDiffB    = sigmaProcessPtr->isDiffB();
  int nFinal = sigmaProcessPtr->nFinal();

  // Pick and create phase space generator. Send pointers where required.
  if (isMinBias) phaseSpacePtr = new PhaseSpace2to2minbias();
  else if (!isResolved) phaseSpacePtr = new PhaseSpace2to2eldiff( isDiffA, 
    isDiffB);
  else if (nFinal == 1) phaseSpacePtr = new PhaseSpace2to1tauy();
  else phaseSpacePtr = new PhaseSpace2to2tauyz();

  // Store common info for PhaseSpace objects.
  double eCM = infoPtr->eCM();
  phaseSpacePtr->initInfo( sigmaProcessPtr, eCM);

  // Reset cross section statistics.
  nTry      = 0;
  nSel      = 0;
  nAcc      = 0;
  nTryStat  = 0;
  sigmaMx   = 0.;
  sigmaSum  = 0.;
  sigma2Sum = 0.;
  sigmaNeg  = 0.;
  sigmaAvg  = 0.;
  sigmaFin  = 0.;
  deltaFin  = 0.;

  // Initialize process. Remove empty inFlux channels and optionally list.
  sigmaProcessPtr->initProc();
  sigmaProcessPtr->initFlux();
  sigmaProcessPtr->checkChannels();

  // Find maximum of differential cross section * phasespace.
  bool physical       = phaseSpacePtr->setupSampling();
  sigmaMx             = phaseSpacePtr->sigmaMax();
  double sigmaHalfWay = sigmaMx;

  // Check maximum by a few events, and extrapolate a further increase.
  if (physical) {
    for (int sample = 0; sample < NSAMPLE; ++sample) 
    while (!phaseSpacePtr->trialKin(false)) { 
      if (sample == NSAMPLE/2) sigmaHalfWay = phaseSpacePtr->sigmaMax();
    }   
    sigmaMx = pow2(phaseSpacePtr->sigmaMax()) / sigmaHalfWay;
    phaseSpacePtr->setSigmaMax(sigmaMx);
  }

  // Done.
  return physical;
}

//*********

// Generate a trial event; selected or not.
 
bool ProcessContainer::trialProcess() { 

  // Update number of tries.
  if (sigmaMx == 0.) return false;
  ++nTry;

  // Generate a trial phase space point, with cross section.
  if (!phaseSpacePtr->trialKin(true)) return false;
  double sigmaNow = phaseSpacePtr->sigmaNow(); 

  // Check that not negative cross section.
  if (sigmaNow < sigmaNeg) {
    ErrorMsg::message("Warning in ProcessContainer::trialProcess:"
      " negative cross section set 0", "for " +  sigmaProcessPtr->name() );
    sigmaNeg = sigmaNow;
  }
  if (sigmaNow < 0.) sigmaNow = 0.;

  // Update statistics and maximum.
  sigmaSum  += sigmaNow;
  sigma2Sum += pow2(sigmaNow);
  sigmaMx    = phaseSpacePtr->sigmaMax();

  // Select or reject trial point.
  bool select = (sigmaNow > Rndm::flat() * sigmaMx);  
  if (select) ++nSel;
  return select;

}

//*********
  
// Give the hard subprocess.

bool ProcessContainer::constructProcess( Event& process, bool isHardest) { 

  // Construct flavour and colours for selected event.
  if (isResolved && !isMinBias) sigmaProcessPtr->pickInState();
  sigmaProcessPtr->setIdColAcol();

  // Construct kinematics from selected phase space point.
  if (!phaseSpacePtr->finalKin()) return false;

  // Basic info on process.
  if (isHardest) infoPtr->setType( name(), code(), nFinal(), isMinBias, 
    isResolved, isDiffA, isDiffB);

  // Let hard process record begin with the event as a whole and
  // the two incoming beam particles.  
  process.append( 90, -11, 0, 0, 1, 2, 0, 0, 
    Vec4(0., 0., 0., infoPtr->eCM()), infoPtr->eCM(), 0. ); 
  process.append( infoPtr->idA(), -12, 0, 0, 3, 0, 0, 0, 
    Vec4(0., 0., infoPtr->pzA(), infoPtr->eA()), infoPtr->mA(), 0. ); 
  process.append( infoPtr->idB(), -12, 0, 0, 4, 0, 0, 0, 
    Vec4(0., 0., infoPtr->pzB(), infoPtr->eB()), infoPtr->mB(), 0. ); 

  // For minbias process no interaction selected so far, so done.
  if (isMinBias) return true;

  // Further info on process.
  int id1        = sigmaProcessPtr->id(1);
  int id2        = sigmaProcessPtr->id(2);
  double pdf1    = sigmaProcessPtr->pdf1();
  double pdf2    = sigmaProcessPtr->pdf2();
  double Q2Fac   = sigmaProcessPtr->Q2Fac();
  double alphaEM = sigmaProcessPtr->alphaEMH();
  double alphaS  = sigmaProcessPtr->alphaSH();
  double Q2Ren   = sigmaProcessPtr->Q2Ren();
  double x1      = phaseSpacePtr->x1();
  double x2      = phaseSpacePtr->x2();
  double sHat    = phaseSpacePtr->sHat();
  double tHat    = phaseSpacePtr->tHat();
  double uHat    = phaseSpacePtr->uHat();
  double pTHat   = phaseSpacePtr->pTHat();
  double m3      = phaseSpacePtr->m(3);
  double m4      = phaseSpacePtr->m(4);
  double theta   = phaseSpacePtr->thetaHat();
  double phi     = phaseSpacePtr->phiHat();
  if (isHardest) {
    infoPtr->setPDFalpha( id1, id2, pdf1, pdf2, Q2Fac, alphaEM, alphaS, Q2Ren);
    infoPtr->setKin( x1, x2, sHat, tHat, uHat, pTHat, m3, m4, theta, phi);
  }
  infoPtr->setTypeMI( code(), pTHat);

  // Insert the subprocess partons - resolved processes.
  if (isResolved) {

    // Set subprocess scale for ISR and MI.
    double scale = sqrt(Q2Fac);
    process.scale( scale );

    // Loop over incoming and outgoing partons.
    int colOffset = process.lastColTag();
    for (int i = 1; i <= 2 + sigmaProcessPtr->nFinal(); ++i) { 
      int id = sigmaProcessPtr->id(i);
      int status = (i <= 2) ? -21 : 23;
      int mother1 = (i <= 2) ? i : 3;
      int mother2 = (i <= 2) ? 0 : 4;
      int daughter1 = (i <= 2) ? 5 : 0;
      int daughter2 = (i <= 2) ? 4 + sigmaProcessPtr->nFinal() : 0;
      int col = sigmaProcessPtr->col(i);
      if (col > 0) col += colOffset;
      int acol = sigmaProcessPtr->acol(i);
      if (acol > 0) acol += colOffset;
      process.append( id, status, mother1, mother2, daughter1, daughter2, 
        col, acol, phaseSpacePtr->p(i), phaseSpacePtr->m(i), scale);
    }

  // Insert the outgoing particles - unresolved processes.
  } else { 
    process.append( sigmaProcessPtr->id(3), 23, 1, 0, 0, 0, 0, 0, 
      phaseSpacePtr->p(3), phaseSpacePtr->m(3));
    process.append( sigmaProcessPtr->id(4), 23, 2, 0, 0, 0, 0, 0, 
      phaseSpacePtr->p(4), phaseSpacePtr->m(4));
  }

  // Done.
  return true;

}

//*********
  
// Handle resonance decays.

bool ProcessContainer::decayResonances( Event& process) {

  // Save current event-record size.
  process.saveSize();
  bool physical    = true;
  bool newFlavours = false;

  // Do sequential chain of uncorrelated isotropic decays.
  do {
    physical = resonanceDecaysPtr->next( process);
    if (!physical) return false;

    // Check whether flavours should be correlated.
    // (Currently only relevant for f fbar -> gamma*/Z0 gamma*/Z0.)
    newFlavours = ( sigmaProcessPtr->weightDecayFlav( process) 
                  < Rndm::flat() ); 

    // Reset the decay chains if have to redo.
    if (newFlavours) {
      process.restoreSize();
      for (int i = 5; i < process.size(); ++i) process[i].statusPos();
    } 

  // Loop back where required to generate new decays with new flavours.    
  } while (newFlavours);

  // Correct to nonisotropic decays.
  phaseSpacePtr->decayKinematics( process); 

  // Done.
  return true;

}

//*********

// Estimate integrated cross section and its uncertainty.

void ProcessContainer::sigmaDelta() {

  // Initial values. No analysis meaningful unless accepted events.
  nTryStat = nTry;
  sigmaAvg = 0.;
  sigmaFin = 0.;
  deltaFin = 0.;
  if (nAcc == 0) return;

  // Average value.
  double nTryInv  = 1. / nTry;
  double nSelInv  = 1. / nSel;
  double nAccInv  = 1. / nAcc;
  sigmaAvg = sigmaSum * nTryInv ;
  double fracAcc  = nAcc * nSelInv;
  sigmaFin        = sigmaAvg * fracAcc;

  // Estimated error. Quadratic sum of cross section term and
  // binomial from accept/reject step.
  deltaFin = sigmaFin;
  if (nAcc == 1) return;
  double delta2Sig   = (sigma2Sum *nTryInv - pow2(sigmaAvg)) * nTryInv
    / pow2(sigmaAvg);
  double delta2Veto  = (nSel - nAcc) * nAccInv * nSelInv;
  double delta2Sum   = delta2Sig + delta2Veto;
  deltaFin           = sqrtpos(delta2Sum) * sigmaFin; 

}
 
//**************************************************************************

// SetupContainer class.
// Turns list of user-desired processes into a vector of containers.

//*********

// Main routine to initialize list of processes.

bool SetupContainers::init(vector<ProcessContainer*>& containerPtrs) {

  // Reset process list, if filled in previous subrun.
  if (containerPtrs.size() > 0) {
    for (int i = 0; i < int(containerPtrs.size()); ++i) 
      delete containerPtrs[i];
    containerPtrs.clear(); 
  }
  SigmaProcess* sigmaPtr;

  // Set up requested objects for soft QCD processes.
  bool softQCD = Settings::flag("SoftQCD:all");
  if (softQCD || Settings::flag("SoftQCD:minBias")) {
    sigmaPtr = new Sigma0minBias;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (softQCD || Settings::flag("SoftQCD:elastic")) {
    sigmaPtr = new Sigma0AB2AB;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (softQCD || Settings::flag("SoftQCD:singleDiffractive")) {
    sigmaPtr = new Sigma0AB2XB;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma0AB2AX;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (softQCD || Settings::flag("SoftQCD:doubleDiffractive")) {
    sigmaPtr = new Sigma0AB2XX;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  
  // Set up requested objects for hard QCD processes.
  bool hardQCD = Settings::flag("HardQCD:all");
  if (hardQCD || Settings::flag("HardQCD:gg2gg")) {
    sigmaPtr = new Sigma2gg2gg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hardQCD || Settings::flag("HardQCD:gg2qqbar")) {
    sigmaPtr = new Sigma2gg2qqbar;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hardQCD || Settings::flag("HardQCD:qg2qg")) {
    sigmaPtr = new Sigma2qg2qg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hardQCD || Settings::flag("HardQCD:qq2qq")) {
    sigmaPtr = new Sigma2qq2qqDiff;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qq2qqSame;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2qqbarSame;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hardQCD || Settings::flag("HardQCD:qqbar2qqbarNew")) {
    sigmaPtr = new Sigma2qqbar2qqbarNew;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hardQCD || Settings::flag("HardQCD:qqbar2gg")) {
    sigmaPtr = new Sigma2qqbar2gg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  
  // Set up requested objects for c cbar and b bbar, also hard QCD.
  if (hardQCD || Settings::flag("HardQCD:gg2ccbar")) {
    sigmaPtr = new Sigma2gg2QQbar(4, 121);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hardQCD || Settings::flag("HardQCD:qqbar2ccbar")) {
    sigmaPtr = new Sigma2qqbar2QQbar(4, 122);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hardQCD || Settings::flag("HardQCD:gg2bbbar")) {
    sigmaPtr = new Sigma2gg2QQbar(5, 123);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (hardQCD || Settings::flag("HardQCD:qqbar2bbbar")) {
    sigmaPtr = new Sigma2qqbar2QQbar(5, 124);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // Set up requested objects for prompt photon processes.
  bool promptPhotons = Settings::flag("PromptPhoton:all");
  if (promptPhotons
    || Settings::flag("PromptPhoton:qg2qgamma")) {
    sigmaPtr = new Sigma2qg2qgamma;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (promptPhotons 
    || Settings::flag("PromptPhoton:qqbar2ggamma")) {
    sigmaPtr = new Sigma2qqbar2ggamma;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (promptPhotons
    || Settings::flag("PromptPhoton:gg2ggamma")) {
    sigmaPtr = new Sigma2gg2ggamma;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (promptPhotons 
    || Settings::flag("PromptPhoton:qqbar2gammagamma")) {
    sigmaPtr = new Sigma2qqbar2gammagamma;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (promptPhotons
    || Settings::flag("PromptPhoton:gg2gammagamma")) {
    sigmaPtr = new Sigma2gg2gammagamma;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // Set up requested objects for weak gauge boson t-channel exchange.
  bool weakBosonExchanges = Settings::flag("WeakBosonExchange:all");
  if (weakBosonExchanges
    || Settings::flag("WeakBosonExchange:ff2ff(t:gmZ)")) {
    sigmaPtr = new Sigma2ff2fftgmZ;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakBosonExchanges
    || Settings::flag("WeakBosonExchange:ff2ff(t:W)")) {
    sigmaPtr = new Sigma2ff2fftW;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // Set up requested objects for weak gauge boson processes.
  bool weakSingleBosons = Settings::flag("WeakSingleBoson:all");
  if (weakSingleBosons
    || Settings::flag("WeakSingleBoson:ffbar2gmZ")) {
    sigmaPtr = new Sigma1ffbar2gmZ;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakSingleBosons
    || Settings::flag("WeakSingleBoson:ffbar2W")) {
    sigmaPtr = new Sigma1ffbar2W;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // Set up requested object for s-channel gamma exchange.
  // Subset of gamma*/Z0 above, intended for Multiple interactions.
  if (Settings::flag("WeakSingleBoson:ffbar2ffbar(s:gm)")) {
    sigmaPtr = new Sigma2ffbar2ffbarsgm;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
   
  // Set up requested objects for weak gauge boson pair processes.
  bool weakDoubleBosons = Settings::flag("WeakDoubleBoson:all");
  if (weakDoubleBosons
    || Settings::flag("WeakDoubleBoson:ffbar2gmZgmZ")) {
    sigmaPtr = new Sigma2ffbar2gmZgmZ;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakDoubleBosons
    || Settings::flag("WeakDoubleBoson:ffbar2ZW")) {
    sigmaPtr = new Sigma2ffbar2ZW;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakDoubleBosons
    || Settings::flag("WeakDoubleBoson:ffbar2WW")) {
    sigmaPtr = new Sigma2ffbar2WW;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  
  // Set up requested objects for weak gauge boson + parton processes.
  bool weakBosonAndPartons = Settings::flag("WeakBosonAndParton:all");
  if (weakBosonAndPartons
    || Settings::flag("WeakBosonAndParton:qqbar2gmZg")) {
    sigmaPtr = new Sigma2qqbar2gmZg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakBosonAndPartons
    || Settings::flag("WeakBosonAndParton:qg2gmZq")) {
    sigmaPtr = new Sigma2qg2gmZq;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakBosonAndPartons
    || Settings::flag("WeakBosonAndParton:ffbar2gmZgm")) {
    sigmaPtr = new Sigma2ffbar2gmZgm;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakBosonAndPartons
    || Settings::flag("WeakBosonAndParton:fgm2gmZf")) {
    sigmaPtr = new Sigma2fgm2gmZf;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakBosonAndPartons
    || Settings::flag("WeakBosonAndParton:qqbar2Wg")) {
    sigmaPtr = new Sigma2qqbar2Wg;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakBosonAndPartons
    || Settings::flag("WeakBosonAndParton:qg2Wq")) {
    sigmaPtr = new Sigma2qg2Wq;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakBosonAndPartons
    || Settings::flag("WeakBosonAndParton:ffbar2Wgm")) {
    sigmaPtr = new Sigma2ffbar2Wgm;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (weakBosonAndPartons
    || Settings::flag("WeakBosonAndParton:fgm2Wf")) {
    sigmaPtr = new Sigma2fgm2Wf;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  
  // Set up requested objects for charmonium production
  bool charmoniums = Settings::flag("Charmonium:all");
  if (charmoniums || Settings::flag("Charmonium:gg2QQbar[3S1(1)]g")) {
    sigmaPtr = new Sigma2gg2QQbar3S11g(4, 401);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || Settings::flag("Charmonium:gg2QQbar[3P0(1)]g")) {
    sigmaPtr = new Sigma2gg2QQbar3PJ1g(4, 0, 402);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || Settings::flag("Charmonium:gg2QQbar[3P1(1)]g")) {
    sigmaPtr = new Sigma2gg2QQbar3PJ1g(4, 1, 403);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || Settings::flag("Charmonium:gg2QQbar[3P2(1)]g")) {
    sigmaPtr = new Sigma2gg2QQbar3PJ1g(4, 2, 404);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || Settings::flag("Charmonium:qg2QQbar[3P0(1)]q")) {
    sigmaPtr = new Sigma2qg2QQbar3PJ1q(4, 0, 405);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || Settings::flag("Charmonium:qg2QQbar[3P1(1)]q")) {
    sigmaPtr = new Sigma2qg2QQbar3PJ1q(4, 1, 406);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || Settings::flag("Charmonium:qg2QQbar[3P2(1)]q")) {
    sigmaPtr = new Sigma2qg2QQbar3PJ1q(4, 2, 407);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || Settings::flag("Charmonium:qqbar2QQbar[3P0(1)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(4, 0, 408);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || Settings::flag("Charmonium:qqbar2QQbar[3P1(1)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(4, 1, 409);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || Settings::flag("Charmonium:qqbar2QQbar[3P2(1)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(4, 2, 410);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || Settings::flag("Charmonium:gg2QQbar[3S1(8)]g")) {
    sigmaPtr = new Sigma2gg2QQbarX8g(4, 0, 411);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || Settings::flag("Charmonium:gg2QQbar[1S0(8)]g")) {
    sigmaPtr = new Sigma2gg2QQbarX8g(4, 1, 412);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || Settings::flag("Charmonium:gg2QQbar[3PJ(8)]g")) {
    sigmaPtr = new Sigma2gg2QQbarX8g(4, 2, 413);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || Settings::flag("Charmonium:qg2QQbar[3S1(8)]q")) {
    sigmaPtr = new Sigma2qg2QQbarX8q(4, 0, 414);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || Settings::flag("Charmonium:qg2QQbar[1S0(8)]q")) {
    sigmaPtr = new Sigma2qg2QQbarX8q(4, 1, 415);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || Settings::flag("Charmonium:qg2QQbar[3PJ(8)]q")) {
    sigmaPtr = new Sigma2qg2QQbarX8q(4, 2, 416);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || Settings::flag("Charmonium:qqbar2QQbar[3S1(8)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbarX8g(4, 0, 417);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || Settings::flag("Charmonium:qqbar2QQbar[1S0(8)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbarX8g(4, 1, 418);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (charmoniums || Settings::flag("Charmonium:qqbar2QQbar[3PJ(8)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbarX8g(4, 2, 419);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
    
  // Set up requested objects for bottomonium production
  bool bottomoniums = Settings::flag("Bottomonium:all");
  if (bottomoniums || Settings::flag("Bottomonium:gg2QQbar[3S1(1)]g")) {
    sigmaPtr = new Sigma2gg2QQbar3S11g(5, 501);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || Settings::flag("Bottomonium:gg2QQbar[3P0(1)]g")) {
    sigmaPtr = new Sigma2gg2QQbar3PJ1g(5, 0, 502);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || Settings::flag("Bottomonium:gg2QQbar[3P1(1)]g")) {
    sigmaPtr = new Sigma2gg2QQbar3PJ1g(5, 1, 503);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || Settings::flag("Bottomonium:gg2QQbar[3P2(1)]g")) {
    sigmaPtr = new Sigma2gg2QQbar3PJ1g(5, 2, 504);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qg2QQbar[3P0(1)]q")) {
    sigmaPtr = new Sigma2qg2QQbar3PJ1q(5, 0, 505);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qg2QQbar[3P1(1)]q")) {
    sigmaPtr = new Sigma2qg2QQbar3PJ1q(5, 1, 506);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qg2QQbar[3P2(1)]q")) {
    sigmaPtr = new Sigma2qg2QQbar3PJ1q(5, 2, 507);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qqbar2QQbar[3P0(1)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(5, 0, 508);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qqbar2QQbar[3P1(1)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(5, 1, 509);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qqbar2QQbar[3P2(1)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(5, 2, 510);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || Settings::flag("Bottomonium:gg2QQbar[3S1(8)]g")) {
    sigmaPtr = new Sigma2gg2QQbarX8g(5, 0, 511);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || Settings::flag("Bottomonium:gg2QQbar[1S0(8)]g")) {
    sigmaPtr = new Sigma2gg2QQbarX8g(5, 1, 512);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || Settings::flag("Bottomonium:gg2QQbar[3PJ(8)]g")) {
    sigmaPtr = new Sigma2gg2QQbarX8g(5, 2, 513);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qg2QQbar[3S1(8)]q")) {
    sigmaPtr = new Sigma2qg2QQbarX8q(5, 0, 514);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qg2QQbar[1S0(8)]q")) {
    sigmaPtr = new Sigma2qg2QQbarX8q(5, 1, 515);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qg2QQbar[3PJ(8)]q")) {
    sigmaPtr = new Sigma2qg2QQbarX8q(5, 2, 516);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qqbar2QQbar[3S1(8)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbarX8g(5, 0, 517);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qqbar2QQbar[1S0(8)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbarX8g(5, 1, 518);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qqbar2QQbar[3PJ(8)]g")) {
    sigmaPtr = new Sigma2qqbar2QQbarX8g(5, 2, 519);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  
  // Set up requested objects for top production
  bool tops = Settings::flag("Top:all");
  if (tops || Settings::flag("Top:gg2ttbar")) {
    sigmaPtr = new Sigma2gg2QQbar(6, 601);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (tops || Settings::flag("Top:qqbar2ttbar")) {
    sigmaPtr = new Sigma2qqbar2QQbar(6, 602);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (tops || Settings::flag("Top:qq2tq(t:W)")) {
    sigmaPtr = new Sigma2qq2QqtW(6, 603, "q q -> t q (t-channel W+-)");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (tops || Settings::flag("Top:ffbar2ttbar(s:gmZ)")) {
    sigmaPtr = new Sigma2ffbar2FFbarsgmZ(6, 604, 
      "f fbar -> t tbar (s-channel gamma*/Z0)");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (tops || Settings::flag("Top:ffbar2tqbar(s:W)")) {
    sigmaPtr = new Sigma2ffbar2FfbarsW(6, 605, 
      "f fbar -> t qbar (s-channel W+-)");
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  
  // Set up requested objects for Standard-Model Higgs production
  bool SMHiggses = Settings::flag("SMHiggs:all");
  if (SMHiggses || Settings::flag("SMHiggs:ffbar2H")) {
    sigmaPtr = new Sigma1ffbar2H;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (SMHiggses || Settings::flag("SMHiggs:gg2H")) {
    sigmaPtr = new Sigma1gg2H;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (SMHiggses || Settings::flag("SMHiggs:gmgm2H")) {
    sigmaPtr = new Sigma1gmgm2H;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (SMHiggses || Settings::flag("SMHiggs:ffbar2HZ")) {
    sigmaPtr = new Sigma2ffbar2HZ;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  if (SMHiggses || Settings::flag("SMHiggs:ffbar2HW")) {
    sigmaPtr = new Sigma2ffbar2HW;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 
  // This process not included in "all" so as not to doublecount.
  if (Settings::flag("SMHiggs:qg2Hq")) {
    sigmaPtr = new Sigma2qg2Hq;
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 


  // Set up requested objects for neutralino pair processes.
  bool SUSYs = Settings::flag("SUSY:all");
  if (SUSYs || Settings::flag("SUSY:qqbar2chi0chi0")) {
    sigmaPtr = new Sigma2qqbar2chi0chi0(1, 1, 1001);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2chi0chi0(1, 2, 1002); 
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2chi0chi0(1, 3, 1003); 
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2chi0chi0(1, 4, 1004);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2chi0chi0(2, 2, 1005);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2chi0chi0(2, 3, 1006);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2chi0chi0(2, 4, 1007);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2chi0chi0(3, 3, 1008); 
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2chi0chi0(3, 4, 1009);
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2chi0chi0(4, 4, 1010); 
    containerPtrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // Done. 
  return true;

}

//*********

// Routine to initialize list of second hard processes.

bool SetupContainers::init2(vector<ProcessContainer*>& container2Ptrs) {

  // Reset process list, if filled in previous subrun.
  if (container2Ptrs.size() > 0) {
    for (int i = 0; i < int(container2Ptrs.size()); ++i) 
      delete container2Ptrs[i];
    container2Ptrs.clear(); 
  }
  SigmaProcess* sigmaPtr;

  // Two hard QCD jets.
  if (Settings::flag("SecondHard:TwoJets")) {
    sigmaPtr = new Sigma2gg2gg;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2qqbar;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qg2qg;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qq2qqDiff;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qq2qqSame;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2qqbarSame;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2qqbarNew;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2gg;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2QQbar(4, 121);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbar(4, 122);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2QQbar(5, 123);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbar(5, 124);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // A prompt photon and a hard jet.
  if (Settings::flag("SecondHard:PhotonAndJet")) {
    sigmaPtr = new Sigma2qg2qgamma;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2ggamma;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2ggamma;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // Two prompt photons.
  if (Settings::flag("SecondHard:TwoPhotons")) {
    sigmaPtr = new Sigma2qqbar2gammagamma;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2gg2gammagamma;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // A single gamma*/Z0.
  if (Settings::flag("SecondHard:SingleGmZ")) {
    sigmaPtr = new Sigma1ffbar2gmZ;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // A single W+-.
  if (Settings::flag("SecondHard:SingleW")) {
    sigmaPtr = new Sigma1ffbar2W;
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  } 

  // Two b jets - already part of TwoJets sample above.
  if (Settings::flag("SecondHard:TwoBJets")) {
    sigmaPtr = new Sigma2gg2QQbar(5, 123);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
    sigmaPtr = new Sigma2qqbar2QQbar(5, 124);
    container2Ptrs.push_back( new ProcessContainer(sigmaPtr) );
  }

  // Done. 
  return true;

}

//**************************************************************************

} // end namespace Pythia8
