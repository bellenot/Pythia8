// Function definitions (not found in the header) for the 
// ProcessContainer and SetupContainers classes.
// Copyright C 2007 Torbjorn Sjostrand

#include "ProcessContainer.h"

// Internal headers for special processes.
#include "SigmaQCD.h"
#include "SigmaEW.h"
#include "SigmaOnia.h"
#include "SigmaSUSY.h"

namespace Pythia8 {

//**************************************************************************

// ProcessContainer class.
// Information allowing the generation of a specific process.

//*********
 
// Definitions of static variables and functions.

// Pointer to the information object.
Info* ProcessContainer::infoPtr;

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
  nAcc      = 0;
  sigmaMx   = 0.;
  sigmaSum  = 0.;
  sigma2Sum = 0.;
  sigmaNeg  = 0.;

  // Initialize process. Remove empty inFlux channels and optionally list.
  sigmaProcessPtr->initProc();
  sigmaProcessPtr->initFlux();
  sigmaProcessPtr->checkChannels();

  // Find maximum of differential cross section * phasespace.
  bool physical = phaseSpacePtr->setupSampling();
  sigmaMx = phaseSpacePtr->sigmaMax();
  double sigmaHalfWay = sigmaMx;

  // Check maximum by a few events, and extrapolate a further increase.
  if (physical) {
    for (int sample = 0; sample < NSAMPLE; ++sample) 
    while (!phaseSpacePtr->trialKin()) { 
      if (sample == NSAMPLE/2) sigmaHalfWay = phaseSpacePtr->sigmaMax();
    }   
    sigmaMx = pow2(phaseSpacePtr->sigmaMax()) / sigmaHalfWay;
    phaseSpacePtr->setSigmaMax(sigmaMx);
  }

  // Done.
  return physical;
}

//*********

// Generate a trial event; accepted or not.
 
bool ProcessContainer::trialProcess() { 

  // Update number of tries.
  if (sigmaMx == 0.) return false;
  ++nTry;

  // Generate a trial phase space point, with cross section.
  if (!phaseSpacePtr->trialKin()) return false;
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

  // Accept or reject trial point.
  bool accept = (sigmaNow > Rndm::flat() * sigmaMx);  
  if (accept) ++nAcc;
  return accept;

}

//*********
  
// Give the hard subprocess.

bool ProcessContainer::constructProcess( Event& process) { 

  // Construct flavour and colours for accepted event.
  if (isResolved && !isMinBias) sigmaProcessPtr->pickInState();
  sigmaProcessPtr->setIdColAcol();

  // Construct kinematics from accepted phase space point.
  if (!phaseSpacePtr->finalKin()) return false;

  // Basic info on process.
  infoPtr->setType( name(), code(), nFinal(), isMinBias, isResolved, 
    isDiffA, isDiffB);

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
  infoPtr->setPDFalpha( id1, id2, pdf1, pdf2, Q2Fac, alphaEM, alphaS, Q2Ren);
  infoPtr->setKin( x1, x2, sHat, tHat, uHat, pTHat, m3, m4, theta, phi);
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

  // Set up requested objects for soft QCD processes.
  bool softQCD = Settings::flag("SoftQCD:all");
  if (softQCD || Settings::flag("SoftQCD:minBias")) {
    SigmaProcess* sigmaPtr = new Sigma0minBias;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (softQCD || Settings::flag("SoftQCD:elastic")) {
    SigmaProcess* sigmaPtr = new Sigma0AB2AB;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (softQCD || Settings::flag("SoftQCD:singleDiffractive")) {
    SigmaProcess* sigmaPtr = new Sigma0AB2XB;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
    sigmaPtr = new Sigma0AB2AX;
    containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (softQCD || Settings::flag("SoftQCD:doubleDiffractive")) {
    SigmaProcess* sigmaPtr = new Sigma0AB2XX;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  
  // Set up requested objects for hard QCD processes.
  bool hardQCD = Settings::flag("HardQCD:all");
  if (hardQCD || Settings::flag("HardQCD:gg2gg")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2gg;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (hardQCD || Settings::flag("HardQCD:gg2qqbar")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2qqbar;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (hardQCD || Settings::flag("HardQCD:qg2qg")) {
    SigmaProcess* sigmaPtr = new Sigma2qg2qg;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (hardQCD || Settings::flag("HardQCD:qq2qq")) {
    SigmaProcess* sigmaPtr = new Sigma2qq2qqDiff;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
    sigmaPtr = new Sigma2qq2qqSame;
    containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
    sigmaPtr = new Sigma2qqbar2qqbarSame;
    containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (hardQCD || Settings::flag("HardQCD:qqbar2qqbarNew")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2qqbarNew;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (hardQCD || Settings::flag("HardQCD:qqbar2gg")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2gg;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  
  // Set up requested objects for c cbar and b bbar, also hard QCD.
  if (hardQCD || Settings::flag("HardQCD:gg2ccbar")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbar(4, 121);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (hardQCD || Settings::flag("HardQCD:qqbar2ccbar")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2QQbar(4, 122);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (hardQCD || Settings::flag("HardQCD:gg2bbbar")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbar(5, 123);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (hardQCD || Settings::flag("HardQCD:qqbar2bbbar")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2QQbar(5, 124);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 

  // Set up requested objects for prompt photon processes.
  bool promptPhotons = Settings::flag("PromptPhoton:all");
  if (promptPhotons
    || Settings::flag("PromptPhoton:qg2qgamma")) {
    SigmaProcess* sigmaPtr = new Sigma2qg2qgamma;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (promptPhotons 
    || Settings::flag("PromptPhoton:qqbar2ggamma")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2ggamma;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (promptPhotons
    || Settings::flag("PromptPhoton:gg2ggamma")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2ggamma;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (promptPhotons 
    || Settings::flag("PromptPhoton:qqbar2gammagamma")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2gammagamma;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (promptPhotons
    || Settings::flag("PromptPhoton:gg2gammagamma")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2gammagamma;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 

  // Set up requested objects for weak gauge boson t-channel exchange.
  bool weakBosonExchanges = Settings::flag("WeakBosonExchange:all");
  if (weakBosonExchanges
    || Settings::flag("WeakBosonExchange:ff2ff(t:gmZ)")) {
    SigmaProcess* sigmaPtr = new Sigma2ff2fftgmZ;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (weakBosonExchanges
    || Settings::flag("WeakBosonExchange:ff2ff(t:W)")) {
    SigmaProcess* sigmaPtr = new Sigma2ff2fftW;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 

  // Set up requested objects for weak gauge boson processes.
  bool weakSingleBosons = Settings::flag("WeakSingleBoson:all");
  if (weakSingleBosons
    || Settings::flag("WeakSingleBoson:ffbar2gmZ")) {
    SigmaProcess* sigmaPtr = new Sigma1ffbar2gmZ;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (weakSingleBosons
    || Settings::flag("WeakSingleBoson:ffbar2W")) {
    SigmaProcess* sigmaPtr = new Sigma1ffbar2W;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 

  // Set up requested object for s-channel gamma exchange.
  // Subset of gamma*/Z0 above, intended for Multiple interactions.
  if (Settings::flag("WeakSingleBoson:ffbar2ffbar(s:gm)")) {
    SigmaProcess* sigmaPtr = new Sigma2ffbar2ffbarsgm;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
   
  // Set up requested objects for weak gauge boson pair processes.
  bool weakDoubleBosons = Settings::flag("WeakDoubleBoson:all");
  if (weakDoubleBosons
    || Settings::flag("WeakDoubleBoson:ffbar2ZW")) {
    SigmaProcess* sigmaPtr = new Sigma2ffbar2ZW;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (weakDoubleBosons
    || Settings::flag("WeakDoubleBoson:ffbar2WW")) {
    SigmaProcess* sigmaPtr = new Sigma2ffbar2WW;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  
  // Set up requested objects for weak gauge boson + parton processes.
  bool weakBosonAndPartons = Settings::flag("WeakBosonAndParton:all");
  if (weakBosonAndPartons
    || Settings::flag("WeakBosonAndParton:qqbar2Wg")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2Wg;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (weakBosonAndPartons
    || Settings::flag("WeakBosonAndParton:qg2Wq")) {
    SigmaProcess* sigmaPtr = new Sigma2qg2Wq;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (weakBosonAndPartons
    || Settings::flag("WeakBosonAndParton:ffbar2Wgm")) {
    SigmaProcess* sigmaPtr = new Sigma2ffbar2Wgm;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  
  // Set up requested objects for charmonium production
  bool charmoniums = Settings::flag("Charmonium:all");
  if (charmoniums || Settings::flag("Charmonium:gg2QQbar[3S1(1)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbar3S11g(4, 401);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (charmoniums || Settings::flag("Charmonium:gg2QQbar[3P0(1)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbar3PJ1g(4, 0, 402);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (charmoniums || Settings::flag("Charmonium:gg2QQbar[3P1(1)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbar3PJ1g(4, 1, 403);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (charmoniums || Settings::flag("Charmonium:gg2QQbar[3P2(1)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbar3PJ1g(4, 2, 404);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (charmoniums || Settings::flag("Charmonium:qg2QQbar[3P0(1)]q")) {
    SigmaProcess* sigmaPtr = new Sigma2qg2QQbar3PJ1q(4, 0, 405);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (charmoniums || Settings::flag("Charmonium:qg2QQbar[3P1(1)]q")) {
    SigmaProcess* sigmaPtr = new Sigma2qg2QQbar3PJ1q(4, 1, 406);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (charmoniums || Settings::flag("Charmonium:qg2QQbar[3P2(1)]q")) {
    SigmaProcess* sigmaPtr = new Sigma2qg2QQbar3PJ1q(4, 2, 407);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (charmoniums || Settings::flag("Charmonium:qqbar2QQbar[3P0(1)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(4, 0, 408);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (charmoniums || Settings::flag("Charmonium:qqbar2QQbar[3P1(1)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(4, 1, 409);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (charmoniums || Settings::flag("Charmonium:qqbar2QQbar[3P2(1)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(4, 2, 410);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (charmoniums || Settings::flag("Charmonium:gg2QQbar[3S1(8)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbarX8g(4, 0, 411);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (charmoniums || Settings::flag("Charmonium:gg2QQbar[1S0(8)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbarX8g(4, 1, 412);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (charmoniums || Settings::flag("Charmonium:gg2QQbar[3PJ(8)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbarX8g(4, 2, 413);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (charmoniums || Settings::flag("Charmonium:qg2QQbar[3S1(8)]q")) {
    SigmaProcess* sigmaPtr = new Sigma2qg2QQbarX8q(4, 0, 414);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (charmoniums || Settings::flag("Charmonium:qg2QQbar[1S0(8)]q")) {
    SigmaProcess* sigmaPtr = new Sigma2qg2QQbarX8q(4, 1, 415);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (charmoniums || Settings::flag("Charmonium:qg2QQbar[3PJ(8)]q")) {
    SigmaProcess* sigmaPtr = new Sigma2qg2QQbarX8q(4, 2, 416);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (charmoniums || Settings::flag("Charmonium:qqbar2QQbar[3S1(8)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2QQbarX8g(4, 0, 417);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (charmoniums || Settings::flag("Charmonium:qqbar2QQbar[1S0(8)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2QQbarX8g(4, 1, 418);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (charmoniums || Settings::flag("Charmonium:qqbar2QQbar[3PJ(8)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2QQbarX8g(4, 2, 419);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
    
  // Set up requested objects for bottomonium production
  bool bottomoniums = Settings::flag("Bottomonium:all");
  if (bottomoniums || Settings::flag("Bottomonium:gg2QQbar[3S1(1)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbar3S11g(5, 501);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (bottomoniums || Settings::flag("Bottomonium:gg2QQbar[3P0(1)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbar3PJ1g(5, 0, 502);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (bottomoniums || Settings::flag("Bottomonium:gg2QQbar[3P1(1)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbar3PJ1g(5, 1, 503);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (bottomoniums || Settings::flag("Bottomonium:gg2QQbar[3P2(1)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbar3PJ1g(5, 2, 504);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qg2QQbar[3P0(1)]q")) {
    SigmaProcess* sigmaPtr = new Sigma2qg2QQbar3PJ1q(5, 0, 505);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qg2QQbar[3P1(1)]q")) {
    SigmaProcess* sigmaPtr = new Sigma2qg2QQbar3PJ1q(5, 1, 506);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qg2QQbar[3P2(1)]q")) {
    SigmaProcess* sigmaPtr = new Sigma2qg2QQbar3PJ1q(5, 2, 507);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qqbar2QQbar[3P0(1)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(5, 0, 508);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qqbar2QQbar[3P1(1)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(5, 1, 509);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qqbar2QQbar[3P2(1)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2QQbar3PJ1g(5, 2, 510);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (bottomoniums || Settings::flag("Bottomonium:gg2QQbar[3S1(8)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbarX8g(5, 0, 511);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (bottomoniums || Settings::flag("Bottomonium:gg2QQbar[1S0(8)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbarX8g(5, 1, 512);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (bottomoniums || Settings::flag("Bottomonium:gg2QQbar[3PJ(8)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbarX8g(5, 2, 513);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qg2QQbar[3S1(8)]q")) {
    SigmaProcess* sigmaPtr = new Sigma2qg2QQbarX8q(5, 0, 514);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qg2QQbar[1S0(8)]q")) {
    SigmaProcess* sigmaPtr = new Sigma2qg2QQbarX8q(5, 1, 515);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qg2QQbar[3PJ(8)]q")) {
    SigmaProcess* sigmaPtr = new Sigma2qg2QQbarX8q(5, 2, 516);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qqbar2QQbar[3S1(8)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2QQbarX8g(5, 0, 517);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qqbar2QQbar[1S0(8)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2QQbarX8g(5, 1, 518);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (bottomoniums || Settings::flag("Bottomonium:qqbar2QQbar[3PJ(8)]g")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2QQbarX8g(5, 2, 519);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  
  // Set up requested objects for top production
  bool tops = Settings::flag("Top:all");
  if (tops || Settings::flag("Top:gg2ttbar")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbar(6, 601);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (tops || Settings::flag("Top:qqbar2ttbar")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2QQbar(6, 602);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (tops || Settings::flag("Top:qq2tq(t:W)")) {
    SigmaProcess* sigmaPtr = new Sigma2qq2QqtW(6, 603, 
      "q q -> t q (t-channel W+-)");
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 

  // Set up requested objects for neutralino pair processes.
  bool SUSYs = Settings::flag("SUSY:all");
  if (SUSYs || Settings::flag("SUSY:qqbar2chi0chi0")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2chi0chi0(1, 1, 1001);
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
    sigmaPtr = new Sigma2qqbar2chi0chi0(1, 2, 1002); 
    containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
    sigmaPtr = new Sigma2qqbar2chi0chi0(1, 3, 1003); 
    containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
    sigmaPtr = new Sigma2qqbar2chi0chi0(1, 4, 1004);
    containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
    sigmaPtr = new Sigma2qqbar2chi0chi0(2, 2, 1005);
    containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
    sigmaPtr = new Sigma2qqbar2chi0chi0(2, 3, 1006);
    containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
    sigmaPtr = new Sigma2qqbar2chi0chi0(2, 4, 1007);
    containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
    sigmaPtr = new Sigma2qqbar2chi0chi0(3, 3, 1008); 
    containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
    sigmaPtr = new Sigma2qqbar2chi0chi0(3, 4, 1009);
    containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
    sigmaPtr = new Sigma2qqbar2chi0chi0(4, 4, 1010); 
    containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 

  // Done. 
  return true;
}

//**************************************************************************

} // end namespace Pythia8
