// Function definitions (not found in the header) for the 
// ProcessContainer and SetupContainers classes.
// Copyright C 2006 Torbjorn Sjostrand

#include "ProcessContainer.h"

namespace Pythia8 {

//**************************************************************************

// ProcessContainer class.
// Information allowing the generation of a specific process.

//*********

// Initialize 

bool ProcessContainer::init(Info* infoPtrIn, BeamParticle* beamAPtr, 
  BeamParticle* beamBPtr) {

  // Store input pointers for future use.
  infoPtr = infoPtrIn;

  // Extract info about current process from SigmaProcess object.
  hasSigmaTot = sigmaProcessPtr->hasSigmaTot();  
  isMinBias = sigmaProcessPtr->isMinBias();
  isResolved = sigmaProcessPtr->isResolved();
  isDiffA = sigmaProcessPtr->isDiffA();
  isDiffB = sigmaProcessPtr->isDiffB();
  int nFinal = sigmaProcessPtr->nFinal();

  // Pick and create phase space generator. Send pointers where required.
  if (isMinBias) phaseSpacePtr = new PhaseSpace2to2minbias();
  else if (!isResolved) phaseSpacePtr = new PhaseSpace2to2eldiff( isDiffA, 
    isDiffB);
  else if (nFinal == 1) phaseSpacePtr = new PhaseSpace2to1tauy();
  else phaseSpacePtr = new PhaseSpace2to2tauyz();

  // For total cross section processes need cross section and beam info.
  if (hasSigmaTot) {
    int idA = beamAPtr->id();
    int idB = beamBPtr->id();
    double mA = beamAPtr->m();
    double mB = beamBPtr->m();
    sigmaProcessPtr->setSigmaTotalPtr(sigmaTotPtr, idA, idB, mA, mB);
    phaseSpacePtr->setSigmaTotalPtr(sigmaTotPtr, idA, idB, mA, mB);
  }

  // Store common info for PhaseSpace objects.
  double eCM = infoPtr->eCM();
  phaseSpacePtr->initInfo( sigmaProcessPtr, eCM);

  // Reset cross section statistics.
  nTry = 0;
  nAcc = 0;
  sigmaMx = 0.;
  sigmaSum = 0.;
  sigma2Sum = 0.;

  // Set up incoming flux (= product of parton densities).
  sigmaProcessPtr->initFlux( beamAPtr->pdfPtr(), beamBPtr->pdfPtr());

  // Find maximum of differential cross section * phasespace.
  phaseSpacePtr->setupSampling();
  sigmaMx = phaseSpacePtr->sigmaMax();

  // Done.
  return true;
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

  // Update statistics.
  sigmaSum += sigmaNow;
  sigma2Sum += pow2(sigmaNow);

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
  int id1 = sigmaProcessPtr->id(1);
  int id2 = sigmaProcessPtr->id(2);
  double pdf1 = sigmaProcessPtr->pdf1();
  double pdf2 = sigmaProcessPtr->pdf2();
  double Q2Fac = sigmaProcessPtr->Q2Fac();
  double alphaEM = sigmaProcessPtr->alphaEM();
  double alphaS = sigmaProcessPtr->alphaS();
  double Q2Ren = sigmaProcessPtr->Q2Ren();
  double x1 = phaseSpacePtr->x1();
  double x2 = phaseSpacePtr->x2();
  double sHat = phaseSpacePtr->sHat();
  double tHat = phaseSpacePtr->tHat();
  double uHat = phaseSpacePtr->uHat();
  double pTHat = phaseSpacePtr->pTHat();
  double m3 = phaseSpacePtr->m(3);
  double m4 = phaseSpacePtr->m(4);
  double theta = phaseSpacePtr->thetaHat();
  double phi = phaseSpacePtr->phiHat();
  infoPtr->setPDFalpha( id1, id2, pdf1, pdf2, Q2Fac, alphaEM, alphaS, Q2Ren);
  infoPtr->setKin( x1, x2, sHat, tHat, uHat, pTHat, m3, m4, theta, phi);

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
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbar(4, 121, "g g -> c cbar");
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (hardQCD || Settings::flag("HardQCD:qqbar2ccbar")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2QQbar(4, 122, "q qbar -> c cbar");
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (hardQCD || Settings::flag("HardQCD:gg2bbbar")) {
    SigmaProcess* sigmaPtr = new Sigma2gg2QQbar(5, 123, "g g -> b bbar");
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (hardQCD || Settings::flag("HardQCD:qqbar2bbbar")) {
    SigmaProcess* sigmaPtr = new Sigma2qqbar2QQbar(5, 124, "q qbar -> b bbar");
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

  // Set up requested objects for weak gauge boson processes.
  bool weakSingleBosons = Settings::flag("WeakSingleBoson:all");
  if (weakSingleBosons
    || Settings::flag("WeakSingleBoson:ffbar2W")) {
    SigmaProcess* sigmaPtr = new Sigma1ffbar2W;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  


  // Done. 
  return true;
}

//**************************************************************************

} // end namespace Pythia8
