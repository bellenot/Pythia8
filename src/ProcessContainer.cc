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

bool ProcessContainer::init(Info& info, PDF* pdfAPtr, PDF* pdfBPtr) {

  // Extract info about current process from SigmaHat object.  
  isResolved = sigmaHatPtr->isResolved();
  isDiffA = sigmaHatPtr->isDiffA();
  isDiffB = sigmaHatPtr->isDiffB();

  // Pick and create phase space generator.
  if (!isResolved) phaseSpacePtr = new PhaseSpace2to2eldiff( isDiffA, 
    isDiffB);
  else phaseSpacePtr = new PhaseSpace2to2tauyz();

  // Send pointers to SigmaTotal object where required.
  if (!isResolved) {
    sigmaHatPtr->setSigmaTotalPtr(sigmaTotPtr);
    phaseSpacePtr->setSigmaTotalPtr(sigmaTotPtr);
  }

  // Store common info for SigmaHat and PhaseSpace objects.
  sigmaHatPtr->initInfo( info);
  phaseSpacePtr->initInfo( sigmaHatPtr, info);

  // Reset cross section statistics.
  nTry = 0;
  nAcc = 0;
  sigmaMx = 0.;
  sigmaSum = 0.;
  sigma2Sum = 0.;

  // Set up incoming flux (= product of parton densities).
  sigmaHatPtr->initFlux( pdfAPtr, pdfBPtr);

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

bool ProcessContainer::constructProcess( Info& info, Event& process) { 

  // Construct flavour and colours for accepted event.
  if (isResolved) sigmaHatPtr->pickInState();
  sigmaHatPtr->setIdColAcol();

  // Construct kinematics from accepted phase space point.
  if (!phaseSpacePtr->finalKin()) return false;

  // Store info on process.
  int id1 = sigmaHatPtr->id(1);
  int id2 = sigmaHatPtr->id(2);
  double pdf1 = sigmaHatPtr->pdf1();
  double pdf2 = sigmaHatPtr->pdf2();
  double Q2pdf = sigmaHatPtr->Q2pdf();
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
  info.setType( name(), code(), nFinal(), isResolved, isDiffA, isDiffB);
  info.setPDF( id1, id2, pdf1, pdf2, Q2pdf);
  info.setKin( x1, x2, sHat, tHat, uHat, pTHat, m3, m4, theta, phi);

  // Let hard process record begin with the event as a whole and
  // the two incoming beam particles.  
  process.append( 90, -11, 0, 0, 1, 2, 0, 0, 
    Vec4(0., 0., 0., info.eCM()), info.eCM(), 0. ); 
  process.append( info.idA(), -12, 0, 0, 3, 0, 0, 0, 
    Vec4(0., 0., info.pzA(), info.eA()), info.mA(), 0. ); 
  process.append( info.idB(), -12, 0, 0, 4, 0, 0, 0, 
    Vec4(0., 0., info.pzB(), info.eB()), info.mB(), 0. ); 

  // Insert the subprocess partons.
  int colOffset = process.lastColTag();
  for (int i = 1; i <= 2 + sigmaHatPtr->nFinal(); ++i) { 
    int id = sigmaHatPtr->id(i);
    int status = (i <= 2) ? -21 : 23;
    int mother1 = (i <= 2) ? i : 3;
    int mother2 = (i <= 2) ? 0 : 4;
    int daughter1 = (i <= 2) ? 5 : 0;
    int daughter2 = (i <= 2) ? 4 + sigmaHatPtr->nFinal() : 0;
    int col = sigmaHatPtr->col(i);
    if (col > 0) col += colOffset;
    int acol = sigmaHatPtr->acol(i);
    if (acol > 0) acol += colOffset;
    process.append( id, status, mother1, mother2, daughter1, daughter2, 
      col, acol, phaseSpacePtr->p(i), phaseSpacePtr->m(i));
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
  if (softQCD || Settings::flag("SoftQCD:elastic")) {
    SigmaHat* sigmaPtr = new SigmaHAB2AB;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (softQCD || Settings::flag("SoftQCD:singleDiffractive")) {
    SigmaHat* sigmaPtr = new SigmaHAB2XB;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
    sigmaPtr = new SigmaHAB2AX;
    containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (softQCD || Settings::flag("SoftQCD:doubleDiffractive")) {
    SigmaHat* sigmaPtr = new SigmaHAB2XX;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  
  // Set up requested objects for hard QCD processes.
  bool hardQCD = Settings::flag("HardQCD:all");
  if (hardQCD || Settings::flag("HardQCD:gg2gg")) {
    SigmaHat* sigmaPtr = new SigmaHgg2gg;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (hardQCD || Settings::flag("HardQCD:gg2qqbar")) {
    SigmaHat* sigmaPtr = new SigmaHgg2qqbar;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (hardQCD || Settings::flag("HardQCD:qg2qg")) {
    SigmaHat* sigmaPtr = new SigmaHqg2qg;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (hardQCD || Settings::flag("HardQCD:qq2qq")) {
    SigmaHat* sigmaPtr = new SigmaHqq2qqDiff;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
    sigmaPtr = new SigmaHqq2qqSame;
    containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
    sigmaPtr = new SigmaHqqbar2qqbarSame;
    containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (hardQCD || Settings::flag("HardQCD:qqbar2qqbarNew")) {
    SigmaHat* sigmaPtr = new SigmaHqqbar2qqbarNew;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (hardQCD || Settings::flag("HardQCD:qqbar2gg")) {
    SigmaHat* sigmaPtr = new SigmaHqqbar2gg;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  
  // Set up requested objects for c cbar and b bbar, also hard QCD.
  // Warning: not finished yet! No masses in phase space! ??
  /*
  if (hardQCD || Settings::flag("HardQCD:qqbar2ccbar")) {
    SigmaHat* sigmaPtr = new SigmaHqqbar2QQbar(4, 122, "q qbar -> c cbar");
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (hardQCD || Settings::flag("HardQCD:qqbar2bbbar")) {
    SigmaHat* sigmaPtr = new SigmaHqqbar2QQbar(5, 124, "q qbar -> b bbar");
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  */

  // Set up requested objects for prompt photon processes.
  bool promptPhotons = Settings::flag("PromptPhoton:all");
  if (promptPhotons
    || Settings::flag("PromptPhoton:qg2qgamma")) {
    SigmaHat* sigmaPtr = new SigmaHqg2qgamma;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (promptPhotons 
    || Settings::flag("PromptPhoton:qqbar2ggamma")) {
    SigmaHat* sigmaPtr = new SigmaHqqbar2ggamma;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (promptPhotons
    || Settings::flag("PromptPhoton:gg2ggamma")) {
    SigmaHat* sigmaPtr = new SigmaHgg2ggamma;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (promptPhotons 
    || Settings::flag("PromptPhoton:qqbar2gammagamma")) {
    SigmaHat* sigmaPtr = new SigmaHqqbar2gammagamma;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 
  if (promptPhotons
    || Settings::flag("PromptPhoton:gg2gammagamma")) {
    SigmaHat* sigmaPtr = new SigmaHgg2gammagamma;
    ProcessContainer* containerPtr = new ProcessContainer(sigmaPtr);
    containerPtrs.push_back( containerPtr);
  } 

  // Done. 
  return true;
}

//**************************************************************************

} // end namespace Pythia8
