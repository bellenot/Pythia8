// Function definitions (not found in the header) for the ProcessLevel class.
// Copyright C 2007 Torbjorn Sjostrand

#include "ProcessLevel.h"

namespace Pythia8 {
 
//**************************************************************************

// Main routine to initialize generation process.

bool ProcessLevel::init( Info* infoPtrIn, BeamParticle* beamAPtrIn, 
  BeamParticle* beamBPtrIn, bool doLHAin, LHAinit* lhaInitPtrIn, 
  LHAevnt* lhaEvntPtrIn, UserHooks* userHooksPtrIn) {

  // Store input pointers for future use. 
  infoPtr       = infoPtrIn;
  beamAPtr      = beamAPtrIn;
  beamBPtr      = beamBPtrIn;
  userHooksPtr  = userHooksPtrIn;

  // Store pointers to Les Houches Accord input if any.
  doLHA         = doLHAin;
  lhaInitPtr    = lhaInitPtrIn;
  lhaEvntPtr    = lhaEvntPtrIn;
  strategyLHA   = (doLHA) ? lhaInitPtr->strategy() : 0;  

  // Initialize information on resonances and their decay treatment.
  if (!initResonances()) return false;

  // Option to allow second hard interaction.
  doSecondHard  = Settings::flag("SecondHard:generate");

  // If not Les Houches then internal machinery.
  doInternal    = !doLHA;
  if (doInternal) return initInternal();

  // Done. (Check return values from other classes??)
  return true;
}

//*********

// Main routine to generate the hard process.
// Currently rather primitive.

  bool ProcessLevel::next( Event& process) {

  // Starting value.
  bool physical = false;  

  // Generate the next internal event with two or one hard interactions. 
  if (doInternal) physical = (doSecondHard) 
    ? next2Internal( process) : nextInternal( process);

  // Read in a simple event in the LHAevnt format. Default info.
  else if (strategyLHA >= 10) {
    infoPtr->setType( "Simple LHA process", 0, 0, false, true, false, false);
    infoPtr->setTypeMI( 0, 0.);
    physical = nextSimpleLHA( process);
  }

  // Read in an event in the LHAevnt format.
  else if (doLHA) physical = nextLHA( process);

  // Check that colour assignments make sense.
  if (physical) physical = checkColours( process);

  // Done.
  return physical;
}

//*********

// Accumulate and update statistics (after possible user veto).
  
void ProcessLevel::accumulate() {

  // Currently does not handle LHA processes ??
  if (doLHA) return;

  // Increase number of accepted events.
  containerPtrs[iContainer]->accumulate();

  // Provide current generated cross section estimate.
  long   nTrySum    = 0; 
  long   nSelSum    = 0; 
  long   nAccSum    = 0;
  double sigmaSum   = 0.;
  double delta2Sum  = 0.;
  double sigSelSum  = 0.;
  for (int i = 0; i < int(containerPtrs.size()); ++i) 
  if (containerPtrs[i]->sigmaMax() != 0.) {
    nTrySum        += containerPtrs[i]->nTried();
    nSelSum        += containerPtrs[i]->nSelected();
    nAccSum        += containerPtrs[i]->nAccepted();
    sigmaSum       += containerPtrs[i]->sigmaMC();
    delta2Sum      += pow2(containerPtrs[i]->deltaMC()); 
    sigSelSum      += containerPtrs[i]->sigmaSelMC();
  }

  // Normally only one hard interaction. Then store info and done.
  if (!doSecondHard) {
    double deltaSum = sqrtpos(delta2Sum);
    infoPtr->setSigma( nTrySum, nSelSum, nAccSum, sigmaSum, deltaSum); 
    return;
  }

  // Increase counter for a second hard interaction.
  container2Ptrs[i2Container]->accumulate();

  // Update statistics on average impact factor.
  ++nImpact;
  sumImpactFac     += infoPtr->enhanceMI();
  sum2ImpactFac    += pow2(infoPtr->enhanceMI());

  // Cross section estimate for second hard process.
  double sigma2Sum  = 0.;
  double sig2SelSum = 0.;
  for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2) 
  if (container2Ptrs[i2]->sigmaMax() != 0.) {
    nTrySum        += container2Ptrs[i2]->nTried();
    sigma2Sum      += container2Ptrs[i2]->sigmaMC();
    sig2SelSum     += container2Ptrs[i2]->sigmaSelMC();
  }

  // Average impact-parameter factor and error.
  double invN       = 1. / max(1, nImpact);
  double impactFac  = max( 1., sumImpactFac * invN);
  double impactErr2 = ( sum2ImpactFac * invN / pow2(impactFac) - 1.) * invN;
     
  // Cross section estimate for combination of first and second process.
  // Combine two possible ways and take average.
  double sigmaComb  = 0.5 * (sigmaSum * sig2SelSum + sigSelSum * sigma2Sum);
  sigmaComb        *= impactFac / sigmaND;
  if (allHardSame) sigmaComb *= 0.5; 
  double deltaComb  = sqrtpos(2. / nAccSum + impactErr2) * sigmaComb;

  // Store info and done.
  infoPtr->setSigma( nTrySum, nSelSum, nAccSum, sigmaComb, deltaComb); 
  
}

//*********

// Print statistics on cross sections and number of events.

void ProcessLevel::statistics(ostream& os) {

  // Internal statistics only.
  if (!doInternal) return;

  // Special processing if two hard interactions selected.
  if (doSecondHard) { 
    statistics2(os);
    return;
  } 
    
  // Header.
  os << "\n *-------  PYTHIA Event and Cross Section Statistics  ------"
     << "--------------------------------------------------*\n"
     << " |                                                            "
     << "                                                |\n" 
     << " | Subprocess                               Code |            "
     << "Number of events       |      sigma +- delta    |\n" 
     << " |                                               |       Tried"
     << "   Selected   Accepted |     (estimated) (mb)   |\n"
     << " |                                               |            "
     << "                       |                        |\n"
     << " |------------------------------------------------------------"
     << "------------------------------------------------|\n"
     << " |                                               |            "
     << "                       |                        |\n";

  // Reset sum counters.
  long   nTrySum   = 0; 
  long   nSelSum   = 0; 
  long   nAccSum   = 0;
  double sigmaSum  = 0.;
  double delta2Sum = 0.;

  // Loop over existing processes.
  for (int i = 0; i < int(containerPtrs.size()); ++i) 
  if (containerPtrs[i]->sigmaMax() != 0.) {

    // Read info for process. Sum counters.
    long   nTry    = containerPtrs[i]->nTried();
    long   nSel    = containerPtrs[i]->nSelected();
    long   nAcc    = containerPtrs[i]->nAccepted();
    double sigma   = containerPtrs[i]->sigmaMC();
    double delta   = containerPtrs[i]->deltaMC(); 
    nTrySum       += nTry;
    nSelSum       += nSel;
    nAccSum       += nAcc; 
    sigmaSum      += sigma;
    delta2Sum     += pow2(delta);    

    // Print individual process info.
    os << " | " << left << setw(40) << containerPtrs[i]->name() 
       << right << setw(5) << containerPtrs[i]->code() << " | " 
       << setw(11) << nTry << " " << setw(10) << nSel << " " 
       << setw(10) << nAcc << " | " << scientific << setprecision(3) 
       << setw(11) << sigma << setw(11) << delta << " |\n";
  }

  // Print summed process info.
  os << " |                                               |            "
     << "                       |                        |\n"
     << " | " << left << setw(45) << "sum" << right << " | " << setw(11) 
     << nTrySum << " " << setw(10) << nSelSum << " " << setw(10) 
     << nAccSum << " | " << scientific << setprecision(3) << setw(11) 
     << sigmaSum << setw(11) << sqrtpos(delta2Sum) << " |\n";

  // Listing finished.
  os << " |                                                            "
     << "                                                |\n"
     << " *-------  End PYTHIA Event and Cross Section Statistics -----"
     << "------------------------------------------------*" << endl;

}

//*********

// Initialize information on resonances.

bool ProcessLevel::initResonances() {

  // Initialize static data members for each resonance.
  ResonanceGmZ::initStatic();
  ResonanceW::initStatic();
  ResonanceH::initStatic();

  // Recalculate widths for current masses.
  ResonanceGmZ::widthInit();
  ResonanceW::widthInit();
  ResonanceH::widthInit();

  // Send ResonanceDecays pointer to ProcessContainer.
  ProcessContainer::setResonanceDecaysPtr( &resonanceDecays);

  // Done.
  return true;

}

//*********

// Initialize the internal event generation machinery.
  
bool ProcessLevel::initInternal( ostream& os) {

  // Set up SigmaTotal. Store sigma_nondiffractive for future use.
  int    idA = infoPtr->idA();
  int    idB = infoPtr->idB();
  double eCM = infoPtr->eCM();
  sigmaTot.init( idA, idB, eCM);
  sigmaND = sigmaTot.sigmaND();

  // Send cross section and beam info to influx, processes and phase space.
  double mA = beamAPtr->m();
  double mB = beamBPtr->m();
  InFlux::setStaticPtrs( beamAPtr, beamBPtr); 
  SigmaProcess::setStaticPtrs( &sigmaTot, idA, idB, mA, mB);
  PhaseSpace::setStaticPtrs( beamAPtr, beamBPtr, &sigmaTot, userHooksPtr);

  // Sets up containers for all the hard processes.
  SetupContainers setupContainers;
  setupContainers.init(containerPtrs);

  // If no processes found then refuse to do anything.
  if ( int(containerPtrs.size()) == 0) {
    ErrorMsg::message("Error in ProcessLevel::initInternal: "
      "no process switched on"); 
    return false;
  }

  // Initialize each process. 
  int numberOn = 0;
  for (int i = 0; i < int(containerPtrs.size()); ++i)
    if (containerPtrs[i]->init()) ++numberOn;

  // Sum maxima for Monte Carlo choice.
  sigmaMaxSum = 0.;
  for (int i = 0; i < int(containerPtrs.size()); ++i)
    sigmaMaxSum += containerPtrs[i]->sigmaMax();

  // Option to pick a second hard interaction: repeat as above.
  int number2On = 0;
  if (doSecondHard) {
    setupContainers.init2(container2Ptrs);
    if ( int(container2Ptrs.size()) == 0) {
      ErrorMsg::message("Error in ProcessLevel::initInternal: "
        "no second hard process switched on"); 
      return false;
    }
    for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2)
      if (container2Ptrs[i2]->init()) ++number2On;
    sigma2MaxSum = 0.;
    for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2)
      sigma2MaxSum += container2Ptrs[i2]->sigmaMax();
  }

  // Construct string with incoming beams and for cm energy.
  string collision = "We collide " + ParticleDataTable::name(idA)
    + " with " + ParticleDataTable::name(idB) + " at a CM energy of "; 
  string pad( 46 - collision.length(), ' ');

  // Print initialization information: header.
  os << "\n *-------  PYTHIA Process Initialization  ---------------------*\n"
     << " |                                                             |\n" 
     << " | " << collision << scientific << setprecision(3)<< setw(9) << eCM 
     << " GeV" << pad << " |\n"
     << " |                                                             |\n"
     << " |-------------------------------------------------------------|\n"
     << " |                                               |             |\n"
     << " | Subprocess                               Code |   Estimated |\n" 
     << " |                                               |    max (mb) |\n"
     << " |                                               |             |\n"
     << " |-------------------------------------------------------------|\n"
     << " |                                               |             |\n";


  // Loop over existing processes: print individual process info.
  for (int i = 0; i < int(containerPtrs.size()); ++i) 
  os << " | " << left << setw(40) << containerPtrs[i]->name() 
     << right << setw(5) << containerPtrs[i]->code() << " | " 
     << scientific << setprecision(3) << setw(11)  
     << containerPtrs[i]->sigmaMax() << " |\n";

  // Loop over second hard processes, if any, and repeat as above.
  if (doSecondHard) {
    os << " |                                               |             |\n"
       << " |-------------------------------------------------------------|\n"
       << " |                                               |             |\n";
    for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2) 
    os << " | " << left << setw(40) << container2Ptrs[i2]->name() 
       << right << setw(5) << container2Ptrs[i2]->code() << " | " 
       << scientific << setprecision(3) << setw(11)  
       << container2Ptrs[i2]->sigmaMax() << " |\n";
  }

  // Listing finished.
  os << " |                                                             |\n" 
     << " *-------  End PYTHIA Process Initialization ------------------*" 
     << endl;

  // If sum of maxima vanishes then refuse to do anything.
  if ( numberOn == 0  || sigmaMaxSum <= 0.) {
    ErrorMsg::message("Error in ProcessLevel::initInternal: "
      "all processes have vanishing cross sections"); 
    return false;
  }
  if ( doSecondHard && (number2On == 0  || sigma2MaxSum <= 0.) ) {
    ErrorMsg::message("Error in ProcessLevel::initInternal: "
      "all second hard processes have vanishing cross sections"); 
    return false;
  }
  
  // If two hard processes then check whether some (but not all) agree.
  allHardSame  = true;
  noneHardSame = true;
  if (doSecondHard) {
    bool foundMatch = false;
    
    // Check for each first process if matched in second.
    for (int i = 0; i < int(containerPtrs.size()); ++i) {
      foundMatch = false;
      for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2) 
        if (container2Ptrs[i2]->code() == containerPtrs[i]->code()) 
          foundMatch = true;
      containerPtrs[i]->isSame( foundMatch );
      if (!foundMatch)  allHardSame = false;
      if ( foundMatch) noneHardSame = false; 
    }

    // Check for each second process if matched in first.
    for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2) {
      foundMatch = false;
      for (int i = 0; i < int(containerPtrs.size()); ++i) 
        if (containerPtrs[i]->code() == container2Ptrs[i2]->code()) 
          foundMatch = true;
      container2Ptrs[i2]->isSame( foundMatch );
      if (!foundMatch)  allHardSame = false;
      if ( foundMatch) noneHardSame = false;   
    }
  }
  someHardSame = !allHardSame && !noneHardSame;

  // Reset counters for average impact-parameter enhancement.
  nImpact       = 0;
  sumImpactFac  = 0.;
  sum2ImpactFac = 0.;

  // Done.
  return true;

}

//*********

// Generate the next internal event.
  
bool ProcessLevel::nextInternal( Event& process) {

  // Loop over tries until trial event succeeds.
  for ( ; ; ) {

    // Pick one of the subprocesses.
    double sigmaMaxNow = sigmaMaxSum * Rndm::flat();
    int iMax = containerPtrs.size() - 1;
    iContainer = -1;
    do sigmaMaxNow -= containerPtrs[++iContainer]->sigmaMax();
    while (sigmaMaxNow > 0. && iContainer < iMax);
    
    // Do a trial event of this subprocess; accept or not.
    if (containerPtrs[iContainer]->trialProcess()) break;
  }

  // Construct kinematics of acceptable process.
  containerPtrs[iContainer]->constructProcess( process);

  // Do all resonance decays.
  if ( !containerPtrs[iContainer]->decayResonances( process) ) return false;

  // Add any junctions to the process event record list.
  findJunctions( process);

  // Done.
  return true;
}

//*********

// Generate the next internal event with two hard interactions.
  
bool ProcessLevel::next2Internal( Event& process) {

  // Loop over both hard processes to find consistent common kinematics.
  for ( ; ; ) {
   
    // Loop internally over tries for hardest process until succeeds.
    for ( ; ; ) {

      // Pick one of the subprocesses.
      double sigmaMaxNow = sigmaMaxSum * Rndm::flat();
      int iMax = containerPtrs.size() - 1;
      iContainer = -1;
      do sigmaMaxNow -= containerPtrs[++iContainer]->sigmaMax();
      while (sigmaMaxNow > 0. && iContainer < iMax);
    
      // Do a trial event of this subprocess; accept or not.
      if (containerPtrs[iContainer]->trialProcess()) break;
    }

    // Loop internally over tries for second hardest process until succeeds.
    for ( ; ; ) {

      // Pick one of the subprocesses.
      double sigma2MaxNow = sigma2MaxSum * Rndm::flat();
      int i2Max = container2Ptrs.size() - 1;
      i2Container = -1;
      do sigma2MaxNow -= container2Ptrs[++i2Container]->sigmaMax();
      while (sigma2MaxNow > 0. && i2Container < i2Max);
    
      // Do a trial event of this subprocess; accept or not.
      if (container2Ptrs[i2Container]->trialProcess()) break;
    }

    // Check whether common set of x values is kinematically possible.
    double xA1      = containerPtrs[iContainer]->x1();
    double xB1      = containerPtrs[iContainer]->x2();
    double xA2      = container2Ptrs[i2Container]->x1();
    double xB2      = container2Ptrs[i2Container]->x2();    
    if (xA1 + xA2 >= 1. || xB1 + xB2 >= 1.) continue;

    // Reset beam contents. Naive parton densities for second interaction.
    // (Subsequent procedure could be symmetrized, but would be overkill.)
    beamAPtr->clear();    
    beamBPtr->clear();    
    int    idA1     = containerPtrs[iContainer]->id1();
    int    idB1     = containerPtrs[iContainer]->id2();
    int    idA2     = container2Ptrs[i2Container]->id1();
    int    idB2     = container2Ptrs[i2Container]->id2();
    double Q2Fac1   = containerPtrs[iContainer]->Q2Fac();
    double Q2Fac2   = container2Ptrs[i2Container]->Q2Fac();
    double pdfA2Raw = beamAPtr->xf( idA2, xA2,Q2Fac2);
    double pdfB2Raw = beamBPtr->xf( idB2, xB2,Q2Fac2);

    // Remove partons in first interaction from beams. 
    beamAPtr->append( 3, idA1, xA1);
    beamAPtr->xfISR( 0, idA1, xA1, Q2Fac1);
    beamAPtr->pickValSeaComp(); 
    beamBPtr->append( 4, idB1, xB1);
    beamBPtr->xfISR( 0, idB1, xB1, Q2Fac1);
    beamBPtr->pickValSeaComp(); 

    // Reevaluate pdf's for second interaction and weight by reduction.
    double pdfA2Mod = beamAPtr->xfMI( idA2, xA2,Q2Fac2);
    double pdfB2Mod = beamBPtr->xfMI( idB2, xB2,Q2Fac2);
    double wtPdfMod = (pdfA2Mod * pdfB2Mod) / (pdfA2Raw * pdfB2Raw); 
    if (wtPdfMod < Rndm::flat()) continue;

    // Reduce by a factor of 2 for identical processes when others not.
    if ( someHardSame && containerPtrs[iContainer]->isSame() 
      && container2Ptrs[i2Container]->isSame() && Rndm::flat() > 0.5) 
      continue;

    // If come this far then acceptable event.
    break;
  }

  // Construct kinematics of acceptable processes.
  Event process2;
  process2.initColTag();
  startColTag2 = process2.lastColTag();
  containerPtrs[iContainer]->constructProcess( process);
  container2Ptrs[i2Container]->constructProcess( process2, false);

  // Do all resonance decays.
  if ( !containerPtrs[iContainer]->decayResonances( process) ) 
    return false;
  if ( !container2Ptrs[i2Container]->decayResonances( process2) ) 
    return false;

  // Append second hard interaction to normal process object.
  combineProcessRecords( process, process2);

  // Add any junctions to the process event record list.
  findJunctions( process);

  // Done.
  return true;
}

//*********

// Append second hard interaction to normal process object.
// Complication: all resonance decay chains must be put at the end.

void ProcessLevel::combineProcessRecords( Event& process, Event& process2) {

  // Find first event record size, excluding resonances.
  int nSize = process.size();
  int nHard = 5;
  while (nHard < nSize && process[nHard].mother1() == 3) ++nHard;

  // Save resonance products temporarily elsewhere. 
  vector<Particle> resProd;
  if (nSize > nHard) {
    for (int i = nHard; i < nSize; ++i) resProd.push_back( process[i] );
    process.popBack(nSize - nHard);
  }

  // Find second event record size, excluding resonances.
  int nSize2 = process2.size();
  int nHard2 = 5;
  while (nHard2 < nSize2 && process2[nHard2].mother1() == 3) ++nHard2;

  // Find amount of necessary position and colour offset for second process.
  int addPos  = nHard  - 3;
  int addCol  = process.lastColTag() - startColTag2;

  // Loop over all particles (except beams) from second process.
  for (int i = 3; i < nSize2; ++i) {

    // Offset mother and daughter pointers and colour tags of particle.
    process2[i].offsetHistory( 2, addPos, 2, addPos);
    process2[i].offsetCol( addCol);

    // Append hard-process particles from process2 to process.
    if (i < nHard2) process.append( process2[i] );
  }

  // Reinsert resonance decay chains of first hard process.
  int addPos2 = nHard2 - 3;
  if (nHard < nSize) {

    // Offset daughter pointers of unmoved mothers.
    for (int i = 5; i < nHard; ++i) 
      process[i].offsetHistory( 0, 0, nHard - 1, addPos2); 
    
    // Modify history of resonance products when restoring. 
    for (int i = 0; i < int(resProd.size()); ++i) {
      resProd[i].offsetHistory( nHard - 1, addPos2, nHard - 1, addPos2);
      process.append( resProd[i] );
    } 
  }   

  // Insert resonance decay chains of second hard process.
  if (nHard2 < nSize2) {
    int nHard3  = nHard + nHard2 - 3;
    int addPos3 = nSize - nHard;

    // Offset daughter pointers of second-process mothers.
    for (int i = nHard + 2; i < nHard3; ++i) 
      process[i].offsetHistory( 0, 0, nHard3 - 1, addPos3); 
    
    // Modify history of second-process resonance products and insert. 
    for (int i = nHard2; i < nSize2; ++i) {
      process2[i].offsetHistory( nHard3 - 1, addPos3, nHard3 - 1, addPos3);
      process.append( process2[i] );
    } 
  }

  // Store PDF scale for second interaction.
  process.scaleSecond( process2.scale() );   

}

//*********

// Read in the hard process from the Les Houches Accord.
// Many more checks to be done for valid input??

bool ProcessLevel::nextLHA( Event& process) {

  // Generate the next Les Houches event.
  if (!lhaEvntPtr->set()) return false;
     
  // Let hard process record begin with the event as a whole and
  // the two incoming beam particles.  
  process.append( 90, -11, 0, 0, 1, 2, 0, 0, 
    Vec4(0., 0., 0., infoPtr->eCM()), infoPtr->eCM(), 0. ); 
  process.append( infoPtr->idA(), -12, 0, 0, 3, 0, 0, 0, 
    Vec4(0., 0., infoPtr->pzA(), infoPtr->eA()), infoPtr->mA(), 0. ); 
  process.append( infoPtr->idB(), -12, 0, 0, 4, 0, 0, 0, 
    Vec4(0., 0., infoPtr->pzB(), infoPtr->eB()), infoPtr->mB(), 0. ); 

  // Since LHA partons may be out of order, determine correct one.
  // (Recall that zeroth particle is empty.) 
  vector<int> newPos;
  newPos.reserve(lhaEvntPtr->size());
  newPos.push_back(0);
  for (int iNew = 0; iNew < lhaEvntPtr->size(); ++iNew) {
    // For iNew == 0 look for the two incoming partons, then for
    // partons having them as mothers, and so on layer by layer.
    for (int i = 1; i < lhaEvntPtr->size(); ++i)
      if (lhaEvntPtr->mother1(i) == newPos[iNew]) newPos.push_back(i);
    if (int(newPos.size()) <= iNew) break;
  } 

  // Find scale from which to begin evolution.
  double scale = lhaEvntPtr->scale();
  process.scale( scale);

  // Copy over info from LHA event to process, in proper order.
  for (int i = 1; i < lhaEvntPtr->size(); ++i) {
    int iOld = newPos[i];
    int id = lhaEvntPtr->id(iOld);

    // Translate from LHA status codes.
    int lhaStatus =  lhaEvntPtr->status(iOld);
    int status = -21;
    if (lhaStatus == 2 || lhaStatus == 3) status = -22;
    if (lhaStatus == 1) status = 23;

    // Find where mothers have been moved by reordering.
    int mother1Old = lhaEvntPtr->mother1(iOld);   
    int mother2Old = lhaEvntPtr->mother2(iOld);   
    int mother1 = 0;
    int mother2 = 0; 
    for (int im = 1; im < i; ++im) {
      if (mother1Old == newPos[im]) mother1 = im + 2; 
      if (mother2Old == newPos[im]) mother2 = im + 2; 
    } 
    if (i <= 2) mother1 = i;

    // Find daughters and where they have been moved by reordering. 
    // (Values shifted two steps to account for inserted beams.)
    int daughter1 = 0;
    int daughter2 = 0;
    for (int im = i + 1; im < lhaEvntPtr->size(); ++im) { 
      if (lhaEvntPtr->mother1(newPos[im]) == iOld
        || lhaEvntPtr->mother2(newPos[im]) == iOld) {
        if (daughter1 == 0 || im + 2 < daughter1) daughter1 = im + 2;
        if (daughter2 == 0 || im + 2 > daughter2) daughter2 = im + 2;
      }
    }
    // For 2 -> 1 hard scatterings reset second daughter to 0.
    if (daughter2 == daughter1) daughter2 = 0;

    // Colour trivial, except reset irrelevant colour indices.
    int colType = ParticleDataTable::colType(id);
    int col1 = (colType == 1 || colType == 2) ? lhaEvntPtr->col1(iOld) : 0;   
    int col2 = (colType == -1 || colType == 2) ?  lhaEvntPtr->col2(iOld) : 0; 

    // Momentum trivial.
    double px = lhaEvntPtr->px(iOld);  
    double py = lhaEvntPtr->py(iOld);  
    double pz = lhaEvntPtr->pz(iOld);  
    double e  = lhaEvntPtr->e(iOld);  
    double m  = lhaEvntPtr->m(iOld);

    // For resonance decay products use resonance mass as scale.
    double scaleNow = scale;
    if (mother1 > 4) scaleNow = process[mother1].m();
    process.append( id, status, mother1, mother2, daughter1, daughter2, 
      col1, col2, Vec4(px, py, pz, e), m, scaleNow);
  }  

  // Add any junctions to the process event record list.
  findJunctions( process);

  // Extract information that is guaranteed available.
  string name = "External LHA process"; 
  int code = lhaEvntPtr->idProc();
  int nFinal = 0;
  for (int i = 5; i < process.size(); ++i) 
    if (process[i].mother1() == 3) ++nFinal;
  int    id1     =  process[3].id(); 
  int    id2     =  process[4].id(); 
  double Q2Fac   = pow2(lhaEvntPtr->scale());
  double alphaEM = lhaEvntPtr->alphaQED();
  double alphaS  = lhaEvntPtr->alphaQCD();
  double Q2Ren   = Q2Fac;
  double x1      = 2. * process[3].e() / infoPtr->eCM();
  double x2      = 2. * process[4].e() / infoPtr->eCM();
  Vec4   pSum    = process[3].p() + process[4].p();
  double sHat    = pSum*pSum;

  // Reset quantities that may or may not be known.
  double pdf1  = 0.;
  double pdf2  = 0.;
  double tHat  = 0.;
  double uHat  = 0.;
  double pTHat = 0.;
  double m3    = 0.;
  double m4    = 0.;
  double theta = 0.;
  double phi   = 0.;

  // Read info on parton densities if provided.
  if (lhaEvntPtr->pdfIsSet()) {
    pdf1  = lhaEvntPtr->xpdf1();
    pdf2  = lhaEvntPtr->xpdf2();
    Q2Fac = pow2(lhaEvntPtr->scalePDF());
    x1    = lhaEvntPtr->x1();
    x2    = lhaEvntPtr->x2();
  }

  // Reconstruct kinematics of 2 -> 2 processes from momenta.
  if (nFinal == 2) {
    Vec4 pDifT = process[3].p() - process[5].p();
    tHat = pDifT*pDifT;    
    Vec4 pDifU = process[3].p() - process[6].p();
    uHat = pDifU*pDifU;
    pTHat = process[5].pT();
    m3 = process[5].m();    
    m4 = process[6].m(); 
    Vec4 p5 = process[5].p();
    p5.bstback(pSum);
    theta = p5.theta();   
    phi = process[5].phi();   
  }

  // Store information in Info object.
  infoPtr->setType( name, code, nFinal, false, true, false, false);
  infoPtr->setPDFalpha( id1, id2, pdf1, pdf2, Q2Fac, alphaEM, alphaS, Q2Ren);
  infoPtr->setKin( x1, x2, sHat, tHat, uHat, pTHat, m3, m4, theta, phi);  
  infoPtr->setTypeMI( code, pTHat);

  // Do all resonance decays. First draft?? Junctions in decays??
  resonanceDecays.next( process);

  // Done.
  return true;

}

//*********

// Read in the hard process, special case if all partons already given.

bool ProcessLevel::nextSimpleLHA( Event& process) {

  // Generate the next Les Houches event.
  if (!lhaEvntPtr->set()) return false;

  // Calculate the total four-momentum of the event.
  // Isolate original system, i.e. those with mother = 0.
  int nTot = lhaEvntPtr->size();
  int nDaughter = 0;
  Vec4 pSum;
  for (int i = 1; i < nTot; ++i) {
    if (lhaEvntPtr->status(i) == 1) pSum += Vec4( lhaEvntPtr->px(i),
      lhaEvntPtr->py(i), lhaEvntPtr->pz(i), lhaEvntPtr->e(i));
    if (lhaEvntPtr->mother1(i) == 0) nDaughter = i;
  }
  double eCM = pSum.mCalc();
    
  // Let hard process record begin with the event as a whole.
  process.append( 90, -11, 0, 0, 1, nDaughter, 0, 0, pSum, eCM, 0. ); 

  // Copy over info from LHA event to process, keeping the given order.
  for (int i = 1; i < nTot; ++i) {

    // Translate from LHA status codes.
    int lhaStatus =  lhaEvntPtr->status(i);
    int status = -21;
    if (lhaStatus == 2 || lhaStatus == 3) status = -22;
    if (lhaStatus == 1) status = 23;

    // Read id and mother information.
    int id = lhaEvntPtr->id(i);
    int mother1 = lhaEvntPtr->mother1(i);   
    int mother2 = lhaEvntPtr->mother2(i);   

    // Find daughters. 
    int daughter1 = 0;
    int daughter2 = 0;
    for (int im = 1; im < nTot; ++im)  
    if (lhaEvntPtr->mother1(im) == i || lhaEvntPtr->mother2(im) == i) {
      if (daughter1 == 0 || im < daughter1) daughter1 = im;
      if (daughter2 == 0 || im > daughter2) daughter2 = im;
    }
    // For 2 -> 1 hard scatterings reset second daughter to 0.
    if (daughter2 == daughter1) daughter2 = 0;

    // Colour trivial, except reset irrelevant colour indices.
    int colType = ParticleDataTable::colType(id);
    int col1 = (colType == 1 || colType == 2) ? lhaEvntPtr->col1(i) : 0;   
    int col2 = (colType == -1 || colType == 2) ?  lhaEvntPtr->col2(i) : 0; 

    // Momentum trivial.
    double px = lhaEvntPtr->px(i);  
    double py = lhaEvntPtr->py(i);  
    double pz = lhaEvntPtr->pz(i);  
    double e  = lhaEvntPtr->e(i);  
    double m  = lhaEvntPtr->m(i);

    // Store the information, particle by particle.
    process.append( id, status, mother1, mother2, daughter1, daughter2, 
      col1, col2, Vec4(px, py, pz, e), m, 0.);
  }  

  // Add any junctions to the process event record list.
  findJunctions( process);

  // Done.
  return true;
}

//*********

// Add any junctions to the process event record list.
// First try, so still incomplete. ??

void ProcessLevel::findJunctions( Event& process) {

  // Loop though event; isolate all uncoloured particles.
  for (int i = 0; i < process.size(); ++i) 
  if ( process[i].col() == 0 && process[i].acol() == 0) {

    // Find all daughters and store daughter colours and anticolours.
    vector<int> daughters = process.daughterList(i);
    vector<int> cols, acols;
    for (int j = 0; j < int(daughters.size()); ++j) {
      int colDau  = process[ daughters[j] ].col();
      int acolDau = process[ daughters[j] ].acol();
      if (colDau > 0)  cols.push_back( colDau);      
      if (acolDau > 0) acols.push_back( acolDau);      
    }

    // Remove all matching colour-anticolour pairs.
    bool foundPair = true;
    while (foundPair && cols.size() > 0 && acols.size() > 0) {
      foundPair = false;
      for (int j = 0; j < int(cols.size()); ++j) {
        for (int k = 0; k < int(acols.size()); ++k) {
	  if (acols[k] == cols[j]) { 
            cols[j]  = cols.back();  cols.pop_back();     
            acols[k] = acols.back(); acols.pop_back();     
            foundPair = true; break;
	  }
	} if (foundPair) break;
      }
    } 

    // Store an (anti)junction when three (anti)coloured daughters.
    if (cols.size() == 3 && acols.size() == 0) 
      process.appendJunction( 1, cols[0], cols[1], cols[2]);
    if (acols.size() == 3 && cols.size() == 0) 
      process.appendJunction( 2, acols[0], acols[1], acols[2]);
  }

  // Done.
}

//*********

// Check that colours match up.

bool ProcessLevel::checkColours( Event& process) {

  // Variables and arrays for common usage.
  bool physical = true;
  bool match;
  int colType, col, acol, iPos, iNow, iNowA;
  vector<int> colTags, colPos, acolPos;

  // Check that each particle has the kind of colours expected of it.
  for (int i = 0; i < process.size(); ++i) {
    colType = process[i].colType();
    col     = process[i].col();
    acol    = process[i].acol();
    if      (colType ==  0 && (col != 0 || acol != 0)) physical = false;
    else if (colType ==  1 && (col <= 0 || acol != 0)) physical = false;
    else if (colType == -1 && (col != 0 || acol <= 0)) physical = false;
    else if (colType ==  2 && (col <= 0 || acol <= 0)) physical = false;
    else if (colType < -1 || colType > 2)              physical = false; 

    // Add to the list of colour tags.
    if (col > 0) {
      match = false;
      for (int ic = 0; ic < int(colTags.size()) ; ++ic)
        if (col == colTags[ic]) match = true;
      if (!match) colTags.push_back(col);
    } else if (acol > 0) {
      match = false;
      for (int ic = 0; ic < int(colTags.size()) ; ++ic)
        if (acol == colTags[ic]) match = true;
      if (!match) colTags.push_back(acol);
    }
  }

  // Warn and give up if particles did not have the expected colours.
  if (!physical) {
    ErrorMsg::message("Error in ProcessLevel::checkColours: "
      "incorrect colour assignment"); 
    return false;
  }

  // Loop through all colour tags and find their positions (by sign). 
  for (int ic = 0; ic < int(colTags.size()); ++ic) {
    col = colTags[ic];
    colPos.resize(0);
    acolPos.resize(0);
    for (int i = 0; i < process.size(); ++i) {
      if (process[i].col() == col) colPos.push_back(i); 
      if (process[i].acol() == col) acolPos.push_back(i); 
    }

    // Trace colours back through decays; remove daughters.
    while (colPos.size() > 1) {
      iPos = colPos.size() - 1; 
      iNow = colPos[iPos]; 
      if ( process[iNow].mother1() == colPos[iPos - 1]
        && process[iNow].mother2() == 0) colPos.pop_back();
      else break;
    }           
    while (acolPos.size() > 1) {
      iPos = acolPos.size() - 1; 
      iNow = acolPos[iPos]; 
      if ( process[iNow].mother1() == acolPos[iPos - 1]
        && process[iNow].mother2() == 0) acolPos.pop_back();
      else break;
    } 

    // Now colour should exist in only 2 copies.
    if (colPos.size() + acolPos.size() != 2) physical = false;

    // If both colours or both anticolours then one mother of the other.
    else if (colPos.size() == 2) {
      iNow = colPos[1];
      if ( process[iNow].mother1() != colPos[0] 
        && process[iNow].mother2() != colPos[0] ) physical = false;
    }
    else if (acolPos.size() == 2) {
      iNowA = acolPos[1];
      if ( process[iNowA].mother1() != acolPos[0] 
        && process[iNowA].mother2() != acolPos[0] ) physical = false;
    }
    
    // If one of each then should have same mother(s), or point to beams.
    else {
      iNow  = colPos[0];
      iNowA = acolPos[0];
      if ( process[iNow].status() == -21 &&  process[iNowA].status() == -21 );
      else if ( (process[iNow].mother1() != process[iNowA].mother1()) 
             || (process[iNow].mother2() != process[iNowA].mother2()) ) 
             physical = false;
    }

  }

  // Error message if problem found. Done.
  if (!physical) ErrorMsg::message("Error in ProcessLevel::checkColours: "
    "unphysical colour flow"); 
  return physical;

}

//*********

// Print statistics when two hard processes allowed.

void ProcessLevel::statistics2(ostream& os) {

  // Average impact-parameter factor and error.
  double invN          = 1. / max(1, nImpact);
  double impactFac     = max( 1., sumImpactFac * invN);
  double impactErr2    = ( sum2ImpactFac * invN / pow2(impactFac) - 1.) * invN;

  // Derive scaling factor to be applied to first set of processes.
  double sigma2SelSum  = 0.;
  int    n2SelSum      = 0;
  for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2) {
    sigma2SelSum      += container2Ptrs[i2]->sigmaSelMC();
    n2SelSum          += container2Ptrs[i2]->nSelected();
  }
  double factor1       = impactFac * sigma2SelSum / sigmaND;  
  double rel1Err       = sqrt(1. / max(1, n2SelSum) + impactErr2);    
  if (allHardSame) factor1 *= 0.5;

  // Derive scaling factor to be applied to second set of processes.
  double sigma1SelSum  = 0.;
  int    n1SelSum      = 0;
  for (int i = 0; i < int(containerPtrs.size()); ++i) {
    sigma1SelSum      += containerPtrs[i]->sigmaSelMC();
    n1SelSum          += containerPtrs[i]->nSelected();
  }
  double factor2       = impactFac * sigma1SelSum / sigmaND;       
  if (allHardSame) factor2 *= 0.5;
  double rel2Err       = sqrt(1. / max(1, n1SelSum) + impactErr2);    
    
  // Header.
  os << "\n *-------  PYTHIA Event and Cross Section Statistics  ------"
     << "--------------------------------------------------*\n"
     << " |                                                            "
     << "                                                |\n" 
     << " | Subprocess                               Code |            "
     << "Number of events       |      sigma +- delta    |\n" 
     << " |                                               |       Tried"
     << "   Selected   Accepted |     (estimated) (mb)   |\n"
     << " |                                               |            "
     << "                       |                        |\n"
     << " |------------------------------------------------------------"
     << "------------------------------------------------|\n"
     << " |                                               |            "
     << "                       |                        |\n"
     << " | First hard process:                           |            "
     << "                       |                        |\n"
     << " |                                               |            "
     << "                       |                        |\n";

  // Reset sum counters.
  long   nTrySum   = 0; 
  long   nSelSum   = 0; 
  long   nAccSum   = 0;
  double sigmaSum  = 0.;
  double delta2Sum = 0.;

  // Loop over existing first processes.
  for (int i = 0; i < int(containerPtrs.size()); ++i) 
  if (containerPtrs[i]->sigmaMax() != 0.) {

    // Read info for process. Sum counters.
    long   nTry    = containerPtrs[i]->nTried();
    long   nSel    = containerPtrs[i]->nSelected();
    long   nAcc    = containerPtrs[i]->nAccepted();
    double sigma   = containerPtrs[i]->sigmaMC() * factor1;
    double delta2  = pow2( containerPtrs[i]->deltaMC() * factor1 );
    nTrySum       += nTry;
    nSelSum       += nSel;
    nAccSum       += nAcc; 
    sigmaSum      += sigma;
    delta2Sum     += delta2; 
    delta2        += pow2( sigma * rel1Err );   

    // Print individual process info.
    os << " | " << left << setw(40) << containerPtrs[i]->name() 
       << right << setw(5) << containerPtrs[i]->code() << " | " 
       << setw(11) << nTry << " " << setw(10) << nSel << " " 
       << setw(10) << nAcc << " | " << scientific << setprecision(3) 
       << setw(11) << sigma << setw(11) << sqrtpos(delta2) << " |\n";
  }

  // Print summed info for first processes.
  delta2Sum       += pow2( sigmaSum * rel1Err ); 
  os << " |                                               |            "
     << "                       |                        |\n"
     << " | " << left << setw(45) << "sum" << right << " | " << setw(11) 
     << nTrySum << " " << setw(10) << nSelSum << " " << setw(10) 
     << nAccSum << " | " << scientific << setprecision(3) << setw(11) 
     << sigmaSum << setw(11) << sqrtpos(delta2Sum) << " |\n";

 
  // Separation lines to second hard processes.
  os << " |                                               |            "
     << "                       |                        |\n"
     << " |------------------------------------------------------------"
     << "------------------------------------------------|\n"
     << " |                                               |            "
     << "                       |                        |\n"
     << " | Second hard process:                          |            "
     << "                       |                        |\n"
     << " |                                               |            "
     << "                       |                        |\n";

  // Reset sum counters.
  nTrySum   = 0; 
  nSelSum   = 0; 
  nAccSum   = 0;
  sigmaSum  = 0.;
  delta2Sum = 0.;

  // Loop over existing second processes.
  for (int i2 = 0; i2 < int(container2Ptrs.size()); ++i2)
  if (container2Ptrs[i2]->sigmaMax() != 0.) {

    // Read info for process. Sum counters.
    long   nTry    = container2Ptrs[i2]->nTried();
    long   nSel    = container2Ptrs[i2]->nSelected();
    long   nAcc    = container2Ptrs[i2]->nAccepted();
    double sigma   = container2Ptrs[i2]->sigmaMC() * factor2;
    double delta2  = pow2( container2Ptrs[i2]->deltaMC() * factor2 );
    nTrySum       += nTry;
    nSelSum       += nSel;
    nAccSum       += nAcc; 
    sigmaSum      += sigma;
    delta2Sum     += delta2;    
    delta2        += pow2( sigma * rel2Err );   

    // Print individual process info.
    os << " | " << left << setw(40) << container2Ptrs[i2]->name() 
       << right << setw(5) << container2Ptrs[i2]->code() << " | " 
       << setw(11) << nTry << " " << setw(10) << nSel << " " 
       << setw(10) << nAcc << " | " << scientific << setprecision(3) 
       << setw(11) << sigma << setw(11) << sqrtpos(delta2) << " |\n";
  }

  // Print summed info for second processes.
  delta2Sum       += pow2( sigmaSum * rel2Err ); 
  os << " |                                               |            "
     << "                       |                        |\n"
     << " | " << left << setw(45) << "sum" << right << " | " << setw(11) 
     << nTrySum << " " << setw(10) << nSelSum << " " << setw(10) 
     << nAccSum << " | " << scientific << setprecision(3) << setw(11) 
     << sigmaSum << setw(11) << sqrtpos(delta2Sum) << " |\n";

  // Print information on how the two processes were combined.
  os << " |                                               |            "
     << "                       |                        |\n"
     << " |------------------------------------------------------------"
     << "------------------------------------------------|\n"
     << " |                                                            "
     << "                                                |\n"
     << " | Uncombined cross sections for the two event sets were " 
     << setw(10) << sigma1SelSum << " and " << sigma2SelSum << " mb, "
     << "respectively, combined  |\n"
     << " | using a sigma(nonDiffractive) of " << setw(10) << sigmaND 
     << " mb and an impact-parameter enhancement factor of "
     << setw(10) << impactFac << ".   |\n";

  // Listing finished.
  os << " |                                                            "
     << "                                                |\n"
     << " *-------  End PYTHIA Event and Cross Section Statistics -----"
     << "------------------------------------------------*" << endl;
     
}

//**************************************************************************

} // end namespace Pythia8
