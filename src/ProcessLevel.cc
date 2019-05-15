// Function definitions (not found in the header) for the ProcessLevel class.
// Copyright C 2007 Torbjorn Sjostrand

#include "ProcessLevel.h"

namespace Pythia8 {
 
//**************************************************************************

// Main routine to initialize generation process.

bool ProcessLevel::init( Info* infoPtrIn, BeamParticle* beamAPtrIn, 
  BeamParticle* beamBPtrIn, bool hasLHAin, LHAinit* lhaInitPtrIn, 
  LHAevnt* lhaEvntPtrIn, UserHooks* userHooksPtrIn) {

  // Store input pointers for future use. 
  infoPtr      = infoPtrIn;
  beamAPtr     = beamAPtrIn;
  beamBPtr     = beamBPtrIn;
  userHooksPtr = userHooksPtrIn;

  // Store pointers to Les Houches Accord input if any.
  hasLHA       = hasLHAin;
  lhaInitPtr   = lhaInitPtrIn;
  lhaEvntPtr   = lhaEvntPtrIn;
  strategyLHA  = (hasLHA) ? lhaInitPtr->strategy() : 0;  

  // If not Les Houches then internal machinery.
  hasInternal  = !hasLHA;
  if (hasInternal) return initInternal();

  // Done. (Check return values from other classes??)
  return true;
}

//*********

// Main routine to generate the hard process.
// Currently rather primitive.

  bool ProcessLevel::next( Event& process) {

  // Starting value.
  bool physical = false;  

  // Generate the next internal event. 
  if (hasInternal) physical = nextInternal( process);

  // Read in a simple event in the LHAevnt format. Default info.
  else if (strategyLHA >= 10) {
    infoPtr->setType( "Simple LHA process", 0, 0, false, true, false, false);
    infoPtr->setTypeMI( 0, 0.);
    physical = nextSimpleLHA( process);
  }

  // Read in an event in the LHAevnt format.
  else if (hasLHA) physical = nextLHA( process);

  // Check that colour assignments make sense.
  if (physical) physical = checkColours( process);

  // Done.
  return physical;
}

//*********

// Accumulate and update statistics (after possible user veto).
  
void ProcessLevel::accumulate() {

  // Currently does not handle LHA processes ??
  if (hasLHA) return;

  // Increase number of accepted events.
  containerPtrs[iNow]->accumulate();

  // Provide current generated cross section estimate.
  long   nTrySum   = 0; 
  long   nSelSum   = 0; 
  long   nAccSum   = 0;
  double sigmaSum  = 0.;
  double delta2Sum = 0.;
  for (int i = 0; i < int(containerPtrs.size()); ++i) 
  if (containerPtrs[i]->sigmaMax() != 0.) {
    nTrySum       += containerPtrs[i]->nTried();
    nSelSum       += containerPtrs[i]->nSelected();
    nAccSum       += containerPtrs[i]->nAccepted();
    sigmaSum      += containerPtrs[i]->sigmaMC();
    delta2Sum     += pow2(containerPtrs[i]->deltaMC()); 
  }
  infoPtr->setSigma( nTrySum, nSelSum, nAccSum, sigmaSum, sqrt(delta2Sum)); 

}

//*********

// Print statistics on cross sections and number of events.

void ProcessLevel::statistics(ostream& os) {

  // Internal statistics.
  if (hasInternal) {
    
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

}

//*********

// Initialize the internal event generation machinery.
  
bool ProcessLevel::initInternal( ostream& os) {

  // Set up SigmaTotal.
  int    idA = infoPtr->idA();
  int    idB = infoPtr->idB();
  double eCM = infoPtr->eCM();
  sigmaTot.init( idA, idB, eCM);

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
      "no processes switched on"); 
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
  for (int i = 0; i < int(containerPtrs.size()); ++i) {
    os << " | " << left << setw(40) << containerPtrs[i]->name() 
       << right << setw(5) << containerPtrs[i]->code() << " | " 
       << scientific << setprecision(3) << setw(11)  
       << containerPtrs[i]->sigmaMax() << " |\n";
  }

  // Listing finished.
  os << " |                                                             |\n" 
     << " *-------  End PYTHIA Process Initialization ------------------*" 
     << endl;

  // If sum of maxima vanishes then refuse to do anything. Done.
  if ( numberOn == 0  || sigmaMaxSum <= 0.) {
    ErrorMsg::message("Error in ProcessLevel::initInternal: "
      "all processes have vanishing cross sections"); 
    return false;
  }
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
    iNow = -1;
    do sigmaMaxNow -= containerPtrs[++iNow]->sigmaMax();
    while (sigmaMaxNow > 0. && iNow < iMax);
    
    // Do a trial event of this subprocess; accept or not.
    if (containerPtrs[iNow]->trialProcess()) break;
  }

  // Construct kinematics of acceptable process.
  containerPtrs[iNow]->constructProcess( process);

  // Do all resonance decays. First draft??
  resonanceDecays.next( process);
   
  // Add any junctions to the process event record list.
  findJunctions( process);

  // Done.
  return true;
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

  // Read info om parton densities if provided.
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
      if ( (iNow == 3 && iNowA == 4) || (iNow == 4 && iNowA == 3) );
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

//**************************************************************************

} // end namespace Pythia8
