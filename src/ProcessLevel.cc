// Function definitions (not found in the header) for the ProcessLevel class.
// Copyright C 2006 Torbjorn Sjostrand

#include "ProcessLevel.h"

namespace Pythia8 {
 
//**************************************************************************

// Main routine to initialize generation process.

bool ProcessLevel::init( Info* infoPtrIn, BeamParticle* beamAPtrIn, 
  BeamParticle* beamBPtrIn, bool hasPythia6In, bool hasLHAin, 
  LHAinit* lhaInitPtrIn, LHAevnt* lhaEvntPtrIn) {

  // Store input pointers for future use. 
  infoPtr = infoPtrIn;
  beamAPtr = beamAPtrIn;
  beamBPtr = beamBPtrIn;

  // Store pointers to Les Houches Accord input if any.
  hasLHA = hasLHAin;
  lhaInitPtr = lhaInitPtrIn;
  lhaEvntPtr = lhaEvntPtrIn;
  strategyLHA = (hasLHA) ? lhaInitPtr->strategy() : 0;  

  // Initialize event generation from Pythia 6.3.
  hasPythia6 = hasPythia6In;
  if (hasPythia6) initPythia6( infoPtr->idA(), infoPtr->idB(), 
    infoPtr->eCM() );

  // If not Les Houches or Pythia 6.3 then internal machinery.
  hasInternal = !hasLHA && !hasPythia6;
  if (hasInternal) return initInternal();

  // Done. (Check return values from other classes??)
  return true;
}

//*********

// Main routine to generate the hard process.
// Currently rather primitive.

  bool ProcessLevel::next( Event& process) {

  // Generate the next internal event. 
  if (hasInternal) return getInternalEvnt( process);

  // Generate the next Pythia 6.4 event. 
  if (hasPythia6) Pythia6::pyupev();

  // Read in a simple event in the LHAevnt format. Default info (to be improved??).
  if (strategyLHA >= 10) {
    infoPtr->setType( "Simple LHA process", 0, 0, false, true, false, false);
    return getSimpleLHAevnt( process);
  }

  // Read in an event in the LHAevnt format. Default info (to be improved??).
  if (hasPythia6 || hasLHA) {
    string name = (hasPythia6) ? "Pythia6 LHA process" : "External LHA process"; 
    infoPtr->setType( name, 0, 0, false, true, false, false);
    return getLHAevnt( process);
  }

  // Done.
  return true;
}

//*********

// Print statistics on cross sections and number of events.

void ProcessLevel::statistics(ostream& os) {

  // Internal statistics.
  if (hasInternal) {
    
    // Header.
    os << "\n *-------  PYTHIA Event and Cross Section Statistics  ------"
       << "---------------------------------------*\n"
       << " |                                                            "
       << "                                     |\n" 
       << " | Subprocess                               Code |       Numbe"
       << "r of events |      sigma +- delta    |\n" 
       << " |                                               |       Tried"
       << "   Accepted |     (estimated) (mb)   |\n"
       << " |                                               |            "
       << "            |                        |\n"
       << " |------------------------------------------------------------"
       << "-------------------------------------|\n"
       << " |                                               |            "
       << "            |                        |\n";

    // Reset sum counters.
    int nTrySum = 0; 
    int nAccSum = 0;
    double sigmaSum = 0.;
    double delta2Sum = 0.;

    // Loop over existing processes.
    for (int i = 0; i < int(containerPtrs.size()); ++i) 
    if (containerPtrs[i]->sigmaMax() != 0.) {

      // Read info for process. Sum counters.
      int nTry = containerPtrs[i]->nTried();
      int nAcc = containerPtrs[i]->nAccepted();
      double sigma = containerPtrs[i]->sigmaMC();
      double delta = containerPtrs[i]->deltaMC(); 
      nTrySum += nTry;
      nAccSum += nAcc; 
      sigmaSum += sigma;
      delta2Sum += pow2(delta);    

      // Print individual process info.
      os << " | " << left << setw(40) << containerPtrs[i]->name() 
         << right << setw(5) << containerPtrs[i]->code() << " | " 
         << setw(11) << nTry << " " << setw(10) << nAcc << " | " 
         << scientific << setprecision(3) << setw(11) << sigma 
         << setw(11) << delta << " |\n";
    }

    // Print summed process info.
    os << " |                                               |            "
       << "            |                        |\n"
       << " | " << left << setw(45) << "sum" << right << " | " << setw(11) 
       << nTrySum << " " << setw(10) << nAccSum << " | " << scientific 
       << setprecision(3) << setw(11) << sigmaSum << setw(11) 
       << sqrt(max(0., delta2Sum)) << " |\n";

    // Listing finished.
    os << " |                                                            "
       << "                                     |\n"
       << " *-------  End PYTHIA Event and Cross Section Statistics -----"
       << "-------------------------------------*" << endl;
     
  }

  // Pythia6 statistics.
  if (hasPythia6) Pythia6::pystat(0);

}

//*********

// Initialize the internal event generation machinery.
  
bool ProcessLevel::initInternal( ostream& os) {

  // Set up SigmaTotal.
  int idA = infoPtr->idA();
  int idB = infoPtr->idB();
  double eCM = infoPtr->eCM();
  sigmaTot.init( idA, idB, eCM);

  // Sets up containers for all the hard processes.
  SetupContainers setupContainers;
  setupContainers.init(containerPtrs);

  // If no processes found then refuse to do anything.
  if ( int(containerPtrs.size()) == 0) {
    ErrorMessages::message("Error in ProcessLevel::initInternal: "
      "no processes switched on"); 
    return false;
  }

  // Initialize each process. Special case for minbias/elastic/diffractive.
  for (int i = 0; i < int(containerPtrs.size()); ++i) {
    if (containerPtrs[i]->needsSigmaTotal()) 
      containerPtrs[i]->setSigmaTotalPtr(&sigmaTot);
    containerPtrs[i]->init(infoPtr, beamAPtr, beamBPtr);
  }

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
  if ( sigmaMaxSum <= 0.) {
    ErrorMessages::message("Error in ProcessLevel::initInternal: "
      "all processes have vanishing cross sections"); 
    return false;
  }
  return true;

}

//*********

// Initialize event generation from Pythia 6.3.

void ProcessLevel::initPythia6( int idA, int idB, double eCM) {

  // Simplified translation into input accessible by Pythia 6.3,
  // assuming we are in the cm frame, and only have p/pbar/e+/e- beams.
  string frame = "cms";
  string beam = "p";
  if (idA == -2212) beam = "pbar";
  if (idA == 11) beam = "e-";
  if (idA == -11) beam = "e+";
  string target = "p";
  if (idB == -2212) target = "pbar";
  if (idB == 11) target = "e-";
  if (idB == -11) target = "e+";

  // Initialize Pythia 6.3, without multiple interactions.
  Pythia6::pygive("mstp(81)=0");
  Pythia6::pyinit(frame, beam, target, eCM);

  // Create pointer to LHAevnt. Warning: delete it somewhere??
  lhaEvntPtr = new LHAevntFortran;  

  // Done.
}

//*********

// Generate the next internal event.
  
bool ProcessLevel::getInternalEvnt( Event& process) {

  // Loop over tries until trial event succeeds.
  int iNow;
  for ( ; ; ) {

    // Pick one of the subprocesses.
    double sigmaMaxNow = sigmaMaxSum * Rndm::flat();
    iNow = -1;
    do sigmaMaxNow -= containerPtrs[++iNow]->sigmaMax();
    while (sigmaMaxNow > 0. && iNow < int(containerPtrs.size()) - 1);
    
    // Do a trial event of this subprocess; accept or not.
    if (containerPtrs[iNow]->trialProcess()) break;
  }

  // Construct kinematics of acceptable process.
  containerPtrs[iNow]->constructProcess( process);

  // Do all resonance decays. First draft??
  resonanceDecays( process);
   
  // Add any junctions to the process event record list.
  findJunctions( process);

  // Done.
  return true;
}

//*********

// Do all resonance decays. First draft??
  
bool ProcessLevel::resonanceDecays( Event& process) {

  // Loop over all entries to find resonances that should decay.
  int iDec = 0;
  do {
    Particle& decayer = process[iDec];
    if (decayer.remains() && decayer.canDecay() && decayer.mayDecay() 
    && decayer.isResonance() ) {

      // Particle data for decaying particle.
      int id0 = decayer.id();
      double m0 = decayer.m();

      // Pick a decay channel; allow up to ten tries.
      int NTRYDECAY = 10;
      double mSafety = 1.; 
      int id1 = 0;
      int id2 = 0;
      double m1 = 0.;
      double m2 = 0.;
      for (int iTryChannel = 0; iTryChannel < NTRYDECAY; ++iTryChannel) {
        DecayChannel& channel = decayer.particleData().decay.pick();
        // int mode = channel.modeME();
        int mult = channel.multiplicity();

        // Consider for now only two-body decay. Check phase space. 
        if (mult != 2) continue;
        id1 = channel.product(0);
        if (id0 < 0 && ParticleDataTable::hasAnti(id1)) id1 = -id1;
        id2 = channel.product(1);
        if (id0 < 0 && ParticleDataTable::hasAnti(id2)) id2 = -id2;
        m1 = ParticleDataTable::mass(id1);          
        m2 = ParticleDataTable::mass(id2); 
        if (m1 + m2 + mSafety > m0) continue;

      // End of loop over tries.
      }

      // Energies and absolute momentum in the rest frame.
      double e1 = 0.5 * (m0*m0 + m1*m1 - m2*m2) / m0;
      double e2 = 0.5 * (m0*m0 + m2*m2 - m1*m1) / m0;
      double pAbs = 0.5 * sqrtpos( (m0 - m1 - m2) * (m0 + m1 + m2)
        * (m0 + m1 - m2) * (m0 - m1 + m2) ) / m0;  

      // Isotropic angles give three-momentum.
      double cosTheta = 2. * Rndm::flat() - 1.;
      double sinTheta = sqrt(1. - cosTheta*cosTheta);
      double phi = 2. * M_PI * Rndm::flat();
      double pX = pAbs * sinTheta * cos(phi);  
      double pY = pAbs * sinTheta * sin(phi);  
      double pZ = pAbs * cosTheta;  

      // Fill four-momenta and boost them away from mother rest frame.
      Vec4 p1( pX, pY, pZ, e1);
      Vec4 p2( -pX, -pY, -pZ, e2);
      p1.bst( decayer.p() );
      p2.bst( decayer.p() );

      // Find colours.
      int col = (id1 > 0 && id1 < 9) ? process.nextColTag() : 0; 

      // Append decay products to the event record.
     process.append( id1, 23, iDec, 0, 0, 0, col, 0, p1, m1, m0);
     process.append( id2, 23, iDec, 0, 0, 0, 0, col, p2, m2, m0);

     // Modify mother status to show it is a decayed resonance.
     decayer.status(-22);
                 
    // End of loop over all entries.
    }
  } while (++iDec < process.size());

  // Done.
  return true;
}

//*********

// Read in the hard process from the Les Houches Accord.
// Many more checks to be done for valid input??

bool ProcessLevel::getLHAevnt( Event& process) {

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
    double e = lhaEvntPtr->e(iOld);  
    double m = lhaEvntPtr->m(iOld);

    // For resonance decay products use resonance mass as scale.
    double scaleNow = scale;
    if (mother1 > 4) scaleNow = process[mother1].m();
    process.append( id, status, mother1, mother2, daughter1, daughter2, 
      col1, col2, Vec4(px, py, pz, e), m, scaleNow);
  }  

  // Add any junctions to the process event record list.
  findJunctions( process);

  // Done.
  return true;
}

//*********

// Read in the hard process, special case if all partons already given.

bool ProcessLevel::getSimpleLHAevnt( Event& process) {

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
    double e = lhaEvntPtr->e(i);  
    double m = lhaEvntPtr->m(i);

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
      int colDau = process[ daughters[j] ].col();
      int acolDau = process[ daughters[j] ].acol();
      if (colDau > 0) cols.push_back( colDau);      
      if (acolDau > 0) acols.push_back( acolDau);      
    }

    // Remove all matching colour-anticolour pairs.
    bool foundPair = true;
    while (foundPair && cols.size() > 0 && acols.size() > 0) {
      foundPair = false;
      for (int j = 0; j < int(cols.size()); ++j) {
        for (int k = 0; k < int(acols.size()); ++k) {
	  if (acols[k] == cols[j]) { 
            cols[j] = cols.back(); cols.pop_back();     
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

//**************************************************************************

} // end namespace Pythia8
