// Function definitions (not found in the header) for the ProcessLevel class.
// Copyright © 2005 Torbjörn Sjöstrand

#include "ProcessLevel.h"

namespace Pythia8 {
 
//**************************************************************************

// Main routine to initialize generation process.

bool ProcessLevel::init( int idA, int idB,  double eCM,
  bool hasLHAin, LHAinit* lhaInitPtrIn, LHAevnt* lhaEvntPtrIn) {

  // Store pointers to Les Houches Accord input, if any.
  hasLHA = hasLHAin;
  lhaInitPtr = lhaInitPtrIn;
  lhaEvntPtr = lhaEvntPtrIn;
  strategyLHA = (hasLHA) ? lhaInitPtr->strategy() : 0;  

  // Initialize event generation from Pythia 6.3.
  if (!hasLHA) initPythia6( idA, idB, eCM);

  // Done. (Check return values from other classes??)
  return true;
}

//*********

// Main routine to generate the hard process.
// Currently rather primitive.

bool ProcessLevel::next( BeamParticle& beamA, BeamParticle& beamB,
  Event& process) {

  // Generate the next Pythia 6.3 event.
  if (!hasLHA) Pythia6::pyupev();

  // Read in a simple event in the LHAevnt format.
  if (strategyLHA >= 10) return getSimpleLHAevnt( process);

  // Read in an event in the LHAevnt format.
  return getLHAevnt( beamA, beamB, process);

  // Done.
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

// Read in the hard process from the Les Houches Accord.
// Many more checks to be done for valid input??

bool ProcessLevel::getLHAevnt( BeamParticle& beamA, BeamParticle& beamB,
  Event& process) {

  // Generate the next Les Houches event.
  if (!lhaEvntPtr->set()) return false;
     
  // Let hard process record begin with the event as a whole and
  // the two incoming beam particles.  
  double eCM = beamA.e() + beamB.e();
  process.append( 90, -11, 0, 0, 1, 2, 0, 0, 
    Vec4(0., 0., 0., eCM), eCM, 0. ); 
  process.append( beamA.id(), -12, 0, 0, 3, 0, 0, 0, 
    Vec4(0., 0., beamA.pz(), beamA.e()), beamA.m(), 0. ); 
  process.append( beamB.id(), -12, 0, 0, 4, 0, 0, 0, 
    Vec4(0., 0., beamB.pz(), beamB.e()), beamB.m(), 0. ); 

  // Since LHA partons may be out of order, determine correct one.
  // (Recall that zeroth particle is empty.) 
  newPos.resize(lhaEvntPtr->size());
  newPos[0] = 0;
  int nNew = 1;
  int iNew = -1;
  do {
    ++iNew;
    // For iNew == 0 look for the two incoming partons, then for
    // partons having them as mothers, and so on layer by layer.
    for (int i = 1; i < lhaEvntPtr->size(); ++i)
      if (lhaEvntPtr->mother1(i) == newPos[iNew]) newPos[nNew++] = i;
  } while (iNew < nNew && iNew < lhaEvntPtr->size()); 

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
