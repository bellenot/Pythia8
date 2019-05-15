// Function definitions (not found in the header) for the Pythia class.
// Copyright C 2006 Torbjorn Sjostrand

#include "Pythia.h"

namespace Pythia8 {
 
//**************************************************************************

// Pythia class.
// This class contains the top-level routines to generate an event.

//*********
 
// Definitions of static variables.
// (Values will be overwritten in initStatic call, so are purely dummy.)

bool Pythia::afterProcessLevel = true;
bool Pythia::afterPartonLevel = true;
bool Pythia::checkEvent = true;
int Pythia::nErrList = 3;
double Pythia::epTolerance = 1e-5;

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Maximum number of tries to produce parton level from given input..
const int Pythia::NTRY = 10; 

//*********

// Read in one update for flag/mode/parameter/Pythia6 from a single line.

bool Pythia::readString(string line, bool warn) {

  // If first character is not a letter/digit, then taken to be a comment.
  int firstChar = line.find_first_not_of(" ");
  if (!isalnum(line[firstChar])) return true; 

  // Identify substring up to colon, i.e. first part of name.
  int firstColon = line.find_first_of(":");
  string name( line, firstChar, firstColon - firstChar);
  string lowername = tolower(name);

  // Send on relevant cases to Pythia 6.
  if (lowername == "pythia6") {
    if (line[firstColon + 1] == ':') ++firstColon; 
    string pythia6command( line, firstColon + 1, 
      line.size() - firstColon - 1);
    // If no = sign in the string, then insert one at a blank.
    if (pythia6command.find("=") == string::npos) {
      int firstBlank = pythia6command.find_first_of(" ");
      pythia6command[firstBlank] = '='; 
    }         
    // No Pythia6 header, since already given Pythia8 one.
    if (!hasPythia6) Pythia6::pygive("mstu(12)=12345");
    Pythia6::pygive(pythia6command); 
    hasPythia6 = true;
    return true; 
  }

  // Send on particle data to the ParticleData database.
  if (isdigit(line[firstChar])) 
    return particleData.readString(line, warn);

  // Everything else sent on to Settings.
  return settings.readString(line, warn);

}

//*********

// Read in updates for flag/mode/parameter/Pythia6 from user-defined file.

 bool Pythia::readFile(string updateFile, bool warn) {

  // Open file with updates.
  const char* cstring = updateFile.c_str();
  ifstream is(cstring);  
  if (!is) {cout << "Error: user update file " << updateFile 
     << " not found \n"; return false;}

  // Read in one line at a time.
  bool accepted = true;
  string line;
  while ( getline(is, line) ) {

    // Process the line.
    if (!readString( line, warn)) accepted = false;

  // Reached end of input file.
  };
  return accepted;
}

//*********

// Routine to pass in pointers to PDF's. Usage optional.

bool Pythia::PDFPtr( PDF* pdfAPtrIn, PDF* pdfBPtrIn) {

  // The routine can have no effect if PDF's already assigned.
  if (pdfAPtr != 0 || pdfBPtr != 0) return false;

  // The two PDF objects cannot be one and the same, or unassigned.
  if (pdfAPtrIn == pdfBPtrIn || pdfAPtrIn == 0 || pdfBPtrIn == 0) return false;

  // Save pointers.  
  pdfAPtr = pdfAPtrIn;
  pdfBPtr = pdfBPtrIn;
  
  return true;
}

//*********

// Routine to initialize with two beam energies specified.

bool Pythia::init( int idAin, int idBin, double eAin, double eBin) {

  // Read in and set values.
  idA = idAin;
  idB = idBin;
  eA = eAin;
  eB = eBin;

  // Send on to common initialization.
  return init(false);
}

//*********

// Routine to initialize with CM energy rather tham beam energies.

bool Pythia::init( int idAin, int idBin, double eCMin) {

  // Read in and set values.
  idA = idAin;
  idB = idBin;
  eCM = eCMin;

  // Send on to common initialization.
  return init(true);
}

//*********

// Routine to initialize when beam info is given in an LHAinit object.

bool Pythia::init( LHAinit* lhaInitPtrIn, LHAevnt* lhaEvntPtrIn) {

  // Save and set flag for subsequent usage of LHAevnt object.
  lhaInitPtr = lhaInitPtrIn;
  lhaEvntPtr = lhaEvntPtrIn;
  hasLHA = true;

  // Set LHAinit information (in some external program).
  lhaInitPtr->set();

  // Extract beams from values set in an LHAinit object. 
  idA = lhaInitPtr->idBeamA();
  idB = lhaInitPtr->idBeamB();
  eA = lhaInitPtr->eBeamA();
  eB = lhaInitPtr->eBeamB();

  // Now do normal initialization.
  return init(false);

}

//*********

// Main routine to initialize the generation process.
// (The alternative init forms end up in this one.)

bool Pythia::init(bool inCMframeIn) {

  // Reset error counter.
  nErrEvent = 0;

  // Initialize all accessible static data members.
  initStatic();

  // Set headers to distinguish the two event listing kinds.
  process.header("(hard process)");
  event.header("(complete event)");

  // Strategy for event generation, Les Houches or not.
  strategyLHA = (hasLHA) ? lhaInitPtr->strategy() : 0;  

  // Do not set up beam kinematics when no parton-level processing.
  if (strategyLHA >= 10) inCMframe = true;
  else {

    // Find masses. Initial guess about CM frame.
    mA = ParticleDataTable::m0(idA);
    mB = ParticleDataTable::m0(idB);
    inCMframe = inCMframeIn;
    betaZ = 0.;
    gammaZ = 1.;

    // When not given: find CM energy and set up boost to rest frame.
    if (!inCMframeIn) {
      eA = max(eA, mA);
      eB = max(eB, mB);
      pzA = sqrt(eA*eA - mA*mA);
      pzB = -sqrt(eB*eB - mB*mB);
      eCM = sqrt( pow2(eA + eB) - pow2(pzA + pzB) );
      betaZ = (pzA + pzB) / (eA + eB);
      gammaZ = (eA + eB) / eCM;
      if (abs(betaZ) < 1e-10) inCMframe = true;
    }

    // Set up kinematics in the rest frame.
    if (eCM < mA + mB) return false;
    pzA = 0.5 * sqrtpos( (eCM + mA + mB) * (eCM - mA - mB) 
      * (eCM - mA + mB) * (eCM + mA - mB) ) / eCM;
    pzB = -pzA;
    eA = sqrt(mA*mA + pzA*pzA);
    eB = sqrt(mB*mB + pzB*pzB);

    // Set up the PDF's, if not already done.
    if (pdfAPtr == 0) {pdfAPtr = setPDFPtr(idA); pdfAnew = true;}
    if (pdfBPtr == 0) {pdfBPtr = setPDFPtr(idB); pdfBnew = true;}
  
    // Set up the two beams and the common remnant system.
    beamA.init( idA, pzA, eA, mA, pdfAPtr);
    beamB.init( idB, pzB, eB, mB, pdfBPtr);
  }

  // Store main info for access in process generation.
  info.setBeamA( idA, pzA, eA, mA);
  info.setBeamB( idB, pzB, eB, mB);
  info.setECM( eCM);

  // Set info in and send pointers to the respective program elements.
  if (!processLevel.init( &info, &beamA, &beamB, hasPythia6, 
    hasLHA, lhaInitPtr, lhaEvntPtr)) return false;
  if (!partonLevel.init( &info, &beamA, &beamB, strategyLHA))
    return false;
  if (!hadronLevel.init( &info)) return false;

  // Succeeded.
  isInit = true; 
  return true;
}

//*********

// Initialization routine for all accessible static data members.

void Pythia::initStatic() {

  // Initialize the random number generator.
  if ( settings.flag("Pythia:setSeed") )  
    Rndm::init( settings.mode("Pythia:seed") );

  // Initialize static data members in the Pythia class.
  afterProcessLevel = settings.flag("Pythia:afterProcessLevel");
  afterPartonLevel = settings.flag("Pythia:afterPartonLevel");
  checkEvent = settings.flag("Pythia:checkEvent");
  
  // Initialize all other accessible static data members.
  ErrorMessages::initStatic();
  ParticleDataEntry::initStatic();
  Event::initStatic();
  BeamParticle::initStatic(); 
  BeamRemnants::initStatic();
  InFlux::initStatic(); 
  SigmaProcess::initStatic(); 
  SigmaTotal::initStatic();
  PhaseSpace::initStatic();
  TimeShower::initStatic();
  SpaceShower::initStatic();
  MultipleInteractions::initStatic();
  HadronLevel::initStatic();
  StringFlav::initStatic();
  StringZ::initStatic();
  StringPT::initStatic();
  ColConfig::initStatic();
  StringRegion::initStatic();
  StringFragmentation::initStatic();
  MiniStringFragmentation::initStatic();
  ParticleDecays::initStatic(); 

}

//*********

// Main routine to generate the next event, using internal machinery.

bool Pythia::next() {

  // Can only generate event if initialization worked.
  if (!isInit) {
    ErrorMessages::message("Error in Pythia::next: "
      "not properly initialized so cannot generate events"); 
    return false;
  }

  // Provide the hard process that starts it off. Only one try.
  info.clear();
  process.clear();
  if ( !processLevel.next( process) ) {
    ErrorMessages::message("Error in Pythia::next: "
      "processLevel failed; giving up"); 
    return false;
  }

  // Possibility to stop the generation at this stage.
  if (!afterProcessLevel) {
    if (!inCMframe) process.bst(0., 0., betaZ, gammaZ);
    return true;
  }
  
  // Allow up to ten tries for parton- and hadron-level processing.
  bool physical = true;
  for (int iTry = 0; iTry < NTRY; ++ iTry) {
    physical = true;

    // Reset event record and (extracted partons from) beam remnants.
    event.clear();
    beamA.clear();
    beamB.clear();
   
    // Parton-level evolution: ISR, FSR, MI.
    if ( !partonLevel.next( process, event) ) {
      ErrorMessages::message("Warning in Pythia::next: "
        "partonLevel failed; try again"); 
      physical = false; 
      continue;
    }

    // When required: boost to lab frame (before decays, for vertices).
    if (!inCMframe) {
      process.bst(0., 0., betaZ, gammaZ);
      event.bst(0., 0., betaZ, gammaZ);      
    }    

    // Possibility to stop the generation at this stage.
    if (!afterPartonLevel) return physical;

    // Hadron-level: hadronization, decays.
    if ( !hadronLevel.next( event) ) {
      ErrorMessages::message("Warning in Pythia::next: "
        "hadronLevel failed; try again"); 
      physical = false; 
      continue;
    }

    // Stop looping, for better or worse.
    break;
  }
  if (!physical) {
    ErrorMessages::message("Error in Pythia::next: "
      "parton+hadronLevel failed; giving up");
    // Debug. To remove??
    info.list();
    process.list();
    event.list(); 
    return false;
  }

  // Optionally check final event for problems. Done
  if (checkEvent) return check();
  return true;

}

//*********

// Print statistics on event generation.

void Pythia::statistics() {

  // Statistics on cross section and number of events. 
  processLevel.statistics();  

  // Summary of which and how many warnings/errors encountered.
  ErrorMessages::statistics();

  // Possibility for statistics from other classes, e.g. for debug.
  partonLevel.statistics();  

}

//*********

// Write the Pythia banner, with symbol and version information.

void Pythia::banner() {

  // Read in version number and last date of change.
  double versionNumber = Settings::parameter("Pythia:versionNumber");
  int versionDate = Settings::mode("Pythia:versionDate");
  string month[12] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};  

  // Get date and time.
  time_t t = time(0);
  char dateNow[12];
  strftime(dateNow,12,"%d %b %Y",localtime(&t));
  char timeNow[9];
  strftime(timeNow,9,"%H:%M:%S",localtime(&t));

  
  cout << "\n"
       << " *-------------------------------------------" 
       << "-----------------------------------------* \n"
       << " |                                           "
       << "                                         | \n"
       << " |  *----------------------------------------" 
       << "--------------------------------------*  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |                *......*                "
       << "  Welcome to the Lund Monte Carlo!    |  | \n" 
       << " |  |           *:::!!:::::::::::*           " 
       << "                                      |  | \n"
       << " |  |        *::::::!!::::::::::::::*        "
       << "  PPP  Y   Y TTTTT H   H III   A      |  | \n"
       << " |  |      *::::::::!!::::::::::::::::*      "
       << "  P  P  Y Y    T   H   H  I   A A     |  | \n"
       << " |  |     *:::::::::!!:::::::::::::::::*     "
       << "  PPP    Y     T   HHHHH  I  AAAAA    |  | \n"
       << " |  |     *:::::::::!!:::::::::::::::::*     "
       << "  P      Y     T   H   H  I  A   A    |  | \n"
       << " |  |      *::::::::!!::::::::::::::::*!     "
       << "  P      Y     T   H   H III A   A    |  | \n"
       << " |  |        *::::::!!::::::::::::::* !!     "
       << "                                      |  | \n"
       << " |  |        !! *:::!!:::::::::::*    !!     "
       << "  This is PYTHIA version " << fixed << setprecision(3) 
       << setw(5) << versionNumber << "        |  | \n"
       << " |  |        !!     !* -><- *         !!     "
       << "  Last date of change: " << setw(2) << versionDate%100 
       << " " << month[ (versionDate/100)%100 - 1 ] 
       << " " << setw(4) << versionDate/10000 <<  "    |  | \n"
       << " |  |        !!     !!                !!     " 
       << "                                      |  | \n"
       << " |  |        !!     !!                !!     "
       << "  Now is " << dateNow << " at " << timeNow << "      |  | \n"
       << " |  |        !!                       !!     " 
       << "                                      |  | \n"
       << " |  |        !!        lh             !!     " 
       << "  Disclaimer: this program comes      |  | \n"
       << " |  |        !!                       !!     "
       << "  without any guarantees. Beware      |  | \n"
       << " |  |        !!                 hh    !!     "
       << "  of errors and use common sense      |  | \n"
       << " |  |        !!    ll                 !!     " 
       << "  when interpreting results.          |  | \n"
       << " |  |        !!                       !!     " 
       << "                                      |  | \n"
       << " |  |        !!                              " 
       << "  Copyright T. Sjostrand (2006)       |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |   An archive of program versions and do" 
       << "cumentation is found on the web:      |  | \n"
       << " |  |   http://www.thep.lu.se/~torbjorn/Pythi" 
       << "a.html                                |  | \n"
       << " |  |                                        " 
       << "                                      |  | \n"
       << " |  |   The current version is intended for t" 
       << "ryout and feedback only,              |  | \n"
       << " |  |   and should not be used for any physic" 
       << "s studies or production runs.         |  | \n"
       << " |  |                                        " 
       << "                                      |  | \n"
       << " |  |   Currently PYTHIA 6.4 can be used to g"
       << "enerate processes, see the manual     |  | \n"
       << " |  |   T. Sjostrand, S. Mrenna and P. Skands"
       << ", JHEP05 (2006) 026 [hep-ph/0603175]. |  | \n"
       << " |  |                                        " 
       << "                                      |  | \n"
       << " |  |   Main author: Torbjorn Sjostrand; CERN" 
       << "/PH, CH-1211 Geneva, Switzerland,     |  | \n"
       << " |  |     and Department of Theoretical Physi"
       << "cs, Lund University, Lund, Sweden;    |  | \n"
       << " |  |     phone: + 41 - 22 - 767 82 27; e-mai"
       << "l: torbjorn@thep.lu.se                |  | \n"
       << " |  |   Author: Stephen Mrenna; Computing Div"
       << "ision, Simulations Group,             |  | \n"
       << " |  |     Fermi National Accelerator Laborato"
       << "ry, MS 234, Batavia, IL 60510, USA;   |  | \n"
       << " |  |     phone: + 1 - 630 - 840 - 2556; e-ma"
       << "il: mrenna@fnal.gov                   |  | \n"
       << " |  |   Author: Peter Skands; Theoretical Phy"
       << "sics Department,                      |  | \n"
       << " |  |     Fermi National Accelerator Laborato"
       << "ry, MS 106, Batavia, IL 60510, USA;   |  | \n"
       << " |  |     phone: + 1 - 630 - 840 - 2270; e-ma"
       << "il: skands@fnal.gov                   |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  |                                        "
       << "                                      |  | \n"
       << " |  *----------------------------------------" 
       << "--------------------------------------*  | \n"
       << " |                                           "
       << "                                         | \n"
       << " *-------------------------------------------" 
       << "-----------------------------------------* " << endl;

}

//*********

// Check that the final event makes sense: no unknown id codes;
// charge and energy-momentum conserved.

bool Pythia::check() {

  // Reset. Incoming beams counted with negative momentum and charge.
  bool physical = true;
  iErrId.resize(0);
  iErrNan.resize(0);
  Vec4 pSum = - (event[1].p() + event[2].p());
  double eLab = abs(pSum.e());
  double chargeSum = - (event[1].charge() + event[2].charge());

  // Loop over particles in the event. 
  for (int i = 0; i < event.size(); ++i) {

    // Look for any unrecognized particle codes.
    int id = event[i].id();
    if (id == 0 || !ParticleDataTable::isParticle(id)) {
      string errCode;
      ostringstream writeCode(errCode);
      writeCode << ", id = " << id;
      ErrorMessages::message("Error in Pythia::check: "
        "unknown particle code", errCode); 
      physical = false;
      iErrId.push_back(i);
    }

    // Look for particle with not-a-number energy/momentum/mass.
    if (abs(event[i].px()) >= 0. && abs(event[i].py()) >= 0. 
      && abs(event[i].pz()) >= 0.  && abs(event[i].e()) >= 0. 
      && abs(event[i].m()) >= 0.) ;
    else {   
      ErrorMessages::message("Error in Pythia::check: "
        "not-a-number energy/momentum/mass"); 
      physical = false;
      iErrNan.push_back(i);
    }

    // Add final-state four-momentum and charge.      
    if (event[i].remains()) {
      pSum += event[i].p();
      chargeSum += event[i].charge();
    }

  // End of particle loop.
  }

  // Check energy-momentum/charge conservation.
  double epDev = abs(pSum.e()) + abs(pSum.px()) + abs(pSum.py())
    + abs(pSum.pz());
  if (epDev > epTolerance * eLab) { 
    ErrorMessages::message("Error in Pythia::check: "
      "energy-momentum badly conserved"); 
    physical = false;
  }
  if (abs(chargeSum) > 0.1) {
    ErrorMessages::message("Error in Pythia::check: "
      "charge not conserved"); 
    physical = false;
  }

  // Done for sensible events.
  if (physical) return true;

  // Print (the first few) flawed events.
  if (nErrEvent < nErrList) {
    cout << " Erroneous event info: \n";
    if (iErrId.size() > 0) {
      cout << " unknown particle codes in lines ";
      for (int i = 0; i < int(iErrId.size()); ++i) 
        cout << iErrId[i] << " ";
      cout << "\n";
    }
    if (iErrNan.size() > 0) {
      cout << " not-a-number energy/momentum/mass in lines ";
      for (int i = 0; i < int(iErrNan.size()); ++i) 
        cout << iErrNan[i] << " ";
      cout << "\n";
    }
    if (epDev > epTolerance * eLab) cout << scientific << setprecision(3)
      << " total energy-momentum non-conservation = " << epDev << "\n";
    if (abs(chargeSum) > 0.1) cout << fixed << setprecision(2) 
      << " total charge non-conservation = " << chargeSum << "\n"; 
    event.list();
  }

  // Update error counter. Done also for flawed event.
  ++nErrEvent;
  return false;

}

//*********

// Routine to set up a PDF pointer.

PDF* Pythia::setPDFPtr(int idIn) {

  PDF* tempPDFPtr;

  // Proton beam.
  if (abs(idIn) == 2212) {
    int pPDFset = settings.mode("Pythia:pPDFset");
    if (pPDFset == 1) tempPDFPtr = new GRV94L(idIn);
    else tempPDFPtr = new CTEQ5L(idIn);

  // Lepton beam.
    } else tempPDFPtr = new Lepton(idIn);
  
  return tempPDFPtr; 
}

//*********

} // end namespace Pythia8

