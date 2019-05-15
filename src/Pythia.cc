
// Function definitions (not found in the header) for the Pythia class.
// Copyright C 2007 Torbjorn Sjostrand

#include "Pythia.h"

// Access time information.
#include <ctime>

// Allow string and character manipulation.
#include <cctype>

namespace Pythia8 {
 
//**************************************************************************

// Pythia class.
// This class contains the top-level routines to generate an event.
//*********
 
// Definitions of static variables.
// (Values will be overwritten in initStatic call, so are purely dummy.)

bool   Pythia::doPartonLevel = true;
bool   Pythia::doHadronLevel = true;
bool   Pythia::checkEvent    = true;
int    Pythia::nErrList      = 3;
double Pythia::epTolerance   = 1e-5;

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Maximum number of tries to produce parton level from given input.
const int Pythia::NTRY       = 10; 

//*********
  
// Constructor. 

Pythia::Pythia(string xmlDir) { 
    
  // Initial values for pointers to PDF's.
  pdfAnew       = false; 
  pdfBnew       = false;
  pdfHardNew    = false;
  pdfAPtr       = 0; 
  pdfBPtr       = 0; 
  pdfHardAPtr   = 0; 
  pdfHardBPtr   = 0; 

  // Initial value for pointer to user hooks.
  userHooksPtr  = 0;

  // Initial values for pointers to timelike and spacelike showers.
  timesNew      = false;
  spaceNew      = false;
  timesDecPtr   = 0;
  timesPtr      = 0;
  spacePtr      = 0;

  // Find path to data files, i.e. xmldoc directory location.
  // Environment variable takes precedence, else use constructor input. 
  string path = "";
  const char* PYTHIA8DATA = "PYTHIA8DATA"; 
  char* envPath = getenv(PYTHIA8DATA);
  if (envPath != 0 && *envPath != '\0') {
    int i = 0;
    while (*(envPath+i) != '\0') path += *(envPath+(i++)); 
  } 
  else path = xmlDir;
  if (path[ path.length() - 1 ] != '/') path += "/";

  // Read in files with all flags, modes, parms and words.
  string initFile = path + "Index.xml";
  isConstructed = settings.init( initFile);
  if (!isConstructed) { 
    ErrorMsg::message("Abort from Pythia::Pythia: "
    "settings unavailable");
    return;
  }

  // Read in files with all particle data.
  string dataFile = path + "ParticleData.xml";
  isConstructed = particleData.init( dataFile);
  if (!isConstructed) {
    ErrorMsg::message("Abort from Pythia::Pythia: "
    "particle data unavailable");
    return;
  }

  // Write the Pythia banner to output. 
  banner();

  // Default for some flags.
  doLHA         = false;
  isInit        = false;
 
} 

//*********
  
// Destructor.

Pythia::~Pythia() { 

  // Delete the PDF's created with new.
  if (pdfHardNew && pdfHardAPtr != pdfAPtr) delete pdfHardAPtr; 
  if (pdfHardNew && pdfHardBPtr != pdfBPtr) delete pdfHardBPtr; 
  if (pdfAnew) delete pdfAPtr; 
  if (pdfBnew) delete pdfBPtr; 

  // Delete the timelike and spacelike showers created with new.
  if (timesNew) delete timesPtr;
  if (spaceNew) delete spacePtr;

} 

//*********

// Read in one update for a setting or particle data from a single line.

bool Pythia::readString(string line, bool warn) {

  // Check that constructor worked.
  if (!isConstructed) return false;

  // If first character is not a letter/digit, then taken to be a comment.
  int firstChar = line.find_first_not_of(" ");
  if (!isalnum(line[firstChar])) return true; 

  // Send on particle data to the ParticleData database.
  if (isdigit(line[firstChar])) 
    return particleData.readString(line, warn);

  // Everything else sent on to Settings.
  return settings.readString(line, warn);

}

//*********

// Read in updates for settings or particle data from user-defined file.

 bool Pythia::readFile(string updateFile, bool warn) {

  // Check that constructor worked.
  if (!isConstructed) return false;

  // Open file with updates.
  const char* cstring = updateFile.c_str();
  ifstream is(cstring);  
  if (!is) {
    ErrorMsg::message("Error in Pythia::readFile: did not find file",
      updateFile);
    return false;
  }

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

bool Pythia::setPDFPtr( PDF* pdfAPtrIn, PDF* pdfBPtrIn, PDF* pdfHardAPtrIn, 
    PDF* pdfHardBPtrIn) {

  // The routine can have no effect if PDF's already assigned.
  if (pdfAPtr != 0 || pdfBPtr != 0) return false;

  // The two PDF objects cannot be one and the same, or unassigned.
  if (pdfAPtrIn == pdfBPtrIn || pdfAPtrIn == 0 || pdfBPtrIn == 0) return false;

  // Save pointers.  
  pdfAPtr = pdfAPtrIn;
  pdfBPtr = pdfBPtrIn;

  // By default same pointers for hard-process PDF's.
  pdfHardAPtr = pdfAPtrIn;
  pdfHardBPtr = pdfBPtrIn;
  
  // Optionally allow separate pointers for hard process.
  if (pdfHardAPtrIn == 0 || pdfHardBPtrIn == 0) return true;
  if (pdfHardAPtrIn == pdfHardBPtrIn) return false;
  pdfHardAPtr = pdfHardAPtrIn;
  pdfHardBPtr = pdfHardBPtrIn;
  
  // Done.
  return true;
}

//*********

// Routine to initialize with two beam energies specified.

bool Pythia::init( int idAin, int idBin, double eAin, double eBin) {

  // Read in and set values.
  idA = idAin;
  idB = idBin;
  eA  = eAin;
  eB  = eBin;
  inCMframe = false;

  // Send on to common initialization. 
  bool status = init();
  if (!status) ErrorMsg::message("Abort from Pythia::init: "
    "initialization failed");
  return status;

}

//*********

// Routine to initialize with CM energy rather tham beam energies.

bool Pythia::init( int idAin, int idBin, double eCMin) {

  // Read in and set values.
  idA = idAin;
  idB = idBin;
  eCM  = eCMin;
  inCMframe = true;

  // Send on to common initialization.
  bool status = init();
  if (!status) ErrorMsg::message("Abort from Pythia::init: "
    "initialization failed");
  return status;

}

//*********

// Routine to initialize when beam info is given in an LHAinit object.

bool Pythia::init( LHAinit* lhaInitPtrIn, LHAevnt* lhaEvntPtrIn) {

  // Save and set flag for subsequent usage of LHAevnt object.
  lhaInitPtr = lhaInitPtrIn;
  lhaEvntPtr = lhaEvntPtrIn;
  doLHA      = true;

  // Set LHAinit information (in some external program).
  bool status = lhaInitPtr->set();
  if (!status) ErrorMsg::message("Abort from Pythia::init: "
    "Les Houches initialization failed");
  if (!status) return false;

  // Extract beams from values set in an LHAinit object. 
  idA = lhaInitPtr->idBeamA();
  idB = lhaInitPtr->idBeamB();
  eA  = lhaInitPtr->eBeamA();
  eB  = lhaInitPtr->eBeamB();
  inCMframe = false;

  // Now do normal initialization. List info if there.
  status = init();
  if (strategyLHA < 10) lhaInitPtr->list();
  if (!status) ErrorMsg::message("Abort from Pythia::init: "
    "initialization failed");
  return status;

}

//*********

// Routine to initialize when all info is given in a Les Houches Event File.

bool Pythia::init( string LesHouchesEventFile) {

  // Create LHAinit and LHAevnt objects. 
  const char* cstring = LesHouchesEventFile.c_str();
  lhaInitPtr = new LHAinitLHEF(cstring);
  lhaEvntPtr = new LHAevntLHEF(cstring);
  doLHA      = true;

  // Set LHAinit information (in some external program).
  bool status = lhaInitPtr->set();
  if (!status) ErrorMsg::message("Abort from Pythia::init: "
    "Les Houches initialization failed");
  if (!status) return false;

  // Extract beams from values set in an LHAinit object. 
  idA = lhaInitPtr->idBeamA();
  idB = lhaInitPtr->idBeamB();
  eA  = lhaInitPtr->eBeamA();
  eB  = lhaInitPtr->eBeamB();
  inCMframe = false;

  // Now do normal initialization. List info if there.
  status = init();
  if (strategyLHA < 10) lhaInitPtr->list(); 
  if (!status) ErrorMsg::message("Abort from Pythia::init: "
    "initialization failed");
  return status;

}

//*********

// Main routine to initialize the generation process.
// (The alternative init forms end up in this one.)

bool Pythia::init() {

  // Check that constructor worked.
  if (!isConstructed) return false;

  // Reset error counter.
  nErrEvent = 0;

  // Strategy for event generation, Les Houches or not.
  strategyLHA = (doLHA) ? lhaInitPtr->strategy() : 0;  

  // Check that beams and beam combination can be handled.
  // Only allow neutrinos as beams when leptons unresolved.
  bool canHandleBeams = false;
  int idAabs = abs(idA);
  int idBabs = abs(idB);
  if (idAabs == 2212 && idBabs == 2212) canHandleBeams = true;
  else if(settings.flag("Pythia:lPDF")) {
    if ( idA + idB == 0 && (idAabs == 11 || idAabs == 13
      || idAabs == 15) ) canHandleBeams = true;
  } else if (idAabs > 10 && idAabs < 17 && idA * idB < 0) {
    if (idA + idB == 0) canHandleBeams = true;
    int idMax  = max(idAabs, idBabs);
    int idMin  = min(idAabs, idBabs);
    if (idMax - idMin == 1 && idMax%2 == 0) canHandleBeams = true; 
  }
  if (!canHandleBeams && strategyLHA < 10) {
    ErrorMsg::message("Error in Pythia::init: "
      "cannot handle this beam combination");
    return false;
  }

  //Set up values related to user hooks.
  hasUserHooks  = (userHooksPtr > 0);
  doVetoProcess = (hasUserHooks) 
                ? userHooksPtr->canVetoProcessLevel() : false;
  doVetoPartons = (hasUserHooks) 
                ? userHooksPtr->canVetoPartonLevel() : false;  

  // Set up objects for timelike and spacelike showers.
  if (timesDecPtr == 0 || timesPtr == 0) {
    TimeShower* timesNow = new TimeShower();
    if (timesDecPtr == 0) timesDecPtr = timesNow;
    if (timesPtr == 0) timesPtr = timesNow; 
    timesNew = true;
  }
  if (spacePtr == 0) {
    spacePtr = new SpaceShower();
    spaceNew = true;
  }

  // Initialize for simple showers in decays. 
  timesDecPtr->init( 0, 0);

  // Initialize all accessible static data members (except resonances).
  initStatic();

  // Initialize SUSY Les Houches Accord data
  if (!initSLHA()) return false;

  // Initialize Breit-Wigner mass selection.
  particleData.initBWmass();

  // Set headers to distinguish the two event listing kinds.
  process.header("(hard process)");
  event.header("(complete event)");

  // Do not set up beam kinematics when no parton-level processing.
  if (strategyLHA >= 10) inCMframe = true;
  else {

    // Find masses. Initial guess about CM frame.
    mA = ParticleDataTable::m0(idA);
    mB = ParticleDataTable::m0(idB);
    betaZ = 0.;
    gammaZ = 1.;

    // When not given: find CM energy and set up boost to rest frame.
    if (!inCMframe) {
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
    eA  = sqrt(mA*mA + pzA*pzA);
    eB  = sqrt(mB*mB + pzB*pzB);

    // Set up the PDF's, if not already done.
    if (pdfAPtr == 0) {
      pdfAPtr     = getPDFPtr(idA); 
      if (!pdfAPtr->isSetup()) return false;
      pdfHardAPtr = pdfAPtr;
      pdfAnew     = true;
    }
    if (pdfBPtr == 0) {
      pdfBPtr     = getPDFPtr(idB); 
      if (!pdfBPtr->isSetup()) return false;
      pdfHardBPtr = pdfBPtr;
      pdfBnew     = true;
    }

    // Optionally set up separate PDF's for hard process.
    if (settings.flag("Pythia:useHardPFD")) {
      pdfHardAPtr = getPDFPtr(idA, 2);      
      if (!pdfHardAPtr->isSetup()) return false;
      pdfHardBPtr = getPDFPtr(idB, 2);      
      if (!pdfHardBPtr->isSetup()) return false;
      pdfHardNew  = true;
    }
  
    // Set up the two beams and the common remnant system.
    bool isUnresolvedA = ( ParticleDataTable::isLepton(idA) 
      && !settings.flag("Pythia:lPDF") );
    bool isUnresolvedB = ( ParticleDataTable::isLepton(idB) 
      && !settings.flag("Pythia:lPDF") );
    beamA.init( idA, pzA, eA, mA, pdfAPtr, pdfHardAPtr, isUnresolvedA);
    beamB.init( idB, pzB, eB, mB, pdfBPtr, pdfHardBPtr, isUnresolvedB);
  }

  // Store main info for access in process generation.
  info.setBeamA( idA, pzA, eA, mA);
  info.setBeamB( idB, pzB, eB, mB);
  info.setECM( eCM);

  // Set info in and send pointers to the respective program elements.
  ProcessContainer::setInfoPtr( &info);
  if (!processLevel.init( &info, &beamA, &beamB, 
    doLHA, lhaInitPtr, lhaEvntPtr, userHooksPtr)) return false;
  if (!partonLevel.init( &info, &beamA, &beamB, timesDecPtr, timesPtr,
    spacePtr, strategyLHA, userHooksPtr)) return false;
  if (!hadronLevel.init( &info, timesDecPtr)) return false;

  // Succeeded.
  isInit = true; 
  return true;
}

//*********

// Initialization routine for all accessible static data members,
// except for individual resonances, which are in ProcesLevel::init(). 

void Pythia::initStatic() {

  // Initialize printing of error messages.
  ErrorMsg::initStatic();

  // Initialize the random number generator.
  if ( settings.flag("Pythia:setSeed") )  
    Rndm::init( settings.mode("Pythia:seed") );

  // Initialize static data members in the Pythia class.
  doPartonLevel = settings.flag("Pythia:partonLevel");
  doHadronLevel = settings.flag("Pythia:hadronLevel");
  checkEvent    = settings.flag("Pythia:checkEvent");
  nErrList      = settings.mode("Pythia:nErrList");
  epTolerance   = settings.parm("Pythia:epTolerance");
  
  // Initialize all other accessible static data members.
  ParticleDataEntry::initStatic();
  AlphaEM::initStatic(); 
  CoupEW::initStatic(); 
  VCKM::initStatic(); 
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
  ResonanceProperties::initStatic();

}

//*********

// Initialize SUSY Les Houches Accord data.

bool Pythia::initSLHA() {

  // Check whether SUSY is on.
  if ( !settings.flag("SUSY") ) return true;      

  // Read SUSY Les Houches Accord File.
  string slhafile = settings.word("SUSY:SusyLesHouchesFile");
  int ifail = slha.readFile(slhafile);

  // In case of problems, print error and fail init.
  if (ifail != 0) {
    ErrorMsg::message("Error from Pythia::initSLHA: "
      "problem reading SLHA file", slhafile);
    return false;
  };

  // Update particle data.
  int id = slha.mass.first();
  for (int i = 1; i <= slha.mass.size() ; i++) {
    double mass = abs(slha.mass(id));
    particleData.m0(id,mass);
    id = slha.mass.next();
  };

  // Check spectrum for consistency. Switch off SUSY if necessary.
  ifail = slha.checkSpectrum();
  if (ifail != 0) {
    ErrorMsg::message("Warning from Pythia::initSLHA: "
      "Problem with SLHA spectrum.", 
      "\n Only using masses and switching off SUSY.");
    settings.flag("SUSY", false);
    return true;
  };

  // Store pointer as static member in SigmaProcess.
  SigmaProcess::setSlhaPtr(&slha);

  return true;

}

//*********

// Main routine to generate the next event, using internal machinery.

bool Pythia::next() {

  // Can only generate event if initialization worked.
  if (!isInit) {
    ErrorMsg::message("Abort from Pythia::next: "
      "not properly initialized so cannot generate events"); 
    return false;
  }

  // Outer loop over hard processes; only relevant for user-set vetoes.
  for ( ; ; ) {
    bool hasVetoed = false;

    // Provide the hard process that starts it off. Only one try.
    info.clear();
    process.clear();
    if ( !processLevel.next( process) ) {
      ErrorMsg::message("Abort from Pythia::next: "
        "processLevel failed; giving up"); 
      return false;
    }

    // Possibility for a user veto of the process-level event.
    if (doVetoProcess) {
      hasVetoed = userHooksPtr->doVetoProcessLevel( process);
      if (hasVetoed) continue;
    }

    // Possibility to stop the generation at this stage.
    if (!doPartonLevel) {
      if (!inCMframe) process.bst(0., 0., betaZ, gammaZ);
      processLevel.accumulate();
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
        // Skip to next hard process for failure owing to deliberate veto.
        hasVetoed = partonLevel.hasVetoed(); 
        if (hasVetoed) break;
        // Else make a new try for other failures.
        ErrorMsg::message("Error in Pythia::next: "
          "partonLevel failed; try again"); 
        physical = false; 
        continue;
      }

      // Possibility for a user veto of the parton-level event.
      if (doVetoPartons) {
        hasVetoed = userHooksPtr->doVetoPartonLevel( event);
        if (hasVetoed) break;
      }

      // When required: boost to lab frame (before decays, for vertices).
      if (!inCMframe) {
        process.bst(0., 0., betaZ, gammaZ);
        event.bst(0., 0., betaZ, gammaZ);      
      }    

      // Possibility to stop the generation at this stage.
      if (!doHadronLevel) {
        processLevel.accumulate();
        partonLevel.accumulate();
        return true;
      }

      // Hadron-level: hadronization, decays.
      if ( !hadronLevel.next( event) ) {
        ErrorMsg::message("Error in Pythia::next: "
          "hadronLevel failed; try again"); 
        physical = false; 
        continue;
      }

      // Stop parton- and hadron-level looping if you got this far.
      break;
    }

    // If event vetoed then to make a new try.
    if (hasVetoed) continue;

    // If event failed any other way (after ten tries) then give up.
    if (!physical) {
      ErrorMsg::message("Abort from Pythia::next: "
        "parton+hadronLevel failed; giving up");
      return false;
    }

    // Process- and parton-level statistics. 
    processLevel.accumulate();
    partonLevel.accumulate();

    // End of outer loop over hard processes.
    break;
  }

  // Optionally check final event for problems.
  if (checkEvent && !check()) {
    ErrorMsg::message("Abort from Pythia::next: "
      "check of event revealed problems");
    return false;
  }

  // Done.
  return true;

}

//*********

// Print statistics on event generation.

void Pythia::statistics(bool all) {

  // Statistics on cross section and number of events. 
  processLevel.statistics();  

  // Statistics from other classes, e.g. multiple interactions.
  if (all) partonLevel.statistics();  

  // Summary of which and how many warnings/errors encountered.
  ErrorMsg::statistics();

}

//*********

// Write the Pythia banner, with symbol and version information.

void Pythia::banner(ostream& os) {

  // Read in version number and last date of change.
  double versionNumber = Settings::parm("Pythia:versionNumber");
  int versionDate = Settings::mode("Pythia:versionDate");
  string month[12] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};  

  // Get date and time.
  time_t t = time(0);
  char dateNow[12];
  strftime(dateNow,12,"%d %b %Y",localtime(&t));
  char timeNow[9];
  strftime(timeNow,9,"%H:%M:%S",localtime(&t));

  
  os << "\n"
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
     << " |  |   PPP   Y   Y  TTTTT  H   H  III    A  "
     << "    Welcome to the Lund Monte Carlo!  |  | \n" 
     << " |  |   P  P   Y Y     T    H   H   I    A A " 
     << "    This is PYTHIA version " << fixed << setprecision(3) 
     << setw(5) << versionNumber << "      |  | \n"
     << " |  |   PPP     Y      T    HHHHH   I   AAAAA"
     << "    Last date of change: " << setw(2) << versionDate%100 
     << " " << month[ (versionDate/100)%100 - 1 ] 
     << " " << setw(4) << versionDate/10000 <<  "  |  | \n"
     << " |  |   P       Y      T    H   H   I   A   A"
     << "                                      |  | \n"
     << " |  |   P       Y      T    H   H  III  A   A"
     << "    Now is " << dateNow << " at " << timeNow << "    |  | \n"
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
     << " |  |   The main physics reference is the 'PY"
     << "THIA 6.4 Physics and Manual',         |  | \n"
     << " |  |   T. Sjostrand, S. Mrenna and P. Skands"
     << ", JHEP05 (2006) 026 [hep-ph/0603175]. |  | \n"
     << " |  |   In addition, for PYTHIA 8.0, also quo"
     << "te the 'Brief Introduction',          |  | \n"
     << " |  |   T. Sjostrand, CERN-LCGAPP-2007-03.   "
     << "                                      |  | \n" 
     << " |  |                                        "
     << "                                      |  | \n"
     << " |  |   An archive of program versions and do" 
     << "cumentation is found on the web:      |  | \n"
     << " |  |   http://www.thep.lu.se/~torbjorn/Pythi" 
     << "a.html                                |  | \n"
     << " |  |                                        "
     << "                                      |  | \n"
     << " |  |   Disclaimer: this program comes withou"
     << "t any guarantees.                     |  | \n"
     << " |  |   Beware of errors and use common sense"
     << " when interpreting results.           |  | \n"
     << " |  |   In addition, the current 8.0 version " 
     << "is intended for tryout and feedback   |  | \n"
     << " |  |   only, and should not be used for any " 
     << "physics studies or production runs.   |  | \n"
     << " |  |                                        "
     << "                                      |  | \n"
     << " |  |   Copyright C 2007 Torbjorn Sjostrand  " 
     << "                                      |  | \n"
     << " |  |                                        "
     << "                                      |  | \n"
     << " |  |                                        "
     << "                                      |  | \n"
     << " |  *----------------------------------------" 
     << "--------------------------------------*  | \n"
     << " |                                           "
     << "                                         | \n"
     << " *-------------------------------------------" 
     << "-----------------------------------------* \n" << endl;

}

//*********

// Check that the final event makes sense: no unknown id codes;
// charge and energy-momentum conserved.

bool Pythia::check(ostream& os) {

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
      ErrorMsg::message("Error in Pythia::check: "
        "unknown particle code", errCode); 
      physical = false;
      iErrId.push_back(i);
    }

    // Look for particle with not-a-number energy/momentum/mass.
    if (abs(event[i].px()) >= 0. && abs(event[i].py()) >= 0. 
      && abs(event[i].pz()) >= 0.  && abs(event[i].e()) >= 0. 
      && abs(event[i].m()) >= 0.) ;
    else {   
      ErrorMsg::message("Error in Pythia::check: "
        "not-a-number energy/momentum/mass"); 
      physical = false;
      iErrNan.push_back(i);
    }

    // Add final-state four-momentum and charge.      
    if (event[i].isFinal()) {
      pSum += event[i].p();
      chargeSum += event[i].charge();
    }

  // End of particle loop.
  }

  // Check energy-momentum/charge conservation.
  double epDev = abs(pSum.e()) + abs(pSum.px()) + abs(pSum.py())
    + abs(pSum.pz());
  if (epDev > epTolerance * eLab) { 
    ErrorMsg::message("Error in Pythia::check: "
      "energy-momentum badly conserved"); 
    physical = false;
  }
  if (abs(chargeSum) > 0.1) {
    ErrorMsg::message("Error in Pythia::check: "
      "charge not conserved"); 
    physical = false;
  }

  // Done for sensible events.
  if (physical) return true;

  // Print (the first few) flawed events.
  if (nErrEvent < nErrList) {
    os << " PYTHIA erroneous event info: \n";
    if (iErrId.size() > 0) {
      os << " unknown particle codes in lines ";
      for (int i = 0; i < int(iErrId.size()); ++i) 
        os << iErrId[i] << " ";
      os << "\n";
    }
    if (iErrNan.size() > 0) {
      os << " not-a-number energy/momentum/mass in lines ";
      for (int i = 0; i < int(iErrNan.size()); ++i) 
        os << iErrNan[i] << " ";
      os << "\n";
    }
    if (epDev > epTolerance * eLab) os << scientific << setprecision(3)
      << " total energy-momentum non-conservation = " << epDev << "\n";
    if (abs(chargeSum) > 0.1) os << fixed << setprecision(2) 
      << " total charge non-conservation = " << chargeSum << "\n"; 
    info.list();
    event.list();
  }

  // Update error counter. Done also for flawed event.
  ++nErrEvent;
  return false;

}

//*********

// Routine to set up a PDF pointer.

PDF* Pythia::getPDFPtr(int idIn, int sequence) {

  PDF* tempPDFPtr;

  // Proton beam, normal choice.
  if (abs(idIn) == 2212 && sequence == 1) {
    int  pPDFset   = settings.mode("Pythia:pPDFset");
    bool useLHAPDF = settings.flag("Pythia:useLHAPDF");

    // Use internal sets.
    if (!useLHAPDF) {
      if (pPDFset == 1) tempPDFPtr = new GRV94L(idIn);
      else tempPDFPtr = new CTEQ5L(idIn);
    }
    
    // Use sets from LHAPDF.
    else {
      string LHAPDFset    = settings.word("Pythia:LHAPDFset");
      int    LHAPDFmember = settings.mode("Pythia:LHAPDFmember");
      tempPDFPtr = new LHAPDF(idIn, LHAPDFset, LHAPDFmember);

      // Set x and Q2 limits.
      if (settings.flag("Pythia:limitLHAPDF")) {
        double xMin  = settings.parm("Pythia:xMinLHAPDF");
        double xMax  = settings.parm("Pythia:xMaxLHAPDF");
        double Q2Min = settings.parm("Pythia:Q2MinLHAPDF");
        double Q2Max = settings.parm("Pythia:Q2MaxLHAPDF");
        tempPDFPtr->setLimits( xMin, xMax, Q2Min, Q2Max);
      }
    }
  }

  // Proton beam, special choice for the hard process..
  else if (abs(idIn) == 2212) {
    int  pPDFset   = settings.mode("Pythia:hardpPDFset");
    bool useLHAPDF = settings.flag("Pythia:useHardLHAPDF");

    // Use internal sets.
    if (!useLHAPDF) {
      if (pPDFset == 1) tempPDFPtr = new GRV94L(idIn);
      else tempPDFPtr = new CTEQ5L(idIn);
    }
    
    // Use sets from LHAPDF.
    else {
      string LHAPDFset    = settings.word("Pythia:hardLHAPDFset");
      int    LHAPDFmember = settings.mode("Pythia:hardLHAPDFmember");
      tempPDFPtr = new LHAPDF(idIn, LHAPDFset, LHAPDFmember, 2);

      // Set x and Q2 limits.
      if (settings.flag("Pythia:limitHardLHAPDF")) {
        double xMin  = settings.parm("Pythia:xMinHardLHAPDF");
        double xMax  = settings.parm("Pythia:xMaxHardLHAPDF");
        double Q2Min = settings.parm("Pythia:Q2MinHardLHAPDF");
        double Q2Max = settings.parm("Pythia:Q2MaxHardLHAPDF");
        tempPDFPtr->setLimits( xMin, xMax, Q2Min, Q2Max);
      }
    }
  }

  // Lepton beam; resolved or not.
  else {
    if (settings.flag("Pythia:lPDF") && abs(idIn)%2 == 1) 
      tempPDFPtr = new Lepton(idIn);
    else tempPDFPtr = new LeptonPoint(idIn);
  }
  
  // Done.
  return tempPDFPtr; 
}

//*********

} // end namespace Pythia8

