// Function definitions (not found in the header) for the Info
// and ErrorMessages classes.
// Copyright C 2006 Torbjorn Sjostrand

#include "Information.h"

namespace Pythia8 {

//**************************************************************************

// Info class.
// This class contains a mixed bag of information on the event generation 
// activity, especially on the current subprocess properties.

//*********

// List (almost) all information currently set.

void Info::list(ostream& os) {

  // Header and beam info.
  os << "\n --------  PYTHIA Info Listing  ------------------------"
     << "---------------- \n \n" 
     << scientific << setprecision(3) 
     << " Beam A: id = " << setw(6) << idAM << ", pz = " << setw(10) 
     << pzAM << ", e = " << setw(10) << eAM << ", m = " << setw(10) 
     << mAM << "\n"
     << " Beam B: id = " << setw(6) << idBM << ", pz = " << setw(10) 
     << pzBM << ", e = " << setw(10) << eBM << ", m = " << setw(10) 
     << mBM << "\n\n";

  // Colliding parton info.
  if (isRes) 
    os << " In 1: id = " << setw(4) << id1H << ", x = " << setw(10)
       << x1H << ", pdf = " << setw(10) << pdf1H << " at Q2 = " 
       << setw(10) << Q2pdfH << "\n"  
       << " In 2: id = " << setw(4) << id2H << ", x = " << setw(10)
       << x2H << ", pdf = " << setw(10) << pdf2H << " at same Q2\n\n";  

  // Process info.
  os << ((isRes) ? " Subprocess " : " Process ") << nameSave 
     << " with code " << codeSave << " is 2 -> " << nFinalSave << "\n";
  if (isRes && nFinalSave == 1) 
    os << " has  sHat = " << setw(10) << sH << " \n";  
  else if ( isRes && nFinalSave == 2)  
    os << " has   sHat = " << setw(10) << sH << ",    tHat = " 
       << setw(10) << tH << ",    uHat = " << setw(10) << uH << ",\n"
       << "      pTHat = " << setw(10) << pTH << ",   m3Hat = " 
       << setw(10) << m3H << ",   m4Hat = " << setw(10) << m4H << ",\n"
       << "   thetaHat = " << setw(10) << thetaH << ",  phiHat = " 
       << setw(10) << phiH << "\n";
  
  else if ( nFinalSave == 2)  
    os << " has   s = " << setw(10) << sH << ",    t = " << setw(10) 
       << tH << ",    u = " << setw(10) << uH << ",\n"
       << "      pT = " << setw(10) << pTH << ",   m3 = " << setw(10) 
       << m3H << ",   m4 = " << setw(10) << m4H << ",\n" 
       << "   theta = " << setw(10) << thetaH << ",  phi = " << setw(10) 
       << phiH << "\n";
  
  // Listing finished.
  os << "\n --------  End PYTHIA Info Listing  --------------------"
     << "----------------" << endl; 

}

//**************************************************************************

// ErrorMessages class.
// This class holds info on all messages received and how many times.

//*********

// Definitions of static variables. 
// (Values will be overwritten in initStatic call, so are purely dummy.)

int ErrorMessages::timesToPrint = 1;
map<string, int> ErrorMessages::messages;

//*********

// Initialize static data members.

void ErrorMessages::initStatic() {

  timesToPrint = Settings::mode("ErrorMessages:timesToPrint"); 

}

//*********
  
// Print a message the first few times. Insert in database.
 
void ErrorMessages::message(string messageIn, string extraIn) {
   
  // Recover number of times message occured. Also inserts new string.
  int times = messages[messageIn];
  ++messages[messageIn];

  // Print message the first few times.
  if (times < timesToPrint) cout << " " << messageIn << " " 
    << extraIn << "\n";

}

//*********

// Print statistics on errors/warnings.

void ErrorMessages::statistics(ostream& os) {

  // Header.
  os << "\n *-------  PYTHIA Error and Warning Messages Statistics  "
     << "---------------------------------------------------------* \n"
     << " |                                                       "
     << "                                                         | \n"
     << " |  times   message                                      "
     << "                                                         | \n" 
     << " |                                                       "
     << "                                                         | \n";

  // Loop over all messages
  map<string, int>::iterator messageEntry = messages.begin();
  if (messageEntry == messages.end()) 
    os << " |      0   no errors or warnings to report              "
       << "                                                         | \n";
  while (messageEntry != messages.end()) {
    // Debug printout.
    string temp = messageEntry->first;
    int len = temp.length();
    temp.insert( len, max(0, 101 - len), ' ');
    os << " | " << setw(6) << messageEntry->second << "   " 
       << temp << " | \n";
    ++messageEntry;
  } 

  // Done. 
  os << " |                                                       "
     << "                                                         | \n"
     << " *-------  End PYTHIA Error and Warning Messages Statistics"
     << "  -----------------------------------------------------* " 
     << endl;

}

//**************************************************************************

} // end namespace Pythia8

