// main32.cc is a part of the PYTHIA event generator.
// Copyright (C) 2012 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Richard Corke (richard.corke@thep.lu.se)
// This is a sample program showing the usage of:
//  AlpgenHooks: a UserHooks derived class provided by 'LHAupAlpgen.h'
//               for reading in ALPGEN format event files
//  MLMhooks:    a UserHooks derived class provided by 'MLMhooks.h'
//               for performing the MLM merging procedure
// Some user supplied options are available, made functional through
// these two UserHooks, and are described in the 'Alpgen and MLM Merging'
// manual page. Further details of the different classes provided
// in these header files is also given on the manual page.

// Includes and namespace
#include "Pythia.h"
#include "LHAupAlpgen.h"
#include "MLMhooks.h"
using namespace Pythia8;

//==========================================================================

// AlpgenAndMLMhooks:
//   A small UserHooks class that gives the functionality of
//   both AlpgenHooks and MLMhooks. These classes only have two
//   overlapping functions, the constructor and 'initAfterBeams()',
//   which are both overridden here such that both are called.

class AlpgenAndMLMhooks : public AlpgenHooks, public MLMhooks {

public:
  AlpgenAndMLMhooks(Pythia &pythia) 
    : AlpgenHooks(pythia), MLMhooks() {}

  virtual bool initAfterBeams() {
    // Call init of both AlpgenHooks and MLMhooks (order important)
    if (!AlpgenHooks::initAfterBeams()) return false;
    if (!MLMhooks::initAfterBeams())    return false;
    return true;
  }
};

//==========================================================================

// Main program: initialise Pythia with AlpgenAndMLMhooks as above

int main() {

  // Generator and read in commands
  Pythia pythia;
  pythia.readFile("main32.cmnd");

  // Extract settings to be used in the main program.
  int nEvent = pythia.mode("Main:numberOfEvents");
  int nAbort = pythia.mode("Main:timesAllowErrors");

  // Create AlpgenAndMLMhooks class and pass to Pythia
  UserHooks *alpgenAndMLMhooks = new AlpgenAndMLMhooks(pythia);
  pythia.setUserHooksPtr(alpgenAndMLMhooks);

  // Initialise Pythia
  if (!pythia.init()) {
    cout << "Error: could not initialise Pythia" << endl;
    return 1;
  };

  // Begin event loop.
  int iAbort = 0, iEvent = 0;
  for (; ; ++iEvent) {
    // Check if iEvent > nEvent
    if (nEvent > 0 && iEvent >= nEvent) break;

    // Generate events. Quit if many failures.
    if (!pythia.next()) {
      if (pythia.info.atEndOfFile()) {
        cout << "Info: end of input file reached" << endl;
        break;
      }
      if (++iAbort < nAbort) continue;
      cout << "Abort: too many errors in generation" << endl;
      break;
    }

    // Event analysis here

  // End of event loop.
  }
   
  // Final statistics.
  pythia.stat();

  // Clean up and done.
  delete alpgenAndMLMhooks;
  return 0;
}
