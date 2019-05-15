// File: main02.cc
// This is a simple test program. 
// It illustrates how Les Houches Accord input can be used in Pythia8.
// It uses the ttsample.init and ttsample.evnt input files, 
// the latter only with 100 events.
// Copyright C 2006 Torbjorn Sjostrand

#include "Pythia.h"
using namespace Pythia8; 
int main() {

  // Number of events to print.
  int nPrint = 2;             

  // Generator           
  Pythia pythia;                            

  // Stick with default values, so do not bother with a separate file
  // for changes. However, do one change, to show readString in action.
   pythia.readString("PartonLevel:MI = off"); 

  // Les Houches initialization and event objects.
  LHAinitPythia6 lhaInit("ttsample.init");  
  LHAevntPythia6 lhaEvnt("ttsample.evnt"); 

  // Initialize with pointers. List initialization information.
  pythia.init(&lhaInit, &lhaEvnt);      
  cout << lhaInit;                       

  // Book histogram.
  Hist nCharged("charged particle multiplicity",100,-0.5,199.5); 

  // Allow for possibility of a few faulty events.
  int nAbort = 10;
  int iAbort = 0;

  // Begin event loop; generate until none left in input file.     
  for (int iEvent = 0; ; ++iEvent) {

    // First few failures write off as potentially error, then quit.
    // (But probably end of file reached already first time.)
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      break;
    }
 
    // List first few events: Les Houches, hard process and complete.
    if (iEvent < nPrint) {     
      cout << lhaEvnt;               
      pythia.process.list();          
      pythia.event.list();           
    }                           

    // Sum up final charged multiplicity and fill in histogram.
    int nChg = 0;                 
    for (int i = 0; i < pythia.event.size(); ++i) 
    if (pythia.event[i].remains() && pythia.event[i].isCharged()) 
      ++nChg;
    nCharged.fill(nChg);               

  // End of event loop.        
  }                                           

  // Give statistics. Print histogram.
  pythia.statistics();
  cout << nCharged;  

  // Done.                           
  return 0;
}
