// This file contains the class for "cluster" fragmentation.
// MiniStringFragmentation: handle the fragmentation of low-mass systems.
// Copyright © 2005 Torbjörn Sjöstrand

#ifndef Pythia8_MiniStringFragmentation_H
#define Pythia8_MiniStringFragmentation_H

#include "Stdlib.h"
#include "Basics.h"
#include "Settings.h"
#include "ParticleData.h"
#include "Event.h"
#include "FragmentationSystems.h"
#include "FragmentationFlavZpT.h"

namespace Pythia8 {

//**************************************************************************

// The MiniStringFragmentation class contains the routines to fragment 
// occasional low-mass colour singlet partonic systems, where the string 
// approach is not directly applicable (for technical reasons).

class MiniStringFragmentation {

public:

  // Constructor. 
  MiniStringFragmentation() {}

  // Initialize static data members.
  static void initStatic();

  // Do the fragmentation: driver routine.
  bool fragment( int iSub, ColConfig& colConfig, Event& event);

private: 

  // Initialization data, normally only set once.
  static int nTryMass;
  static double sigma, sigma2Had, bLund;

  // Constants: could only be changed in the code itself.
  static const int NTRYLASTRESORT, NTRYFLAV;
  static const double SIGMAMIN;

  // Attempt to produce two particles from a cluster.
  bool ministring2two( int nTry, Event& event);

  // Attempt to produce one particle from a cluster.
  bool ministring2one( int iSub, ColConfig& colConfig, Event& event);

  // Class for flavour generation.
  StringFlav flavSel;

  // Data members.
  int id1, id2;
  Vec4 pSum;
  double mSum, m2Sum;
  bool isClosed;
  vector<int> iParton;

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_MiniStringFragmentation_H
