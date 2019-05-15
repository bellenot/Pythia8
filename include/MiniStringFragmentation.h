// This file contains the class for "cluster" fragmentation.
// MiniStringFragmentation: handle the fragmentation of low-mass systems.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_MiniStringFragmentation_H
#define Pythia8_MiniStringFragmentation_H

#include "Basics.h"
#include "Event.h"
#include "FragmentationFlavZpT.h"
#include "FragmentationSystems.h"
#include "ParticleData.h"
#include "PythiaStdlib.h"
#include "Settings.h"

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
  bool fragment( int iSub, ColConfig& colConfig, Event& event, 
    bool isDiff = false);

private: 

  // Initialization data, normally only set once.
  static int nTryMass;
  static double sigma, sigma2Had, bLund;

  // Constants: could only be changed in the code itself.
  static const int NTRYDIFFRACTIVE, NTRYLASTRESORT, NTRYFLAV;
  static const double SIGMAMIN;

  // Attempt to produce two particles from a cluster.
  bool ministring2two( int nTry, Event& event);

  // Attempt to produce one particle from a cluster.
  bool ministring2one( int iSub, ColConfig& colConfig, Event& event);

  // Class for flavour generation.
  StringFlav flavSel;

  // Data members.
  vector<int> iParton;
  int id1, id2;
  Vec4 pSum;
  double mSum, m2Sum;
  bool isClosed;

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_MiniStringFragmentation_H
