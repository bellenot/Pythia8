// Function definitions (not found in the header) for the LHAinit,
// LHAevnt, LHAinitFortran, LHAevntFortran, LHAinitPythia6 and 
// LHAevntPythia6 classes.
// Copyright C 2006 Torbjorn Sjostrand

#include "LesHouches.h"

//**************************************************************************

// LHAinit class.

//*********

// Print the info; to check it worked.

ostream& operator<<(ostream& os, const LHAinit& lha) {

  // Header.
  os << "\n --------  LHA initialization information  ------------ \n"; 

  // Beam info.
  os << fixed << setprecision(3) 
     << "\n  beam    kind      energy  pdfgrp  pdfset \n" 
     << "     A  " << setw(6) << lha.idBeamA() <<  setw(12) 
     << lha.eBeamA() << setw(8) << lha.pdfGroupBeamA() << setw(8) 
     << lha.pdfSetBeamA() << "\n" 
     << "     B  " << setw(6) << lha.idBeamB() <<  setw(12) 
     << lha.eBeamB() << setw(8) << lha.pdfGroupBeamB() << setw(8) 
     << lha.pdfSetBeamB() << "\n"; 

  // Event weighting strategy.
  os << "\n  Event weighting strategy = " << setw(2) 
     << lha.strategy() << "\n" ;

  // Process list.
  os << scientific << setprecision(4) 
     << "\n  Processes, with strategy-dependent cross section info \n" 
     << "  number      xsec (pb)      xerr (pb)      xmax (pb) \n" ;
  for (int ip = 0; ip < lha.size(); ++ip) {
    os << setw(8) << lha.idProcess(ip) << setw(15) << lha.xSec(ip) 
       << setw(15) << lha.xErr(ip) << setw(15) << lha.xMax(ip) 
       << "\n";
  }

  // Finished.
  os << "\n --------  End LHA initialization information  -------- \n"; 
  return os;

}
 
//**************************************************************************

// LHAevnt class.

//*********

// Print the info; to check it worked.

ostream& operator<<(ostream& os, const LHAevnt& lha) {

  // Header.
  os << "\n --------  LHA event information and listing  -------------"
     << "--------------------------------------------------------- \n"; 

  // Basic event info.
  os << scientific << setprecision(4) 
     << "\n    process = " << setw(8) << lha.idProc() 
     << "    weight = " << setw(12) << lha.weight() 
     << "     scale = " << setw(12) << lha.scale() << " (GeV) \n"
     << "                   "
     << "     alpha_em = " << setw(12) << lha.alphaQED() 
     << "    alpha_strong = " << setw(12) << lha.alphaQCD() << "\n";

  // Particle list
  os << fixed << setprecision(3) 
     << "\n    Participating Particles \n" 
     << "    no        id stat     mothers     colours      p_x        "
     << "p_y        p_z         e          m        tau    spin \n" ;
  for (int ip = 0; ip < lha.size(); ++ip) {
    os << setw(6) << ip << setw(10) << lha.id(ip) << setw(5)
       << lha.status(ip) << setw(6) << lha.mother1(ip) << setw(6) 
       << lha.mother2(ip) << setw(6) << lha.col1(ip) << setw(6) 
       << lha.col2(ip) << setw(11) << lha.px(ip) << setw(11) 
       << lha.py(ip) << setw(11) << lha.pz(ip) << setw(11) << lha.e(ip)
       << setw(11) << lha.m(ip) << setw(8) << lha.tau(ip) << setw(8) 
       << lha.spin(ip) << "\n";
  }

  // Finished.
  os << "\n --------  End LHA event information and listing  ---------"
     << "--------------------------------------------------------- \n"; 
  return os;

}
 
//**************************************************************************

// LHAinitFortran class.

//*********

// Give access to the HEPRUP Fortran commonblock.

extern "C" {

  extern struct {
    int idbmup[2];
    double ebmup[2];
    int pdfgup[2], pdfsup[2], idwtup, nprup;
    double xsecup[100], xerrup[100], xmaxup[100];
    int lprup[100];
  } heprup_;

}

//*********

// Read in information stored in the HEPRUP Fortran commonblock.

bool LHAinitFortran::set() {

  // Store beam and strategy info. 
  beamA(heprup_.idbmup[0], heprup_.ebmup[0], heprup_.pdfgup[0], 
    heprup_.pdfsup[0]);
  beamB(heprup_.idbmup[1], heprup_.ebmup[1], heprup_.pdfgup[1], 
    heprup_.pdfsup[1]);
  strategy(heprup_.idwtup);

  // Store process info.
  for (int ip = 0; ip < heprup_.nprup; ++ip) process(heprup_.lprup[ip], 
    heprup_.xsecup[ip], heprup_.xerrup[ip], heprup_.xmaxup[ip]) ;

  // Done.
  return true;

}
 
//**************************************************************************

// LHAevntFortran class.

//*********

// Give access to the HEPEUP Fortran commonblock.

extern "C" {

  extern struct {
    int nup, idprup;
    double xwgtup, scalup, aqedup, aqcdup;
    int idup[500], istup[500], mothup[500][2], icolup[500][2];
    double pup[500][5], vtimup[500],spinup[500];
  } hepeup_;

}

//*********

// Read in information stored in the HEPEUP Fortran commonblock.
bool LHAevntFortran::set() {

  // Store process info.
  process(hepeup_.idprup, hepeup_.xwgtup, hepeup_.scalup, hepeup_.aqedup, 
    hepeup_.aqcdup);
  
  // Store particle info.
  for (int ip = 0; ip < hepeup_.nup; ++ip) particle(hepeup_.idup[ip], 
    hepeup_.istup[ip], hepeup_.mothup[ip][0], hepeup_.mothup[ip][1], 
    hepeup_.icolup[ip][0], hepeup_.icolup[ip][1], hepeup_.pup[ip][0], 
    hepeup_.pup[ip][1], hepeup_.pup[ip][2], hepeup_.pup[ip][3], 
    hepeup_.pup[ip][4], hepeup_.vtimup[ip], hepeup_.spinup[ip]) ;

  // Done.
  return true;

}

//**************************************************************************

// LHAinitPythia6 class.

//*********

// Read in a file previously written from Pythia 6.3 with PYUPIN.
// It should be transparent how to apply it also to other programs.

bool LHAinitPythia6::set() {

  // Reset infile error flag.
  is.clear(); 
  
  // Read in beam and strategy info, and store it. 
  int idbmupA, idbmupB;
  double ebmupA, ebmupB;
  int pdfgupA, pdfgupB, pdfsupA, pdfsupB, idwtup, nprup;
  is >> idbmupA >> idbmupB >> ebmupA >> ebmupB >> pdfgupA 
     >> pdfgupB >> pdfsupA >> pdfsupB >> idwtup >> nprup;
  if (!is) return false;
  beamA(idbmupA, ebmupA, pdfgupA, pdfsupA);
  beamB(idbmupB, ebmupB, pdfgupB, pdfsupB);
  strategy(idwtup);

  // Read in process info, one process at a time, and store it.
  double xsecup, xerrup, xmaxup;
  int lprup; 
  for (int ip = 0; ip < nprup; ++ip) { 
    is >> xsecup >> xerrup >> xmaxup >> lprup ;
    if (!is) return false;
    process(lprup, xsecup, xerrup, xmaxup);
  }

  // Reading worked.
  return true;

}
 
//**************************************************************************

// LHAevntPythia6 class.

//*********

// Read in a file previously written from Pythia 6.3 with PYUPEV.
// It should be transparent how to apply it also to other programs.

bool LHAevntPythia6::set() {

  // Reset infile error flag.
  is.clear(); 
  
  // Read in process info and store it.
  int nup, idprup;
  double xwgtup, scalup, aqedup, aqcdup;
  is >> nup >> idprup >> xwgtup >> scalup >> aqedup >> aqcdup;
  if (!is) return false;
  process(idprup, xwgtup, scalup, aqedup, aqcdup);

  // Read in particle info one by one, and store it.
  // Note unusual C++ loop range, to better reflect LHA/Fortran standard.
  // (Recall that process(...) above added empty particle at index 0.) 
  int idup, istup, mothup1, mothup2, icolup1, icolup2; 
  double pup1, pup2, pup3, pup4, pup5, vtimup, spinup;
  for (int ip = 1; ip <= nup; ++ip) { 
    is >> idup >> istup >> mothup1 >> mothup2 >> icolup1 >> icolup2 
       >> pup1 >> pup2 >> pup3 >> pup4 >> pup5 >> vtimup >> spinup;
    if (!is) return false;   
    particle(idup, istup, mothup1, mothup2, icolup1, icolup2,
      pup1, pup2, pup3, pup4, pup5, vtimup, spinup) ;
  }

  // Reading worked.
  return true;
}
 
//**************************************************************************
