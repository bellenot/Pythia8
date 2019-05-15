// LesHouches.cc is a part of the PYTHIA event generator.
// Copyright (C) 2007 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the LHAinit,
// LHAevnt, LHAinitLHEF and LHAevntLHEF classes.

#include "LesHouches.h"

namespace Pythia8 {

//**************************************************************************

// LHAinit class.

//*********

// Print the info; to check it worked.

void LHAinit::list(ostream& os) {

  // Header.
  os << "\n --------  LHA initialization information  ------------ \n"; 

  // Beam info.
  os << fixed << setprecision(3) 
     << "\n  beam    kind      energy  pdfgrp  pdfset \n" 
     << "     A  " << setw(6) << idBeamASave 
     <<  setw(12) << eBeamASave 
     << setw(8) << pdfGroupBeamASave 
     << setw(8) << pdfSetBeamASave << "\n" 
     << "     B  " << setw(6) << idBeamBSave 
     <<  setw(12) << eBeamBSave 
     << setw(8) << pdfGroupBeamBSave 
     << setw(8) << pdfSetBeamBSave << "\n"; 

  // Event weighting strategy.
  os << "\n  Event weighting strategy = " << setw(2) 
     << strategySave << "\n" ;

  // Process list.
  os << scientific << setprecision(4) 
     << "\n  Processes, with strategy-dependent cross section info \n" 
     << "  number      xsec (pb)      xerr (pb)      xmax (pb) \n" ;
  for (int ip = 0; ip < int(processes.size()); ++ip) {
    os << setw(8) << processes[ip].idProc 
       << setw(15) << processes[ip].xSecProc 
       << setw(15) << processes[ip].xErrProc 
       << setw(15) << processes[ip].xMaxProc << "\n";
  }

  // Finished.
  os << "\n --------  End LHA initialization information  -------- \n"; 

}
 
//**************************************************************************

// LHAevnt class.

//*********

// Print the info; to check it worked.

void LHAevnt::list(ostream& os) {

  // Header.
  os << "\n --------  LHA event information and listing  -------------"
     << "--------------------------------------------------------- \n"; 

  // Basic event info.
  os << scientific << setprecision(4) 
     << "\n    process = " << setw(8) << idProc 
     << "    weight = " << setw(12) << weightProc 
     << "     scale = " << setw(12) << scaleProc << " (GeV) \n"
     << "                   "
     << "     alpha_em = " << setw(12) << alphaQEDProc 
     << "    alpha_strong = " << setw(12) << alphaQCDProc << "\n";

  // Particle list
  os << fixed << setprecision(3) 
     << "\n    Participating Particles \n" 
     << "    no        id stat     mothers     colours      p_x        "
     << "p_y        p_z         e          m        tau    spin \n" ;
  for (int ip = 0; ip < int(particles.size()); ++ip) {
    os << setw(6) << ip 
       << setw(10) << particles[ip].idPart 
       << setw(5) << particles[ip].statusPart 
       << setw(6) << particles[ip].mother1Part
       << setw(6) << particles[ip].mother2Part 
       << setw(6) << particles[ip].col1Part
       << setw(6) << particles[ip].col2Part 
       << setw(11) << particles[ip].pxPart
       << setw(11) << particles[ip].pyPart
       << setw(11) << particles[ip].pzPart 
       << setw(11) << particles[ip].ePart 
       << setw(11) <<  particles[ip].mPart 
       << setw(8) <<  particles[ip].tauPart 
       << setw(8) <<  particles[ip].spinPart << "\n";
  }

  // PDF info - optional.
  if (pdfIsSetSave) os << "\n   pdf: id1 =" << setw(5) << id1Save  
    << " id2 =" << setw(5) << id2Save 
    << " x1 ="  << scientific << setw(10) << x1Save    
    << " x2 =" << setw(10) << x2Save 
    << " scalePDF =" << setw(10) << scalePDFSave 
    << " xpdf1 =" << setw(10) << xpdf1Save    
    << " xpdf2 =" << setw(10) << xpdf2Save << "\n";    

  // Finished.
  os << "\n --------  End LHA event information and listing  ---------"
     << "--------------------------------------------------------- \n"; 

}
 
//**************************************************************************

// LHAinitLHEF class.

//*********

// Read in a Les Houches Event File.

bool LHAinitLHEF::set() {

  // Check that first line is consistent with proper LHEF file.
  string line;
  if (!getline(is, line)) return false;
  if (line.find("<LesHouchesEvents") == string::npos) return false;  
  if (line.find("version=\"1.0\"" ) == string::npos ) return false;
 
  // Loop over lines until an <init tag is found first on a line.
  string tag = " ";
  do { 
    if (!getline(is, line)) return false;
    if (line.find_first_not_of(" ") != string::npos) {
      istringstream getfirst(line);
      getfirst >> tag;
      if (!getfirst) return false;
    }
  } while (tag != "<init>" && tag != "<init"); 
  
  // Read in beam and strategy info, and store it. 
  int idbmupA, idbmupB;
  double ebmupA, ebmupB;
  int pdfgupA, pdfgupB, pdfsupA, pdfsupB, idwtup, nprup;
  if (!getline(is, line)) return false;
  istringstream getbms(line);
  getbms >> idbmupA >> idbmupB >> ebmupA >> ebmupB >> pdfgupA 
     >> pdfgupB >> pdfsupA >> pdfsupB >> idwtup >> nprup;
  if (!getbms) return false;
  beamA(idbmupA, ebmupA, pdfgupA, pdfsupA);
  beamB(idbmupB, ebmupB, pdfgupB, pdfsupB);
  strategy(idwtup);

  // Read in process info, one process at a time, and store it.
  double xsecup, xerrup, xmaxup;
  int lprup; 
  for (int ip = 0; ip < nprup; ++ip) { 
    if (!getline(is, line)) return false;
    istringstream getpro(line);
    getpro >> xsecup >> xerrup >> xmaxup >> lprup ;
    if (!getpro) return false;
    process(lprup, xsecup, xerrup, xmaxup);
  }

  // Reading worked.
  return true;

}
 
//**************************************************************************

// LHAevntLHEF class.

//*********

// Read in a Les Houches Event File.

bool LHAevntLHEF::set( int ) {
  
  // Loop over lines until an <event tag is found first on a line.
  string line, tag;
  do { 
    if (!getline(is, line)) return false;
    if (line.find_first_not_of(" ") != string::npos) {
      istringstream getfirst(line);
      getfirst >> tag;
      if (!getfirst) return false;
    }
  } while (tag != "<event>" && tag != "<event"); 

  // Read in process info and store it.
  int nup, idprup;
  double xwgtup, scalup, aqedup, aqcdup;
  if (!getline(is, line)) return false;
  istringstream getpro(line);
  getpro >> nup >> idprup >> xwgtup >> scalup >> aqedup >> aqcdup;
  if (!getpro) return false;
  process(idprup, xwgtup, scalup, aqedup, aqcdup);

  // Read in particle info one by one, and store it.
  // Note unusual C++ loop range, to better reflect LHA/Fortran standard.
  // (Recall that process(...) above added empty particle at index 0.) 
  int idup, istup, mothup1, mothup2, icolup1, icolup2; 
  double pup1, pup2, pup3, pup4, pup5, vtimup, spinup;
  for (int ip = 1; ip <= nup; ++ip) { 
    if (!getline(is, line)) return false;
    istringstream getall(line);
    getall >> idup >> istup >> mothup1 >> mothup2 >> icolup1 >> icolup2 
      >> pup1 >> pup2 >> pup3 >> pup4 >> pup5 >> vtimup >> spinup;
    if (!getall) return false;   
    particle(idup, istup, mothup1, mothup2, icolup1, icolup2,
      pup1, pup2, pup3, pup4, pup5, vtimup, spinup) ;
  }

  // Continue parsing till </event>. Extract pdf info if present.
  do { 
    if (!getline(is, line)) return false;
    istringstream getpdf(line);
    getpdf >> tag;
    if (!getpdf) return false;
    if (tag == "#pdf") {
      int id1, id2;
      double x1, x2, scalePDF, xpdf1, xpdf2;
      getpdf >> id1 >> id2 >>  x1 >> x2 >> scalePDF >> xpdf1 >> xpdf2;
      if (!getpdf) return false;
      pdf(id1, id2, x1, x2, scalePDF, xpdf1, xpdf2);  
    }
  } while (tag != "</event>" && tag != "</event"); 
  

  // Reading worked.
  return true;
}

//**************************************************************************

} // end namespace Pythia8
