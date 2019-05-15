// Function definitions (not found in the header) for the LHAinit,
// LHAevnt, LHAinitFortran, LHAevntFortran, LHAinitPythia6 and 
// LHAevntPythia6 classes.
// Copyright C 2006 Torbjorn Sjostrand

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
     << "     A  " << setw(6) << idBeamAsave 
     <<  setw(12) << eBeamAsave 
     << setw(8) << pdfGroupBeamAsave 
     << setw(8) << pdfSetBeamAsave << "\n" 
     << "     B  " << setw(6) << idBeamBsave 
     <<  setw(12) << eBeamBsave 
     << setw(8) << pdfGroupBeamBsave 
     << setw(8) << pdfSetBeamBsave << "\n"; 

  // Event weighting strategy.
  os << "\n  Event weighting strategy = " << setw(2) 
     << strategySave << "\n" ;

  // Process list.
  os << scientific << setprecision(4) 
     << "\n  Processes, with strategy-dependent cross section info \n" 
     << "  number      xsec (pb)      xerr (pb)      xmax (pb) \n" ;
  for (int ip = 0; ip < int(processes.size()); ++ip) {
    os << setw(8) << processes[ip].idPr 
       << setw(15) << processes[ip].xSecPr 
       << setw(15) << processes[ip].xErrPr 
       << setw(15) << processes[ip].xMaxPr << "\n";
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
     << "\n    process = " << setw(8) << idPr 
     << "    weight = " << setw(12) << weightPr 
     << "     scale = " << setw(12) << scalePr << " (GeV) \n"
     << "                   "
     << "     alpha_em = " << setw(12) << alphaQEDPr 
     << "    alpha_strong = " << setw(12) << alphaQCDPr << "\n";

  // Particle list
  os << fixed << setprecision(3) 
     << "\n    Participating Particles \n" 
     << "    no        id stat     mothers     colours      p_x        "
     << "p_y        p_z         e          m        tau    spin \n" ;
  for (int ip = 0; ip < int(particles.size()); ++ip) {
    os << setw(6) << ip 
       << setw(10) << particles[ip].idPa 
       << setw(5) << particles[ip].statusPa 
       << setw(6) << particles[ip].mother1Pa
       << setw(6) << particles[ip].mother2Pa 
       << setw(6) << particles[ip].col1Pa
       << setw(6) << particles[ip].col2Pa 
       << setw(11) << particles[ip].pxPa
       << setw(11) << particles[ip].pyPa
       << setw(11) << particles[ip].pzPa 
       << setw(11) << particles[ip].ePa 
       << setw(11) <<  particles[ip].mPa 
       << setw(8) <<  particles[ip].tauPa 
       << setw(8) <<  particles[ip].spinPa << "\n";
  }

  // PDF info - optional.
  if (pdfIsSetSv) os << "\n   pdf: id1 =" << setw(5) << id1Sv  
    << " id2 =" << setw(5) << id2Sv 
    << " x1 ="  << scientific << setw(10) << x1Sv    
    << " x2 =" << setw(10) << x2Sv 
    << " scalePDF =" << setw(10) << scalePDFSv 
    << " xpdf1 =" << setw(10) << xpdf1Sv    
    << " xpdf2 =" << setw(10) << xpdf2Sv << "\n";    

  // Finished.
  os << "\n --------  End LHA event information and listing  ---------"
     << "--------------------------------------------------------- \n"; 

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

  // The following is used to transfer pdf info from Pythia 6.4.
  // For another program it would be of no use, but also of no real harm.
  // Is optional, so can be commented out, if desired.
  extern struct {
    int mstp[200];
    double parp[200];
    int msti[200];
    double pari[200];
  } pypars_;    

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

  // Store pdf info. This part works only for Pythia 6.4 
  // (cf. extern pypars_ above) but should do no real harm for others.
  // Is optional, so can be commented out, if desired.
  pdf(pypars_.msti[14], pypars_.msti[15], 
      pypars_.pari[32], pypars_.pari[33], pypars_.pari[22], 
      pypars_.pari[28], pypars_.pari[29]);

  // Done.
  return true;

}

//**************************************************************************

// LHAinitLHEF class.

//*********

// Read in a Les Houches Event File.

bool LHAinitLHEF::set() {

  // Check that first line is consistent with proper LHEF file.
  string line, tag;
  getline(is, line);
  if (line.find("<LesHouchesEvents") == string::npos) return false;  
  if (line.find("version=\"1.0\"" ) == string::npos ) return false;
 
  // Loop over lines until an <init tag is found first on a line.
  do { 
    getline(is, line);
    istringstream getfirst(line);
    getfirst >> tag;
    if (!getfirst) return false;
  } while (tag != "<init>" && tag != "<init"); 
  
  // Read in beam and strategy info, and store it. 
  int idbmupA, idbmupB;
  double ebmupA, ebmupB;
  int pdfgupA, pdfgupB, pdfsupA, pdfsupB, idwtup, nprup;
  getline(is, line);
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
    getline(is, line);
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

bool LHAevntLHEF::set() {
  
  // Loop over lines until an <event tag is found first on a line.
  string line, tag;
  do { 
    getline(is, line);
    istringstream getfirst(line);
    getfirst >> tag;
    if (!getfirst) return false;
  } while (tag != "<event>" && tag != "<event"); 
  
  // Read in process info and store it.
  int nup, idprup;
  double xwgtup, scalup, aqedup, aqcdup;
  getline(is, line);
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
    getline(is, line);
    istringstream getall(line);
    getall >> idup >> istup >> mothup1 >> mothup2 >> icolup1 >> icolup2 
      >> pup1 >> pup2 >> pup3 >> pup4 >> pup5 >> vtimup >> spinup;
    if (!getall) return false;   
    particle(idup, istup, mothup1, mothup2, icolup1, icolup2,
      pup1, pup2, pup3, pup4, pup5, vtimup, spinup) ;
  }

  // Continue parsing till </event>. Extract pdf info.
  do { 
    getline(is, line);
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
