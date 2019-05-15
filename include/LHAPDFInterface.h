// Header file for the LHAPDF f77 external linkage to C++.
// All required code is contained here, i.e. there is no matching .cc file.
// Copyright C 2007 Torbjorn Sjostrand

#ifndef Pythia8_LHAPDFInterface_H
#define Pythia8_LHAPDFInterface_H

namespace Pythia8 {
 
//**************************************************************************

// Declare the LHAPDF f77 subroutines that are needed.

extern "C" {

  extern void initpdfsetm_(int&, const char*, int);

  extern void initpdfsetbynamem_(int&, const char*, int);

  extern void initpdfm_(int&, int&);

  extern void evolvepdfm_(int&, double&, double&, double*);
    
}

//**************************************************************************

// Interfaces to the above routines, to make the C++ calls similar to f77.

class LHAPDFInterface {

public:

  // Initialize set with full parthname, allowing multiple sets.
  static void initPDFsetM( int nSet, string name) {
    const char* cName = name.c_str(); int lenName = name.length();
    initpdfsetm_( nSet, cName, lenName);
  }

  // Initialize set with simple name, allowing multiple sets.
  static void initPDFsetByNameM( int nSet, string name) {
    const char* cName = name.c_str(); int lenName = name.length();
    initpdfsetbynamem_( nSet, cName, lenName);
  }

  // Initialize member of set.
  static void initPDFM(int nSet, int member) {
    initpdfm_(nSet, member);
  }

  // Evaluate x f_i(x, Q).
  static void evolvePDFM( int nSet, double x, double Q, double* xfArray) {
    evolvepdfm_( nSet, x, Q, xfArray);
  }

};
 
//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_LHAPDFInterface_H
