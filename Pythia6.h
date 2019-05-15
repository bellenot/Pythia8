// Header file for the Pythia 6.3 f77 external linkage to C++.
// Copyright © 2005 Torbjörn Sjöstrand

#ifndef Pythia6_H
#define Pythia6_H

namespace Pythia8 {

//**************************************************************************

// Declare the f77 subroutines that will be used.

extern "C" {

  extern void pygive_(const char*, int);

  extern void pyinit_(const char*, const char*, const char*, double&,
    int, int, int);

  extern void pyupev_();

  extern void pylist_(int&);

  extern void pystat_(int&);

}

//**************************************************************************

// Interfaces to the above routines, to make the C++ calls similar to f77.

class Pythia6 {

public:

  static void pygive(const string cmnd) { const char* cstring = cmnd.c_str();
    int len = cmnd.length(); pygive_(cstring, len);}

  static void pyinit(const string frame, const string beam, 
    const string target, double wIn) { 
    const char* cframe = frame.c_str(); int lenframe = frame.length();
    const char* cbeam = beam.c_str(); int lenbeam = beam.length();
    const char* ctarget = target.c_str(); int lentarget = target.length();
    pyinit_(cframe, cbeam, ctarget, wIn, lenframe, lenbeam, lentarget); }

  static void pyupev() {pyupev_();}

  static void pylist(int mode) {pylist_(mode);}

  static void pystat(int mode) {pystat_(mode);}

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia6_H
