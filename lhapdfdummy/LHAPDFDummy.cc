// Dummy routines to link when LHAPDF not linked.
// Copyright C 2007 Torbjorn Sjostrand

extern "C" {

  void initpdfsetm_(int& nSet, const char*, int) {nSet = -1;}

  void initpdfsetbynamem_(int& nSet, const char*, int) {nSet = -1;}

  void initpdfm_(int& nSet, int&) {nSet = -1;}

  void evolvepdfm_(int& nSet, double&, double&, double*) {nSet = -1;}

}
