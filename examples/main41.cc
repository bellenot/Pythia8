// main41.cc is a part of the PYTHIA event generator.
// Copyright (C) 2007 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Test of LHAPDF interface and whether PDF's behave sensibly.
// Histogram F(x, Q2) = (9/4) xg(x, Q2) + sum_{i = q, qbar} xf_i(x, Q2)
// for range 10^{-8} < x < 1 and for Q2 = 4 and 100.

#include "Pythia.h"

using namespace Pythia8; 

int main() {

  // Pointer to PDF sets: old reference and new to try.
  PDF* oldPDF = new CTEQ5L(2212);
  //PDF* newPDF = new LHAPDF(2212, "cteq5l.LHgrid", 0);
  PDF* newPDF = new LHAPDF(2212, "cteq61.LHgrid", 0);
  //PDF* newPDF = new LHAPDF(2212, "cteq61.LHpdf", 0);
  //PDF* newPDF = new LHAPDF(2212, "MRST2004nlo.LHgrid", 0);
  //PDF* newPDF = new LHAPDF(2212, "Alekhin_100.LHpdf", 0);

  // Set x and Q2 limits. May be required to get sensible small-x shape.
  //newPDF->setLimits( 1e-5, 1., 1., 1e10);

  // Histogram.
  Hist oldF4("F( x, Q2 = 4) old", 80 , -8., 0.);
  Hist newF4("F( x, Q2 = 4) new", 80 , -8., 0.);
  Hist ratF4("F( x, Q2 = 4) new/old", 80 , -8., 0.);
  Hist oldF100("F( x, Q2 = 100) old", 80 , -8., 0.);
  Hist newF100("F( x, Q2 = 100) new", 80 , -8., 0.);
  Hist ratF100("F( x, Q2 = 100) new/old", 80 , -8., 0.);

  // Loop over two Q2 values.
  for (int iQ = 0; iQ < 2; ++iQ) {
    double Q2 = (iQ == 0) ? 4. : 100;
    
    // Loop over x values, in a logarithmic scale
    for (int iX = 0; iX < 80; ++iX) {
      double xLog = -(0.1 * iX + 0.05);
      double x = pow( 10., xLog); 

      // Evaluate old summed PDF.
      double oldSum = 2.25 * oldPDF->xf( 21, x, Q2);   
      for (int i = 1; i < 6; ++i) 
        oldSum += oldPDF->xf( i, x, Q2) + oldPDF->xf( -i, x, Q2);  
      if (iQ == 0) oldF4.fill ( xLog, oldSum ); 
      else       oldF100.fill ( xLog, oldSum ); 

      // Evaluate new summed PDF.
      double newSum = 2.25 * newPDF->xf( 21, x, Q2);   
      for (int i = 1; i < 6; ++i) 
        newSum += newPDF->xf( i, x, Q2) + newPDF->xf( -i, x, Q2);  
      if (iQ == 0) newF4.fill ( xLog, newSum ); 
      else       newF100.fill ( xLog, newSum ); 

    //End loops over x and Q2 values.
    }
  } 

  // Done.
  ratF4 = newF4 / oldF4;
  ratF100 = newF100 / oldF100;
  cout << oldF4 << newF4 << ratF4 << oldF100 << newF100 << ratF100;

  return 0;
}
