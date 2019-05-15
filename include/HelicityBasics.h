// HelicityBasics.h is a part of the PYTHIA event generator.
// Copyright (C) 2011 Philip Ilten, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for a number of helper classes used in tau decays.

#ifndef Pythia8_HelicityBasics_H
#define Pythia8_HelicityBasics_H

#include "Basics.h"
#include "Event.h"
#include "PythiaComplex.h"
#include "PythiaStdlib.h"

namespace Pythia8 {

//==========================================================================
  
// The Wave4 class provides a class for complex four-vector wave functions.
// The Wave4 class can be multiplied with the GammaMatrix class to allow
// for the writing of helicity matrix elements.
  
class Wave4 {
  
  protected:

    complex values[4];

  public:

    // Constructors and destructors.
    Wave4() {};
    Wave4(complex v0, complex v1, complex v2, complex v3) {values[0] = v0;
      values[1] = v1; values[2] = v2; values[3] = v3;}
    Wave4(Vec4 v) {values[0] = v.e(); values[1] = v.px(); values[2] = v.py();
      values[3] = v.pz();}
    ~Wave4() {};

    // Access an element of the wave vector.
    complex& operator() (int i) {return values[i];}
    // Wave4 + Wave4
    Wave4 operator+(Wave4 w) {return Wave4(values[0] + w.values[0],
					   values[1] + w.values[1],
					   values[2] + w.values[2],
					   values[3] + w.values[3]);}
    // Wave4 - Wave4
    Wave4 operator-(Wave4 w) {return Wave4(values[0] - w.values[0],
					   values[1] - w.values[1],
					   values[2] - w.values[2],
					   values[3] - w.values[3]);}
    // - Wave4
    Wave4 operator-() {return Wave4(-values[0], -values[1], -values[2], 
				    -values[3]);}
    // Wave4 * Wave4
    complex operator*(Wave4 w) {return values[0]*w.values[0]
	+ values[1]*w.values[1] + values[2]*w.values[2]
	+ values[3]*w.values[3];}
    // Wave4 * complex
    Wave4 operator*(complex s) {return Wave4(values[0]*s, values[1]*s,
					     values[2]*s, values[3]*s);}
    // complex * Wave4
    friend Wave4 operator*(complex s, const Wave4& w); 
    // Wave4 / complex
    Wave4 operator/(complex s) {return Wave4(values[0]/s, values[1]/s,
					     values[2]/s, values[3]/s);}
    // Wave4 / double
    Wave4 operator/(double s) {return Wave4(values[0]/s, values[1]/s,
					    values[2]/s, values[3]/s);}
    // Complex conjugate.
    friend Wave4 conj(Wave4 w) ;
    // Invariant squared mass for REAL Wave4 (to save time).
    friend double m2(Wave4 w1, Wave4 w2);
    
    // Wave4 * GammaMatrix multiplication is defined in the GammaMatrix class.
    // Print a Wave4 vector.
    friend ostream& operator<<(ostream& output, Wave4 w);

};

//--------------------------------------------------------------------------

// Namespace function declarations; friends of Wave4 class.
Wave4 operator*(complex s, const Wave4& w);
Wave4 conj(Wave4 w);
double m2(Wave4 w1, Wave4 w2);
ostream& operator<< (ostream& output, Wave4 w);

//==========================================================================

// The GammaMatrix class is a special sparse matrix class used to write
// helicity matrix elements in conjuction with the Wave4 class. Note that
// only left to right multplication of Wave4 vectors with the GammaMatrix
// class is allowed. Additionally, subtracting a scalar from a GammaMatrix
// (or subtracting a GammaMatrix from a scalar) subtracts the scalar from 
//each non-zero element of the GammaMatrix. This is designed specifically 
// with the (1 - gamma^5) structure of matrix elements in mind.
 
class GammaMatrix {
  
  protected:

    complex values[4];
    int     index[4];
    complex zero;

  public:

    // Constructors and destructors.
    GammaMatrix() {};
    GammaMatrix(int mu);
    ~GammaMatrix() {};

    // Access an element of the matrix.
    complex& operator() (int I, int J) {
      if (index[J] == I) return values[J];
      else return zero;
    }

    // Wave4 * GammaMatrix
    friend Wave4 operator*(Wave4 w, GammaMatrix g);

    // GammaMatrix * Scalar
    GammaMatrix operator*(complex s) {values[0] = s*values[0];
      values[1] = s*values[1]; values[2] = s*values[2];
      values[3] = s*values[3]; return *this;}

    // Scalar * GammaMatrix
    friend GammaMatrix operator*(complex s, GammaMatrix g); 
    
    // Gamma5 - I*Scalar
    GammaMatrix operator-(complex s) {values[0] = values[0]-s;
      values[1] = values[1]-s; values[2] = values[2]-s;
      values[3] = values[3]-s; return *this;}

    // I*Scalar - Gamma5
    friend GammaMatrix operator-(complex s, GammaMatrix g);

    // Gamma5 + I*Scalar
    GammaMatrix operator+(complex s) {values[0] = values[0]+s;
      values[1] = values[1]+s; values[2] = values[2]+s;
      values[3] = values[3]+s; return *this;}

    // I*Scalar + Gamma5
    friend GammaMatrix operator+(complex s, GammaMatrix g);

    // << GammaMatrix
    friend ostream& operator<< (ostream& output, GammaMatrix g);

};

//--------------------------------------------------------------------------

// Namespace function declarations; friends of GammaMatrix class.
Wave4 operator*(Wave4 w, GammaMatrix g);
GammaMatrix operator*(complex s, GammaMatrix g); 
GammaMatrix operator-(complex s, GammaMatrix g); 
GammaMatrix operator+(complex s, GammaMatrix g); 
ostream& operator<< (ostream& output, GammaMatrix g);

//==========================================================================

// Helicity particle class containing helicity information, derived from
// particle base class.
 
class HelicityParticle : public Particle {
    
  public:
    
    // Event record position.
    int idx;
    
    // Flag for whether particle is incoming (-1) or outgoing (1).
    int direction;
    
    // Helicity density matrix.
    vector< vector<complex> > rho;
    
    // Decay matrix.
    vector< vector<complex> > D;
    
    // Constructors.
  HelicityParticle() : Particle() {
      direction = 1;
    };
    
  HelicityParticle(int idIn, int statusIn = 0, int mother1In = 0,
		   int mother2In = 0, int daughter1In = 0, int daughter2In = 0,
		   int colIn = 0, int acolIn = 0, double pxIn = 0.,
		   double pyIn = 0., double pzIn = 0., double eIn = 0.,
		   double mIn = 0., double scaleIn = 0., ParticleData* ptr = 0)
    : Particle(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In,
	       colIn, acolIn, pxIn, pyIn, pzIn, eIn, mIn, scaleIn) {
      if (ptr) {
	setPDTPtr(ptr);
	setPDEPtr();
      }
      rho = vector< vector<complex> >(spinType(),
				      vector<complex>(spinType(), 0));
      D   = vector< vector<complex> >(spinType(),
				      vector<complex>(spinType(), 0));
      for (int i = 0; i < spinType(); i++) {
	rho[i][i] = 0.5; D[i][i] = 1;
      }
      direction = 1;
    };
    
  HelicityParticle(int idIn, int statusIn, int mother1In, int mother2In,
		   int daughter1In, int daughter2In, int colIn, int acolIn,
		   Vec4 pIn, double mIn = 0., double scaleIn = 0.,
		   ParticleData* ptr = 0)
    : Particle(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In,
	       colIn, acolIn, pIn, mIn, scaleIn) {
      if (ptr) {
	setPDTPtr(ptr);
	setPDEPtr();
      }
      rho = vector< vector<complex> >(spinType(),
				      vector<complex>(spinType(), 0));
      D   = vector< vector<complex> >(spinType(),
				      vector<complex>(spinType(), 0));
      for (int i = 0; i < spinType(); i++) {
	rho[i][i] = 0.5; D[i][i] = 1;
      }
      direction = 1;
    };
    
  HelicityParticle(const Particle& ptIn,
		   ParticleData* ptr = 0) : Particle(ptIn) {
      if (ptr) {
	setPDTPtr(ptr);
	setPDEPtr();
      }
      rho = vector< vector<complex> >(spinType(),
				      vector<complex>(spinType(), 0));
      D   = vector< vector<complex> >(spinType(),
				      vector<complex>(spinType(), 0));
      for (int i = 0; i < spinType(); i++) {
	rho[i][i] = 0.5; D[i][i] = 1;
      }
      direction = 1;
    };
    
    // Methods.
    Wave4 wave(int h);
    Wave4 waveBar(int h);
    void normalize(vector< vector<complex> >& m);

};

//==========================================================================
  
} // end namespace Pythia8

#endif // end Pythia8_HelicityBasics_H
