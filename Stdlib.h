// Header file for functionality pulled in from stdlib,
// plus powers of small integers and strings to lowercase. 
// Copyright © 2005 Torbjörn Sjöstrand

#ifndef Pythia8_Stdlib_H
#define Pythia8_Stdlib_H

// Stdlib header files for character manipulation.
#include <cctype>

// Stdlib header files for math and time.
#include <cmath>
#include <ctime>

// Stdlib header files for containers.
#include <vector>
#include <map>

// Stdlib header files for input/output.
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

// Generic utilities.
using std::tolower; 
using std::swap;

// Mathematical functions.
using std::max;
using std::min; 
using std::abs; 

// Standard containers.
using std::string; 
using std::vector; 
using std::map; 

// Input/output.
using std::cin; 
using std::cout; 
using std::cerr; 
using std::endl; 
using std::fixed; 
using std::scientific; 
using std::setw; 
using std::setprecision; 
using std::ios; 
using std::istream; 
using std::ostream; 
using std::ifstream; 
using std::ofstream; 
using std::istringstream; 
using std::ostringstream; 

// Powers of small integers - for balance speed/code clarity.
inline double pow2(const double& x) {return x*x;}
inline double pow3(const double& x) {return x*x*x;}
inline double pow4(const double& x) {return x*x*x*x;}
inline double pow5(const double& x) {return x*x*x*x*x;}

// Avoid problem with negative square root argument (from roundoff).
inline double sqrtpos(const double& x) {return sqrt( max( 0., x));}

// Convert string to lowercase for case-insensitive comparisons.
inline string tolower(const string& name) { string temp(name);
  for (int i = 0; i < int(temp.length()); ++i) temp[i] = tolower(temp[i]); 
  return temp; }

#endif // Pythia8_Stdlib_H
