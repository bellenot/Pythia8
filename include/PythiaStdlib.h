// Header file for functionality pulled in from Stdlib,
// plus a few useful untilities (small powers, string manipulation).
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_PythiaStdlib_H
#define Pythia8_PythiaStdlib_H

// Stdlib header files for string and character manipulation.
#include <string>
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

// Stdlib complex numbers specialized to double precision.
#include <complex>
typedef std::complex<double> complex;

// Define pi if not yet done.
#ifndef M_PI
#define M_PI 3.1415926535897932385
#endif

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

// Input/output streams.
using std::cin; 
using std::cout; 
using std::cerr; 
using std::ios; 
using std::istream; 
using std::ostream; 
using std::ifstream; 
using std::ofstream; 
using std::istringstream; 
using std::ostringstream; 

// Input/output formatting.
using std::endl; 
using std::fixed; 
using std::scientific; 
using std::left; 
using std::right; 
using std::setw; 
using std::setprecision; 

// Powers of small integers - for balance speed/code clarity.
inline double pow2(const double& x) {return x*x;}
inline double pow3(const double& x) {return x*x*x;}
inline double pow4(const double& x) {return x*x*x*x;}
inline double pow5(const double& x) {return x*x*x*x*x;}

// Avoid problem with negative square root argument (from roundoff).
inline double sqrtpos(const double& x) {return sqrt( max( 0., x));}

// Convert string to lowercase for case-insensitive comparisons.
inline string tolower(const string& name) { 
  string temp(name);
  for (int i = 0; i < int(temp.length()); ++i) temp[i] = tolower(temp[i]); 
  return temp; 
}

// Allow several alternative inputs for true/false.
inline bool boolString(string tag) {
  string tagLow = tolower(tag);
  return (tagLow == "true" || tagLow == "1" || tagLow == "on" 
  || tagLow == "yes" || tagLow == "ok" ) ? true : false ; 
}  

// Extract XML value string following XML attribute.
inline string attributeValue(string line, string attribute) {
  if (line.find(attribute) == string::npos) return "";
  int iBegAttri = line.find(attribute); 
  int iBegQuote = line.find("\"", iBegAttri + 1);
  int iEndQuote = line.find("\"", iBegQuote + 1);
  return line.substr(iBegQuote + 1, iEndQuote - iBegQuote - 1);
}

// Extract XML bool value following XML attribute.
inline bool boolAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return false;
  return boolString(valString);   
}

// Extract XML int value following XML attribute.
inline int intAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return 0; istringstream valStream(valString);
  int intVal; valStream >> intVal; return intVal;     
}

// Extract XML double value following XML attribute.
inline double doubleAttributeValue(string line, string attribute) {
  string valString = attributeValue(line, attribute);
  if (valString == "") return 0.; istringstream valStream(valString);
  double doubleVal; valStream >> doubleVal; return doubleVal;     
}

#endif // Pythia8_PythiaStdlib_H
