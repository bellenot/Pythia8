// Header file for the settings database; and for error statistics.
// Flag: helper class with bool flags.
// Mode: helper class with int modes.
// Parameter: helper class with double parameters.
// Settings: maps of flags, modes and parameters with input/output.
// Copyright C 2006 Torbjorn Sjostrand

#ifndef Pythia8_Settings_H
#define Pythia8_Settings_H

#include "Stdlib.h"

namespace Pythia8 {

//**************************************************************************

// Class for bool flags.

class Flag {

public:

  // Constructor
  Flag(string nameIn = " ", bool defaultIn = false) : name(nameIn), 
    valNow(defaultIn) , valDefault(defaultIn) { }

  // Data members.
  string name;
  bool valNow, valDefault;

};

//**************************************************************************

// Class for integer modes.

class Mode {

public:

  // Constructor
  Mode(string nameIn = " ", int defaultIn = 0, bool hasMinIn = false,
    bool hasMaxIn = false, int minIn = 0,  int maxIn = 0) :  name(nameIn), 
    valNow(defaultIn), valDefault(defaultIn), hasMin(hasMinIn),
    hasMax(hasMaxIn), valMin(minIn), valMax(maxIn) { }

  // Data members.
  string name;
  int valNow, valDefault;
  bool hasMin, hasMax;
  int valMin, valMax;

};

//**************************************************************************

// Class for double parameters.

class Parameter {

public:

  // Constructor
  Parameter(string nameIn = " ", double defaultIn = 0., 
    bool hasMinIn = false, bool hasMaxIn = false, double minIn = 0., 
    double maxIn = 0.) :  name(nameIn), valNow(defaultIn), 
    valDefault(defaultIn), hasMin(hasMinIn), hasMax(hasMaxIn), 
    valMin(minIn), valMax(maxIn) { }

  // Data members.
  string name;
  double valNow, valDefault;
  bool hasMin, hasMax;
  double valMin, valMax;

};

//**************************************************************************

// This class holds info on flags (bool), modes (int) and 
// parameters (double).

class Settings {

public:

  // Constructor.
  Settings() {}
 
  // Read in database from specific file.
  static bool init(string startFile = "../doc/Index.xml") ;

  // Overwrite existing database by reading from specific file.
  static bool reInit(string startFile = "../doc/Index.xml") ;

  // Read in one update from a single line.
  static bool readString(string line, bool warn = true) ; 
 
  // Read in updates from user-defined file.
  static bool readFile(string updateFile, bool warn = true) ;

  // Print out table of database, either all or only changed ones,
  // or ones containing a given string.
  static void listAll(ostream& os = cout) { 
    list( true, false, " ", os); } 
  static void listChanged(ostream& os = cout) { 
    list (false, false, " ", os); } 
  static void list(string match, ostream& os = cout) { 
    list (false, true, match, os); } 

  // Reset all values to their defaults.
  static void resetAll() ;

  // Query existence of an entry.
  static bool isFlag(string keyIn) {
    return (flags.find(tolower(keyIn)) == flags.end()) ? false : true ; }
  static bool isMode(string keyIn) { 
    return (modes.find(tolower(keyIn)) == modes.end()) ? false : true ; }
  static bool isParameter(string keyIn) {
    return (parameters.find(tolower(keyIn)) == parameters.end()) ? false 
    : true; }
 
  // Add new entry.
  static void addFlag(string keyIn, bool defaultIn) {
    flags[tolower(keyIn)] = Flag(keyIn, defaultIn); }  
  static void addMode(string keyIn, int defaultIn, bool hasMinIn, 
    bool hasMaxIn, int minIn, int maxIn) { modes[tolower(keyIn)] 
    = Mode(keyIn, defaultIn, hasMinIn, hasMaxIn, minIn, maxIn); }      
  static void addParameter(string keyIn, double defaultIn, bool hasMinIn, 
    bool hasMaxIn, double minIn, double maxIn) { parameters[tolower(keyIn)] 
    = Parameter(keyIn, defaultIn, hasMinIn, hasMaxIn, minIn, maxIn); }  

  // Give back current value. 
  static bool flag(string keyIn) {
    return isFlag(keyIn) ? flags[tolower(keyIn)].valNow : false ; } 
  static int mode(string keyIn) { 
    return isMode(keyIn) ? modes[tolower(keyIn)].valNow : 0 ; }
  static double parameter(string keyIn) {
    return isParameter(keyIn) ? parameters[tolower(keyIn)].valNow : 0 ; }
  
  // Change current value, respecting limits.
  static void flag(string keyIn, bool nowIn); 
  static void mode(string keyIn, int nowIn);
  static void parameter(string keyIn, double nowIn); 

  // Change current value, disregarding limits.
  static void forceMode(string keyIn, int nowIn) { 
    if (isMode(keyIn)) modes[tolower(keyIn)].valNow = nowIn; }
  static void forceParameter(string keyIn, double nowIn) { 
    if (isParameter(keyIn)) parameters[tolower(keyIn)].valNow = nowIn; }
     
  // Restore current value to default. 
  static void resetFlag(string keyIn) {
    if (isFlag(keyIn)) flags[tolower(keyIn)].valNow 
      = flags[tolower(keyIn)].valDefault ; }
  static void resetMode(string keyIn) {
    if (isMode(keyIn)) modes[tolower(keyIn)].valNow 
      = modes[tolower(keyIn)].valDefault ; }
  static void resetParameter(string keyIn) {
    if (isParameter(keyIn)) parameters[tolower(keyIn)].valNow 
      = parameters[tolower(keyIn)].valDefault ; }

private:

  // Map for bool flags.
  static map<string, Flag> flags;

  // Map for integer modes.
  static map<string, Mode> modes;

  // Map for double parameters.
  static map<string, Parameter> parameters;

  // Flag that initialization has been performed.
  static bool isInit;

  // Print out table of database, called from listAll and listChanged.
  static void list(bool listAll, bool listString, string match,
    ostream& os = cout) ; 

};

//**************************************************************************

} // end namespace Pythia8

#endif // Pythia8_Settings_H
