// Function definitions (not found in the header) for the Settings 
// and ErrorMessages classes.
// Copyright C 2006 Torbjorn Sjostrand

#include "Settings.h"

namespace Pythia8 {

//**************************************************************************

// Settings class.
// This class contains flags, modes and parameters used in generation.

//*********

// Definitions of static variables. 
map<string, Flag> Settings::flags;
map<string, Mode> Settings::modes;
map<string, Parameter> Settings::parameters;
bool Settings::isInit = false;

//*********

// Read in database from specific file.

bool Settings::init(string startFile) {

  // Don't initialize if it has already been done.
  if (isInit) return true;

  // List of files to be checked.
  vector<string> files;
  files.push_back(startFile);

  // Loop over files. Open them for read.
  for (int i = 0; i < int(files.size()); ++i) {
    const char* cstring = files[i].c_str();
    ifstream is(cstring);  

    // Check that instream is OK.
    if (!is) {cout << "\n Error: settings file " << files[i] 
       << " not found \n"; return false;}

    // Read in one line at a time.
    string line;
    while ( getline(is, line) ) {

      // Get first word of a line.
      istringstream getfirst(line);
      string word1;
      getfirst >> word1;
    
      // Check for occurence of a bool and add to flag map.
      if (word1 == "<flag>") {
        for ( int j = 0; j < int(line.size()); ++ j) 
          if (line[j] == '=') line[j] = ' ';
        istringstream flagData(line);
        string dummy, name, valueString;
        flagData >> dummy >> name >> valueString;
        if (!flagData) cout << " Error: incomplete bool " << name 
          << " in file " << files[i] << "\n";
        else if (isFlag(name)) cout << " Error: bool " << name 
          << " is already defined\n";
        else { 
          string valueLower = tolower(valueString);
          bool value = (valueLower == "true" || valueLower == "on"
            || valueLower == "yes" || valueLower == "ok" 
            || valueLower == "1" ) ? true : false ;
          addFlag( name, value);
        }
    
      // Check for occurence of an int and add to mode map.
      } else if (word1 == "<mode>") {
        for ( int j = 0; j < int(line.size()); ++ j) 
          if (line[j] == '=' || line[j] == ',' || line[j] == ';') 
            line[j] = ' ';
        istringstream modeData(line);
        string dummy, name, limit;
        int value;
        bool hasMin = false;
        bool hasMax = false;
        int minVal = 0;
        int maxVal = 0;
        modeData >> dummy >> name >> value;
        if (!modeData) cout << " Error: incomplete int " << name  
          << " in file " << files[i] << "\n";
        else if (isMode(name)) cout << " Error: int " << name 
          << " is already defined\n";
	else {
          do {
            modeData >> limit;
            if (modeData) {
              if (tolower(limit) == "min") 
                {hasMin = true; modeData >> minVal;}
	      else if (tolower(limit) == "max") 
                {hasMax = true; modeData >> maxVal;}
	    }
          } while (modeData); 
          addMode( name, value, hasMin, hasMax, minVal, maxVal);
	}	  
    
      // Check for occurence of a double and add to parameters map.
      } else if (word1 == "<parameter>") {
        for ( int j = 0; j < int(line.size()); ++ j) 
          if (line[j] == '=' || line[j] == ',' || line[j] == ';') 
            line[j] = ' ';
        istringstream paramData(line);
        string dummy, name, limit;
        double value;
        bool hasMin = false;
        bool hasMax = false;
        double minVal = 0;
        double maxVal = 0;
        paramData >> dummy >> name >> value;
        if (!paramData) cout << " Error: incomplete double " << name
          << " in file " << files[i] << "\n";
        else if (isParameter(name)) cout << " Error: int " << name 
          << " is already defined\n";
        else {
          do {
            paramData >> limit;
            if (paramData) {
              if (tolower(limit) == "min") 
                {hasMin = true; paramData >> minVal;}
	      else if (tolower(limit) == "max") 
                {hasMax = true; paramData >> maxVal;}
	    }
          } while (paramData); 
          addParameter( name, value, hasMin, hasMax, minVal, maxVal);
	}

      // Check for occurence of a file also to be read.
      } else if (word1 == "<file>") {
        istringstream fileData(line);
        string dummy, file;
        fileData >> dummy >> file;
        if (!fileData) cout << " Error: incomplete filename " << file
          << " in file " << files[i] << "\n";
        else files.push_back(file);
      }

    // End of loop over lines in input file and loop over files.
    };
  };

  // Done.
  isInit = true;
  return true;

}

//*********

// Overwrite existing database by reading from specific file.

bool Settings::reInit(string startFile) {

  // Reset maps to empty.
  flags.clear();
  modes.clear();
  parameters.clear();

  // Then let normal init do the rest.
  isInit = false;
  return init(startFile);

} 

//*********

// Read in updates from a character string, like a line of a file. 
// Is used by readFile, and by readString (and readFile) in Pythia.

bool Settings::readString(string line, bool warn) {

  // If first character is not a letter, then taken to be a comment line.
  int firstChar = line.find_first_not_of(" ");
  if (!isalpha(line[firstChar])) return true; 

  // Replace an equal sign by a blank to make parsing simpler.
  if (line.find("=") != string::npos) {
    int firstEqual = line.find_first_of("=");
    line.replace(firstEqual, 1, " ");   
  }

  // Get first word of a line.
  istringstream getWord(line);
  string name;
  getWord >> name;

  // Replace two colons by one (:: -> :) to allow for such mistakes.
  while (name.find("::") != string::npos) {
    int firstColonColon = name.find_first_of("::");
    name.replace(firstColonColon, 2, ":");   
  }
     
  // Check whether this is in the database. Done if not.
  int inDataBase = 0;
  if (isFlag(name)) inDataBase = 1;   
  else if (isMode(name)) inDataBase = 2;   
  else if (isParameter(name)) inDataBase = 3; 
  if (inDataBase == 0) {
    if (warn) cout << "\n Warning: input string not found in settings"
      << " databases; skip:\n   " << line << "\n";
    return false;  
  }  

  // Find value.
  string valueString;
  getWord >> valueString;

  // Update flag map; allow many ways to say yes.
  if (getWord && inDataBase == 1) {
    string valueLower = tolower(valueString);
    bool value = (valueLower == "true" || valueLower == "on"
      || valueLower == "yes" || valueLower == "ok" 
      || valueLower == "1" ) ? true : false ;
    flag(name, value);

  // Update mode map.
  } else if (getWord && inDataBase == 2) {
    istringstream modeData(valueString);
    int value;
    modeData >> value;
    if (modeData) mode(name, value);
        
  // Update parameter map.
  } else if (getWord && inDataBase == 3) {
    istringstream paramData(valueString);
    double value;
    paramData >> value;
    if (paramData) parameter(name, value);
  }

  // Done.
  return true;
}

//*********

// Read in updates from user-defined file.

 bool Settings::readFile(string updateFile, bool warn) {

  // Open file with updates.
  const char* cstring = updateFile.c_str();
  ifstream is(cstring);  
  if (!is) {cout << " Error: user update file " << updateFile 
     << " not found \n"; return false;}

  // Read in one line at a time.
  string line;
  while ( getline(is, line) ) {

    // Process the line.
    readString( line, warn);

  // Reached end of input file.
  };
  return true;
}

//*********
 
// Print out table of database in lexigraphical order.

void Settings::list(bool listAll, ostream& os) {

  // Table header; output for bool as off/on.
  if (listAll) 
    os << "\n --------  Pythia Flag + Mode + Parameter Settings (all)  -"
       << "------------------------------ \n";
  else
    os << "\n --------  Pythia Flag + Mode + Parameter Settings (changes" 
       << " only)  ---------------------- \n" ;
  os << " \n   Kind  Name                                           Now"
     << "   Default       Min       Max \n";
 
  // Iterators for the flag, mode and param tables.
  map<string, Flag>::iterator flagEntry = flags.begin();
  map<string, Mode>::iterator modeEntry = modes.begin();
  map<string, Parameter>::iterator paramEntry = parameters.begin();

  // Loop while there is something left to do.
  while (flagEntry != flags.end() || modeEntry != modes.end() 
    || paramEntry != parameters.end()) {

    // Check if a flag is next in lexigraphical order; if so print it.
    if ( flagEntry != flags.end() && ( modeEntry == modes.end() 
      || flagEntry->first < modeEntry->first ) && ( paramEntry 
      == parameters.end() || flagEntry->first < paramEntry->first ) ) {
      string state[2] = {"off", "on"};
      bool valNow = flagEntry->second.valNow;
      bool valDefault = flagEntry->second.valDefault;
      if ( listAll || valNow != valDefault )
        os << "   bool  " << setw(40) << left 
           << flagEntry->second.name << setw(10) << right
           << state[valNow] << setw(10) << state[valDefault] << "\n";
      ++flagEntry;
     
    // Else check if mode is next, and if so print it.
    } else if ( modeEntry != modes.end() && ( paramEntry  
      == parameters.end() || modeEntry->first < paramEntry->first ) ) {
      int valNow = modeEntry->second.valNow;
      int valDefault = modeEntry->second.valDefault;
      if ( listAll || valNow != valDefault ) {
        os << "    int  " << setw(40) << left 
           << modeEntry->second.name << setw(10) << right 
           << valNow << setw(10) << valDefault; 
        if (modeEntry->second.hasMin) 
          os << setw(10) << modeEntry->second.valMin; 
        else os << "          ";
        if (modeEntry->second.hasMax) 
          os << setw(10) << modeEntry->second.valMax; 
        os << "\n";
      }
      ++modeEntry;
      
    // Else print param; fixed or scientific depending on value.
    } else {
      double valNow = paramEntry->second.valNow;
      double valDefault = paramEntry->second.valDefault;      
      if ( listAll || valNow != valDefault ) {
        os << " double  " << setw(40) << left 
           << paramEntry->second.name << right;
	for (int i = 0; i < 4; ++i) { 
          if (i == 1) valNow = valDefault;  
          if (i == 2) valNow = paramEntry->second.valMin;  
          if (i == 3) valNow = paramEntry->second.valMax;  
          if ( (i == 2 && !paramEntry->second.hasMin)
	    || (i == 3 && !paramEntry->second.hasMax) )
            os << "          ";
          else if ( valNow == 0. 
            || ( abs(valNow) > 0.01 && abs(valNow) < 10000. ) )
            os << fixed << setprecision(4) << setw(10) << valNow; 
          else os << scientific << setprecision(2) << setw(10) 
            << valNow; 
	}  
        cout << "\n";
      }
      ++paramEntry;
    }
  } ;

  // End of loop over database contents.
  os << "\n --------  End Pythia Flag + Mode + Parameter Settings  ---"
     << "------------------------------ \n" << endl;

}

//*********
 
// Reset all values to their defaults.

void Settings::resetAll() {

  // Loop through the flags table, resetting all entries.
  for (map<string, Flag>::iterator flagEntry = flags.begin();
    flagEntry != flags.end(); ++flagEntry) {
    string name = flagEntry->first;
    resetFlag(name);
  }

  // Loop through the modes table, resetting all entries.
  for (map<string, Mode>::iterator modeEntry = modes.begin();
    modeEntry != modes.end(); ++modeEntry) {
    string name = modeEntry->first;
    resetMode(name);
  }

  // Loop through the parameters table, resetting all entries.
  for (map<string, Parameter>::iterator paramEntry = parameters.begin(); 
    paramEntry != parameters.end(); ++paramEntry) {
    string name = paramEntry->first;
    resetParameter(name);
  }

}

//*********
 
// Change current value, respecting limits.

void Settings::flag(string keyIn, bool nowIn) { 
    if (isFlag(keyIn)) flags[tolower(keyIn)].valNow = nowIn; 
}

void Settings:: mode(string keyIn, int nowIn) { 
  if (isMode(keyIn)) { 
    Mode& modeNow = modes[tolower(keyIn)];
    if (modeNow.hasMin && nowIn < modeNow.valMin) 
      modeNow.valNow = modeNow.valMin; 
    else if (modeNow.hasMax && nowIn > modeNow.valMax) 
      modeNow.valNow = modeNow.valMax;
    else modeNow.valNow = nowIn; 
  } 
} 

void Settings::parameter(string keyIn, double nowIn) { 
  if (isParameter(keyIn)) {
    Parameter& paramNow = parameters[tolower(keyIn)];
    if (paramNow.hasMin && nowIn < paramNow.valMin) 
      paramNow.valNow = paramNow.valMin; 
    else if (paramNow.hasMax && nowIn > paramNow.valMax) 
      paramNow.valNow = paramNow.valMax;
    else paramNow.valNow = nowIn; 
  } 
}  

//**************************************************************************

// ErrorMessages class.
// This class holds info on all messages received and how many times.

//*********

// Definitions of static variables. 
// (Values will be overwritten in initStatic call, so are purely dummy.)

int ErrorMessages::timesToPrint = 1;
map<string, int> ErrorMessages::messages;

//*********

// Initialize static data members.

void ErrorMessages::initStatic() {

  timesToPrint = Settings::mode("ErrorMessages:timesToPrint"); 

}

//*********
  
// Print a message the first few times. Insert in database.
 
void ErrorMessages::message(string messageIn, string extraIn) {
   
  // Recover number of times message occured. Also inserts new string.
  int times = messages[messageIn];
  ++messages[messageIn];

  // Print message the first few times.
  if (times < timesToPrint) cout << " " << messageIn << " " 
    << extraIn << "\n";

}

//*********

// Print statistics on errors/warnings.

void ErrorMessages::statistics(ostream& os) {

  // Header.
  os << "\n --------  Pythia warning and error messages statistics  "
     << "------------------------------ \n\n times   message \n";

  // Loop over all messages
  map<string, int>::iterator messageEntry = messages.begin();
  while (messageEntry != messages.end()) {
    os << setw(6) << messageEntry->second << "   " 
       << messageEntry->first << "\n";
    ++messageEntry;
  } 

  // Done. 
  os << "\n --------  End Pythia warning and error messages statistics"
     << "  -------------------------- \n";


}

//**************************************************************************

} // end namespace Pythia8
