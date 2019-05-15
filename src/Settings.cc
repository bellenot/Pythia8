// Function definitions (not found in the header) for the Settings class.
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

  // List of files to be checked. Start with input file. 
  vector<string> files;
  files.push_back(startFile);

  // If nontrivial startfile path, then use that for other files as well.
  string pathName = "";
  if (startFile.rfind("/") != string::npos) 
    pathName = startFile.substr(0, startFile.rfind("/") + 1);

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

      // Skip ahead if not interesting. Only look for new files in startfile.
      if (word1 != "<flag" && word1 != "<mode" && word1 != "<parameter"
        && (word1 != "<a" || i > 0) ) continue;

      // Read and append continuation line if line does not contain >.
      if (line.find('>') == string::npos) {
        string line2;
        getline(is, line2);
        line += " " + line2;
      }
      
      // Remove extra blanks before an = sign.
      while (line.find(" =") != string::npos) line.erase( line.find(" ="), 1);

      // Find name token.
      string nameToken = (word1 != "<a") ? "name=" : "href=";
      if (line.find(nameToken) == string::npos) {
        cout << " Error: failed to find name token in line " << line;
        continue;
      }        
      int iName = line.find(nameToken);
      int iBegQuote = line.find("\"", iName + 1);
      int iEndQuote = line.find("\"", iBegQuote + 1);
      string name = line.substr(iBegQuote + 1, iEndQuote - iBegQuote - 1);

      // Add file also to be read.
      if (word1 == "<a") {
        files.push_back(pathName + name);
        continue;
      }

      // Find default value token, without/with blank before =.
      if (line.find("default=") == string::npos) {
        cout << " Error: failed to find default value token in line " << line;
        continue;
      }        
      int iDef = line.find("default=");
      iBegQuote = line.find("\"", iDef + 1);
      iEndQuote = line.find("\"", iBegQuote + 1);
      string stringDef = line.substr(iBegQuote + 1, iEndQuote - iBegQuote - 1);

      // Find min and max value tokens, if any.
      bool hasMin = (line.find("min=") != string::npos) ? true : false;
      bool hasMax = (line.find("max=") != string::npos) ? true : false;
      int iMin, iMax;
      string stringMin, stringMax;
      if (hasMin) {
        iMin = line.find("min=");
        iBegQuote = line.find("\"", iMin + 1);
        iEndQuote = line.find("\"", iBegQuote + 1);
        stringMin = line.substr(iBegQuote + 1, iEndQuote - iBegQuote - 1);
      }
      if (hasMax) {
        iMax = line.find("max=");
        iBegQuote = line.find("\"", iMax + 1);
        iEndQuote = line.find("\"", iBegQuote + 1);
        stringMax = line.substr(iBegQuote + 1, iEndQuote - iBegQuote - 1);
      }
    
      // Check for occurence of a bool and add to flag map.
      if (word1 == "<flag") {
        string valueLower = tolower(stringDef);
        bool value = (valueLower == "true" || valueLower == "on"
          || valueLower == "yes" || valueLower == "ok" 
          || valueLower == "1" ) ? true : false ;
        addFlag( name, value);
    
      // Check for occurence of an int and add to mode map.
      } else if (word1 == "<mode") {
        istringstream modeDef(stringDef);
        int value;
        modeDef >> value; 
        int minVal = 0;
        if (hasMin) {
          istringstream modeMin(stringMin);
          modeMin >> minVal;
	}
        int maxVal = 0;
        if (hasMax) {
          istringstream modeMax(stringMax);
          modeMax >> maxVal;
	}
        addMode( name, value, hasMin, hasMax, minVal, maxVal);	  
    
      // Check for occurence of a double and add to parameters map.
      } else if (word1 == "<parameter") {
        istringstream parameterDef(stringDef);
        double value;
        parameterDef >> value; 
        double minVal = 0.;
        if (hasMin) {
          istringstream parameterMin(stringMin);
          parameterMin >> minVal;
	}
        double maxVal = 0.;
        if (hasMax) {
          istringstream parameterMax(stringMax);
          parameterMax >> maxVal;
	}
        addParameter( name, value, hasMin, hasMax, minVal, maxVal);
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

void Settings::list(bool listAll,  bool listString, string match,
  ostream& os) {

  // Table header; output for bool as off/on.
  if (listAll) 
    os << "\n *-------  PYTHIA Flag + Mode + Parameter Settings (all)  -"
       << "---------------------------------------* \n";
  else if (!listString) 
    os << "\n *-------  PYTHIA Flag + Mode + Parameter Settings (changes" 
       << " only)  -------------------------------* \n" ;
  else
    os << "\n *-------  PYTHIA Flag + Mode + Parameter Settings (with re" 
       << "quested string)  ----------------------* \n" ;
  os << " |                                                           "
     << "                                     | \n"
     << " | Name                                     |          Now | "
     << "     Default         Min         Max | \n"
     << " |                                          |              | "
     << "                                     | \n";
 
  // Convert input string to lowercase for match.
  match = tolower(match);

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
      if ( listAll || (!listString && valNow != valDefault)
        || (listString && flagEntry->first.find(match) != string::npos) )
        os << " | " << setw(40) << left 
           << flagEntry->second.name << " | " << setw(12) << right
           << state[valNow] << " | " << setw(12) << state[valDefault] 
           << "                         | \n";
      ++flagEntry;
     
    // Else check if mode is next, and if so print it.
    } else if ( modeEntry != modes.end() && ( paramEntry  
      == parameters.end() || modeEntry->first < paramEntry->first ) ) {
      int valNow = modeEntry->second.valNow;
      int valDefault = modeEntry->second.valDefault;
      if ( listAll || (!listString && valNow != valDefault)
        || (listString && modeEntry->first.find(match) != string::npos) ) {
        os << " | " << setw(40) << left 
           << modeEntry->second.name << " | " << setw(12) << right 
           << valNow << " | " << setw(12) << valDefault; 
        if (modeEntry->second.hasMin) 
          os << setw(12) << modeEntry->second.valMin; 
        else os << "            ";
        if (modeEntry->second.hasMax) 
          os << setw(12) << modeEntry->second.valMax; 
        else os << "            ";
        os << " | \n";
      }
      ++modeEntry;
      
    // Else print param; fixed or scientific depending on value.
    } else {
      double valNow = paramEntry->second.valNow;
      double valDefault = paramEntry->second.valDefault;      
      if ( listAll || (!listString && valNow != valDefault ) 
        || (listString && paramEntry->first.find(match) != string::npos) ) {
        os << " | " << setw(40) << left 
           << paramEntry->second.name << right << " | ";
	for (int i = 0; i < 4; ++i) { 
          if (i == 1) valNow = valDefault;  
          if (i == 2) valNow = paramEntry->second.valMin;  
          if (i == 3) valNow = paramEntry->second.valMax;  
          if ( (i == 2 && !paramEntry->second.hasMin)
	    || (i == 3 && !paramEntry->second.hasMax) )
            os << "            ";
          else if ( valNow == 0. ) 
            os << fixed << setprecision(1) << setw(12) << valNow; 
          else if ( abs(valNow) < 0.001 ) 
            os << scientific << setprecision(4) << setw(12) << valNow;  
          else if ( abs(valNow) < 0.1 )
            os << fixed << setprecision(7) << setw(12) << valNow; 
          else if ( abs(valNow) < 1000. )
            os << fixed << setprecision(5) << setw(12) << valNow; 
          else if ( abs(valNow) < 1000000. )
            os << fixed << setprecision(3) << setw(12) << valNow; 
          else 
            os << scientific << setprecision(4) << setw(12) << valNow; 
          if (i == 0) os << " | ";
	}  
        os << " | \n";
      }
      ++paramEntry;
    }
  } ;

  // End of loop over database contents.
  os << " |                                                           "
     << "                                     | \n"
     << " *-------  End PYTHIA Flag + Mode + Parameter Settings  -----"
     << "-------------------------------------* " << endl;

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

} // end namespace Pythia8
