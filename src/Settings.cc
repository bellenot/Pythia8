// Function definitions (not found in the header) for the Settings class.
// Copyright C 2006 Torbjorn Sjostrand

#include "Settings.h"

namespace Pythia8 {

//**************************************************************************

// Settings class.
// This class contains flags, modes, parms and words used in generation.

//*********

// Definitions of static variables. 
map<string, Flag> Settings::flags;
map<string, Mode> Settings::modes;
map<string, Parm> Settings::parms;
map<string, Word> Settings::words;
bool Settings::isInit = false;

//*********

// Read in database from specific file.

bool Settings::init(string startFile, ostream& os) {

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
    if (!is) {os << "\n Error: settings file " << files[i] 
       << " not found \n"; return false;}

    // Read in one line at a time.
    string line;
    while ( getline(is, line) ) {

      // Get first word of a line, to interpret it as tag.
      istringstream getfirst(line);
      string tag;
      getfirst >> tag;

      // Skip ahead if not interesting. Only look for new files in startfile.
      if (tag != "<flag" && tag != "<mode" && tag != "<parm"
         && tag != "<word" && (tag != "<a" || i > 0) ) continue;

      // Read and append continuation line(s) if line does not contain >.
      while (line.find(">") == string::npos) {   
        string addLine;
        getline(is, addLine);
        line += " " + addLine;
      } 
      
      // Remove extra blanks before an = sign.
      while (line.find(" =") != string::npos) line.erase( line.find(" ="), 1);

      // Find name attribute.
      string nameAttribute = (tag != "<a") ? "name=" : "href=";
      string name = attributeValue( line, nameAttribute);
      if (name == "") {
        os << " Error: failed to find name attribute in line " << line;
        continue;
      }        

      // Add file also to be read.
      if (tag == "<a") {
        files.push_back(pathName + name);
        continue;
      }

      // Check that default value attribute present, and whether max and min.
      if (line.find("default=") == string::npos) {
        os << " Error: failed to find default value token in line " << line;
        continue;
      }        
      bool hasMin = (line.find("min=") != string::npos) ? true : false;
      bool hasMax = (line.find("max=") != string::npos) ? true : false;
    
      // Check for occurence of a bool and add to flag map.
      if (tag == "<flag") {
        bool value = boolAttributeValue( line, "default=");
        addFlag( name, value);
    
      // Check for occurence of an int and add to mode map.
      } else if (tag == "<mode") {
        int value = intAttributeValue( line, "default="); 
        int minVal = intAttributeValue( line, "min=");
        int maxVal = intAttributeValue( line, "max=");
        addMode( name, value, hasMin, hasMax, minVal, maxVal);	  
    
      // Check for occurence of a double and add to parm map.
      } else if (tag == "<parm") {
        double value = doubleAttributeValue( line, "default="); 
        double minVal = doubleAttributeValue( line, "min=");
        double maxVal = doubleAttributeValue( line, "max=");
        addParm( name, value, hasMin, hasMax, minVal, maxVal);
    
      // Check for occurence of a string and add to word map.
      } else if (tag == "<word") {
        string value = attributeValue( line, "default="); 
        addWord( name, value);
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
  parms.clear();
  words.clear();

  // Then let normal init do the rest.
  isInit = false;
  return init(startFile);

} 

//*********

// Read in updates from a character string, like a line of a file. 
// Is used by readString (and readFile) in Pythia.

bool Settings::readString(string line, bool warn, ostream& os) {

  // If first character is not a letter, then taken to be a comment line.
  int firstChar = line.find_first_not_of(" ");
  if (!isalpha(line[firstChar])) return true; 

  // Replace an equal sign by a blank to make parsing simpler.
  while (line.find("=") != string::npos) {
    int firstEqual = line.find_first_of("=");
    line.replace(firstEqual, 1, " ");   
  }

  // Get first word of a line.
  istringstream splitLine(line);
  string name;
  splitLine >> name;

  // Replace two colons by one (:: -> :) to allow for such mistakes.
  while (name.find("::") != string::npos) {
    int firstColonColon = name.find_first_of("::");
    name.replace(firstColonColon, 2, ":");   
  }
     
  // Check whether this is in the database. Done if not.
  int inDataBase = 0;
  if (isFlag(name)) inDataBase = 1;   
  else if (isMode(name)) inDataBase = 2;   
  else if (isParm(name)) inDataBase = 3; 
  else if (isWord(name)) inDataBase = 4; 
  if (inDataBase == 0) {
    if (warn) os << "\n Warning: input string not found in settings"
      << " databases; skip:\n   " << line << "\n";
    return false;  
  }  

  // Find value.
  string valueString;
  splitLine >> valueString;

  // Update flag map; allow many ways to say yes.
  if (splitLine && inDataBase == 1) {
    bool value = boolString(valueString);
    flag(name, value);

  // Update mode map.
  } else if (splitLine && inDataBase == 2) {
    istringstream modeData(valueString);
    int value;
    modeData >> value;
    if (modeData) mode(name, value);
        
  // Update parm map.
  } else if (splitLine && inDataBase == 3) {
    istringstream parmData(valueString);
    double value;
    parmData >> value;
    if (parmData) parm(name, value);
        
  // Update word map.
  } else if (splitLine && inDataBase == 4) {
    word(name, valueString);
  }

  // Done.
  return true;
}

//*********
 
// Write updates or everything to user-defined file.

bool Settings::writeFile(string toFile, bool writeAll) {

  // Open file for writing.
  const char* cstring = toFile.c_str();
  ofstream os(cstring);  
  if (!os) {
    ErrorMessages::message("Error in Settings::writeFile:"
      " could not open file", toFile);
    return false;
  }
  return writeFile( os, writeAll);

}

//*********
 
// Write updates or everything to user-defined file.

bool Settings::writeFile(ostream& os, bool writeAll) {

  // Write simple header as comment.
  if (writeAll) os << "! List of all current Pythia ";
  else          os << "! List of all modified Pythia ";
  os << fixed << setprecision(3) << parm("Pythia:versionNumber")
     << " settings.\n";

  // Iterators for the flag, mode and parm tables.
  map<string, Flag>::iterator flagEntry = flags.begin();
  map<string, Mode>::iterator modeEntry = modes.begin();
  map<string, Parm>::iterator parmEntry = parms.begin();
  map<string, Word>::iterator wordEntry = words.begin();

  // Loop while there is something left to do.
  while (flagEntry != flags.end() || modeEntry != modes.end() 
    || parmEntry != parms.end() || wordEntry != words.end()) {

    // Check if a flag is next in lexigraphical order; if so print it.
    if ( flagEntry != flags.end() 
      && ( modeEntry == modes.end() || flagEntry->first < modeEntry->first ) 
      && ( parmEntry == parms.end() || flagEntry->first < parmEntry->first )
      && ( wordEntry == words.end() || flagEntry->first < wordEntry->first ) 
      ) {
      string state[2] = {"off", "on"};
      bool valNow = flagEntry->second.valNow;
      bool valDefault = flagEntry->second.valDefault;
      if ( writeAll || valNow != valDefault )
        os << flagEntry->second.name << " = " << state[valNow] << "\n";
      ++flagEntry;
     
    // Else check if mode is next, and if so print it.
    } else if ( modeEntry != modes.end() 
      && ( parmEntry == parms.end() || modeEntry->first < parmEntry->first ) 
      && ( wordEntry == words.end() || modeEntry->first < wordEntry->first ) 
      ) {
      int valNow = modeEntry->second.valNow;
      int valDefault = modeEntry->second.valDefault;
      if ( writeAll || valNow != valDefault ) 
        os << modeEntry->second.name << " = " << valNow << "\n"; 
      ++modeEntry;
      
    // Else check if parm is next, and if so print it; 
    // fixed or scientific depending on value.
    } else if ( parmEntry != parms.end()
      && ( wordEntry == words.end() || parmEntry->first < wordEntry->first ) 
      ) {
      double valNow = parmEntry->second.valNow;
      double valDefault = parmEntry->second.valDefault;      
      if ( writeAll || valNow != valDefault ) {
        os  << parmEntry->second.name << " = ";
        if ( valNow == 0. ) os << fixed << setprecision(1); 
        else if ( abs(valNow) < 0.001 ) os << scientific << setprecision(4); 
        else if ( abs(valNow) < 0.1 ) os << fixed << setprecision(7); 
        else if ( abs(valNow) < 1000. ) os << fixed << setprecision(5); 
        else if ( abs(valNow) < 1000000. ) os << fixed << setprecision(3); 
        else os << scientific << setprecision(4); 
        os << valNow << "\n";
      }
      ++parmEntry;

    // Else print word. 
    } else {
      string valNow = wordEntry->second.valNow;
      string valDefault = wordEntry->second.valDefault; 
      if ( writeAll || valNow != valDefault ) 
        os << wordEntry->second.name << " = " << valNow << "\n";
      ++wordEntry;
    }
  } ;

  // Done.
  return true;
}

//*********
 
// Print out table of database in lexigraphical order.

void Settings::list(bool listAll,  bool listString, string match,
  ostream& os) {

  // Table header; output for bool as off/on.
  if (listAll) 
    os << "\n *-------  PYTHIA Flag + Mode + Parm + Word Settings (all) "
       << " --------------------------------------------------* \n";
  else if (!listString) 
    os << "\n *-------  PYTHIA Flag + Mode + Parm + Word Settings (chang" 
       << "es only)  -----------------------------------------* \n" ;
  else
    os << "\n *-------  PYTHIA Flag + Mode + Parm + Word Settings (with " 
       << "requested string)  --------------------------------* \n" ;
  os << " |                                                           "
     << "                                                 | \n"
     << " | Name                                     |                "
     << "      Now |      Default         Min         Max | \n"
     << " |                                          |                "
     << "          |                                      | \n";
 
  // Convert input string to lowercase for match.
  match = tolower(match);

  // Iterators for the flag, mode and parm tables.
  map<string, Flag>::iterator flagEntry = flags.begin();
  map<string, Mode>::iterator modeEntry = modes.begin();
  map<string, Parm>::iterator parmEntry = parms.begin();
  map<string, Word>::iterator wordEntry = words.begin();

  // Loop while there is something left to do.
  while (flagEntry != flags.end() || modeEntry != modes.end() 
    || parmEntry != parms.end() || wordEntry != words.end()) {

    // Check if a flag is next in lexigraphical order; if so print it.
    if ( flagEntry != flags.end() 
      && ( modeEntry == modes.end() || flagEntry->first < modeEntry->first ) 
      && ( parmEntry == parms.end() || flagEntry->first < parmEntry->first )
      && ( wordEntry == words.end() || flagEntry->first < wordEntry->first ) 
      ) {
      string state[2] = {"off", "on"};
      bool valNow = flagEntry->second.valNow;
      bool valDefault = flagEntry->second.valDefault;
      if ( listAll || (!listString && valNow != valDefault)
        || (listString && flagEntry->first.find(match) != string::npos) )
        os << " | " << setw(40) << left 
           << flagEntry->second.name << " | " << setw(24) << right
           << state[valNow] << " | " << setw(12) << state[valDefault] 
           << "                         | \n";
      ++flagEntry;
     
    // Else check if mode is next, and if so print it.
    } else if ( modeEntry != modes.end() 
      && ( parmEntry == parms.end() || modeEntry->first < parmEntry->first ) 
      && ( wordEntry == words.end() || modeEntry->first < wordEntry->first ) 
      ) {
      int valNow = modeEntry->second.valNow;
      int valDefault = modeEntry->second.valDefault;
      if ( listAll || (!listString && valNow != valDefault)
        || (listString && modeEntry->first.find(match) != string::npos) ) {
        os << " | " << setw(40) << left 
           << modeEntry->second.name << " | " << setw(24) << right 
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
      
    // Else check if parm is next, and if so print it; 
    // fixed or scientific depending on value.
    } else if ( parmEntry != parms.end()
      && ( wordEntry == words.end() || parmEntry->first < wordEntry->first ) 
      ) {
      double valNow = parmEntry->second.valNow;
      double valDefault = parmEntry->second.valDefault;      
      if ( listAll || (!listString && valNow != valDefault ) 
        || (listString && parmEntry->first.find(match) != string::npos) ) {
        os << " | " << setw(40) << left 
           << parmEntry->second.name << right << " |             ";
	for (int i = 0; i < 4; ++i) { 
          if (i == 1) valNow = valDefault;  
          if (i == 2) valNow = parmEntry->second.valMin;  
          if (i == 3) valNow = parmEntry->second.valMax;  
          if ( (i == 2 && !parmEntry->second.hasMin)
	    || (i == 3 && !parmEntry->second.hasMax) )
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
      ++parmEntry;

    // Else print word. 
    } else {
      string valNow = wordEntry->second.valNow;
      string valDefault = wordEntry->second.valDefault; 
      int blankLeft = max(0, 60 - max(24, int(valNow.length()) ) 
        - max(12, int(valDefault.length()) ) );  
      string blankPad( blankLeft, ' '); 
      if ( listAll || (!listString && valNow != valDefault)
        || (listString && wordEntry->first.find(match) != string::npos) )
        os << " | " << setw(40) << left 
           << wordEntry->second.name << " | " << setw(24) << right
           << valNow << " | " << setw(12) << valDefault << blankPad 
           << " | \n";
      ++wordEntry;
    }
  } ;

  // End of loop over database contents.
  os << " |                                                           "
     << "                                                 | \n"
     << " *-------  End PYTHIA Flag + Mode + Parm + Word Settings  ---"
     << "-------------------------------------------------* " << endl;

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

  // Loop through the parms table, resetting all entries.
  for (map<string, Parm>::iterator parmEntry = parms.begin(); 
    parmEntry != parms.end(); ++parmEntry) {
    string name = parmEntry->first;
    resetParm(name);
  }

  // Loop through the words table, resetting all entries.
  for (map<string, Word>::iterator wordEntry = words.begin();
    wordEntry != words.end(); ++wordEntry) {
    string name = wordEntry->first;
    resetWord(name);
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

void Settings::parm(string keyIn, double nowIn) { 
  if (isParm(keyIn)) {
    Parm& parmNow = parms[tolower(keyIn)];
    if (parmNow.hasMin && nowIn < parmNow.valMin) 
      parmNow.valNow = parmNow.valMin; 
    else if (parmNow.hasMax && nowIn > parmNow.valMax) 
      parmNow.valNow = parmNow.valMax;
    else parmNow.valNow = nowIn; 
  } 
}  

void Settings::word(string keyIn, string nowIn) { 
    if (isWord(keyIn)) words[tolower(keyIn)].valNow = nowIn; 
}

//**************************************************************************

} // end namespace Pythia8
