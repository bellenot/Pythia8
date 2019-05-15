// Function definitions (not found in the header) for the
// DecayChannel, ParticleDataEntry and ParticleDataTable classes.
// Copyright C 2006 Torbjorn Sjostrand

#include "ParticleData.h"

namespace Pythia8 {

//**************************************************************************

// DecayChannel class.
// This class holds info on all decay channels of a particle. 

//*********

// Pick a decay channel according to branching ratios.

DecayChannel& DecayTable::pick() {

  // Find sum of branching ratios, in case it is not unity.
  double sumBR = 0.;
  for ( int i = 0; i < size(); ++ i) sumBR += channel[i].branchingRatio(); 

  // Find channel in table, assuming normalization to unity.
  double rand = Rndm::flat() * sumBR;
  int i = -1;
  do rand -= channel[++i].branchingRatio(); 
  while (rand > 0.);
  return channel[i];

}

//*********

// Rescale all branching ratios to assure normalization to unity.

void DecayTable::rescaleBR(double newSumBR) {

  // Sum up branching ratios. Find rescaling factor. Rescale.
  double oldSumBR = 0.;
  for ( int i = 0; i < size(); ++ i) 
    oldSumBR += channel[i].branchingRatio(); 
  double rescaleFactor = newSumBR / oldSumBR;
  for ( int i = 0; i < size(); ++ i) channel[i].rescaleBR(rescaleFactor); 

}

//**************************************************************************

// ParticleDataEntry class.
// This class holds info on a single particle species.

//*********
 
// Definitions of static variables.
// (Values will be overwritten in initStatic call, so are purely dummy.)

int ParticleDataEntry::modeBreitWigner = 1;

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Narrow states are assigned nominal mass.
const double ParticleDataEntry::NARROWMASS = 1e-6;

//*********

// Initialize static data members.

void ParticleDataEntry::initStatic() {

  // Mass generation: fixed mass or linear/quadratic Breit-Wigner.
  modeBreitWigner = Settings::mode("ParticleData:modeBreitWigner");

}

//*********

// Function to give mass of a particle, either at the nominal value
// or picked according to a (linear or quadratic) Breit-Wigner. 

double ParticleDataEntry::mass() {

  // Nominal value.
  if (modeBreitWigner == 0 || widthSave < NARROWMASS 
    || rangeSave < NARROWMASS) return m0Save;

  // Mass according to a Breit-Wigner linear in m.
  if (modeBreitWigner == 1) {  
    double truncatedBW = tan( (2. * Rndm::flat() - 1.) 
      * atan(2. * rangeSave / widthSave) );
    return m0Save + 0.5 * widthSave * truncatedBW;
  }

  // Mass according to a Breit-Wigner quadratic in m.
  double atanM2Low = atan( (pow2(max(0., m0Save - rangeSave)) - pow2(m0Save))
    / (m0Save * widthSave) );
  double atanM2High = atan( (pow2(m0Save + rangeSave) - pow2(m0Save))
    / (m0Save * widthSave) );
  double m2 = pow2(m0Save) + m0Save * widthSave 
    * tan(atanM2Low + (atanM2High - atanM2Low) * Rndm::flat());
  return sqrt(max(0., m2));

}

//*********

// Constituent masses for (d, u, s, c, b) quarks and diquarks.
// Hardcoded here so that they are not overwritten by mistake,
// and separated from the "normal" masses. 
  
void ParticleDataEntry::constituentMassCalc() {

  // Equate with the normal masses as default guess.
  constituentMassSave = m0Save;

  // The constituent masses.
  const double constituentMassTable[6] 
    = {0., 0.325, 0.325, 0.50, 1.60, 5.00};

  // Quark masses trivial.
  if (idAbs < 6) constituentMassSave = constituentMassTable[idAbs];
 
  // Diquarks as simple sum of constituent quarks.  
  if (idAbs > 1000 && idAbs < 10000 && (idAbs/10)%10 == 0) {
    int id1 = idAbs/1000;
    int id2 = (idAbs/100)%10;
    if (id1 <6 && id2 < 6) constituentMassSave = constituentMassTable[id1] 
      + constituentMassTable[id2];
  }

}

//*********

// Function to tell which particles are vsible in a detector.
// May not cover all particles ??.

bool ParticleDataEntry::isVisible() const {

  const int SUSY = 1000000;
  if (idAbs == 12 || idAbs == 14 || idAbs == 16 || idAbs == 18) 
    return false;
  if (idAbs == 39 || idAbs == SUSY + 22 || idAbs == SUSY + 39)
    return false;
  return true;

}

//*********

// Function to give spin of particle as 2 * s + 1.
// May not cover all particles ??.

int ParticleDataEntry::spinType() const {
  
  // Preliminaries and default value (= unknown).
  const int SUSY = 1000000;
  int spin = 0;

  // Go through various known cases; could be expanded.
       if (idAbs == 130 || idAbs == 310) spin = 1; 
  else if (idAbs > 100 && idAbs <= SUSY) spin = idAbs%10; 
  else if (idAbs > 0 && idAbs <= 18) spin = 2; 
  else if (idAbs == 21) spin = 3; 
  else if (idAbs >= 22 && idAbs <= 24) spin = 3; 
  else if (idAbs >= 32 && idAbs <= 34) spin = 3; 
  else if (idAbs == 25 || (idAbs >= 35 && idAbs <= 37) ) spin = 1;
  else if (idAbs > SUSY && idAbs <= SUSY+18) spin = 1;  
  else if (idAbs > 2*SUSY && idAbs <= 2*SUSY+18) spin = 1; 
  else if (idAbs == SUSY+21) spin = 2; 
  else if (idAbs >= SUSY+22 && idAbs <= SUSY+37) spin = 2;

  return spin;
}

//**************************************************************************

// ParticleDataTable class.
// This class holds a map of all ParticleDataEntries,
// each entry containing info on a particle species.

//*********

// Definitions of static variables. 

map<int, ParticleDataEntry> ParticleDataTable::pdt;
bool ParticleDataTable::isInit = false;
ParticleDataEntry* ParticleDataTable::particlePtr = 0;

//*********

// Read in database from specific file.

bool ParticleDataTable::init(string startFile) {

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
    if (!is) {cout << "Error: particle data file " << files[i] 
       << "not found \n"; return false;}

    // Read in one line at a time.
    particlePtr = 0;
    string line;
    while ( getline(is, line) ) {

      // Get first word of a line.
      istringstream getfirst(line);
      string word1;
      getfirst >> word1;
    
      // Check for occurence of a particle and add to particle data table.
      if (word1 == "<particle>") readParticle(line);
    
      // Check for occurence of a decay channel and add to decay table. 
      else if (word1 == "<channel>") readChannel(line);

      // Check for occurence of a file also to be read.
      else if (word1 == "<file>") {
        istringstream fileData(line);
        string dummy, file;
        fileData >> dummy >> file;
        if (!fileData) cout << "Error: incomplete file " << file << "\n";
        else files.push_back(file);
      }

    // End of loop over lines in input file and loop over files.
    };
  };

  // All particle data at this stage defines baseline original.
  for (map<int, ParticleDataEntry>::iterator pdtEntry 
    = pdt.begin(); pdtEntry != pdt.end(); ++pdtEntry) {
    ParticleDataEntry* particlePtr = &pdtEntry->second;
    particlePtr->setHasChangedAll(false);
  }

  // Done.
  isInit = true;
  return true;
}

//*********

// Overwrite existing database by reading from specific file.

bool ParticleDataTable::reInit(string startFile) {

  // Reset map to empty and then let normal init do the rest.
  pdt.clear();
  isInit = false;
  return init(startFile);

} 

//*********

// Read in updates from a character string, like a line of a file. 
// Is used by readFile, and by readString (and readFile) in Pythia.

bool ParticleDataTable::readString(string lineIn, bool warn) {

  // Take copy that will be modified.
  string line = lineIn;

  // If first character is not a letter/digit, then taken to be a comment.
  int firstChar = line.find_first_not_of(" ");
  if (!isalnum(line[firstChar])) return true; 

  // Replace colons and equal signs by blanks to make parsing simpler.
  for ( int j = 0; j < int(line.size()); ++ j) 
     if (line[j] == ':' || line[j] == '=') line[j] = ' ';

  // Get particle id and property name.
  int id;
  string property;
  istringstream getWord(line);
  getWord >> id >> property;
  property = tolower(property);
  
  // Check that valid particle.
  if ( (!isParticle(id) && property  != "all" && property  != "new") 
  || id <= 0) {
    if (warn) cout << "\n Warning: input particle not found in Particle"
      << " Data Table; skip:\n   " << lineIn << "\n";
    return false;
  }

  // Identify particle property and read + set its value, case by case.
  if (property == "name") {
    string name;
    getWord >> name;    
    pdt[id].setName(name);
    return true; 
  }
  if (property == "antiname") {
    string antiName;
    getWord >> antiName;    
    pdt[id].setAntiName(antiName);
    return true; 
  }
  if (property == "names") {
    string name, antiName; 
    getWord >> name >> antiName;
    pdt[id].setNames(name, antiName); 
    return true; 
  }
  if (property.substr(0,6) == "charge") {
    int charge3;
    getWord >> charge3;
    pdt[id].setCharge3(charge3); 
    return true; 
  }
  if (property.substr(0,7) == "coltype") {
    int colType;
    getWord >> colType;
    pdt[id].setColType(colType); 
    return true; 
  }
  if (property == "m0") {
    double m0;
    getWord >> m0;
    pdt[id].setM0(m0); 
    return true; 
  }
  if (property == "width") {
    double width;
    getWord >> width;
    pdt[id].setWidth(width); 
    return true; 
  }
  if (property == "range") {
    double range;
    getWord >> range;
    pdt[id].setRange(range); 
    return true; 
  }
  if (property == "tau0") {
    double tau0;
    getWord >> tau0;
    pdt[id].setTau0(tau0); 
    return true; 
  }
  if (property == "maydecay") {
    string may;
    getWord >> may;
    may = tolower(may);
    bool mayDecay = (may == "true" || may == "on" || may == "yes" 
      || may == "ok" || may == "1" ) ? true : false ;
    pdt[id].setMayDecay(mayDecay);
    return true; 
  }  
  if (property == "isresonance") {
    string isres;
    getWord >> isres;
    isres = tolower(isres);
    bool isResonance = (isres == "true" || isres == "on" || isres == "yes" 
      || isres == "ok" || isres == "1" ) ? true : false ;
    pdt[id].setIsResonance(isResonance);
    return true; 
  }

    // Rescale all branching ratios by common factor.
  if (property == "rescalebr") {
    double factor;
    getWord >> factor;
    pdt[id].rescaleBR(factor); 
    return true; 
  }        
   
  // Addition or complete replacement of a particle.
  if (property == "all" || property == "new") {

    // Properties to be read. 
    string name, antiName, may, isres;
    int charge3, col;
    double m0, width, range, tau0;

    // Read in data from stream.
    getWord >> name >> antiName >> charge3 >> col >> m0 >> width 
            >> range >> tau0 >> may >> isres;   
    may = tolower(may);
    bool mayDecay = (may == "true" || may == "on" || may == "yes" 
      || may == "ok" || may == "1" ) ? true : false ;
    isres = tolower(isres);
    bool isResonance = (isres == "true" || isres == "on" || isres == "yes" 
      || isres == "ok" || isres == "1" ) ? true : false ;

    // Error printout if something went wrong.
    if (!getWord) {
      cout << "Error: incomplete particle " << id << "\n";
      return false;
    }

    // To keep existing decay channels, only overwrite particle data.
    if (property == "all" && isParticle(id)) {
      allParticle( id, name, antiName, charge3, col, m0, width, 
        range, tau0, mayDecay, isResonance);   

    // Else start over completely from scratch.
    } else {
      if (isParticle(id)) pdt.erase(id);
      addParticle( id, name, antiName, charge3, col, m0, width, 
        range, tau0, mayDecay, isResonance);
    }
    return true;
  }

  // Add or change a decay channel: get channel number and new property.
  if (property.substr(0,3) == "add" || isdigit(property[0])) {
    int channel;
    if (property.substr(0,3) == "add") {
      pdt[id].decay.addChannel(); 
      channel = pdt[id].decay.size() - 1; 
      property = "all"; 
    } else{ 
      istringstream getChannel(property);
      getChannel >> channel;
      getWord >> property;
      property = tolower(property);
    }

    // Check that channel exists.
    if (channel < 0 || channel >= pdt[id].decay.size()) return false;   
  
    // Find decay channel property and value, case by case.
    // At same time also do case where all should be replaced.
    if (property.substr(0,2) == "br" || property == "all") {
      double bRat;
      getWord >> bRat;
      pdt[id].decay[channel].branchingRatio(bRat); 
      if (property.substr(0,2) == "br") return true; 
    }
    if (property == "modeme" || property == "all") {
      int mode;
      getWord >> mode;
      pdt[id].decay[channel].modeME(mode); 
      if (property == "modeme") return true; 
    }    

    // Scan for products until end of line.  
    if (property == "products" || property == "all") {
      int nProd = 0;
      for (int iProd = 0; iProd < 8; ++iProd) {
        int idProd;
        getWord >> idProd;
        if (!getWord) break;
        pdt[id].decay[channel].product(iProd, idProd);
        ++nProd;
      }   
      for (int iProd = nProd; iProd < 8; ++iProd) 
        pdt[id].decay[channel].product(iProd, 0);
      pdt[id].decay[channel].multiplicity(nProd);
      return true; 
    }

    // Rescale an existing branching ratio.
    if (property == "rescalebr") {
      double factor;
      getWord >> factor;
      pdt[id].decay[channel].rescaleBR(factor); 
      return true; 
    }        
  }

  // Return false if failed to recognize property.
  if (warn) cout << "\n Warning: input property not found in Particle"
    << " Data Table; skip:\n   " << lineIn << "\n";
  return false;

}

//*********

// Read in updates from user-defined file.

bool ParticleDataTable::readFile(string updateFile, bool warn) {

  // Open file with updates.
  const char* cstring = updateFile.c_str();
  ifstream is(cstring);  
  if (!is) {cout << "Error: user update file " << updateFile 
     << "not found \n"; return false;}

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

// Read in a line with <particle> format information.
  
bool ParticleDataTable::readParticle(string line) {

  // Properties to be read. 
  string dummy, may, isres;
  int id;
  string name, antiName;
  int charge3, col;
  double m0, width, range, tau0;

  // Read in data from stream.
  istringstream particleDataLine(line);
  particleDataLine >> dummy >> id >> name >> antiName >> charge3 
    >> col >> m0 >> width >> range >> tau0 >> may >> isres;   
  may = tolower(may);
  bool mayDecay = (may == "true" || may == "on" || may == "yes" 
    || may == "ok" || may == "1" ) ? true : false ;
  isres = tolower(isres);
  bool isResonance = (isres == "true" || isres == "on" || isres == "yes" 
    || isres == "ok" || isres == "1" ) ? true : false ;

  // Error printout if something went wrong.
  if (!particleDataLine) {
    cout << "Error: incomplete particle " << id << "\n";
    return false;
  }

  // Erase if particle already exists.
  if (isParticle(id)) pdt.erase(id);

  // Store new particle. Save pointer, to be used for decay channels.
  addParticle( id, name, antiName, charge3, col, m0, width, range, 
    tau0, mayDecay, isResonance);
  particlePtr = particleDataPtr(id);

  // Done.
  return true;

}

//*********

// Read in a line with <channel> format information.
  
bool ParticleDataTable::readChannel(string line) {

  // Properties to be read. 
  string dummy;
  double brat;
  int mode = 0;
  int prod0 = 0;
  int prod1 = 0;
  int prod2 = 0;
  int prod3 = 0;
  int prod4 = 0;
  int prod5 = 0;
  int prod6 = 0;
  int prod7 = 0;
  
  // Read in data from stream. Need at least one decay product.
  istringstream channelDataLine(line);
  channelDataLine >> dummy >> brat >> mode >> prod0;
  if (!channelDataLine) {
    cout << "Error: incomplete decay channel\n";
    return false;
  }
  channelDataLine >> prod1 >> prod2 >> prod3 >> prod4 >> prod5 
    >> prod6  >> prod7;   

  // Store new channel.
  if (particlePtr == 0) {
    cout << "Error: no decay particle\n";
    return false;
  }
  particlePtr->decay.addChannel(brat, mode, prod0, prod1, prod2, 
    prod3, prod4, prod5, prod6, prod7);  

  // Done.
  return true;

}

//*********
 
// Print out complete or changed table of database in numerical order.

void ParticleDataTable::list(bool changedOnly, ostream& os) {

  // Table header; output for bool as off/on.
  if (!changedOnly) {
    os << "\n --------  Pythia Particle Data Table (complete)  --------"
       << "----------------------------------------------------- \n \n";
  } else { 
    os << "\n --------  Pythia Particle Data Table (changed only)  ----"
       << "----------------------------------------------------- \n \n";
       }
  os << "      id   name            antiname    3*charge colour      m0 "
     << "       width       range        tau0 decay reson\n"
     << "             no branchratio mode        products \n" ;

  // Iterate through the particle data table. Option to skip unchanged.
  for (map<int, ParticleDataEntry>::const_iterator pdtEntry 
    = pdt.begin(); pdtEntry != pdt.end(); ++pdtEntry) {
    const ParticleDataEntry* particlePtr = &pdtEntry->second;
    if (changedOnly && !particlePtr->hasChangedAny()) continue;

    // Print particle properties.
    string antiName = (particlePtr->name(-1) == "void") ?
      "                " : particlePtr->name(-1) ;
    os << "\n" << setw(8) << particlePtr->id() << "   " << setw(16) 
       << left << particlePtr->name() << setw(16) << antiName 
       << right << setw(4) << particlePtr->charge3() << setw(4) 
       << particlePtr->colType() << fixed << setprecision(4) << setw(12) 
       << particlePtr->m0() << setw(12) << particlePtr->width() 
       << setw(12) << particlePtr->range() << scientific << setw(12) 
       << particlePtr->tau0() << setw(6) << particlePtr->mayDecay() 
       << setw(6) << particlePtr->isResonance() << "\n";

    // Loop through the decay channel table for each particle.
    if (particlePtr->decay.size() > 0) {
      for (int i = 0; i < int(particlePtr->decay.size()); ++i) {
        const DecayChannel& channel = particlePtr->decay[i];
        os << "          " << setw(5) << i << fixed << setprecision(7) 
           << setw(12) << channel.branchingRatio() << setw(5) 
           << channel.modeME();
        for (int j = 0; j < channel.multiplicity(); ++j) 
          os << setw(9) << channel.product(j);
        os << "\n";  
      }
    }  

  } ;

  // End of loop over database contents.
  os << "\n --------  End Particle Data Table  ------------------------"
     << "--------------------------------------------------- \n" << endl;

}

//*********
 
// Print out partial table of database in input order.

void ParticleDataTable::list(vector<int> idList, ostream& os) {

  // Table header; output for bool as off/on.
  os << "\n --------  Pythia Particle Data Table (partial)  ---------"
     << "----------------------------------------------------- \n \n"
     << "      id   name            antiname    3*charge colour      m0 "
     << "       width       range        tau0 decay reson\n"
     << "             no branchratio mode        products \n" ;

  // Iterate through the given list of input particles.
  for (int i = 0; i < int(idList.size()); ++i) {
    const ParticleDataEntry* particlePtr = particleDataPtr(idList[i]);
    string antiName = (particlePtr->name(-1) == "void") ?
      "                " : particlePtr->name(-1) ;
    os << "\n" << setw(8) << particlePtr->id() << "   " << setw(16) 
       << left << particlePtr->name() << setw(16) << antiName 
       << right << setw(4) << particlePtr->charge3() << setw(4) 
       << particlePtr->colType() << fixed << setprecision(4) << setw(12) 
       << particlePtr->m0() << setw(12) << particlePtr->width() 
       << setw(12) << particlePtr->range() << scientific << setw(12) 
       << particlePtr->tau0() << setw(6) << particlePtr->mayDecay() 
       << setw(6) << particlePtr->isResonance() << "\n";

    // Loop through the decay channel table for each particle.
    if (particlePtr->decay.size() > 0) {
      for (int i = 0; i < int(particlePtr->decay.size()); ++i) {
        const DecayChannel& channel = particlePtr->decay[i];
        os << "          " << setw(5) << i << fixed << setprecision(7) 
           << setw(12) << channel.branchingRatio() << setw(5) 
           << channel.modeME();
        for (int j = 0; j < channel.multiplicity(); ++j) 
          os << setw(9) << channel.product(j);
        os << "\n";  
      }
    }  

  } ;

  // End of loop over database contents.
  os << "\n --------  End Particle Data Table  ------------------------"
     << "--------------------------------------------------- \n" << endl;

}

//**************************************************************************

} // end namespace Pythia8
