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

  // Find channel in table, assuming normalization to unity.
  double rand = Rndm::flat();
  int i = -1;
  do rand -= channel[++i].branchingRatio(); 
  while (rand > 0.);
  return channel[i];

}

//*********

// Rescale all branching ratios to assure normalization to unity.

void DecayTable::rescaleBR() {

  // Sum up branching ratios. Find rescaling factor. Rescale.
  double sumBR = 0.;
  for ( int i = 0; i < size(); ++ i) sumBR += channel[i].branchingRatio(); 
  double rescaleFactor = 1. / sumBR;
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

  // Special treatment for lines beginning with <particle> or <channel>.
  int firstChar = line.find_first_not_of(" ");
  if (tolower( line.substr(firstChar, 10) ) == "<particle>")
    readParticle(line);
  if (tolower( line.substr(firstChar, 9) ) == "<channel>")
    readChannel(line);

  // If first character is not a letter/digit, then taken to be a comment.
  if (!isalnum(line[firstChar])) return true; 

  // Replace colons and equal signs by blanks to make parsing simpler.
  for ( int j = 0; j < int(line.size()); ++ j) 
     if (line[j] == ':' || line[j] == '=') line[j] = ' ';

  // Peel off if first word is particle.
  if (tolower( line.substr(firstChar, 8) ) == "particle")
    line.erase(firstChar,8);

  // Get particle id, variable name and value string.
  int id;
  string property, valueString;
  istringstream getWord(line);
  getWord >> id >> property >> valueString;
  property = tolower(property);
  
  // Check that valid particle.
  if (!isParticle(id) || id <= 0) {
    if (warn) cout << "\n Warning: input particle not found in Particle"
      << " Data Table; skip:\n   " << lineIn << "\n";
    return false;
  }

  // Find quantity and value, case by case.
  if (property == "name") {
    pdt[id].setName(valueString);
    return true; 
  }
  if (property == "antiname") {
    pdt[id].setAntiName(valueString);
    return true; 
  }
  if (property == "names") {
    string secondValue; 
    getWord >> secondValue;
    pdt[id].setNames(valueString, secondValue); 
    return true; 
  }
  if (property.substr(0,6) == "charge") {
    istringstream getInt(valueString); 
    int charge3;
    getInt >> charge3;
    pdt[id].setCharge3(charge3); 
    return true; 
  }
  if (property.substr(0,7) == "coltype") {
    istringstream getInt(valueString); 
    int colType;
    getInt >> colType;
    pdt[id].setColType(colType); 
    return true; 
  }
  if (property == "m0") {
    istringstream getDouble(valueString); 
    double m0;
    getDouble >> m0;
    pdt[id].setM0(m0); 
    return true; 
  }
  if (property == "width") {
    istringstream getDouble(valueString); 
    double width;
    getDouble >> width;
    pdt[id].setWidth(width); 
    return true; 
  }
  if (property == "range") {
    istringstream getDouble(valueString); 
    double range;
    getDouble >> range;
    pdt[id].setRange(range); 
    return true; 
  }
  if (property == "tau0") {
    istringstream getDouble(valueString); 
    double tau0;
    getDouble >> tau0;
    pdt[id].setTau0(tau0); 
    return true; 
  }
  if (property == "maydecay") {
    bool mayDecay = (valueString == "true" || valueString == "on"
      || valueString == "yes" || valueString == "ok" 
      || valueString == "1" ) ? true : false ;
    pdt[id].setMayDecay(mayDecay);
    return true; 
  }

  // If is a number then refers to a decay channel.
  if (isdigit(property[0])) {
    int idDummy, channel;
    istringstream getChannel(line);
    getChannel >> idDummy >> channel >> property >> valueString;
    property = tolower(property);
  
    // Find decay channel quantity and value, case by case.
    if (property == "branchingratio") {
      istringstream getDouble(valueString); 
      double bRat;
      getDouble >> bRat;
      pdt[id].decay[channel].branchingRatio(bRat); 
      return true; 
    }
    if (property == "rescalebr") {
      istringstream getDouble(valueString); 
      double factor;
      getDouble >> factor;
      pdt[id].decay[channel].rescaleBR(factor); 
      return true; 
    }
    if (property == "modeme") {
      istringstream getInt(valueString); 
      int mode;
      getInt >> mode;
      pdt[id].decay[channel].modeME(mode); 
      return true; 
    }

    // Scan for products until end of line or encounter non-digit.  
    if (property == "products") {
      istringstream getInt(valueString); 
      int idProd;
      getInt >> idProd;
      pdt[id].decay[channel].product(0, idProd); 
      int nProd = 1;
      for (int iProd = 1; iProd < 8; ++iProd) {
        getChannel >> valueString;
        if (!getChannel || !isdigit(valueString[0])) break;
        istringstream getInt(valueString); 
        getInt >> idProd;
        pdt[id].decay[channel].product(iProd, idProd);
        ++nProd;
      }   
      for (int iProd = nProd; iProd < 8; ++iProd) 
        pdt[id].decay[channel].product(iProd, 0);
      pdt[id].decay[channel].multiplicity(nProd);
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
  string dummy;
  int id;
  string name, antiName;
  int charge3, col;
  double m0, width, range, tau0;
  bool mayDecay;

  // Read in data from stream.
  istringstream particleDataLine(line);
  particleDataLine >> dummy >> id >> name >> antiName >> charge3 
    >> col >> m0 >> width >> range >> tau0 >> mayDecay;   

  // Error printout if something went wrong.
  if (!particleDataLine) {
    cout << "Error: incomplete particle " << id << "\n";
    return false;
  }

  // Erase if particle already exists.
  if (isParticle(id)) pdt.erase(id);

  // Store new particle. Save pointer, to be used for decay channels.
  addParticle( id, name, antiName, charge3, col, m0, width, range, 
    tau0, mayDecay);
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
 
// Print out complete table of database in numerical order.

void ParticleDataTable::list(ostream& os) {

  // Table header; output for bool as off/on.
  os << "\n --------  Pythia Particle Data Table (complete)  ----------"
     << "--------------------------------------------- \n \n"
     << "      id   name            antiname    3*charge colour      m0 "
     << "       width       range        tau0 decay\n"
     << "             no branchratio mode        products \n" ;

  // Iterate through the particle data table.
  for (map<int, ParticleDataEntry>::const_iterator pdtEntry 
    = pdt.begin(); pdtEntry != pdt.end(); ++pdtEntry) {
    const ParticleDataEntry* particlePtr = &pdtEntry->second;
    string antiName = (particlePtr->name(-1) == "void") ?
      "                " : particlePtr->name(-1) ;
    os << "\n" << setw(8) << particlePtr->id() << "   " << setw(16) 
       << left << particlePtr->name() << setw(16) << antiName 
       << right << setw(4) << particlePtr->charge3() << setw(4) 
       << particlePtr->colType() << fixed << setprecision(4) << setw(12) 
       << particlePtr->m0() << setw(12) << particlePtr->width() 
       << setw(12) << particlePtr->range() << scientific << setw(12) 
       << particlePtr->tau0() << setw(6) << particlePtr->mayDecay() 
       << "\n";

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
     << "--------------------------------------------- \n" << endl;

}

//*********
 
// Print out one single particle in the database.

void ParticleDataTable::list(int idList, ostream& os) {

  // Fill one single particle in vector and send it on to other method.
  vector<int> idListTemp;
  idListTemp.push_back(idList);
  list( idListTemp, os); 
 
}

//*********
 
// Print out partial table of database in input order.

void ParticleDataTable::list(vector<int> idList, ostream& os) {

  // Table header; output for bool as off/on.
  os << "\n --------  Pythia Particle Data Table (partial)  ----------"
     << "--------------------------------------------- \n \n"
     << "      id   name            antiname    3*charge colour      m0 "
     << "       width       range        tau0 decay\n"
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
       << "\n";

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
     << "--------------------------------------------- \n" << endl;

}

//**************************************************************************

} // end namespace Pythia8
