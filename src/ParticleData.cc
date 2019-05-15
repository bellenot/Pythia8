// Function definitions (not found in the header) for the
// DecayChannel, ParticleDataEntry and ParticleDataTable classes.
// Copyright C 2006 Torbjorn Sjostrand

#include "ParticleData.h"

namespace Pythia8 {

//**************************************************************************

// DecayChannel class.
// This class holds info on all decay channels of a particle. 

//*********

// Rescale all branching ratios to assure normalization to unity.

void DecayTable::rescaleBR(double newSumBR) {

  // Sum up branching ratios. Find rescaling factor. Rescale.
  double oldSumBR = 0.;
  for ( int i = 0; i < size(); ++ i) 
    oldSumBR += channel[i].bRatio(); 
  double rescaleFactor = newSumBR / oldSumBR;
  for ( int i = 0; i < size(); ++ i) channel[i].rescaleBR(rescaleFactor); 

}

//*********

// Pick a decay channel according to branching ratios.

DecayChannel& DecayTable::pick() {

  // Find sum of allowed branching ratios, in case it is not unity.
  double sumBR = 0.;
  for ( int i = 0; i < size(); ++ i) if (channel[i].onMode() > 0) 
    sumBR += channel[i].bRatio(); 

  // Find channel in table, assuming normalization to unity.
  double rand = Rndm::flat() * sumBR;
  int i = 0;
  do {
    if (channel[i].onMode() > 0) rand -= channel[i].bRatio(); 
    ++i;
  } while (rand > 0.);
  return channel[i - 1];

}

//*********

// Pick a decay channel according to dynamically calculated branching ratios.

DecayChannel& DecayTable::dynamicPick() {

  // Find sum of branching ratios, in case it is not unity.
  double sumBR = 0.;
  for ( int i = 0; i < size(); ++ i) if (channel[i].onMode() > 0) 
    sumBR += channel[i].dynamicBR(); 

  // If vanishing then default back to normal branching ratios.
  if (sumBR <= 0.) return pick();

  // Find channel in table, assuming normalization to unity.
  double rand = Rndm::flat() * sumBR;
  int i = 0;
  do {
    if (channel[i].onMode() > 0) rand -= channel[i].dynamicBR(); 
    ++i;
  } while (rand > 0.);
  return channel[i - 1];

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

// Particles with a read-in tau0 (in mm/c) below this mayDecay by default.
const double ParticleDataEntry::MAXTAU0FORDECAY = 1000.;

// Particles with a read-in m0 above this isResonance by default.
const double ParticleDataEntry::MINMASSRESONANCE = 20.;

// Narrow states are assigned nominal mass.
const double ParticleDataEntry::NARROWMASS = 1e-6;

//*********

// Initialize static data members.

void ParticleDataEntry::initStatic() {

  // Mass generation: fixed mass or linear/quadratic Breit-Wigner.
  modeBreitWigner = Settings::mode("ParticleData:modeBreitWigner");

}

//*********

// Set initial default values for some quantities. 

void ParticleDataEntry::setDefaults() {

  // A particle may decay if it is shortlived enough.
  mayDecaySave = (tau0Save < MAXTAU0FORDECAY); 

  // A particle is a resonance if it is heavy enough.
  isResonanceSave = (m0Save > MINMASSRESONANCE);

  // A particle is invisible if it has neither strong nor electric charge,
  // and is not made up of constituents that have it. Only relevant for
  // long-lived particles. This list may need to be extended.
  isVisibleSave = true;
  const int invisibleTable[25] = { 12, 14, 16, 18, 23, 25, 32, 33, 35, 
    36, 39, 41, 1000012, 1000014, 1000016, 1000018, 1000022, 1000023, 
    1000025, 1000035, 1000039, 2000012, 2000014, 2000016, 2000018};     
  for (int i = 0; i < 25; ++i) 
    if (idSave == invisibleTable[i]) isVisibleSave = false;   

  // A particle by default has no external decays.
  externalDecaySave = false;

  // Set up constituent masses.
  constituentMassCalc();

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
  if (idSave < 6) constituentMassSave = constituentMassTable[idSave];
 
  // Diquarks as simple sum of constituent quarks.  
  if (idSave > 1000 && idSave < 10000 && (idSave/10)%10 == 0) {
    int id1 = idSave/1000;
    int id2 = (idSave/100)%10;
    if (id1 <6 && id2 < 6) constituentMassSave = constituentMassTable[id1] 
      + constituentMassTable[id2];
  }

}

//*********

// Function to give mass of a particle, either at the nominal value
// or picked according to a (linear or quadratic) Breit-Wigner. 

double ParticleDataEntry::mass() {

  // Nominal value.
  if (modeBreitWigner == 0 || mWidthSave < NARROWMASS 
    || (mMaxSave > mMinSave && mMaxSave - mMinSave < NARROWMASS) ) 
    return m0Save;

  // Mass according to a Breit-Wigner linear in m.
  if (modeBreitWigner == 1) {  
    double atanMLow = atan( 2. * (mMinSave - m0Save) / mWidthSave );
    double atanMHigh = (mMaxSave > mMinSave) 
      ? atan( 2. * (mMaxSave - m0Save) / mWidthSave ) : 0.5 * M_PI;
    double atanMRand = atanMLow + (atanMHigh - atanMLow) * Rndm::flat();
    return m0Save + 0.5 * mWidthSave * tan(atanMRand);
  }

  // Mass according to a Breit-Wigner quadratic in m.
  double atanM2Low = atan( (pow2(mMinSave) - pow2(m0Save))
    / (m0Save * mWidthSave) );
  double atanM2High = (mMaxSave > mMinSave)
    ? atan( (pow2(mMaxSave) - pow2(m0Save)) / (m0Save * mWidthSave) )
    : 0.5 * M_PI;
  double atanM2Rand = atanM2Low + (atanM2High - atanM2Low) * Rndm::flat();   
  double m2 = pow2(m0Save) + m0Save * mWidthSave * tan(atanM2Rand);
  return sqrt(max(0., m2));

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

// Read in database from specific XML file (which may refer to others).

bool ParticleDataTable::readXML(string inFile, bool reset) {

  // Normally reset whole database before beginning.
  if (reset) {pdt.clear(); isInit = false;}
 
  // List of files to be checked.
  vector<string> files;
  files.push_back(inFile);

  // Loop over files. Open them for read.
  for (int i = 0; i < int(files.size()); ++i) {
    const char* cstring = files[i].c_str();
    ifstream is(cstring);  

    // Check that instream is OK.
    if (!is) {
      ErrorMessages::message("Error in ParticleDataTable::readXML:"
        " did not find file", files[i]);
      return false;
    }

    // Read in one line at a time.
    particlePtr = 0;
    string line;
    while ( getline(is, line) ) {

      // Get first word of a line.
      istringstream getfirst(line);
      string word1;
      getfirst >> word1;
    
      // Check for occurence of a particle. Add any continuation lines.
      if (word1 == "<particle") {
        while (line.find(">") == string::npos) {   
          string addLine;
          getline(is, addLine);
          line += addLine;
        } 

        // Read in particle properties.
        int id = intAttributeValue( line, "id");
        string name = attributeValue( line, "name");
        string antiName = attributeValue( line, "antiName");
        if (antiName == "") antiName = "void";
        int spinType = intAttributeValue( line, "spinType");
        int chargeType = intAttributeValue( line, "chargeType");
        int colType = intAttributeValue( line, "colType");
        double m0 = doubleAttributeValue( line, "m0");
        double mWidth = doubleAttributeValue( line, "mWidth");
        double mMin = doubleAttributeValue( line, "mMin");
        double mMax = doubleAttributeValue( line, "mMax");
        double tau0 = doubleAttributeValue( line, "tau0");

        // Erase if particle already exists.
        if (isParticle(id)) pdt.erase(id);

        // Store new particle. Save pointer, to be used for decay channels.
        addParticle( id, name, antiName, spinType, chargeType, colType, 
          m0, mWidth, mMin, mMax, tau0);
        particlePtr = particleDataPtr(id);
 
      // Check for occurence of a decay channel. Add any continuation lines. 
      } else if (word1 == "<channel") {
        while (line.find(">") == string::npos) {   
          string addLine;
          getline(is, addLine);
          line += addLine;
        } 

        // Read in channel properties - products so far only as a string.
        int onMode = intAttributeValue( line, "onMode");
        double bRatio = doubleAttributeValue( line, "bRatio");
        int meMode = intAttributeValue( line, "meMode");
        string products = attributeValue( line, "products");
 
        // Read in decay products from stream. Must have at least one.
        istringstream prodStream(products);
        int prod0 = 0; int prod1 = 0; int prod2 = 0; int prod3 = 0;
        int prod4 = 0; int prod5 = 0; int prod6 = 0; int prod7 = 0;
        prodStream >> prod0 >> prod1 >> prod2 >> prod3 >> prod4 >> prod5 
                   >> prod6 >> prod7;   
        if (prod0 == 0) {
          ErrorMessages::message("Error in ParticleDataTable::readXML:"
            " incomplete decay channel");
          return false;
        }

        // Store new channel (if particle already known).
        if (particlePtr == 0) {
          ErrorMessages::message("Error in ParticleDataTable::readXML:"
            " orphan decay channel");
          return false;
        }
        particlePtr->decay.addChannel(onMode, bRatio, meMode, prod0, prod1, 
          prod2, prod3, prod4, prod5, prod6, prod7);  
    
      // Check for occurence of a file also to be read.
      } else if (word1 == "<file") {
        string file = attributeValue(line, "name");
        if (file == "") {
          ErrorMessages::message("Warning in ParticleDataTable::readXML:"
            " skip unrecognized file name");
        } else files.push_back(file);
      }

    // End of loop over lines in input file and loop over files.
    };
  };

  // All particle data at this stage defines baseline original.
  if (reset) for (map<int, ParticleDataEntry>::iterator pdtEntry 
    = pdt.begin(); pdtEntry != pdt.end(); ++pdtEntry) {
    ParticleDataEntry* particlePtr = &pdtEntry->second;
    particlePtr->setHasChanged(false);
  }

  // Done.
  isInit = true;
  return true;

}

//*********
 
// Print out complete database in numerical order as an XML file.

void ParticleDataTable::listXML(string outFile) {

  // Convert file name to ofstream.
    const char* cstring = outFile.c_str();
    ofstream os(cstring);  

  // Iterate through the particle data table.
  for (map<int, ParticleDataEntry>::const_iterator pdtEntry 
    = pdt.begin(); pdtEntry != pdt.end(); ++pdtEntry) {
    const ParticleDataEntry* particlePtr = &pdtEntry->second;

    // Print particle properties.
    os << "<particle id=\"" << particlePtr->id() << "\""
       << " name=\"" << particlePtr->name() << "\"";
    if (particlePtr->hasAnti()) 
      os << " antiName=\"" << particlePtr->name(-1) << "\"";
    os << " spinType=\"" << particlePtr->spinType() << "\"" 
       << " chargeType=\"" << particlePtr->chargeType() << "\""
       << " colType=\"" << particlePtr->colType() << "\"\n" 
       << fixed << setprecision(5) 
       << "          m0=\"" << particlePtr->m0() << "\"";
    if (particlePtr->mWidth() > 0.) 
      os << " mWidth=\"" << particlePtr->mWidth() << "\""
         << " mMin=\"" << particlePtr->mMin() << "\""
         << " mMax=\"" << particlePtr->mMax() << "\"";
    if (particlePtr->tau0() > 0.) os << scientific
         << " tau0=\"" << particlePtr->tau0() << "\"";
    os << ">\n";

    // Loop through the decay channel table for each particle.
    if (particlePtr->decay.size() > 0) {
      for (int i = 0; i < int(particlePtr->decay.size()); ++i) {
        const DecayChannel& channel = particlePtr->decay[i];
        int mult = channel.multiplicity();

        // Print decay channel properties.
        os << " <channel onMode=\"" << channel.onMode() << "\"" 
           << fixed << setprecision(7)
           << " bRatio=\"" << channel.bRatio() << "\"";
        if (channel.meMode() > 0) 
          os << " meMode=\"" << channel.meMode() << "\"";  
        os << " products=\"";
	for (int j = 0; j < mult; ++j) {
          os << channel.product(j);
          if (j < mult - 1) os << " ";
	}

        // Finish off line and loop over allowed decay channels.
        os  << "\"/>\n";
      }
    }
        
    // Finish off existing particle.
    os << "</particle>\n\n";   

  } 

}

//*********

// Read in database from specific free format file.

bool ParticleDataTable::readFF(string inFile, bool reset) {

  // Normally reset whole database before beginning.
  if (reset) {pdt.clear(); isInit = false;}

  // Open file for read and check that instream is OK. 
  const char* cstring = inFile.c_str();
  ifstream is(cstring);  
  if (!is) {
    ErrorMessages::message("Error in ParticleDataTable::readFF:"
      " did not find file", inFile);
    return false;
  }

  // Read in one line at a time.
  particlePtr = 0;
  string line;
  bool readParticle = false;
  while ( getline(is, line) ) {

    // Empty lines begins new particle. 
    if (line.find_first_not_of(" ") == string::npos) {
      readParticle = true;
      continue;
    } 
 
    // Prepare to use standard read from line.
    istringstream readLine(line);

    // Read in a line with particle information.
    if (readParticle) {

      // Properties to be read. 
      int id;
      string name, antiName;
      int spinType, chargeType, colType;
      double m0, mWidth, mMin, mMax, tau0;
      string may;

      // Do the reading.
      readLine >> id >> name >> antiName >> spinType >> chargeType 
        >> colType >> m0 >> mWidth >> mMin >> mMax >> tau0;          

      // Error printout if something went wrong.
      if (!readLine) {
        ErrorMessages::message("Error in ParticleDataTable::readFF:"
          " incomplete particle");
        return false;
      }

      // Erase if particle already exists.
      if (isParticle(id)) pdt.erase(id);

      // Store new particle. Save pointer, to be used for decay channels.
      addParticle( id, name, antiName, spinType, chargeType, colType, 
        m0, mWidth, mMin, mMax, tau0);
      particlePtr = particleDataPtr(id);
      readParticle = false;

    // Read in a line with decay channel information.
    } else {
       
      // Properties to be read. 
      int onMode = 0;
      double bRatio = 0.;
      int meMode = 0;
      int prod0 = 0; int prod1 = 0; int prod2 = 0; int prod3 = 0;
      int prod4 = 0; int prod5 = 0; int prod6 = 0; int prod7 = 0;
  
      // Read in data from stream. Need at least one decay product.
      readLine >> onMode >> bRatio >> meMode >> prod0;
      if (!readLine) {
        ErrorMessages::message("Error in ParticleDataTable::readFF:"
          " incomplete decay channel");
        return false;
      }
      readLine >> prod1 >> prod2 >> prod3 >> prod4 >> prod5 
        >> prod6  >> prod7;   

      // Store new channel.
      if (particlePtr == 0) {
        ErrorMessages::message("Error in ParticleDataTable::readFF:"
          " orphan decay channel");
        return false;
      }
      particlePtr->decay.addChannel(onMode, bRatio, meMode, prod0, prod1, 
        prod2, prod3, prod4, prod5, prod6, prod7);  
    }

  // End of while loop over lines in the file.
  }


  // Done.
  isInit = true;
  return true;
   
}  

//*********
  
// Print out complete database in numerical order as a free format file.

void ParticleDataTable::listFF(string outFile) {

  // Convert file name to ofstream.
    const char* cstring = outFile.c_str();
    ofstream os(cstring);  

  // Iterate through the particle data table. 
  for (map<int, ParticleDataEntry>::const_iterator pdtEntry 
    = pdt.begin(); pdtEntry != pdt.end(); ++pdtEntry) {
    const ParticleDataEntry* particlePtr = &pdtEntry->second;

    // Print particle properties.
    os << "\n" << setw(8) << particlePtr->id() << "   " 
       << left << setw(16) << particlePtr->name() << "  " 
       << setw(16) << particlePtr->name(-1) << "  " 
       << right << setw(2) << particlePtr->spinType() << "  " 
       << setw(2) << particlePtr->chargeType() << "  " 
       << setw(2) << particlePtr->colType() << "  " 
       << fixed << setprecision(5) << setw(10) << particlePtr->m0() 
       << "  " << setw(10) << particlePtr->mWidth() << "  " 
       << setw(10) << particlePtr->mMin() << "  "
       << setw(10) << particlePtr->mMax() << "  " 
       << scientific << setw(12) << particlePtr->tau0() << "\n";

    // Loop through the decay channel table for each particle.
    if (particlePtr->decay.size() > 0) {
      for (int i = 0; i < int(particlePtr->decay.size()); ++i) {
        const DecayChannel& channel = particlePtr->decay[i];
        os << "          " << setw(6) << channel.onMode() << "  "
           << fixed << setprecision(7) << setw(10) 
           << channel.bRatio() << "  " 
           << setw(3) << channel.meMode() << "  ";
        for (int j = 0; j < channel.multiplicity(); ++j) 
          os << setw(8) << channel.product(j) << "  ";
        os << "\n";  
      }
    }  

  } 

}

//*********

// Read in updates from a character string, like a line of a file. 
// Is used by readString (and readFile) in Pythia.

bool ParticleDataTable::readString(string lineIn, bool warn,
  ostream& os) {

  // Take copy that will be modified.
  string line = lineIn;

  // If first character is not a digit then taken to be a comment.
  int firstChar = line.find_first_not_of(" ");
  if (!isdigit(line[firstChar])) return true; 

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
    if (warn) os << "\n Warning: input particle not found in Particle"
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
  if (property == "spintype") {
    int spinType;
    getWord >> spinType;
    pdt[id].setSpinType(spinType); 
    return true; 
  }
  if (property == "chargetype") {
    int chargeType;
    getWord >> chargeType;
    pdt[id].setChargeType(chargeType); 
    return true; 
  }
  if (property == "coltype") {
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
  if (property == "mwidth") {
    double mWidth;
    getWord >> mWidth;
    pdt[id].setMWidth(mWidth); 
    return true; 
  }
  if (property == "mmin") {
    double mMin;
    getWord >> mMin;
    pdt[id].setMMin(mMin); 
    return true; 
  }
  if (property == "mmax") {
    double mMax;
    getWord >> mMax;
    pdt[id].setMMax(mMax); 
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
    bool mayDecay = boolString(may);
    pdt[id].setMayDecay(mayDecay);
    return true; 
  }  
  if (property == "isresonance") {
    string isres;
    getWord >> isres;
    bool isResonance = boolString(isres);
    pdt[id].setIsResonance(isResonance);
    return true; 
  }
  if (property == "isvisible") {
    string isvis;
    getWord >> isvis;
    bool isVisible = boolString(isvis);
    pdt[id].setIsVisible(isVisible);
    return true; 
  }
  if (property == "externaldecay") {
    string extdec;
    getWord >> extdec;
    bool externalDecay = boolString(extdec);
    pdt[id].setExternalDecay(externalDecay);
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
    string name, antiName;
    int spinType, chargeType, colType;
    double m0, mWidth, mMin, mMax, tau0;

    // Read in data from stream.
    getWord >> name >> antiName >> spinType >> chargeType >> colType 
            >> m0 >> mWidth >> mMin >> mMax >> tau0;   

    // Error printout if something went wrong.
    if (!getWord) {
      os << "Error: incomplete particle " << id << "\n";
      return false;
    }

    // To keep existing decay channels, only overwrite particle data.
    if (property == "all" && isParticle(id)) {
      setAll( id, name, antiName, spinType, chargeType, colType, 
        m0, mWidth, mMin, mMax, tau0);   

    // Else start over completely from scratch.
    } else {
      if (isParticle(id)) pdt.erase(id);
      addParticle( id, name, antiName, spinType, chargeType, colType, 
        m0, mWidth, mMin, mMax, tau0);   
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
    if (property == "onmode" || property == "all") {
      int onMode = 0;
      string onModeIn;
      getWord >> onModeIn;
      // For onMode allow the optional possibility of Bool input.
      if (isdigit(onModeIn[0])) {
        istringstream getOnMode(onModeIn);
        getOnMode >> onMode;
      } else onMode = (boolString(onModeIn)) ? 1 : 0;
      pdt[id].decay[channel].onMode(onMode); 
      if (property == "onmode") return true; 
    }
    if (property == "bratio" || property == "all") {
      double bRatio;
      getWord >> bRatio;
      pdt[id].decay[channel].bRatio(bRatio); 
      if (property == "bratio") return true; 
    }
    if (property == "memode" || property == "all") {
      int meMode;
      getWord >> meMode;
      pdt[id].decay[channel].meMode(meMode); 
      if (property == "memode") return true; 
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
  if (warn) os << "\n Warning: input property not found in Particle"
    << " Data Table; skip:\n   " << lineIn << "\n";
  return false;

}

//*********
 
// Print out complete or changed table of database in numerical order.

void ParticleDataTable::list(bool changedOnly, ostream& os) {

  // Table header; output for bool as off/on.
  if (!changedOnly) {
    os << "\n --------  PYTHIA Particle Data Table (complete)  --------"
       << "------------------------------------------------------------"
       << "---------\n \n";

  } else { 
    os << "\n --------  PYTHIA Particle Data Table (changed only)  ----"
       << "------------------------------------------------------------" 
       << "---------\n \n";
  }
  os << "      id   name            antiName         spn chg col     m0"
     << "        mWidth      mMin       mMax       tau0    dec res vis "
     << "ext\n             no onMode   bRatio   meMode     products \n";

  // Iterate through the particle data table. Option to skip unchanged.
  int nList = 0;
  for (map<int, ParticleDataEntry>::const_iterator pdtEntry 
    = pdt.begin(); pdtEntry != pdt.end(); ++pdtEntry) {
    const ParticleDataEntry* particlePtr = &pdtEntry->second;
    if (changedOnly && !particlePtr->hasChanged()) continue;

    // Print particle properties.
    ++nList;
    string antiName = particlePtr->name(-1);
    if (antiName == "void") antiName = " ";
    os << "\n"  << setprecision(4)
       << setw(8) << particlePtr->id() << "   " 
       << left << setw(16) << particlePtr->name() 
       << setw(16) << antiName 
       << right << setw(4) << particlePtr->spinType() 
       << setw(4) << particlePtr->chargeType() 
       << setw(4) << particlePtr->colType() 
       << fixed << setw(10) << particlePtr->m0() << " " 
       << setw(10) << particlePtr->mWidth() << " " 
       << setw(10) << particlePtr->mMin() << " " 
       << setw(10) << particlePtr->mMax() << " "
       << scientific << setw(12) << particlePtr->tau0() << setw(4)
       << (particlePtr->mayDecay() && particlePtr->decay.size() > 0) 
       << setw(4) << particlePtr->isResonance() 
       << setw(4) << particlePtr->isVisible() 
       << setw(4) << particlePtr->externalDecay() << "\n";

    // Loop through the decay channel table for each particle.
    if (particlePtr->decay.size() > 0) {
      for (int i = 0; i < int(particlePtr->decay.size()); ++i) {
        const DecayChannel& channel = particlePtr->decay[i];
        os << "          "  << setprecision(7) 
           << setw(5) << i 
           << setw(6) << channel.onMode() 
           << fixed<< setw(12) << channel.bRatio() 
           << setw(5) << channel.meMode() << " ";
        for (int j = 0; j < channel.multiplicity(); ++j) 
          os << setw(8) << channel.product(j) << " ";
        os << "\n";  
      }
    }  

  } 

  // End of loop over database contents.
  if (changedOnly && nList == 0) os << "\n no particle data has been "
    << "changed from its default value \n";
  os << "\n --------  End Particle Data Table  ------------------------"
     << "--------------------------------------------------------------"
     << "-----\n" << endl;

}

//*********
 
// Print out partial table of database in input order.

void ParticleDataTable::list(vector<int> idList, ostream& os) {

  // Table header; output for bool as off/on.
  os << "\n --------  PYTHIA Particle Data Table (partial)  ---------"
     << "------------------------------------------------------------"
     << "---------\n \n";
  os << "      id   name            antiName         spn chg col     m0"
     << "        mWidth      mMin       mMax       tau0    dec res vis "
     << "ext\n             no onMode   bRatio   meMode     products \n";

  // Iterate through the given list of input particles.
  for (int i = 0; i < int(idList.size()); ++i) {
    const ParticleDataEntry* particlePtr = particleDataPtr(idList[i]);

    // Print particle properties.
    string antiName = particlePtr->name(-1);
    if (antiName == "void") antiName = " ";
    os << "\n"  << setprecision(4)
       << setw(8) << particlePtr->id() << "   " 
       << left << setw(16) << particlePtr->name() 
       << setw(16) << antiName 
       << right << setw(4) << particlePtr->spinType() 
       << setw(4) << particlePtr->chargeType() 
       << setw(4) << particlePtr->colType() 
       << fixed << setw(10) << particlePtr->m0() << " " 
       << setw(10) << particlePtr->mWidth() << " " 
       << setw(10) << particlePtr->mMin() << " " 
       << setw(10) << particlePtr->mMax() << " "
       << scientific << setw(12) << particlePtr->tau0() << setw(4)
       << (particlePtr->mayDecay() && particlePtr->decay.size() > 0) 
       << setw(4) << particlePtr->isResonance() 
       << setw(4) << particlePtr->isVisible() 
       << setw(4) << particlePtr->externalDecay() << "\n";

    // Loop through the decay channel table for each particle.
    if (particlePtr->decay.size() > 0) {
      for (int i = 0; i < int(particlePtr->decay.size()); ++i) {
        const DecayChannel& channel = particlePtr->decay[i];
        os << "          "  << setprecision(7) 
           << setw(5) << i 
           << setw(6) << channel.onMode() 
           << fixed<< setw(12) << channel.bRatio() 
           << setw(5) << channel.meMode() << " ";
        for (int j = 0; j < channel.multiplicity(); ++j) 
          os << setw(8) << channel.product(j) << " ";
        os << "\n";  
      }
    }  

  } 

  // End of loop over database contents.
  os << "\n --------  End Particle Data Table  ------------------------"
     << "--------------------------------------------------------------"
     << "-----\n" << endl;

}

//**************************************************************************

} // end namespace Pythia8
