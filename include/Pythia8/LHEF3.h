// LHEF3.h is a part of the PYTHIA event generator.
// Copyright (C) 2014 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file is written by Stefan Prestel. The code evolved from
// a LHEF 2.0 reader supplied by Leif Lonnblad.
// LHEF3.h contains the main class for LHEF 3.0 functionalities.
// Header file.

#ifndef Pythia8_LHEF3_H
#define Pythia8_LHEF3_H

#include "Pythia8/PythiaStdlib.h"
#include <stdexcept>

namespace Pythia8 {

//==========================================================================

// The XMLTag struct is used to represent all information within an XML tag. 
// It contains the attributes as a map, any sub-tags as a vector of pointers 
// to other XMLTag objects, and any other information as a single string.
// The XMLTag struct written by Leif Lonnblad.

struct XMLTag {

  // Convenient typdef.
  typedef std::string::size_type pos_t;

  // Convenient alias for npos.
  static const pos_t end = std::string::npos;

  // The destructor also destroys any sub-tags.
  ~XMLTag() {
    for ( int i = 0, N = tags.size(); i < N; ++i ) delete tags[i];
  }

  // The name of this tag.
  std::string name;

  // The attributes of this tag.
  std::map<std::string,std::string> attr;

  // A vector of sub-tags.
  std::vector<XMLTag*> tags;

  // The contents of this tag.
  std::string contents;

  // Find an attribute named n and set the double variable v to
  // the corresponding value. Return false if no attribute was found.
  bool getattr(std::string n, double & v) const {
    std::map<std::string,std::string>::const_iterator it = attr.find(n);
    if ( it == attr.end() ) return false;
    v = std::atof(it->second.c_str());
    return true;
  }

  // Find an attribute named n and set the bool variable v to true if the
  // corresponding value is "yes". Return false if no attribute was found.
  bool getattr(std::string n, bool & v) const {
    std::map<std::string,std::string>::const_iterator it = attr.find(n);
    if ( it == attr.end() ) return false;
    if ( it->second == "yes" ) v = true;
    return true;
  }

  // Find an attribute named n and set the long variable v to the
  // corresponding value. Return false if no attribute was found.
  bool getattr(std::string n, long & v) const {
    std::map<std::string,std::string>::const_iterator it = attr.find(n);
    if ( it == attr.end() ) return false;
    v = std::atoi(it->second.c_str());
    return true;
  }

  // Find an attribute named n and set the long variable v to the
  // corresponding value. Return false if no attribute was found.
  bool getattr(std::string n, int & v) const {
    std::map<std::string,std::string>::const_iterator it = attr.find(n);
    if ( it == attr.end() ) return false;
    v = int(std::atoi(it->second.c_str()));
    return true;
  }

  // Find an attribute named n and set the string variable v to the
  // corresponding value. Return false if no attribute was found.
  bool getattr(std::string n, std::string & v) const {
    std::map<std::string,std::string>::const_iterator it = attr.find(n);
    if ( it == attr.end() ) return false;
    v = it->second;
    return true;
  }

  // Scan the given string and return all XML tags found as a vector
  // of pointers to XMLTag objects.
  static std::vector<XMLTag*> findXMLTags(std::string str,
    std::string * leftover = 0) {
    std::vector<XMLTag*> tags;
    pos_t curr = 0;

    while ( curr != end ) {

      // Find the first tag.
      pos_t begin = str.find("<", curr);

      // Skip comments.
      if ( str.find("<!--", curr) == begin ) {
        pos_t endcom = str.find("-->", begin);
        if ( endcom == end ) {
          if ( leftover ) *leftover += str.substr(curr);
          return tags;
        }
        if ( leftover ) *leftover += str.substr(curr, endcom - curr);
        curr = endcom;
        continue;
      }

      if ( leftover ) *leftover += str.substr(curr, begin - curr);
      if ( begin == end || begin > str.length() - 3 || str[begin + 1] == '/' )
        return tags; 

      pos_t close = str.find(">", curr);
      if ( close == end ) return tags;

      // Find the tag name.
      curr = str.find_first_of(" \t\n/>", begin);
      tags.push_back(new XMLTag());
      tags.back()->name = str.substr(begin + 1, curr - begin - 1);

      while ( true ) {

        // Now skip some white space to see if we can find an attribute.
        curr = str.find_first_not_of(" \t\n", curr);
        if ( curr == end || curr >= close ) break;

        pos_t tend = str.find_first_of("= \t\n", curr);
        if ( tend == end || tend >= close ) break;

        std::string name = str.substr(curr, tend - curr);
        curr = str.find("=", curr) + 1;

        // OK now find the beginning and end of the atribute.
        curr = str.find("\"", curr);
        if ( curr == end || curr >= close ) break;
        pos_t bega = ++curr;
        curr = str.find("\"", curr);
        while ( curr != end && str[curr - 1] == '\\' )
          curr = str.find("\"", curr + 1);

        std::string value = str.substr(bega, curr == end? end: curr - bega);

        tags.back()->attr[name] = value;

        ++curr;

      }

      curr = close + 1;
      if ( str[close - 1] == '/' ) continue;

      pos_t endtag = str.find("</" + tags.back()->name + ">", curr);
      if ( endtag == end ) {
        tags.back()->contents = str.substr(curr);
        curr = endtag;
      } else {
        tags.back()->contents = str.substr(curr, endtag - curr);
        curr = endtag + tags.back()->name.length() + 3;
      }

      std::string leftovers;
      tags.back()->tags = findXMLTags(tags.back()->contents, &leftovers);
      if ( leftovers.find_first_not_of(" \t\n") == end ) leftovers="";
      tags.back()->contents = leftovers;

    }
    return tags;

  }

  // Print out this tag to a stream.
  void print(std::ostream & os) const {
    os << "<" << name;
    for ( std::map<std::string,std::string>::const_iterator it = attr.begin();
          it != attr.end(); ++it )
      os << " " << it->first << "=\"" << it->second << "\"";
    if ( contents.empty() && tags.empty() ) {
      os << "/>" << std::endl;
      return;
    }
    os << ">" << std::endl;
    for ( int i = 0, N = tags.size(); i < N; ++i )
      tags[i]->print(os);

    os << "````" << contents << "''''</" << name << ">" << std::endl;
  }

};

//==========================================================================

// The LHAweights struct represents the information in a weights tag.

struct LHAweights {

  // Initialize default values.
  LHAweights() {}

  // Construct from XML tag
  LHAweights(const XMLTag & tag);

  // Print out an XML tag.
  void print(ostream & file) const;

  // Function to reset this object.
  void clear() {
    contents="";
    weights.clear();
    attributes.clear();
  }

  // The weights of this event.
  vector<double> weights;

  // Any other attributes.
  map<string,string> attributes;

  // The contents of the tag.
  string contents;

};

//==========================================================================

// Collect different scales relevant for an event.

struct LHAscales {

  // Empty constructor.
  LHAscales(double defscale = -1.0)
  : muf(defscale), mur(defscale), mups(defscale), SCALUP(defscale) {}

  // Construct from an XML-tag
  LHAscales(const XMLTag & tag, double defscale = -1.0);

  // Print out the corresponding XML-tag.
  void print(ostream & file) const;

  // Function to reset this object.
  void clear() {
    contents="";
    muf=mur=mups=SCALUP;
    attributes.clear();
  }

  // The factorization scale used for this event.
  double muf;

  // The renormalization scale used for this event.
  double mur;

  // The starting scale for the parton shower as suggested by the
  // matrix element generator.
  double mups;

  // Any other scales reported by the matrix element generator.
  map<string,double> attributes;

  // The default scale in this event.
  double SCALUP;

  // The contents of the tag.
  string contents;

};

//==========================================================================

// Collect generator information for an event file.

struct LHAgenerator {

  // Empty constructor.
  LHAgenerator()
  : name(""), version(""), contents("") {}

  // Construct from an XML-tag
  LHAgenerator(const XMLTag & tag, string defname = "");

  // Print out the corresponding XML-tag.
  void print(ostream & file) const;

  // Function to reset this object.
  void clear() {
    contents="";
    name="";
    version="";
    attributes.clear();
  }

  // The generator name used for this file.
  string name;

  // The generator version used for this file.
  string version;

  // Any other attributes.
  map<string,string> attributes;

  // The contents of the tag.
  string contents;

};

//==========================================================================

// Collect the wgt information.

struct LHAwgt {

  // Empty constructor.
  LHAwgt(double defwgt = 1.0)
  : id(""), contents(defwgt) {}

  // Construct from an XML-tag
  LHAwgt(const XMLTag & tag, double defwgt = 1.0);

  // Print out the corresponding XML-tag.
  void print(ostream & file) const;

  // Function to reset this object.
  void clear() {
    contents=0.0;
    id="";
    attributes.clear();
  }

  // The identification number of this wgt tag.
  string id;

  // Any other attributes.
  map<string,string> attributes;

  // The weight associated to this tag.
  double contents;

};

//==========================================================================

// Collect the wgt information.

struct LHAweight {

  // Empty constructor.
  LHAweight(string defname = "")
  : id(defname), contents(defname) {}

  // Construct from an XML-tag
  LHAweight(const XMLTag & tag, string defname = "");

  // Print out the corresponding XML-tag.
  void print(ostream & file) const;

  // Function to reset this object.
  void clear() {
    contents="";
    id="";
    attributes.clear();
  }

  // The identification number of this weight tag.
  string id;

  // Any other attributes.
  map<string,string> attributes;

  // The weight description associated to this tag.
  string contents;

};

//==========================================================================

// The LHAweightgroup assigns a group-name to a set of LHAweight objects.

struct LHAweightgroup {

  // Default constructor;
  LHAweightgroup() {}

  // Construct a group of LHAweight objects from an XML tag and
  // insert them in the given vector.
  LHAweightgroup(const XMLTag & tag);

  // Print out the corresponding XML-tag.
  void printweights(ostream & file) const;

  // Function to reset this object.
  void clear() {
    contents="";
    name="";
    weights.clear();
    attributes.clear();
  }

  // The contents of the tag.
  string contents;

  // The name.
  string name;

  // The vector of weights.
  map<string, LHAweight> weights;

  // Any other attributes.
  map<string,string> attributes;

};

//==========================================================================

// The LHArwgt assigns a group-name to a set of LHAwgt objects.

struct LHArwgt {

  // Default constructor;
  LHArwgt() {}

  // Construct a group of LHAwgt objects from an XML tag and
  // insert them in the given vector.
  LHArwgt(const XMLTag & tag);

  // Function to reset this object.
  void clear() {
    contents="";
    wgts.clear();
    attributes.clear();
  }

  // The contents of the tag.
  string contents;

  // The map of weights.
  map<string, LHAwgt> wgts;

  // Any other attributes.
  map<string,string> attributes;

};

//==========================================================================

// The LHAinitrwgt assigns a group-name to a set of LHAweightgroup objects.

struct LHAinitrwgt {

  // Default constructor;
  LHAinitrwgt() {}

  // Construct a group of LHAweightgroup objects from an XML tag and
  // insert them in the given vector.
  LHAinitrwgt(const XMLTag & tag);

  // Function to reset this object.
  void clear() {
    contents="";
    weights.clear();
    weightgroups.clear();
    attributes.clear();
  }

  // The contents of the tag.
  string contents;

  // The vector of weight's.
  map<string, LHAweight> weights;

  // The vector of weightgroup's.
  map<string, LHAweightgroup> weightgroups;

  // Any other attributes.
  map<string,string> attributes;

};

//==========================================================================

// The HEPRUP class is a simple container corresponding to the Les Houches
// accord (<A HREF="http://arxiv.org/abs/hep-ph/0109068">hep-ph/0109068</A>)
// common block with the same name. The members are named in the same
// way as in the common block. However, fortran arrays are represented
// by vectors, except for the arrays of length two which are
// represented by pair objects.

class HEPRUP {

public:

  // Default constructor.
  HEPRUP() : IDWTUP(0), NPRUP(0) {}

  // Assignment operator.
  HEPRUP & operator=(const HEPRUP & x) {
    IDBMUP = x.IDBMUP;
    EBMUP = x.EBMUP;
    PDFGUP = x.PDFGUP;
    PDFSUP = x.PDFSUP;
    IDWTUP = x.IDWTUP;
    NPRUP = x.NPRUP;
    XSECUP = x.XSECUP;
    XERRUP = x.XERRUP;
    XMAXUP = x.XMAXUP;
    LPRUP = x.LPRUP;
    initrwgt = x.initrwgt;
    generators = x.generators;
    weightgroups = x.weightgroups;
    weights = x.weights;
    return *this;
  }

  // Destructor.
  ~HEPRUP() {}

  // Set the NPRUP variable, corresponding to the number of
  // sub-processes, to \a nrup, and resize all relevant vectors
  // accordingly.
  void resize(int nrup) {
    NPRUP = nrup;
    resize();
  }

  // Assuming the NPRUP variable, corresponding to the number of
  // sub-processes, is correctly set, resize the relevant vectors
  // accordingly.
  void resize() {
    XSECUP.resize(NPRUP);
    XERRUP.resize(NPRUP);
    XMAXUP.resize(NPRUP);
    LPRUP.resize(NPRUP);
  }

  // PDG id's of beam particles. (first/second is in +/-z direction).
  pair<long,long> IDBMUP;

  // Energy of beam particles given in GeV.
  pair<double,double> EBMUP;

  // The author group for the PDF used for the beams according to the
  // PDFLib specification.
  pair<int,int> PDFGUP;

  // The id number the PDF used for the beams according to the
  // PDFLib specification.
  pair<int,int> PDFSUP;

  // Master switch indicating how the ME generator envisages the
  // events weights should be interpreted according to the Les Houches
  // accord.
  int IDWTUP;

  // The number of different subprocesses in this file.
  int NPRUP;

  // The cross sections for the different subprocesses in pb.
  vector<double> XSECUP;

  // The statistical error in the cross sections for the different
  // subprocesses in pb.
  vector<double> XERRUP;

  // The maximum event weights (in HEPEUP::XWGTUP) for different
  // subprocesses.
  vector<double> XMAXUP;

  // The subprocess code for the different subprocesses.
  vector<int> LPRUP;

  // Contents of the LHAinitrwgt tag
  LHAinitrwgt initrwgt;

  // Contents of the LHAgenerator tags.
  vector<LHAgenerator> generators;

  // A map of the LHAweightgroup tags, indexed by name.
  map<string,LHAweightgroup> weightgroups;

  // A map of the LHAweight tags, indexed by name.
  map<string,LHAweight> weights;

};

//==========================================================================

// The HEPEUP class is a simple container corresponding to the Les Houches 
// accord (<A HREF="http://arxiv.org/abs/hep-ph/0109068">hep-ph/0109068</A>)
// common block with the same name. The members are named in the same
// way as in the common block. However, fortran arrays are represented
// by vectors, except for the arrays of length two which are
// represented by pair objects.

class HEPEUP {

public:

  // Default constructor.
  HEPEUP()
    : NUP(0), IDPRUP(0), XWGTUP(0.0), XPDWUP(0.0, 0.0),
      SCALUP(0.0), AQEDUP(0.0), AQCDUP(0.0), heprup(0) {}

  // Copy constructor
  HEPEUP(const HEPEUP & x) {
    operator=(x);
  }

  // Copy information from the given HEPEUP.
  HEPEUP & setEvent(const HEPEUP & x);

  // Assignment operator.
  HEPEUP & operator=(const HEPEUP & x) {
    clear();
    setEvent(x);
    return *this;
  }

  // Destructor.
  ~HEPEUP() {
    clear();
  };

  // Reset the HEPEUP object.
  void reset();

  // Clear the HEPEUP object.
  void clear() {
    reset();
  }

  // Set the NUP variable, corresponding to the number of particles in
  // the current event, to \a nup, and resize all relevant vectors
  // accordingly.
  void resize(int nup) {
    NUP = nup;
    resize();
  }

  // Return the main weight for this event.
  double weight() const {
      return XWGTUP;
  }

  // Assuming the NUP variable, corresponding to the number of
  // particles in the current event, is correctly set, resize the
  // relevant vectors accordingly.
  void resize();

  // The number of particle entries in the current event.
  int NUP;

  // The subprocess code for this event (as given in LPRUP).
  int IDPRUP;

  // The weight for this event.
  double XWGTUP;

  // The PDF weights for the two incoming partons. Note that this
  // variable is not present in the current LesHouches accord
  // (<A HREF="http://arxiv.org/abs/hep-ph/0109068">hep-ph/0109068</A>),
  // hopefully it will be present in a future accord.
  pair<double,double> XPDWUP;

  // The scale in GeV used in the calculation of the PDF's in this
  // event.
  double SCALUP;

  // The value of the QED coupling used in this event.
  double AQEDUP;

  // The value of the QCD coupling used in this event.
  double AQCDUP;

  // The PDG id's for the particle entries in this event.
  vector<long> IDUP;

  // The status codes for the particle entries in this event.
  vector<int> ISTUP;

  // Indices for the first and last mother for the particle entries in
  // this event.
  vector< pair<int,int> > MOTHUP;

  // The colour-line indices (first(second) is (anti)colour) for the
  // particle entries in this event.
  vector< pair<int,int> > ICOLUP;

  // Lab frame momentum (Px, Py, Pz, E and M in GeV) for the particle
  // entries in this event.
  vector< vector<double> > PUP;

  // Invariant lifetime (c*tau, distance from production to decay in
  // mm) for the particle entries in this event.
  vector<double> VTIMUP;

  // Spin info for the particle entries in this event given as the
  // cosine of the angle between the spin vector of a particle and the
  // 3-momentum of the decaying particle, specified in the lab frame.
  vector<double> SPINUP;

  // A pointer to the current HEPRUP object.
  HEPRUP * heprup;

  // The weights associated with this event, as given by the LHAwgt tags.
  map<string,double> weights_detailed;

  // The weights associated with this event, as given by the LHAweights tags.
  vector<double> weights_compressed;

  // Contents of the LHAscales tag
  LHAscales scales;

  // Contents of the LHAweights tag (compressed format)
  LHAweights weights;

  // Contents of the LHArwgt tag (detailed format)
  LHArwgt rwgt;

  // Any other attributes.
  map<string,string> attributes;

};

//==========================================================================

// The Reader class is initialized with a stream from which to read a
// version 1/3 Les Houches Accord event file. In the constructor of
// the Reader object the optional header information is read and then
// the mandatory init is read. After this the whole header block
// including the enclosing lines with tags are available in the public
// headerBlock member variable. Also the information from the init
// block is available in the heprup member variable and any additional
// comment lines are available in initComments. After each successful
// call to the readEvent() function the standard Les Houches Accord
// information about the event is available in the hepeup member
// variable and any additional comments in the eventComments
// variable. A typical reading sequence would look as follows:

class Reader {

public:

  // Initialize the Reader with a stream from which to read an event
  // file. After the constructor is called the whole header block
  // including the enclosing lines with tags are available in the
  // public headerBlock member variable. Also the information from the
  // init block is available in the heprup member variable and any
  // additional comment lines are available in initComments.
  //
  // is: the stream to read from.
  //
  Reader(istream & is)
    : file(is) {
    isGood = init();
  }

  // Initialize the Reader with a filename from which to read an event
  // file. After the constructor is called the whole header block
  // including the enclosing lines with tags are available in the
  // public headerBlock member variable. Also the information from the
  // init block is available in the heprup member variable and any
  // additional comment lines are available in initComments.
  //
  // filename: the name of the file to read from.
  //
  Reader(string filename)
    : intstream(filename.c_str()), file(intstream) {
    init();
  }

private:

  // Used internally in the constructors to read header and init blocks.
  bool init();

public:

  // Read an event from the file and store it in the hepeup
  // object. Optional comment lines are stored in the eventComments
  // member variable.
  bool readEvent(HEPEUP * peup = 0);

protected:

  // Used internally to read a single line from the stream.
  bool getLine() {
    currentLine = "";
    if(!getline(file, currentLine)) return false;
    // Replace single by double quotes
    replace(currentLine.begin(),currentLine.end(),'\'','\"');
    return true;
  }

protected:

  // A local stream which is unused if a stream is supplied from the
  // outside.
  ifstream intstream;

  // The stream we are reading from. This may be a reference to an
  // external stream or the internal intstream.
  istream & file;

  // The last line read in from the stream in getline().
  string currentLine;

public:

  // Save if the initialisation worked.
  bool isGood;

  // XML file version
  int version;

  // All lines (since the last readEvent()) outside the header, init
  // and event tags.
  string outsideBlock;

  // All lines from the header block.
  string headerBlock;
  string headerComments;

  // The standard init information.
  HEPRUP heprup;

  // Additional comments found in the init block.
  string initComments;

  // The standard information about the last read event.
  HEPEUP hepeup;

  // Additional comments found with the last read event.
  string eventComments;

private:

  // The default constructor should never be used.
  Reader();

  // The copy constructor should never be used.
  Reader(const Reader &);

  // The Reader cannot be assigned to.
  Reader & operator=(const Reader &);

};

} // end namespace Pythia8

#endif // end Pythia8_LHEF3_H
