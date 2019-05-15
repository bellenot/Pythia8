#ifndef SLHA_H
#define SLHA_H

// Stdlib header files for string and character manipulation.
#include <string>
#include <cctype>
// Stdlib header files for containers.
#include <vector>
#include <map>
// Stdlib header files for input/output.
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

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
// Generic utilities.
using std::tolower; 
using std::swap;

class SusyLesHouches {

public:

  //Constructor, with and without filename.
  SusyLesHouches(int verboseIn=1) : verbose(verboseIn) {};
  SusyLesHouches(string filename, int verboseIn=1) : verbose(verboseIn) {
    readFile(filename);};

  //***************************** SLHA FILE I/O *****************************//
  //readFile(string filename) : read SLHA file from filename
  //writeFile(string filename): write SLHA file on filename
  //checkSpectrum() : check current blocks for consistency of model
  int readFile(string="slha.spc");
  int checkSpectrum();

  //***************************** SLHA CLASSES *****************************//
  //class block: the generic SLHA block (see below for matrices)
  //Explicit typing required, e.g. block<double> minpar;
  template <class T> class block {    

  public: 

    //Constructor. 
    block<T>() : idnow(0) { } ;    

    //Does block exist?
    bool exists() { return entry.size() == 0 ? false : true ; };
    //Clear block
    void clear() { entry.clear(); };

    //set: set block entry values.
    //Possible return values from set:
    // 0: normal return. Entry did not previously exist and has been created.
    // 1: normal return. Entry did previously exist and has been overwritten.
    //-1: failure. 
    int set(int i,T val) { 
      int alreadyexisting=exists(i)?1:0;
      entry[i]=val; 
      return alreadyexisting;
    };
    // Read index and value from SLHA data line
    int set(istringstream& linestream) {
      linestream >> i >> val;
      return linestream ? set(i,val) : -1;
    };
    // With i already given, read value from remaining SLHA data line
    int set(int i,istringstream& linestream) {
      linestream >> val;
      return linestream ? set(i,val) : -1;
    };
    // Shorthand for entry[0]. Used e.g. for block ALPHA.
    void set(T val) { entry[0]=val; };

    // Does entry i already exist in this block?
    bool exists(int i=0) {return entry.find(i) != entry.end() ? true : false;};

    // Indexing with (). Output only.
    T operator()(int i=0) {
      if (exists(i)) {return entry[i];} else {T dummy(0); return dummy;};
    };

    // Size of map
    int size() {return entry.size();};

    // First and next key code
    int first() { idnow = entry.begin()->first; return idnow; };
    int next() { 
      typename map<int,T>::iterator itnow;
      itnow = ++entry.find(idnow);
      if ( itnow == entry.end() ) itnow=entry.begin();
      return idnow = itnow->first;
    };

    // Simple print utility
    void print() {
      bool finished=false;
      int ibegin=first();
      i=ibegin;
      while (not finished) {
	cout << "  "<< i << " " << entry[i] <<endl;
	i=next();
	if (i == ibegin) finished=true;
      };       
    };

    // Special for DRbar running blocks.
    void setq(double qIn) { qDRbar=qIn; }
    double q() { return qDRbar; }
 
  private:
    map<int,T> entry;    
    int idnow;
    double qDRbar;
    //Auxiliary vars
    int i; 
    T val;
  };

  // class matrixblock: the generic SLHA matrix 
  // Explicit sizing required, e.g. matrixblock<4> nmix;
  template <int size> class matrixblock {    
  public: 
    //Constructor. Set uninitialized and explicitly zero.
    matrixblock<size>() { 
      initialized=false; 
      for (i=1;i<=size;i++) {
	for (j=1;j<=size;j++) {
	  entry[i][j]=0.0;
	};
      };
    };    

    // Assignment
    matrixblock& operator=(const matrixblock& m) { 
      if (this != &m) { 
	for (i=0;i<size;i++) for (j=0;j<=size;j++) entry[i][j] = m(i,j);
	qDRbar = m.qDRbar; 
	initialized = m.initialized; 
      } 
      return *this; };

    // Does this matrix contain any entries?
    bool exists() { return initialized; };
    // Clear initialized flag
    void clear() { initialized=false; };

    // Set matrix entry
    int set(int i,int j, double val) { 
      if (i>0 and j>0 and i<=size and j<=size) {
	entry[i][j]=val;
	initialized=true;
	return 0;
      } else {
	return -1;
      };
    };

    // Set entry from linestream (used during file read)
    int set(istringstream& linestream) {
      linestream >> i >> j >> val;
      return linestream ? set(i,j,val) : -1;
    };

    // () Overloading: Get entry
    double operator()(int i, int j) {
      return (i <= size and j <= size and i > 0 and j > 0) ? 
	entry[i][j] : 0.0;
    };

    // Set and get scale for DRbar running blocks.
    void setq(double qIn) { qDRbar=qIn; }
    double q() { return qDRbar; }

    // Simple print utility, to be elaborated on.
    void print() {
      for (i=1;i<=size;i++) {
	cout << "   "<<i << " " ;
	for (j=1;j<=size;j++) cout << entry[i][j] << " ";
	cout << endl;
      };
    };

  private:
    bool initialized;
    double entry[size+1][size+1];
    double qDRbar;
    //Auxiliary vars
    int i,j; 
    double val;
  };

  //*************************** THE SLHA1 BLOCKS ***************************//
  //blocks for model definition:
  block<int> modsel;
  block<int> modsel21;
  block<double> modsel12;
  block<double> minpar;
  block<double> extpar;
  block<double> sminputs;
  //blocks for RGE program specific output
  block<string> spinfo;
  block<string> spinfo3;
  block<string> spinfo4;
  //blocks for mass and coupling spectrum
  block<double> mass;
  matrixblock<4> nmix;
  matrixblock<2> umix;
  matrixblock<2> vmix;
  matrixblock<2> stopmix;
  matrixblock<2> sbotmix;
  matrixblock<2> staumix;
  block<double> alpha;
  block<double> hmix;
  block<double> gauge;
  block<double> msoft;
  matrixblock<3> au;
  matrixblock<3> ad;
  matrixblock<3> ae;
  matrixblock<3> yu;
  matrixblock<3> yd;
  matrixblock<3> ye;

  //*************************** THE SLHA2 BLOCKS ***************************//
  //Input blocks:
  matrixblock<3> vckmin; // The input CKM matrix in PDG parm.
  matrixblock<3> msq2in; // The input upper off-diagonal msq2
  matrixblock<3> msu2in; // The input upper off-diagonal msu2
  matrixblock<3> msd2in; // The input upper off-diagonal msd2
  //Output blocks:
  matrixblock<3> vckm;   // The output DRbar running VCKM at Q
  matrixblock<3> msq2;   // The output DRbar running msq2 at Q
  matrixblock<3> msu2;   // The output DRbar running msu2 at Q
  matrixblock<3> msd2;   // The output DRbar running msd2 at Q
  matrixblock<6> usqmix; // The up squark mixing matrix
  matrixblock<6> dsqmix; // The down squark mixing matrix
  matrixblock<6> imusqmix; // The up squark mixing matrix
  matrixblock<6> imdsqmix; // The down squark mixing matrix
  matrixblock<4> imnmix; // The up squark mixing matrix

  //*************************** SET BLOCK VALUE ****************************//
  template <class T> int set(string,T);
  template <class T> int set(string,int,T);
  template <class T> int set(string,int,int,T);

  //***************************** SLHA PRIVATE *****************************//
private:
  //SLHA I/O
  string spectrumfile;
  void message(int, string,string ,int);
  int verbose;

};

#endif
