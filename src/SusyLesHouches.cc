#include "SusyLesHouches.h"

//****************************** SLHA FILE I/O ******************************//
int SusyLesHouches::readFile(string slhafile) {

  // Check that input file is OK.
  int ifailfile=0;
  spectrumfile=slhafile;
  const char* cstring = slhafile.c_str();
  ifstream file(cstring);  
  if (!file) {
    message(2,"readFile","slha file "+slhafile+" not found",0);
    return -1;}

  if (verbose >= 1) message(0,"readFile","parsing SLHA file "+slhafile,0);

  // Read in one line at a time.
  string line="";
  string block="";
  string decay="";
  string blockname="";
  int nline=0;

  while ( getline(file, line) ) {
    nline++;

    //Ignore comment lines with # as first character
    if (line.find("#") == 0) continue;

    //Rewrite string in lowercase
    for (unsigned int i=0;i<line.length();i++) line[i]=tolower(line[i]);

    // Remove extra blanks before and after an = sign.
    while (line.find(" =") != string::npos) line.erase( line.find(" ="), 1);
    while (line.find("= ") != string::npos) line.erase( line.find("= ")+1, 1);

    //New block. 
    if (line.substr(0,1) != " ") {
      if (line.find("block") == 0) { 
	block=line ; 
	decay="";
	int NameBegin=6 ;
	int NameEnd=block.find(" ",7);
	blockname=line.substr(NameBegin,NameEnd-NameBegin);

	//Find Q=... for DRbar running blocks
	if (block.find("q=") != string::npos) {
	  int qbegin=block.find("q=")+2;
	  istringstream qstream(block.substr(qbegin,block.length()));
	  double q=0.0;
	  qstream >> q;
	  if (qstream) {
	    if (blockname=="hmix") hmix.setq(q);
	    if (blockname=="yu") yu.setq(q);
	    if (blockname=="yd") yd.setq(q);
	    if (blockname=="ye") ye.setq(q);
	    if (blockname=="au") au.setq(q);
	    if (blockname=="ad") ad.setq(q);
	    if (blockname=="ae") ae.setq(q);
	    if (blockname=="msoft") msoft.setq(q);
	    if (blockname=="gauge") gauge.setq(q);
	  };
	};

	//Skip to next line.
	continue ; 

      } else {
	//Unrecognized first character. Skip until next block.
	block="";
      };
    };

    //Skip not currently reading block data lines.
    if (block == "") continue;

    // Replace an equal sign by a blank to make parsing simpler.
    while (line.find("=") != string::npos) {
      int firstEqual = line.find_first_of("=");
      line.replace(firstEqual, 1, " ");   
    };
    
    //Parse data lines within given block
    //Constructed explicitly so that each block can have its own types and
    //own rules defined. For extra user blocks, just add more recognized 
    //blocknames at the end and implement user defined rules accordingly.
    //string comment = line.substr(line.find("#"),line.length());    
    int ifail=-2;
    istringstream linestream(line);
    //MODEL
    if (blockname == "modsel") {
      int i;
      linestream >> i; 
      if (linestream) {
	if (i == 12) {ifail=modsel12.set(0,linestream);} 
	else if (i == 21) {ifail=modsel21.set(0,linestream);}
	else {ifail=modsel.set(i,linestream);};}
      else {
	ifail = -1;}
    };
    if (blockname == "minpar") ifail=minpar.set(linestream); 
    if (blockname == "sminputs") ifail=sminputs.set(linestream);
    if (blockname == "extpar") ifail=extpar.set(linestream);
    if (blockname == "spinfo") {
      int i;
      string entry;
      linestream >> i >> entry;
      if (linestream) {
	if ( i == 3 ) {
	  string warning=line.substr(line.find("3")+1,line.length());
	  message(1,"readFile","(from RGE program): "+warning,0);
	  spinfo3.set(warning);
	} else if ( i == 4 ) {
	  string error=line.substr(line.find("4")+1,line.length());
	  message(2,"readFile","(from RGE program): "+error,0);
	  spinfo4.set(error);
	} else {
	  ifail=spinfo.set(i,entry);
	};
      } else {
	ifail=-1;
      };
    };
    //SPECTRUM
    //Pole masses
    if (blockname == "mass") ifail=mass.set(linestream);
    //Mixing
    if (blockname == "alpha") ifail=alpha.set(linestream);
    if (blockname == "stopmix") ifail=stopmix.set(linestream);
    if (blockname == "sbotmix") ifail=sbotmix.set(linestream);
    if (blockname == "staumix") ifail=staumix.set(linestream);
    if (blockname == "nmix") ifail=nmix.set(linestream);
    if (blockname == "umix") ifail=umix.set(linestream);
    if (blockname == "vmix") ifail=vmix.set(linestream);
    //DRbar Lagrangian parameters
    if (blockname == "gauge") ifail=gauge.set(linestream);      
    if (blockname == "yu") ifail=yu.set(linestream);
    if (blockname == "yd") ifail=yd.set(linestream);
    if (blockname == "ye") ifail=ye.set(linestream);
    if (blockname == "au") ifail=au.set(linestream);
    if (blockname == "ad") ifail=ad.set(linestream);
    if (blockname == "ae") ifail=ae.set(linestream);
    if (blockname == "hmix") ifail=hmix.set(linestream);
    if (blockname == "msoft") ifail=msoft.set(linestream);

    //Diagnostics
    if (ifail != 0) { 
      if (ifail == -2) {
	message(1,"readFile","Ignoring unknown block: "+blockname,nline);
	block="";
      };
      if (ifail == -1) {
	message(2,"readFile","Error on line.",nline);	
      };
      if (ifail == 1) {
	message(0,"readFile",blockname+" existing entry overwritten.",nline);
      };
    };

  };

  return ifailfile;
    
}

int SusyLesHouches::checkSpectrum() {
  int ifail=0;
  //1) Check MODSEL. Assign default values where applicable.
  if (!modsel.exists(1)) {
    message(1,"checkSpectrum","MODSEL(1) undefined. Assuming =0.",0);
    modsel.set(1,0);
  }
  if (!modsel.exists(3)) modsel.set(3,0);
  if (!modsel.exists(4)) modsel.set(4,0);
  if (!modsel.exists(5)) modsel.set(5,0);
  if (!modsel.exists(6)) modsel.set(6,0);
  if (!modsel.exists(11)) modsel.set(11,1);
  
  //2) Check for existence / duplication of mixing matrices
  if (!stopmix.exists() && !usqmix.exists()) {
    message(1,"checkSpectrum","No mixing matrix for sup sector.",0); 
    //Ensure sensible defaults if no mixing specified.
    stopmix.set(1,1, 0.0) ; stopmix.set(1,2, 1.0);
    stopmix.set(2,1,-1.0) ; stopmix.set(2,2, 0.0);
  };
  if (stopmix.exists() && usqmix.exists()) {
    message(1,"checkSpectrum","STOPMIX and USQMIX given. Using USQMIX.",0); 
    stopmix.clear();
  };
  if (!sbotmix.exists() && !dsqmix.exists()) {
    message(1,"checkSpectrum","No mixing matrix for sdown sector.",0); 
    //Ensure sensible defaults if no mixing.
    sbotmix.set(1,1, 0.0) ; sbotmix.set(1,2, 1.0);
    sbotmix.set(2,1,-1.0) ; sbotmix.set(2,2, 0.0);
  };
  if (sbotmix.exists() && dsqmix.exists()) {
    message(1,"checkSpectrum","SBOTMIX and DSQMIX given. Using DSQMIX.",0);
    sbotmix.clear();
  };  

  //3) SLHA1 --> SLHA2 interoperability
  //Note: the mass basis is NOT mass-ordered in SLHA1, so be careful!
  //Here, the mass basis is hence by PDG code, not by mass-ordered value.
  if (stopmix.exists() && modsel(6) == 0) {
    //1000002 = ~uL, 1000004 = ~cL, 2000002 = ~uR, 1000004 = ~cR 
    usqmix.set(1,1, 1.0);
    usqmix.set(2,2, 1.0); 
    usqmix.set(4,4, 1.0);
    usqmix.set(5,5, 1.0);
    //Fill (1000006,2000006) sector from stopmix
    usqmix.set(3,3, stopmix(1,1));
    usqmix.set(3,6, stopmix(1,2));
    usqmix.set(6,3, stopmix(2,1));
    usqmix.set(6,6, stopmix(2,2));
  };
  if (sbotmix.exists() && modsel(6) == 0) {
    //1000001 = ~dR, 1000003 = ~sR, 2000001 = ~dL, 1000003 = ~sL 
    dsqmix.set(1,1, 1.0);
    dsqmix.set(2,2, 1.0); 
    dsqmix.set(4,4, 1.0);
    dsqmix.set(5,5, 1.0);
    //Fill (1000005,2000005) sector from sbotmix
    dsqmix.set(3,3, sbotmix(1,1));
    dsqmix.set(3,6, sbotmix(1,2));
    dsqmix.set(6,3, sbotmix(2,1));
    dsqmix.set(6,6, sbotmix(2,2));
  };
  return ifail;
}

void SusyLesHouches::message(int level, string place,string themessage,int line) {
  //Send normal messages and warnings to stdout, errors to stderr.
  ostream* outstream = &cerr;
  if (level <= 1) outstream = &cout;
  if (level == 2) { *outstream<<endl; }
  *outstream << "(SLHA::"+place+") ";
  if (level == 0) *outstream<< ": ";
  if (level == 1) *outstream<< "Warning: "; 
  if (level == 2) { *outstream <<"Error: "; } 
  if (line != 0) *outstream<< "line "<<line<<" - ";
  *outstream << themessage << endl;
  if (level == 2) *outstream <<endl;
  return;
}


