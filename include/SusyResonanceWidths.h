// ResonanceWidths.h is a part of the PYTHIA event generator.
// Copyright (C) 2010 
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file for resonance properties: dynamical widths etc. 
// SusyResonanceWidths: base class for all SUSY resonances.


#ifndef Pythia8_SusyResonanceWidths_H
#define Pythia8_SusyResonanceWidths_H

#include "ResonanceWidths.h"
#include "SusyCouplings.h"

namespace Pythia8 {

class ParticleData;

//==========================================================================

class SUSYResonanceWidths : public ResonanceWidths{

public:
  SUSYResonanceWidths(){}

  //Return particle type
  int typeNeut(int idPDG);
  int typeChar(int idPDG); 

  double lam(double x, double y, double z);

protected:

  virtual bool init(Info* infoPtrIn, Settings* settingsPtrIn,
	    ParticleData* particleDataPtrIn, Couplings* couplingsPtrIn);
  
  //SUSY couplings
  CoupSUSY* coupSUSYPtr;
  static const bool DEBUG;


};

//==========================================================================

// The ResonanceSquark class handles the Squark resonances.

class ResonanceSquark : public SUSYResonanceWidths {

//--------------------------------------------------------------------------
public:

  // Constructor. 
  ResonanceSquark(int idResIn) {initBasic(idResIn);} 


private: 

  // Locally stored properties and couplings.
  double lambda;
 
  // Initialize constants.
  virtual void initConstants(); 
 
  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool calledFromInit = false);

  double s2W;

};
  
//==========================================================================

// The ResonanceGluino class handles the Gluino resonances.

//--------------------------------------------------------------------------

class ResonanceGluino : public SUSYResonanceWidths {

public:

  // Constructor. 
  ResonanceGluino(int idResIn) {initBasic(idResIn);} 


private: 

  // Locally stored properties and couplings.
  double lambda;
 
  // Initialize constants.
  virtual void initConstants(); 
 
  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool calledFromInit = false);

  
};
  
//==========================================================================

// The ResonanceNeut class handles the Neutralino resonances.

//--------------------------------------------------------------------------

class ResonanceNeut : public SUSYResonanceWidths {

public:

  // Constructor. 
  ResonanceNeut(int idResIn) {initBasic(idResIn);} 


private: 

  // Locally stored properties and couplings.
  double lambda;
  double kinFac2;

  // Initialize constants.
  virtual void initConstants(); 
 
  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool calledFromInit = false);

  double s2W;

};
  
//==========================================================================

// The ResonanceChar class handles the Chargino resonances.

//--------------------------------------------------------------------------

class ResonanceChar : public SUSYResonanceWidths {

public:

  // Constructor. 
  ResonanceChar(int idResIn) {initBasic(idResIn);} 


private: 

  // Locally stored properties and couplings.
  double lambda;
  double kinFac2;

  // Initialize constants.
  virtual void initConstants(); 
 
  // Calculate various common prefactors for the current mass.
  virtual void calcPreFac(bool = false);

  // Caclulate width for currently considered channel.
  virtual void calcWidth(bool calledFromInit = false);

  double s2W;

};
  
//==========================================================================

} // end namespace Pythia8

#endif
