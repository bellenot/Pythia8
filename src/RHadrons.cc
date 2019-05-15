// RHadrons.cc is a part of the PYTHIA event generator.
// Copyright (C) 2011 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the RHadrons class.

#include "RHadrons.h"

namespace Pythia8 {
 
//==========================================================================

// The RHadrons class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

const int RHadrons::IDRHADSB[14] = {  1000512, 1000522, 1000532,
  1000542, 1000552, 1005113, 1005211, 1005213, 1005223, 1005311,
  1005313, 1005321, 1005323, 1005333 };

const int RHadrons::IDRHADST[14] = {  1000612, 1000622, 1000632,
  1000642, 1000652, 1006113, 1006211, 1006213, 1006223, 1006311,
  1006313, 1006321, 1006323, 1006333 };

// Gluino code and list of gluino R-hadron codes.
const int RHadrons::IDRHADGO[38] = {  1000993, 1009113, 1009213, 
  1009223, 1009313, 1009323, 1009333, 1009413, 1009423, 1009433, 
  1009443, 1009513, 1009523, 1009533, 1009543, 1009553, 1091114, 
  1092114, 1092214, 1092224, 1093114, 1093214, 1093224, 1093314, 
  1093324, 1093334, 1094114, 1094214, 1094224, 1094314, 1094324, 
  1094334, 1095114, 1095214, 1095224, 1095314, 1095324, 1095334 };

// Allow a few tries for flavour selection since failure is allowed.
const int RHadrons::NTRYMAX = 10;

// Safety margin (in GeV) when constructing kinematics of system.
const double RHadrons::MSAFETY = 0.1;

//--------------------------------------------------------------------------

// Main routine to initialize R-hadron handling.

bool RHadrons::init( Info* infoPtrIn, Settings& settings, 
  ParticleData* particleDataPtrIn, Rndm* rndmPtrIn) {

  // Store input pointers for future use. 
  infoPtr          = infoPtrIn;
  particleDataPtr  = particleDataPtrIn;  
  rndmPtr          = rndmPtrIn;

  // Flags and parameters related to R-hadron formation and decay.
  allowRH          = settings.flag("RHadrons:allow");
  maxWidthRH       = settings.parm("RHadrons:maxWidth");
  idRSb            = settings.mode("RHadrons:idSbottom");
  idRSt            = settings.mode("RHadrons:idStop");
  idRGo            = settings.mode("RHadrons:idGluino");
  setMassesRH      = settings.flag("RHadrons:setMasses");
  probGluinoballRH = settings.parm("RHadrons:probGluinoball");
  mOffsetCloudRH   = settings.parm("RHadrons:mOffsetCloud");
  mCollapseRH      = settings.parm("RHadrons:mCollapse");
  diquarkSpin1RH   = settings.parm("RHadrons:diquarkSpin1"); 

  // Check whether sbottom, stop or gluino should form R-hadrons. 
  allowRSb         = allowRH && idRSb > 0 
    && (particleDataPtr->mWidth(idRSb) < maxWidthRH);
  allowRSt         = allowRH && idRSt > 0 
    && (particleDataPtr->mWidth(idRSt) < maxWidthRH);
  allowRGo         = allowRH && idRGo > 0 
    && (particleDataPtr->mWidth(idRGo) < maxWidthRH);
  allowSomeR       = allowRSb || allowRSt || allowRGo;

  // Set nominal masses of sbottom R-mesons and R-baryons.
  if (allowRSb) {
    m0Sb = particleDataPtr->m0(idRSb);
    if (setMassesRH) {
      for (int i = 0; i < 14; ++i) {
        int idR = IDRHADSB[i]; 
        double m0RHad = m0Sb + mOffsetCloudRH;
        m0RHad += particleDataPtr->constituentMass( (idR%100)/10);
        if (i > 4) 
        m0RHad += particleDataPtr->constituentMass( (idR%1000)/100 );
        particleDataPtr->m0( idR, m0RHad);    
      }
    }

    // Set widths and lifetimes of sbottom R-hadrons.
    double mWidthRHad = particleDataPtr->mWidth(idRSb);
    double tau0RHad   = particleDataPtr->tau0(  idRSb);
    for (int i = 0; i < 14; ++i) {
      particleDataPtr->mWidth( IDRHADSB[i], mWidthRHad);
      particleDataPtr->tau0(   IDRHADSB[i],   tau0RHad);
    }
  }

  // Set nominal masses of stop R-mesons and R-baryons.
  if (allowRSt) {
    m0St = particleDataPtr->m0(idRSt);
    if (setMassesRH) {
      for (int i = 0; i < 14; ++i) {
        int idR = IDRHADST[i]; 
        double m0RHad = m0St + mOffsetCloudRH;
        m0RHad += particleDataPtr->constituentMass( (idR%100)/10);
        if (i > 4) 
        m0RHad += particleDataPtr->constituentMass( (idR%1000)/100 );
        particleDataPtr->m0( idR, m0RHad);    
      }
    }

    // Set widths and lifetimes of stop R-hadrons.
    double mWidthRHad = particleDataPtr->mWidth(idRSt);
    double tau0RHad   = particleDataPtr->tau0(  idRSt);
    for (int i = 0; i < 14; ++i) {
      particleDataPtr->mWidth( IDRHADST[i], mWidthRHad);
      particleDataPtr->tau0(   IDRHADST[i],   tau0RHad);
    }
  }

  // Set nominal masses of gluino R-glueballs, R-mesons and R-baryons.
  if (allowRGo) {
    m0Go = particleDataPtr->m0(idRGo);
    if (setMassesRH) {
      particleDataPtr->m0( IDRHADGO[0], m0Go + 2. * mOffsetCloudRH 
        + particleDataPtr->constituentMass(21) );
      for (int i = 1; i < 38; ++i) {
        int idR = IDRHADGO[i]; 
        double m0RHad = m0Go + 2. * mOffsetCloudRH;
        m0RHad += particleDataPtr->constituentMass( (idR%1000)/100 );
        m0RHad += particleDataPtr->constituentMass( (idR%100)/10);
        if (i > 15) 
        m0RHad += particleDataPtr->constituentMass( (idR%10000)/1000 );
        particleDataPtr->m0( idR, m0RHad);    
      }
    }

    // Set widths and lifetimes of gluino R-hadrons.
    double mWidthRHad = particleDataPtr->mWidth(idRGo);
    double tau0RHad   = particleDataPtr->tau0(  idRGo);
    for (int i = 0; i < 38; ++i) {
      particleDataPtr->mWidth( IDRHADGO[i], mWidthRHad);
      particleDataPtr->tau0(   IDRHADGO[i],   tau0RHad);
    }
  }   

  // Done. 
  return true;

}
//--------------------------------------------------------------------------

// Check if a given particle can produce R-hadrons.

bool RHadrons::givesRHadron( int id) {
  if (allowRSb && abs(id) == idRSb) return true;
  if (allowRSt && abs(id) == idRSt) return true;
  if (allowRGo && id == idRGo) return true;
  return false;

}

//--------------------------------------------------------------------------

// Produce R-hadrons by fragmenting them off from existing strings.

bool RHadrons::produce( ColConfig& colConfig, Event& event) {

  // Check whether some sparticles are allowed to hadronize.
  if (!allowSomeR) return true;

  // Reset arrays and values.
  iBefRHad.resize(0);
  iCreRHad.resize(0);
  iRHadron.resize(0);
  iAftRHad.resize(0);
  isTriplet.resize(0);
  nRHad = 0;

  // Loop over event and identify hadronizing sparticles.
  for (int i = 0; i < event.size(); ++i) 
   if (event[i].isFinal() && givesRHadron(event[i].id())) { 
    iBefRHad.push_back(i);
    iCreRHad.push_back(i);
    iRHadron.push_back(0);
    iAftRHad.push_back(0);
    isTriplet.push_back(true);
  } 
  nRHad = iRHadron.size();
  
  // Done if no hadronizing sparticles.
  if (nRHad == 0) return true;

  // Max two R-hadrons. Randomize order of processing.
  if (nRHad > 2) {
     infoPtr->errorMsg("Error in RHadrons::produce: "
       "cannot handle more than two R-hadrons");
     return false;
  }
  if (nRHad > 1 && rndmPtr->flat() > 0.5) swap( iBefRHad[0], iBefRHad[1]);

  // Split up a colour colour singlet system that should give two R-hadrons.
  // For now don't handle systems involving junctions or loops.
  if (nRHad == 2) {
    int iSys1 = colConfig.findSinglet( iBefRHad[0]);
    int iSys2 = colConfig.findSinglet( iBefRHad[1]);
    if (iSys2 == iSys1) { 
      iSys = iSys1;
      systemPtr = &colConfig[iSys];
      if (systemPtr->hasJunction) {
        infoPtr->errorMsg("Error in RHadrons::produce: "
          "cannot handle system with junction");
        return false;
      }
      if (systemPtr->isClosed) {
        infoPtr->errorMsg("Error in RHadrons::produce: "
          "cannot handle closed colour loop");
        return false;
      }
      if ( !splitSystem( colConfig, event) ) { 
        infoPtr->errorMsg("Error in RHadrons::produce: "
          "failed to handle two sparticles in same system");
        return false;
      }
    } 
  }
    
  // Loop over R-hadrons to be formed. Find its colour singlet system.
  for (iRHad = 0; iRHad < nRHad; ++iRHad) {
    iBef = iBefRHad[iRHad];  
    iSys = colConfig.findSinglet( iBef);
    if (iSys < 0) {
      infoPtr->errorMsg("Error in RHadrons::produce: "
        "sparticle not in any colour singlet");
      return false;
    }
    systemPtr = &colConfig[iSys];

    // For now don't handle systems involving junctions or loops.
    if (systemPtr->hasJunction) {
      infoPtr->errorMsg("Error in RHadrons::produce: "
        "cannot handle system with junction");
      return false;
    }
    if (systemPtr->isClosed) {
      infoPtr->errorMsg("Error in RHadrons::produce: "
        "cannot handle closed colour loop");
      return false;
    }

    // Handle formation of R-hadron separately for gluino and squark.
    if (event[iBef].id() == idRGo) isTriplet[iRHad] = false;
    bool formed = (isTriplet[iRHad]) ? produceSquark( colConfig, event)
                                     : produceGluino( colConfig, event);
    if (!formed) return false;

  // End of loop over R-hadrons. Done.
  } 
  return true;

}

//--------------------------------------------------------------------------

// Decay R-hadrons by resolving them into string systems and letting the
// heavy unstable particle decay as normal.

bool RHadrons::decay( Event& event) {

  // Loop over R-hadrons to decay.
  for (iRHad = 0; iRHad < nRHad; ++iRHad) {
    int    iRNow  = iRHadron[iRHad]; 
    int    iRBef  = iBefRHad[iRHad];
    int    idRHad = event[iRNow].id();
    double mRHad  = event[iRNow].m();
    double mRBef  = event[iRBef].m();
    int    iR0    = 0;
    int    iR2    = 0; 

    // Find flavour content of squark or gluino R-hadron.
    pair<int,int> idPair = (isTriplet[iRHad]) 
      ? fromIdWithSquark( idRHad) : fromIdWithGluino( idRHad);
    int id1 = idPair.first;
    int id2 = idPair.second;

    // Sharing of momentum: ideally the squark/gluino should be restored
    // to original mass, but not at expense of negative-mass spectators.
    double fracR = mRBef / mRHad;

    // Check and, if necessary, fix up mass sharing for squark.
    if (isTriplet[iRHad]) {
      double m2Min = 0.5 * particleDataPtr->constituentMass(id2); 
      if (mRBef + m2Min > mRHad) {
        infoPtr->errorMsg("Warning in RHadrons::decay: "
          "unexpectedly low R-hadron mass (squark case)");
        fracR = mRBef / (mRBef + m2Min); 
      }

      // New colour needed in breakup.
      int colNew = event.nextColTag();
      int col    = (event[iRBef].col() != 0) ? colNew : 0;
      int acol   = (col == 0) ? colNew : 0; 
      
      // Store the constituents of a squark R-hadron.
      iR0 = event.append( id1, 106, iRNow, 0, 0, 0, col, acol,
        fracR * event[iRNow].p(), fracR * mRHad, 0.);
      iR2 = event.append( id2, 106, iRNow, 0, 0, 0, acol, col, 
        (1. - fracR) * event[iRNow].p(), (1. - fracR) * mRHad, 0.);

    // Check and, if necessary, fix up mass sharing for gluino.
    } else {
      double m1Min = 0.5 * particleDataPtr->constituentMass(id1);
      double m2Min = 0.5 * particleDataPtr->constituentMass(id2);
      if (mRBef + m1Min + m2Min > mRHad) {
        infoPtr->errorMsg("Warning in RHadrons::decay: "
          "unexpectedly low R-hadron mass (gluino case)");
        fracR = mRBef / (mRBef + m1Min + m2Min); 
      }
  
      // Set mass sharing between two spectators.
      double m1Eff  = 2. * m1Min + mOffsetCloudRH;  
      double m2Eff  = 2. * m2Min + mOffsetCloudRH;   
      double frac1 = (1. - fracR) * m1Eff / ( m1Eff + m2Eff); 
      double frac2 = (1. - fracR) * m2Eff / ( m1Eff + m2Eff); 
   
      // Two new colours in the breakups.
      int col1 = event.nextColTag();
      int col2 = event.nextColTag();

      // Store the constituents of a gluino R-hadron.
      iR0 = event.append( idRGo, 106, iRNow, 0, 0, 0, col2, col1,
        fracR * event[iRNow].p(), fracR * mRHad, 0.);
      event.append( id1, 106, iRNow, 0, 0, 0, col1, 0, 
        frac1 * event[iRNow].p(), frac1 * mRHad, 0.);
      iR2 = event.append( id2, 106, iRNow, 0, 0, 0, 0, col2, 
        frac2 * event[iRNow].p(), frac2 * mRHad, 0.);
    }

    // Mark R-hadron as decayed and update history.
    event[iRNow].statusNeg();
    event[iRNow].daughters( iR0, iR2);
    iAftRHad[iRHad] = iR0;

  // End loop over R-hadron decays.
  }

  // Done.
  return true;

}

//--------------------------------------------------------------------------

// Split up a colour colour singlet system that should give two R-hadrons.
// To fix : if nLeg >= 3 && mMin large handle as in nLeg == 1??

bool RHadrons::splitSystem( ColConfig& colConfig, Event& event) {

  // First and second R-hadron mother. Number of legs in between.
  int iFirst  = -1;
  int iSecond = -1;
  for (int i = 0; i < int(systemPtr->size()); ++i) {
    int  iTmp   = systemPtr->iParton[i];
    if ( givesRHadron( event[iTmp].id()) ) { 
      if (iFirst == -1) iFirst  = i;
      else              iSecond = i;
    }
  }
  int nLeg = iSecond - iFirst;

  // New flavour pair for breaking the string, and its mass.
  int    idNewQ = flavSelPtr->pickLightQ();
  double mNewQ  = particleDataPtr->constituentMass( idNewQ);
  vector<int> iNewSys1, iNewSys2;

  // If sparticles are neighbours then need new q-qbar pair in between.
  if (nLeg == 1) {

    // The location of the two sparticles.
    int i1Old = systemPtr->iParton[iFirst];
    int i2Old = systemPtr->iParton[iSecond];

    // Take a fraction of momentum to give breakup quark pair.
    double mHat = (event[i1Old].p() + event[i2Old].p()).mCalc();
    double mMax = mHat - event[i1Old].m() - event[i2Old].m(); 
    if (mMax < 2. * (mNewQ + MSAFETY)) return false;
    double mEff = min( 2. * (mNewQ + mOffsetCloudRH), mMax - 2. * MSAFETY);
    double frac = mEff / mHat;
    Vec4   pEff = frac * (event[i1Old].p() + event[i2Old].p());
  
    // New kinematics by (1) go to same mHat with bigger masses, and 
    // (2) rescale back down to original masses and reduced mHat.
    Vec4 p1New, p2New;
    if ( !newKin( event[i1Old].p(), event[i2Old].p(), 
      event[i1Old].m() / (1. - frac), event[i2Old].m() / (1. - frac), 
      p1New, p2New) ) return false; 
    p1New *= 1. - frac;
    p2New *= 1. - frac;

    // Fill in new partons after branching, with correct colour/flavour sign.
    int i1New, i2New, i3New, i4New;
    int newCol = event.nextColTag();
    i1New = event.copy( i1Old, 101);
    if (event[i2Old].acol() == event[i1Old].col()) {
      i3New = event.append( -idNewQ, 101, i1Old, 0, 0, 0, 
        0, event[i2Old].acol(), 0.5 * pEff, 0.5 * mEff, 0.);
      i2New = event.copy( i2Old, 101);
      event[i2New].acol( newCol);
      i4New = event.append(  idNewQ, 101, i2Old, 0, 0, 0, 
        newCol, 0, 0.5 * pEff, 0.5 * mEff, 0.); 
    } else {
      i3New = event.append(  idNewQ, 101, i1Old, 0, 0, 0, 
        event[i2Old].col(), 0, 0.5 * pEff, 0.5 * mEff, 0.);
      i2New = event.copy( i2Old, 101);
      event[i2New].col( newCol);
      i4New = event.append( -idNewQ, 101, i2Old, 0, 0, 0, 
        0, newCol, 0.5 * pEff, 0.5 * mEff, 0.); 
    }

    // Modify momenta and history. For iBotCopyId tracing asymmetric 
    // bookkeeping where one sparticle is radiator and the other recoiler.
    event[i1New].p( p1New);
    event[i2New].p( p2New);
    event[i1Old].daughters( i1New, i3New);
    event[i1New].mother2( 0);
    event[i2Old].daughters( i2New, i4New);
    event[i2New].mother2( 0);
    iBefRHad[0] = i1New;
    iBefRHad[1] = i2New;
 
    // Partons in the two new systems.
    for (int i = 0; i < iFirst; ++i) 
      iNewSys1.push_back( systemPtr->iParton[i]);
    iNewSys1.push_back( i1New);
    iNewSys1.push_back( i3New);
    iNewSys2.push_back( i4New);
    iNewSys2.push_back( i2New);
    for (int i = iSecond + 1; i < int(systemPtr->size()); ++i) 
      iNewSys2.push_back( systemPtr->iParton[i]);

  // If one gluon between sparticles then split it into a
  // collinear q-qbar pair.
  } else if (nLeg == 2) {

    // Gluon to be split and its two daughters.
    int iOld  = systemPtr->iParton[iFirst + 1];
    int i1New = event.append(  idNewQ, 101, iOld, 0, 0, 0, 
      event[iOld].col(), 0, 0.5 * event[iOld].p(), 
      0.5 * event[iOld].m(), 0.);
    int i2New = event.append( -idNewQ, 101, iOld, 0, 0, 0, 
      0, event[iOld].acol(), 0.5 * event[iOld].p(), 
      0.5 * event[iOld].m(), 0.);
    event[iOld].statusNeg();
    event[iOld].daughters( i1New, i2New);
     
    // Partons in the two new systems.
    if (event[systemPtr->iParton[iFirst]].col() == event[i2New].acol()) 
      swap( i1New, i2New);
    for (int i = 0; i <= iFirst; ++i) 
      iNewSys1.push_back( systemPtr->iParton[i]);
    iNewSys1.push_back( i1New);
    iNewSys2.push_back( i2New);
    for (int i = iSecond; i < int(systemPtr->size()); ++i) 
      iNewSys2.push_back( systemPtr->iParton[i]);

  // If several gluons between sparticles then find lowest-mass gluon pair
  // and replace it by a q-qbar pair.
  } else {

    // Find lowest-mass gluon pair and adjust effective quark mass.
    int    iMin  = 0;
    int    i1Old = 0;
    int    i2Old = 0;
    double mMin  = 1e20;
    for (int i = iFirst + 1; i < iSecond - 1; ++i) { 
      int    i1Tmp = systemPtr->iParton[i];
      int    i2Tmp = systemPtr->iParton[i + 1];
      double mTmp  = (event[i1Tmp].p() + event[i2Tmp].p()).mCalc();
      if (mTmp < mMin) {
        iMin  = i;
        i1Old = i1Tmp;
        i2Old = i2Tmp;
        mMin  = mTmp;
      }
    }
    double mEff = min( mNewQ + mOffsetCloudRH, 0.4 * mMin);

    // New kinematics  by sharing mHat appropriately.
    Vec4 p1New, p2New;
    if ( !newKin( event[i1Old].p(), event[i2Old].p(), 
      mEff, mEff, p1New, p2New, false) ) return false; 

    // Insert new partons and update history.
    int i1New, i2New;
    if (event[systemPtr->iParton[0]].acol() == 0) {
      i1New = event.append( -idNewQ, 101, i1Old, 0, 0, 0, 
        0, event[i1Old].acol(), p1New, mEff, 0.);
      i2New = event.append(  idNewQ, 101, i2Old, 0, 0, 0, 
        event[i2Old].col(), 0, p2New, mEff, 0.);
    } else {
      i1New = event.append(  idNewQ, 101, i1Old, 0, 0, 0, 
        event[i1Old].col(), 0, p1New, mEff, 0.);
      i2New = event.append( -idNewQ, 101, i2Old, 0, 0, 0, 
        0, event[i2Old].acol(), p2New, mEff, 0.);
    } 
    event[i1Old].statusNeg();
    event[i2Old].statusNeg();
    event[i1Old].daughters( i1New, 0);
    event[i2Old].daughters( i2New, 0);
     
    // Partons in the two new systems.
    for (int i = 0; i < iMin; ++i) 
      iNewSys1.push_back( systemPtr->iParton[i]);
    iNewSys1.push_back( i1New);
    iNewSys2.push_back( i2New);
    for (int i = iMin + 2; i < int(systemPtr->size()); ++i) 
      iNewSys2.push_back( systemPtr->iParton[i]);
  }

  // Erase the old system and insert the two new ones instead.
  colConfig.erase(iSys);
  colConfig.insert( iNewSys1, event);
  colConfig.insert( iNewSys2, event);   

  // Done. 
  return true;

}

//--------------------------------------------------------------------------

// Produce a R-hadron from a squark and another string end.

bool RHadrons::produceSquark( ColConfig& colConfig, Event& event) {

  // Initial values.
  int    nBody  = 0;
  int    iRNow  = 0;
  int    iNewQ  = 0;
  int    iNewL  = 0;

  // Check at which end of the string the squark is located.
  int    idAbsTop = event[ systemPtr->iParton[0] ].idAbs(); 
  bool   sqAtTop  = (allowRSb && idAbsTop == idRSb) 
                 || (allowRSt && idAbsTop == idRSt);

  // Copy down system consecutively from squark end.
  int    iBeg = event.size();
  iCreRHad[iRHad] = iBeg;
  if (sqAtTop) for (int i = 0; i < systemPtr->size(); ++i) 
    event.copy( systemPtr->iParton[i], 102);
  else         for (int i = systemPtr->size() - 1; i >= 0; --i) 
    event.copy( systemPtr->iParton[i], 102);
  int    iEnd = event.size() - 1; 

  // Input flavours. From now on H = heavy and L = light string ends. 
  int    idOldH = event[iBeg].id(); 
  int    idOldL = event[iEnd].id();

  // Pick new flavour to form R-hadron. For now exclude baryons.
  int    idNewQ = flavSelPtr->pickLightQ();
  if (idOldH > 0) idNewQ = -idNewQ;  
  int    idRHad = toIdWithSquark( idOldH, idNewQ);

  // Target mass of R-hadron and z value of fragmentation function.
  double mRHad  = particleDataPtr->m0(idRHad) + event[iBeg].m() 
                - ( (abs(idOldH) == idRSb) ? m0Sb : m0St );
  double z      = zSelPtr->zFrag( idOldH, idNewQ, mRHad*mRHad);

  // Basic kinematics of string piece where break is to occur.
  Vec4   pOldH  = event[iBeg].p();
  int    iOldL  = iBeg + 1;
  Vec4   pOldL  = event[iOldL].p();
  double mOldL  = event[iOldL].m();
  double mNewH  = mRHad / z;
  double sSys   = (pOldH + pOldL).m2Calc();
  double sRem   = (1. - z) * (sSys - mNewH*mNewH);
  double sMin   = pow2(mOldL + mCollapseRH); 

  // If too low remaining mass in system then add one more parton to it. 
  while ( ( sRem < sMin || sSys < pow2(mNewH + mOldL + MSAFETY) )
    && iOldL < iEnd ) {    
    ++iOldL;
    pOldL      += event[iOldL].p();
    mOldL       = event[iOldL].m();
    sSys        = (pOldH + pOldL).m2Calc();
    sRem        = (1. - z) * (sSys - mNewH*mNewH);
    sMin        = pow2(mOldL + mCollapseRH); 
  }

  // If enough mass then split off R-hadron and reduced system.
  if ( sRem > sMin && sSys > pow2(mNewH + mOldL + MSAFETY) ) {
    Vec4 pNewH, pNewL;
    if ( !newKin( pOldH, pOldL, mNewH, mOldL, pNewH, pNewL) ) {
      infoPtr->errorMsg("Error in RHadrons::produceSquark: "
       "failed to construct kinematics with reduced system");
      return false;
    }

    // Insert R-hadron with new momentum. 
    iRNow       = event.append( idRHad, 104, iBeg, iOldL, 0, 0, 0, 0,
      z * pNewH, mRHad, 0.);
 
    // Reduced system with new string endpoint and modified recoiler. 
    idNewQ      = -idNewQ;
    bool hasCol = (idNewQ > 0 && idNewQ < 10) || idNewQ < -10;
    int  col    = (hasCol) ? event[iOldL].acol() : 0;
    int  acol   = (hasCol) ? 0 : event[iOldL].col(); 
    iNewQ       = event.append( idNewQ, 105, iBeg, iOldL, 0, 0, col, acol,
      (1. - z) * pNewH, (1. - z) * mNewH, 0.);
    iNewL       = event.copy( iOldL, 105);
    event[iNewL].mothers( iBeg, iOldL);
    event[iNewL].p( pNewL);

    // Done with processing of split to R-hadron and reduced system.
    nBody = 3;
  }

  // Two-body final state: form light hadron from endpoint and new flavour.
  if (nBody == 0) {
    FlavContainer flav1( idOldL);
    FlavContainer flav2( -idNewQ);
    int iTry   = 0;
    int idNewL = flavSelPtr->combine( flav1, flav2);
    while (++iTry < NTRYMAX && idNewL == 0) 
      idNewL = flavSelPtr->combine( flav1, flav2);    
    if (idNewL == 0) {
       infoPtr->errorMsg("Error in RHadrons::produceSquark: "
         "cannot form light hadron code");
       return false;
    }
    double mNewL = particleDataPtr->mass( idNewL); 
    
    // Check that energy enough for light hadron and R-hadron.
    if ( sSys > pow2(mRHad + mNewL + MSAFETY) ) { 
      Vec4 pRHad, pNewL;
      if ( !newKin( pOldH, pOldL, mRHad, mNewL, pRHad, pNewL) ) {
        infoPtr->errorMsg("Error in RHadrons::produceSquark: "
         "failed to construct kinematics for two-hadron decay");
        return false;
      }

      // Insert R-hadron and light hadron. 
      iRNow = event.append( idRHad, 104, iBeg, iOldL, 0, 0, 0, 0,
        pRHad, mRHad, 0.);
      event.append( idNewL, 105, iBeg, iOldL, 0, 0, 0, 0, 
        pNewL, mNewL, 0.);
   
      // Done for two-body case.
      nBody = 2;
    }
  }

  // Final case: for very low mass collapse to one single R-hadron.  
  if (nBody == 0) { 
    idRHad = toIdWithSquark( idOldH, idOldL);
    if (idRHad == 0) {
       infoPtr->errorMsg("Error in RHadrons::produceSquark: "
         "cannot form R-hadron code");
       return false;
    }

    // Insert R-hadron with new momentum. 
    iRNow = event.append( idRHad, 104, iBeg, iOldL, 0, 0, 0, 0,
      systemPtr->pSum, systemPtr->mass, 0.);

    // Done with one-body case.
    nBody = 1;
  }
      
  // History bookkeeping: the R-hadron and processed partons. 
  iRHadron[iRHad] = iRNow;
  int iLast = event.size() - 1;
  for (int i = iBeg; i <= iOldL; ++i) {
    event[i].statusNeg(); 
    event[i].daughters( iRNow, iLast);
  }  

  // Remove old string system and, if needed, insert a new one.
  colConfig.erase(iSys);
  if (nBody == 3) {
    vector<int> iNewSys;
    iNewSys.push_back( iNewQ);
    iNewSys.push_back( iNewL);
    for ( int i = iOldL + 1; i <= iEnd; ++i) iNewSys.push_back( i);
    colConfig.insert( iNewSys, event);
  }     
 
  // Done with production of a R-hadron from a squark.  
  return true;

}

//--------------------------------------------------------------------------

// Produce a R-hadron from a gluino and two string ends (or a gluon).

bool RHadrons::produceGluino( ColConfig& colConfig, Event& event) {
         
  // Extract one string system on either side of the gluino.
  int    iGlui  = 0; 
  int    idSave = 0; 
  vector<int> iSide1, iSide2, iSideTmp, iNewSys1, iNewSys2;
  Vec4   pGlui, pSide1, pSide2;
  for (int i = 0; i < systemPtr->size(); ++i) {
    int  iTmp   = systemPtr->iParton[i];
    if (iGlui == 0 && event[ iTmp ].id() == idRGo) {
      iGlui     = iTmp;
      pGlui     = event[ iTmp ].p();
    } else if (iGlui == 0) {
      iSideTmp.push_back( iTmp);
      pSide1   += event[ iTmp ].p();
    } else {
      iSide2.push_back( iTmp);
      pSide2   += event[ iTmp ].p();
    }
  }
 
  // Order sides from gluino outwards and with lowest relative mass first.
  for (int i = int(iSideTmp.size()) - 1; i >= 0; --i) 
    iSide1.push_back( iSideTmp[i]);
  double m1H    = (pGlui + pSide1).mCalc() - event[ iSide1[0] ].m();
  double m2H    = (pGlui + pSide2).mCalc() - event[ iSide2[0] ].m();
  if (m2H < m1H) {
    swap( iSide1, iSide2);
    swap( pSide1, pSide2);
  }

  // Begin loop over two sides of gluino, with lowest relative mass first.
  for (int iSide = 1; iSide < 3; ++iSide) {

    // Begin processing of lower-mass string side.
    int    idRHad = 0;
    double mRHad  = 0.;
    int    nBody  = 0;
    int    iRNow  = 0;
    int    iNewQ  = 0;
    int    iNewL  = 0;
    int    statusRHad = 0;

    // Copy down system consecutively from gluino end.
    int    iBeg   = event.size();
    if (iSide == 1) {
      iCreRHad[iRHad] = iBeg;
      event.copy( iGlui, 102);
      for (int i = 0; i < int(iSide1.size()); ++i) 
        event.copy( iSide1[i], 102);
    } else {
      event.copy( iGlui, 102);
      for (int i = 0; i < int(iSide2.size()); ++i) 
        event.copy( iSide2[i], 102);
    }
    int    iEnd   = event.size() - 1; 

    // Pick new flavour to help form R-hadron. For now exclude baryons.
    int    idOldL = event[iEnd].id();
    int    idNewQ = flavSelPtr->pickLightQ();
    if (event[iEnd].col() == 0) idNewQ = -idNewQ;  
    int    colR   = 0;
    int    acolR  = 0;

    // Target intermediary mass or R-hadron mass.
    if (iSide == 1) {
      idSave      = idNewQ;
      idRHad      = (idNewQ > 0) ? 1009002 : -1009002;
      if (idRHad > 0) colR  = event[iBeg].col();
      if (idRHad < 0) acolR = event[iBeg].acol();
      statusRHad  = 103;
      mRHad       = event[iBeg].m() + mOffsetCloudRH
                  + particleDataPtr->constituentMass( idNewQ);
    } else {
      idRHad      = toIdWithGluino( idSave, idNewQ);
      statusRHad  = 104;
      mRHad       = particleDataPtr->m0(idRHad) + event[iBeg].m() - m0Go;
    }
      
    // z value of fragmentation function.
    double z      = zSelPtr->zFrag( idRGo, idNewQ, mRHad*mRHad);

    // Basic kinematics of string piece where break is to occur.
    Vec4   pOldH  = event[iBeg].p();
    int    iOldL  = iBeg + 1;
    Vec4   pOldL  = event[iOldL].p();
    double mOldL  = event[iOldL].m();
    double mNewH  = mRHad / z;
    double sSys   = (pOldH + pOldL).m2Calc();
    double sRem   = (1. - z) * (sSys - mNewH*mNewH);
    double sMin   = pow2(mOldL + mCollapseRH); 

    // If too low remaining mass in system then add one more parton to it. 
    while ( ( sRem < sMin || sSys < pow2(mNewH + mOldL + MSAFETY) )
      && iOldL < iEnd ) {    
      ++iOldL;
      pOldL      += event[iOldL].p();
      mOldL       = event[iOldL].m();
      sSys        = (pOldH + pOldL).m2Calc();
      sRem        = (1. - z) * (sSys - mNewH*mNewH);
      sMin        = pow2(mOldL + mCollapseRH); 
    }

    // If enough mass then split off R-hadron and reduced system.
    if ( sRem > sMin && sSys > pow2(mNewH + mOldL + MSAFETY) ) {
      Vec4 pNewH, pNewL;
      if ( !newKin( pOldH, pOldL, mNewH, mOldL, pNewH, pNewL) ) {
        infoPtr->errorMsg("Error in RHadrons::produceGluino: "
         "failed to construct kinematics with reduced system");
        return false;
      }

      // Insert intermediary or R-hadron with new momentum, less colour.
      iRNow       = event.append( idRHad, statusRHad, iBeg, iOldL, 
        0, 0, colR, acolR, z * pNewH, mRHad, 0.);
 
      // Reduced system with new string endpoint and modified recoiler. 
      idNewQ      = -idNewQ;
      bool hasCol = (idNewQ > 0 && idNewQ < 10) || idNewQ < -10;
      int  col    = (hasCol) ? event[iOldL].acol() : 0;
      int  acol   = (hasCol) ? 0 : event[iOldL].col(); 
      iNewQ       = event.append( idNewQ, 105, iBeg, iOldL, 0, 0, 
        col, acol, (1. - z) * pNewH, (1. - z) * mNewH, 0.);
      iNewL       = event.copy( iOldL, 105);
      event[iNewL].mothers( iBeg, iOldL);
      event[iNewL].p( pNewL);

      // Collect partons of new string system.
      if (iSide == 1) {
        iNewSys1.push_back( iNewQ);
        iNewSys1.push_back( iNewL);
        for ( int i = iOldL + 1; i <= iEnd; ++i) iNewSys1.push_back( i);
      } else {
        iNewSys2.push_back( iNewQ);
        iNewSys2.push_back( iNewL);
        for ( int i = iOldL + 1; i <= iEnd; ++i) iNewSys2.push_back( i);
      }     

      // Done with processing of split to R-hadron and reduced system.
      nBody = 3;
    }

    // Two-body final state: form light hadron from endpoint and new flavour.
    if (nBody == 0) {
      FlavContainer flav1( idOldL);
      FlavContainer flav2( -idNewQ);
      int iTry   = 0;
      int idNewL = flavSelPtr->combine( flav1, flav2);
      while (++iTry < NTRYMAX && idNewL == 0) 
        idNewL = flavSelPtr->combine( flav1, flav2);    
      if (idNewL == 0) {
         infoPtr->errorMsg("Error in RHadrons::produceGluino: "
           "cannot form light hadron code");
         return false;
      }
      double mNewL = particleDataPtr->mass( idNewL); 
    
      // Check that energy enough for light hadron and R-hadron.
      if ( sSys > pow2(mRHad + mNewL + MSAFETY) ) { 
        Vec4 pRHad, pNewL;
        if ( !newKin( pOldH, pOldL, mRHad, mNewL, pRHad, pNewL) ) {
          infoPtr->errorMsg("Error in RHadrons::produceGluino: "
           "failed to construct kinematics for two-hadron decay");
          return false;
        }

        // Insert intermediary or R-hadron and light hadron. 
        iRNow = event.append( idRHad, statusRHad, iBeg, iOldL, 0, 0,
          colR, acolR, pRHad, mRHad, 0.);
        event.append( idNewL, 105, iBeg, iOldL, 0, 0, 0, 0, 
          pNewL, mNewL, 0.);
   
        // Done for two-body case.
        nBody = 2;
      }
    }

    // Final case: for very low mass collapse to one single R-hadron.  
    if (nBody == 0) { 
      if (iSide == 1) idSave = idOldL;
      else            idRHad = toIdWithGluino( idSave, idOldL);
      if (idRHad == 0) {
         infoPtr->errorMsg("Error in RHadrons::produceGluino: "
           "cannot form R-hadron code");
         return false;
      }

      // Insert R-hadron with new momentum. 
      iRNow = event.append( idRHad, statusRHad, iBeg, iOldL, 0, 0, 
        colR, acolR, pOldH + pOldL, (pOldH + pOldL).mCalc(), 0.);

      // Done with one-body case.
      nBody = 1;
    }
      
    // History bookkeeping: the processed partons. 
    int iLast = event.size() - 1;
    for (int i = iBeg; i <= iOldL; ++i) {
      event[i].statusNeg(); 
      event[i].daughters( iRNow, iLast);
    }  

    // End of loop over two sides of the gluino.
    iGlui   = iRNow;
  }

  // History bookkeeping: insert R-hadron and replace old string system. 
  iRHadron[iRHad] = iGlui;
  colConfig.erase(iSys);
  if (iNewSys1.size() > 0) colConfig.insert( iNewSys1, event);
  if (iNewSys2.size() > 0) colConfig.insert( iNewSys2, event);
 
  // Done with production of a R-hadron from a gluino.  
  return true;

}

//--------------------------------------------------------------------------

// Form a R-hadron code from a squark and a (di)quark code.
// First argument is the (anti)squark, second the (anti)(di)quark.

int RHadrons::toIdWithSquark( int id1, int id2) {

  // Check that physical combination; return 0 if not.
  int id1Abs = abs(id1);
  int id2Abs = abs(id2);
  if (id2Abs < 10 && id1 * id2 > 0) return 0;
  if (id2Abs > 10 && id1 * id2 < 0) return 0;

  // Form R-hadron code. Flip sign for antisquark. 
  bool isSt = (id1Abs == idRSt);
  int idRHad = 1000000;
  if (id2Abs < 10) idRHad += ((isSt) ? 600 : 500) + 10 * id2Abs + 2;
  else idRHad += ((isSt) ? 6000 : 5000) + 10 * (id2Abs/100) + id2Abs%10;
  if (id1 < 0) idRHad = -idRHad;

  // Done.
  return idRHad;

}

//--------------------------------------------------------------------------

// Resolve a R-hadron code into a squark and a (di)quark.
// Return a pair, where first is the squark and the second the (di)quark.

pair<int,int> RHadrons::fromIdWithSquark( int idRHad) {

  // Find squark flavour content. 
  int idLight = (abs(idRHad) - 1000000) / 10;
  int idSq    = (idLight < 100) ? idLight/10 : idLight/100;
  int id1     = (idSq == 6) ? idRSt : idRSb;
  if (idRHad < 0) id1 = -id1;

  // Find light (di)quark flavour content. 
  int id2     =  (idLight < 100) ? idLight%10 : idLight%100;
  if (id2 > 10) id2 = 10 * id2 + abs(idRHad)%10;
  if ((id2 < 10 && idRHad > 0) || (id2 > 10 && idRHad < 0)) id2 = -id2;

  // Done.
  return make_pair( id1, id2);

}   
 
//--------------------------------------------------------------------------

// Form a R-hadron code from two quark/diquark endpoints and a gluino.

int RHadrons::toIdWithGluino( int id1, int id2) {

  // Check that physical combination; return 0 if not.
  int id1Abs = abs(id1);
  int id2Abs = abs(id2);
  int idMax  = max( id1Abs, id2Abs);
  int idMin  = min( id1Abs, id2Abs);
  if (idMin > 10) return 0;
  if (idMax > 10 && id1 * id2 < 0) return 0;
  if (idMax < 10 && id1 * id2 > 0) return 0;

  // Form R-hadron code. Flip sign for antiparticle.
  int idRHad = (idMax > 10) ? 1090004 + 100 * (idMax/100) + 10 * idMin
                            : 1009003 + 100 * idMax + 10 * idMin;
  if (idMax > 10 && id1 < 0) idRHad = -idRHad;
  if (idMax < 10 && idMin != idMax && idMax%2 == 1) {
    if (id1Abs == idMax && id1 > 0) idRHad = -idRHad;
    if (id2Abs == idMax && id2 > 0) idRHad = -idRHad;
  }
  if (idMax < 10 && idMin != idMax && idMax%2 == 0) {
    if (id1Abs == idMax && id1 < 0) idRHad = -idRHad;
    if (id2Abs == idMax && id2 < 0) idRHad = -idRHad;
  }

  // Done.
  return idRHad;

}

//--------------------------------------------------------------------------

// Resolve a R-hadron code into two quark/diquark endpoints (and a gluino).
// Return a pair, where first carries colour and second anticolour.

pair<int,int> RHadrons::fromIdWithGluino( int idRHad) {

  // Find light flavour content of R-hadron.
  int idLight = (abs(idRHad) - 1000000) / 10; 
  int id1, id2, idTmp, idA, idB, idC; 

  // Gluinoballs: split g into d dbar or u ubar.
  if (idLight < 100) {
    id1 = 1.5 + rndmPtr->flat();
    id2 = -id1;

  // Gluino-meson: split into q + qbar.
  } else if (idLight < 1000) {
    id1 = (idLight / 10) % 10;  
    id2 = -(idLight % 10);
    // Flip signs when first quark of down-type.
    if (id1%2 == 1) {
      idTmp = id1;
      id1   = -id2;
      id2   = -idTmp;
    }

  // Gluino-baryon: split to q + qq (diquark). 
  // Pick diquark at random, except if c or b involved.
  } else {
    idA = (idLight / 100) % 10;
    idB = (idLight / 10) % 10;
    idC = idLight % 10;
    double rndmQ = 3. * rndmPtr->flat();
    if (idA > 3) rndmQ = 0.5;
    if (rndmQ < 1.) {
      id1 = idA;
      id2 = 1000 * idB + 100 * idC + 3;
      if (idB != idC && rndmPtr->flat() > diquarkSpin1RH) id2 -= 2; 
    } else if (rndmQ < 2.) {
      id1 = idB;
      id2 = 1000 * idA + 100 * idC + 3;
      if (idA != idC && rndmPtr->flat() > diquarkSpin1RH) id2 -= 2; 
    } else {
      id1 = idC;
      id2 = 1000 * idA + 100 * idB +3;
      if (idA != idB && rndmPtr->flat() > diquarkSpin1RH) id2 -= 2; 
    }
  }

  // Flip signs for anti-R-hadron.
  if (idRHad < 0) {
    idTmp = id1;
    id1   = -id2;
    id2   = -idTmp;
  }

  // Done.
  return make_pair( id1, id2);

}   
 
//--------------------------------------------------------------------------

// Construct modified four-vectors to match modified masses:
// minimal reshuffling of momentum along common axis. 
// Note that last two arguments are used to provide output!

bool RHadrons::newKin( Vec4 pOld1, Vec4 pOld2, double mNew1, double mNew2,
  Vec4& pNew1, Vec4& pNew2, bool checkMargin) {

  // Squared masses in initial and final kinematics.
  double sSum   = (pOld1 + pOld2).m2Calc();
  double sOld1  = pOld1.m2Calc();
  double sOld2  = pOld2.m2Calc();
  double sNew1  = mNew1 * mNew1;
  double sNew2  = mNew2 * mNew2;

  // Check that kinematically possible.
  if (checkMargin && pow2(mNew1 + mNew2 + MSAFETY) > sSum) return false;

  // Transfer coefficients to give four-vectors with the new masses.
  double lamOld = sqrt( pow2(sSum - sOld1 - sOld2) - 4. * sOld1 * sOld2 );
  double lamNew = sqrt( pow2(sSum - sNew1 - sNew2) - 4. * sNew1 * sNew2 );   
  double move1  = (lamNew * (sSum - sOld1 + sOld2) 
                -  lamOld * (sSum - sNew1 + sNew2)) / (2. * sSum * lamOld);
  double move2  = (lamNew * (sSum + sOld1 - sOld2) 
                -  lamOld * (sSum + sNew1 - sNew2)) / (2. * sSum * lamOld);
  
  // Construct final vectors. Done.
  pNew1 = (1. + move1) * pOld1 - move2 * pOld2;
  pNew2 = (1. + move2) * pOld2 - move1 * pOld1;
  return true;

}   
 
//==========================================================================

} // end namespace Pythia8
