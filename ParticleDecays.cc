// Function definitions (not found in the header) for the 
// ParticleDecays class.
// Copyright © 2005 Torbjörn Sjöstrand

#include "ParticleDecays.h"

namespace Pythia8 {

//**************************************************************************

// The ParticleDecays class.

//*********
 
// Definitions of static variables.
// (Values will be overwritten in initStatic call, so are purely dummy.)

double ParticleDecays::mSafety = 0.002;
bool ParticleDecays::limitTau0 = false;
double ParticleDecays::tau0Max = 10.0;
bool ParticleDecays::limitTau = false;
double ParticleDecays::tauMax = 10.0;
bool ParticleDecays::limitRadius = false;
double ParticleDecays::rMax = 10.0;
bool ParticleDecays::limitCylinder = false;
double ParticleDecays::xyMax = 10.0;
double ParticleDecays::zMax = 10.0;
bool ParticleDecays::limitDecay = false;
bool ParticleDecays::mixB = true;
double ParticleDecays::xBdMix = 0.771;
double ParticleDecays::xBsMix = 25.0;
double ParticleDecays::multIncrease = 4.5;
double ParticleDecays::multRefMass = 0.7;
double ParticleDecays::multGoffset = 0.0;
double ParticleDecays::colRearrange = 0.5;
double ParticleDecays::probStoU = 0.3;
double ParticleDecays::probQandS = 2.3;
double ParticleDecays::stopMass = 1.0;
double ParticleDecays::sRhoDal = 0.5;
double ParticleDecays::wRhoDal = 0.02;
bool ParticleDecays::FSRinDecays = true;

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of times one tries to let decay happen (for 2 nested loops).
const int ParticleDecays::NTRYDECAY = 10;

// These numbers are hardwired empirical parameters, 
// intended to speed up the M-generator.
const double ParticleDecays::WTCORRECTION[11] = { 1., 1., 1., 
  2., 5., 15., 60., 250., 1250., 7000., 50000. };

//*********

// Method to pass in pointer and particle list for external decay processing.

bool ParticleDecays::decayPtr( DecayHandler* decayHandlePtrIn, 
  vector<int> handledParticles) {

  // Save pointer.
  if (decayHandlePtrIn == 0) return false; 
  decayHandlePtr = decayHandlePtrIn;
  
  // Set which particles should be handled externally.
  for (int i = 0; i < int(handledParticles.size()); ++i) 
    ParticleDataTable::externalDecay(handledParticles[i], true);

  // Done.
  return true;

}

//*********

// Initialize parameters of particle decays.

void ParticleDecays::initStatic() {

  // Safety margin in mass to avoid troubles.
  mSafety = Settings::parameter("ParticleDecays:mSafety");

  // Lifetime and vertex rules for determining whether decay allowed.
  limitTau0 = Settings::flag("ParticleDecays:limitTau0");
  tau0Max = Settings::parameter("ParticleDecays:tau0Max");
  limitTau = Settings::flag("ParticleDecays:limitTau");
  tauMax = Settings::parameter("ParticleDecays:tauMax");
  limitRadius = Settings::flag("ParticleDecays:limitRadius");
  rMax = Settings::parameter("ParticleDecays:rMax");
  limitCylinder = Settings::flag("ParticleDecays:limitCylinder");
  xyMax = Settings::parameter("ParticleDecays:xyMax");
  zMax = Settings::parameter("ParticleDecays:zMax");
  limitDecay = limitTau0 || limitTau || limitRadius || limitCylinder;

  // B-Bbar mixing parameters.
  mixB = Settings::flag("ParticleDecays:mixB");
  xBdMix = Settings::parameter("ParticleDecays:xBdMix");
  xBsMix = Settings::parameter("ParticleDecays:xBsMix");

  // Selection of multiplicity and colours in "phase space" model.
  multIncrease = Settings::parameter("ParticleDecays:multIncrease");
  multRefMass = Settings::parameter("ParticleDecays:multRefMass");
  multGoffset = Settings::parameter("ParticleDecays:multGoffset");
  colRearrange = Settings::parameter("ParticleDecays:colRearrange");

  // Flavour parameter from StringFlav and StringFragmentation.
  probStoU = Settings::parameter("StringFlav:probStoU");
  probQandS = 2. + probStoU;
  stopMass = Settings::parameter("StringFragmentation:stopMass");

  // Parameters for Dalitz decay virtual gamma mass spectrum.
  sRhoDal = pow2(ParticleDataTable::m0(113)); 
  wRhoDal = pow2(ParticleDataTable::width(113));  

  // Allow showers in decays to qqbar/gg/ggg/gammagg.
  FSRinDecays = Settings::mode("ParticleDecays:FSRinDecays"); 

}

//*********

// Decay a particle; main method.

bool ParticleDecays::decay( int iDec, Event& event) {

  // Check whether a decay is allowed, given the upcoming decay vertex.
  Particle& decayer = event[iDec];
  if (limitDecay && !checkVertex(decayer)) return true; 

  // Fill the decaying particle in slot 0 of arrays.  
  int idDec = decayer.id();
  iProd.resize(0);
  idProd.resize(0);
  mProd.resize(0);
  iProd.push_back( iDec );
  idProd.push_back( idDec );
  mProd.push_back( decayer.m() );
  bool foundChannel = false;

  // Check for oscillations B0 <-> B0bar or B_s0 <-> B_s0bar.
  bool hasOscillated = (abs(idDec) == 511 || abs(idDec) == 531) 
    ? oscillateB(decayer) : false;
  if (hasOscillated) {idDec = - idDec; idProd[0] = idDec;} 

  // Particle data for decaying particle.
  ParticleDataEntry& particleData = decayer.particleData();

  // Optionally send on to external decay program.
  bool doneExternally = false;
  if (particleData.externalDecay()) {
    pProd.resize(0);
    pProd.push_back(decayer.p());
    doneExternally = decayHandlePtr->decay(idProd, mProd, pProd, 
      iDec, event);

    // If it worked, then store the decay products in the event record.
    if (doneExternally) {
      mult = idProd.size() - 1;
      int status = (hasOscillated) ? 94 : 93;
      for (int i = 1; i <= mult; ++i) {
        int iPos = event.append( idProd[i], status, iDec, 0, 0, 0, 
        0, 0, pProd[i], mProd[i]); 
        iProd.push_back( iPos);
      }

      // Also mark mother decayed and store daughters.
      event[iDec].statusNeg(); 
      event[iDec].daughters( iProd[1], iProd[mult]);
    }
  }
    
  // Now begin normal internal decay treatment.
  if (!doneExternally) {

    // Pick a decay channel; allow up to ten tries.
    for (int iTryChannel = 0; iTryChannel < NTRYDECAY; ++iTryChannel) {
      DecayChannel& channel = particleData.decay.pick();
      mode = channel.modeME();
      mult = channel.multiplicity();

      // Allow up to ten tries for each channel (e.g with different masses).
      for (int iTryMode = 0; iTryMode < NTRYDECAY; ++iTryMode) {
        idProd.resize(1);
        mProd.resize(1);
      
        // Extract and store the decay products.
        bool toPartons = false;
        for (int i = 0; i < mult; ++i) {
          int idNow = channel.product(i);
          int idAbs = abs(idNow);
          if (idAbs < 10 || idAbs == 21 || idAbs == 81 || idAbs == 82) 
            toPartons = true;
          if (idDec < 0 && ParticleDataTable::hasAnti(idNow)) idNow = - idNow;
          double mNow = ParticleDataTable::mass(idNow);
          idProd.push_back( idNow);
          mProd.push_back( mNow);
        }  

        // Decays into partons usually translate into hadrons.
        if (toPartons && (mode < 31 || mode > 40) 
          && !pickHadrons(event) ) continue;

        // Need to set colour flow if explicit decay to partons.
        cols.resize(0);
        acols.resize(0);
        for (int i = 0; i <= mult; ++i) {
          cols.push_back(0);
          acols.push_back(0);
        }
        if (toPartons && (mode > 30 && mode < 41) 
          && !setColours(event) ) continue; 

        // Check that enough phase space for decay.
        if (mult > 1) {
          double mDiff = mProd[0];
          for (int i = 1; i <= mult; ++i) mDiff -= mProd[i];
          if (mDiff < mSafety) continue;
        }
  
        // End of two trial loops. Check if succeeded or not.
        foundChannel = true;
        break;
      }
      if (foundChannel) break;
    }
    if (!foundChannel) return false;
    
    // Store decay products in the event record.
    int status = (hasOscillated) ? 92 : 91;
    for (int i = 1; i <= mult; ++i) {
      int iPos = event.append( idProd[i], status, iDec, 0, 0, 0, 
        cols[i], acols[i], Vec4(0., 0., 0., 0.), mProd[i]); 
      iProd.push_back( iPos);
    }
    
    // Do a decay, split by multiplicity.
    bool decayed = false;
    if (mult == 1) decayed = oneBody(event);
    else if (mult == 2) decayed = twoBody(event);
    else if (mult == 3) decayed = threeBody(event);
    else decayed = mGenerator(event);

    // If the decay worked, then mark mother decayed and store daughters.
    if (decayed) {
      event[iDec].statusNeg(); 
      event[iDec].daughters( iProd[1], iProd[mult]);
  
    // Else remove unused daughters and return failure.
    } else {
      event.popBack(mult);
      return false;
    }

  // Now finished normal internal decay treatment. 
  }

  // Set decay vertex when this is displaced.
  if (event[iDec].hasVertex() || event[iDec].tau() > 0.) {
    Vec4 vDec = event[iDec].vDec();
    for (int i = 1; i <= mult; ++i) event[iProd[i]].vProd( vDec );
  }

  // Set lifetime of hadrons.
  for (int i = 1; i <= mult; ++i) 
    event[iProd[i]].tau( event[iProd[i]].tau0() * Rndm::exp() );
  
  // In a decay explicitly to partons then optionally do a shower,
  // and always flag that partonic system should be fragmented. 
  if ( (mode > 30 && mode < 41) && FSRinDecays ) 
    times.shower( event, iProd[1], iProd.back(), mProd[0]);
  moreHadronization = (mode > 30 && mode < 41) ? true : false;

  // Done.
  return true;

}

//*********

// Check whether a decay is allowed, given the upcoming decay vertex.

bool ParticleDecays::checkVertex(Particle& decayer) {

  // Check whether any of the conditions are not fulfilled.
  if (limitTau0 && decayer.tau0() > tau0Max) return false;
  if (limitTau && decayer.tau() > tauMax) return false;
  if (limitRadius && pow2(decayer.xDec()) + pow2(decayer.yDec())
    + pow2(decayer.zDec()) > pow2(rMax)) return false;
  if (limitCylinder && (pow2(decayer.xDec()) + pow2(decayer.yDec())
    > pow2(xyMax) || abs(decayer.zDec()) > zMax) ) return false;

  // Done.
  return true;

}

//*********

// Check for oscillations B0 <-> B0bar or B_s0 <-> B_s0bar.

bool ParticleDecays::oscillateB(Particle& decayer) {

  // Extract relevant information and decide.
  double xBmix = (abs(decayer.id()) == 511) ? xBdMix : xBsMix;
  double tau = decayer.tau();
  double tau0 = decayer.tau0();
  double probosc = pow2(sin(0.5 * xBmix * tau / tau0));
  return (probosc > Rndm::flat()) ? true : false;  

}

//*********

// Do a one-body decay. (Rare; e.g. for K0 -> K0_short.)

bool ParticleDecays::oneBody(Event& event) {

  // References to the particles involved.
  Particle& decayer = event[iProd[0]];
  Particle& prod = event[iProd[1]];
   
  // Set momentum and expand mother information.
  prod.p( decayer.p() );
  prod.m( decayer.m() );
  prod.mother2( iProd[0] );

  // Done.
  return true;

}

//*********

// Do a two-body decay.

bool ParticleDecays::twoBody(Event& event) {

  // References to the particle involved.
  Particle& decayer = event[iProd[0]];
  Particle& prod1 = event[iProd[1]]; 
  Particle& prod2 = event[iProd[2]]; 

  // Masses. 
  double m0 = mProd[0];
  double m1 = mProd[1];    
  double m2 = mProd[2];    

  // Eenergies and absolute momentum in the rest frame.
  if (m1 + m2 + mSafety > m0) return false;
  double e1 = 0.5 * (m0*m0 + m1*m1 - m2*m2) / m0;
  double e2 = 0.5 * (m0*m0 + m2*m2 - m1*m1) / m0;
  double pAbs = 0.5 * sqrtpos( (m0 - m1 - m2) * (m0 + m1 + m2)
     * (m0 + m1 - m2) * (m0 - m1 + m2) ) / m0;  

  // When mode = 3, for V -> PS2 + PS3 (V = vector, pseudoscalar),
  // need to check if production is PS0 -> PS1/gamma + V.
  int iMother = event[iProd[0]].mother1();
  int idSister = 0;
  if (mode == 3) {
    if (iMother <= 0 || iMother >= iProd[0]) mode = 0;
    else { 
      int iDaughter1 = event[iMother].daughter1();
      int iDaughter2 = event[iMother].daughter2();
      if (iDaughter2 != iDaughter1 + 1) mode = 0;
      else {
        int idMother = abs( event[iMother].id() );
        if (idMother <= 100 || idMother%10 !=1 
          || (idMother/1000)%10 != 0) mode = 0; 
        else {
          int iSister = (iProd[0] == iDaughter1) ? iDaughter2 : iDaughter1;
          idSister = abs( event[iSister].id() );
          if ( (idSister <= 100 || idSister%10 !=1 
            || (idSister/1000)%10 != 0) && idSister != 22) mode = 0; 
        } 
      } 
    } 
  }

  // Begin loop over matrix element corrections.
  double wtME, wtMEmax;
  do {
    wtME = 1.;
    wtMEmax = 1.;

    // Isotropic angles give three-momentum.
    double cosTheta = 2. * Rndm::flat() - 1.;
    double sinTheta = sqrt(1. - cosTheta*cosTheta);
    double phi = 2. * M_PI * Rndm::flat();
    double pX = pAbs * sinTheta * cos(phi);  
    double pY = pAbs * sinTheta * sin(phi);  
    double pZ = pAbs * cosTheta;  

    // Fill four-momenta and boost them away from mother rest frame.
    prod1.p( pX, pY, pZ, e1);
    prod2.p( -pX, -pY, -pZ, e2);
    prod1.bst( decayer.p() );
    prod2.bst( decayer.p() );

    // Matrix element for PS0 -> PS1 + V1 -> PS1 + PS2 + PS3 of form 
    // cos**2(theta02) in V1 rest frame, and for PS0 -> gamma + V1 
    // -> gamma + PS2 + PS3 of form sin**2(theta02).
    if (mode == 3) {
      double p10 = decayer.p() * event[iMother].p();
      double p12 = decayer.p() * prod1.p();
      double p02 = event[iMother].p() * prod1.p();
      double s0 = pow2(event[iMother].m());
      double s1 = pow2(decayer.m());
      double s2 =  pow2(prod1.m());
      if (idSister != 22) wtME = pow2(p10 * p12 - s1 * p02);
      else wtME = s1 * (2. * p10 * p12 * p02 - s1 * p02*p02 
        - s0 * p12*p12 - s2 * p10*p10 + s1 * s0 * s2);
      wtME = max( wtME, 1e-6 * s1*s1 * s0 * s2);
      wtMEmax = (p10*p10 - s1 * s0) * (p12*p12 - s1 * s2);
    } 

  // If rejected, try again with new invariant masses.
  } while ( wtME < Rndm::flat() * wtMEmax ); 

  // Done.
  return true;

}

//*********

// Do a three-body decay.

bool ParticleDecays::threeBody(Event& event) {

  // References to the particle involved.
  Particle& decayer = event[iProd[0]];
  Particle& prod1 = event[iProd[1]]; 
  Particle& prod2 = event[iProd[2]]; 
  Particle& prod3 = event[iProd[3]]; 

  // Mother and sum daughter masses. Fail if too close. 
  double m0 = mProd[0];
  double m1 = mProd[1];    
  double m2 = mProd[2];    
  double m3 = mProd[3]; 
  double mSum = m1 + m2 + m3;
  double mDiff = m0 - mSum;   
  if (mDiff < mSafety) return false; 

  // Kinematical limits for 2+3 mass. Variables inside loops.
  double m23Min = m2 + m3;
  double m23Max = m0 - m1;
  double wtPS, wtME, wtMEmax, wtDal, p1Abs, p23Abs;
  double wtPSmax = 1.;
  double m23 = 0.;
  double s23 = 0.; 
  double sMinDal = 1.000001 * m23Min*m23Min; 
  double sMaxDal = m0*m0;

  // Calculate the maximum phase space weight for normal decays.
  if (mode != 2) {
    double p1Max = 0.5 * sqrtpos( (m0 - m1 - m23Min) * (m0 + m1 + m23Min)
      * (m0 + m1 - m23Min) * (m0 - m1 + m23Min) ) / m0; 
    double p23Max = 0.5 * sqrtpos( (m23Max - m2 - m3) * (m23Max + m2 + m3)
      * (m23Max + m2 - m3) * (m23Max - m2 + m3) ) / m23Max;
    wtPSmax = 0.5 * p1Max * p23Max;

  // Alternatively select virtual gamma mass for Dalitz decay.
  } else { 
    do {
      s23 = sMinDal * pow( sMaxDal / sMinDal, Rndm::flat() );
      wtDal = (1. + 0.5 * sMinDal / s23) *  sqrt(1. - sMinDal / s23) 
        * pow3(1. - s23 / sMaxDal) * sRhoDal * (sRhoDal + wRhoDal) 
        / ( pow2(s23 - sRhoDal) + sRhoDal * wRhoDal ); 
    } while ( wtDal < Rndm::flat() );   
    m23 = sqrt(s23);
  }

  // Begin loop over matrix element corrections.
  do {
    wtME = 1.;
    wtMEmax = 1.;

    // Begin loop to find the intermediate m23 invariant masses.
    do {
      wtPS = 1.;

      // Pick an intermediate mass flat in the allowed range.
      // (Except for Dalitz, where a mass has already been picked.)
      if (mode != 2) m23 = m23Min + Rndm::flat() * mDiff;

      // Translate into relative momenta and find weight.
      p1Abs = 0.5 * sqrtpos( (m0 - m1 - m23) * (m0 + m1 + m23)
        * (m0 + m1 - m23) * (m0 - m1 + m23) ) / m0; 
      p23Abs = 0.5 * sqrtpos( (m23 - m2 - m3) * (m23 + m2 + m3)
        * (m23 + m2 - m3) * (m23 - m2 + m3) ) / m23;
      wtPS = p1Abs * p23Abs;

    // If rejected, try again with new invariant masses.
    } while ( mode !=2 && wtPS < Rndm::flat() * wtPSmax ); 

    // Set up m23 -> m2 + m3 isotropic in its rest frame.
    double cosTheta = 2. * Rndm::flat() - 1.;
    double sinTheta = sqrt(1. - cosTheta*cosTheta);
    double phi = 2. * M_PI * Rndm::flat();
    double pX = p23Abs * sinTheta * cos(phi);  
    double pY = p23Abs * sinTheta * sin(phi);  
    double pZ = p23Abs * cosTheta;  
    double e2 = sqrt( m2*m2 + p23Abs*p23Abs);
    double e3 = sqrt( m3*m3 + p23Abs*p23Abs);
    prod2.p( pX, pY, pZ, e2);
    prod3.p( -pX, -pY, -pZ, e3);

    // Set up m0 -> m1 + m23 isotropic in its rest frame.
    cosTheta = 2. * Rndm::flat() - 1.;
    sinTheta = sqrt(1. - cosTheta*cosTheta);
    phi = 2. * M_PI * Rndm::flat();
    pX = p1Abs * sinTheta * cos(phi);  
    pY = p1Abs * sinTheta * sin(phi);  
    pZ = p1Abs * cosTheta;  
    double e1 = sqrt( m1*m1 + p1Abs*p1Abs);
    double e23 = sqrt( m23*m23 + p1Abs*p1Abs);
    prod1.p( pX, pY, pZ, e1);

    // Boost 2 + 3 to the m0 rest frame.
    Vec4 p23( -pX, -pY, -pZ, e23);
    prod2.bst( p23 );
    prod3.bst( p23 );

    // Matrix element weight for omega/phi -> pi+ pi- pi0.
    if (mode == 1) {
      double p12 = prod1.p() * prod2.p(); 
      double p13 = prod1.p() * prod3.p(); 
      double p23 = prod2.p() * prod3.p(); 
      wtME = pow2(m1 * m2 * m3) - pow2(m1 * p23) - pow2(m2 * p13) 
        - pow2(m3 * p12) + 2. * p12 * p13 * p23;
      wtMEmax = pow3(m0 * m0) / 150.;

    // Matrix element weight for Dalitz decay pi0/eta -> gamma e+ e-..
    } else if (mode == 2) {
      double p12 = prod1.p() * prod2.p(); 
      double p13 = prod1.p() * prod3.p(); 
      wtME = (s23 - 0.5 * sMinDal) * (p12*p12 + p13*p13) 
        + sMinDal * (p12*p12 + p13*p13 + p12 * p13);
      wtMEmax = 0.25 * s23 * pow2(sMaxDal - s23);

    // Matrix element weight for "onium" -> g + g + g or gamma + g + g.
    } else if (mode == 33) {
      double x1 = 2. * prod1.e() / m0;
      double x2 = 2. * prod2.e() / m0;
      double x3 = 2. * prod3.e() / m0;
      wtME = pow2( (1. - x1) / (x2 * x3) ) + pow2( (1. - x2) / (x1 * x3) )
        + pow2( (1. - x3) / (x1 * x2) );
      wtMEmax = 2.;
      // For gamma + g + g require minimum mass for g + g system.
      if (prod1.id() == 22 && sqrt(1. - x1) * m0 < 2. * stopMass) wtME = 0.;
      if (prod2.id() == 22 && sqrt(1. - x2) * m0 < 2. * stopMass) wtME = 0.;
      if (prod3.id() == 22 && sqrt(1. - x3) * m0 < 2. * stopMass) wtME = 0.;

    // Effective matrix element for nu spectrum in tau -> nu + hadrons.
    } else if (mode == 41) {
      double x1 = 2. *  prod1.e() / m0;
      wtME = x1 * (3. - 2. * x1);
      double xMax = min( 0.75, 2. * (1. - mSum / m0) ); 
      wtMEmax = xMax * (3. - 2. * xMax); 
    } 

  // If rejected, try again with new invariant masses.
  } while ( wtME < Rndm::flat() * wtMEmax ); 

  // Boost 1 + 2 + 3 to the current frame. 
  prod1.bst( decayer.p() ); 
  prod2.bst( decayer.p() ); 
  prod3.bst( decayer.p() ); 

  // Done.
  return true;

}

//*********

// Do a multibody decay using the M-generator algorithm.

bool ParticleDecays::mGenerator(Event& event) {

  // Mother and sum daughter masses. Fail if too close.
  double m0 = mProd[0];
  double mSum = mProd[1];
  for (int i = 2; i <= mult; ++i) mSum += mProd[i]; 
  double mDiff = m0 - mSum;
  if (mDiff < mSafety) return false; 
   
  // Begin setup of intermediate invariant masses.
  mInv.resize(0);
  for (int i = 0; i <= mult; ++i) mInv.push_back( mProd[i]);

  // Calculate the maximum weight in the decay.
  double wtPS, wtME, wtMEmax;
  double wtPSmax = 1. / WTCORRECTION[mult];
  double mMax = mDiff + mProd[mult];
  double mMin = 0.; 
  for (int i = mult - 1; i > 0; --i) {
    mMax += mProd[i];
    mMin += mProd[i+1];
    double mNow = mProd[i];
    wtPSmax *= 0.5 * sqrtpos( (mMax - mMin - mNow) * (mMax + mMin + mNow)
    * (mMax + mMin - mNow) * (mMax - mMin + mNow) ) / mMax;  
  }

  // Begin loop over matrix element corrections.
  do {
    wtME = 1.;
    wtMEmax = 1.;

    // Begin loop to find the set of intermediate invariant masses.
    do {
      wtPS = 1.;

      // Find and order random numbers in descending order.
      rndmOrd.resize(0);
      rndmOrd.push_back(1.);
      for (int i = 1; i < mult - 1; ++i) { 
        double rndm = Rndm::flat();
        rndmOrd.push_back(rndm);
        for (int j = i - 1; j > 0; --j) {
          if (rndm > rndmOrd[j]) swap( rndmOrd[j], rndmOrd[j+1] );
          else break;
        } 
      }
      rndmOrd.push_back(0.);
  
      // Translate into intermediate masses and find weight.
      for (int i = mult - 1; i > 0; --i) {
        mInv[i] = mInv[i+1] + mProd[i] + (rndmOrd[i-1] - rndmOrd[i]) * mDiff; 
        wtPS *= 0.5 * sqrtpos( (mInv[i] - mInv[i+1] - mProd[i]) 
          * (mInv[i] + mInv[i+1] + mProd[i]) * (mInv[i] + mInv[i+1] - mProd[i]) 
          * (mInv[i] - mInv[i+1] + mProd[i]) ) / mInv[i];  
      }

    // If rejected, try again with new invariant masses.
    } while ( wtPS < Rndm::flat() * wtPSmax ); 

    // Perform two-particle decays in the respective rest frame.
    pInv.resize(mult + 1);
    for (int i = 1; i < mult; ++i) {
      double pAbs = 0.5 * sqrtpos( (mInv[i] - mInv[i+1] - mProd[i]) 
        * (mInv[i] + mInv[i+1] + mProd[i]) * (mInv[i] + mInv[i+1] - mProd[i])
        * (mInv[i] - mInv[i+1] + mProd[i]) ) / mInv[i]; 

      // Isotropic angles give three-momentum.
      double cosTheta = 2. * Rndm::flat() - 1.;
      double sinTheta = sqrt(1. - cosTheta*cosTheta);
      double phi = 2. * M_PI * Rndm::flat();
      double pX = pAbs * sinTheta * cos(phi);  
      double pY = pAbs * sinTheta * sin(phi);  
      double pZ = pAbs * cosTheta;  

      // Calculate energies, fill four-momentum.
      double eHad = sqrt( mProd[i]*mProd[i] + pAbs*pAbs);
      double eInv = sqrt( mInv[i+1]*mInv[i+1] + pAbs*pAbs);
      event[iProd[i]].p( pX, pY, pZ, eHad);
      pInv[i+1].p( -pX, -pY, -pZ, eInv);
    }       
  
    // Boost decay products to the mother rest frame.
    event[iProd[mult]].p( pInv[mult] );
    for (int iFrame = mult - 1; iFrame > 1; --iFrame) 
      for (int i = iFrame; i <= mult; ++i) event[iProd[i]].bst(pInv[iFrame]);

    // Effective matrix element for nu spectrum in tau -> nu + hadrons.
    if (mode == 41) {
      double x1 = 2. * event[iProd[1]].e() / m0;
      wtME = x1 * (3. - 2. * x1);
      double xMax = min( 0.75, 2. * (1. - mSum / m0) ); 
      wtMEmax = xMax * (3. - 2. * xMax); 
    } 

  // If rejected, try again with new invariant masses.
  } while ( wtME < Rndm::flat() * wtMEmax ); 

  // Boost decay products to the current frame. 
  pInv[1].p( event[iProd[0]].p() );
  for (int i = 1; i <= mult; ++i) event[iProd[i]].bst(pInv[1]);

  // Done.
  return true;

}

//*********

// Translate a partonic content into a set of actual hadrons.

bool ParticleDecays::pickHadrons(Event& event) {

  // Find partonic decay products, assumed at the end of the list.
  int nPartons = 0;
  idPartons.resize(0);
  for (int i = 1; i <= mult; ++i) {
    int idAbs = abs(idProd[i]);
    if ( idAbs < 9 || (idAbs > 1000 && idAbs < 10000 && (idAbs/10)%10 == 0)
      || idAbs == 81 || idAbs == 82) {
      ++nPartons;
      idPartons.push_back(idProd[i]);   
    } else if (nPartons > 0) {
      ErrorMessages::message("Error in ParticleDecays::pickHadrons: "
        "hadrons and partons out of order"); 
      return false;
    } 
  }

  // Number of known decay products.
  int nKnown = mult - nPartons;
  int nFilled = nKnown + 1;

  // Replace generic spectator flavour code by the actual one.
  for (int i = 0; i < nPartons; ++i) {
    int idProd = idPartons[i];
    int idNew = idProd;
    if (idProd == 81) { 
      int idDec = event[iProd[0]].id();
      int idAbs = abs(idDec);
      // Is this right sign for B mesons, recall 521 = bbar u ! ??
      if ( (idAbs/1000)%10 == 0 ) idNew = - (idAbs/10)%10;  
      else if ( (idAbs/100)%10 >= (idAbs/10)%10 ) 
        idNew = 100 * ((idAbs/10)%100) + 3;
      else idNew = 1000 * ((idAbs/10)%10) + 100 * ((idAbs/100)%10) + 1;
      if (idDec < 0) idNew = -idNew;

    // Replace generic random flavour by a randomly selected one.
    } else if (idProd == 82) {
      double mFlav;
      do {
        double rndmFlav = probQandS * Rndm::flat();
        int idDummy = -1;
        if (rndmFlav > 1.) idDummy = -2;
        if (rndmFlav > 2.) idDummy = -3;
        do idNew = flavSel.pick(idDummy); 
        while (idNew == 0);  
        mFlav = ParticleDataTable::constituentMass(idNew);
      } while (2. * mFlav + stopMass > mProd[0]);
    } else if (idProd == -82) {
      idNew = -idPartons[i-1];
    } 
    idPartons[i] = idNew;
  }

  // Determine whether fixed multiplicity or to be selected at random.
  int nMin = (mode > 10 && mode <= 15) ? max( 2, mode - 10) : 2;
  if (mode == 43) nMin = 3;
  int nFix = (mode > 20 && mode <= 30) ? max( 2, mode - 20) : 0; 

  // Now begin loop to set new hadronic content.
  int nTot, nNew, nSpec, nLeft;
  double mDiff;
  bool diquarkClash = false;
  do {
    nFilled = nKnown + 1;
    idProd.resize(nFilled);
    mProd.resize(nFilled);      
    nTot = nKnown;
    nNew = 0;
    nSpec = 0;
    nLeft = nPartons;
    mDiff = mProd[0]; 
    for (int i = 1; i < nFilled; ++i) mDiff -= mProd[i];
    diquarkClash = false;

    // For weak decays collapse spectators to one particle.
    if ( (mode == 42 || mode == 43) && nLeft > 1) {
      int id1 = idPartons[nPartons - 2];
      int id2 = idPartons[nPartons - 1];
      int idHad; 
      do idHad = flavSel.combine( id1, id2); 
      while (idHad == 0);
      double mHad = ParticleDataTable::mass(idHad);
      mDiff -= mHad;
      idProd.push_back( idHad);
      mProd.push_back( mHad);
      ++nFilled;
      nSpec = 1;
      nLeft -= 2;
    } 

    // Find mass excess and from there expected multiplicity.
    if (nLeft > 0) {
      double mDiffPS = mDiff;
      for (int i = 0; i <= nLeft; ++i) 
        mDiffPS -= ParticleDataTable::constituentMass( idPartons[i] );
      double coeff = multIncrease * log( max( 1.1, mDiffPS / multRefMass ) );
      if (mode == 12) coeff += multGoffset;
      if (nFix == 0) {
        do { 
          double gauss = sqrt(coeff) * Rndm::gauss();
          nTot = int( 0.5 + 0.5 * (nKnown + nSpec) + 0.25 * nPartons 
            + coeff + gauss); 
	  nNew = nTot - nKnown - nSpec; 
        } while (nNew < nLeft/2 || nTot < nMin || nTot > 10); 

      // Alternatively predetermined multiplicity.
      } else {
        nTot = nFix;
        nNew = nTot - nKnown - nSpec;
        if (nNew < nLeft/2) return false;
      } 

      // Set up ends of fragmented system, as copy of idPartons.
      idEnds.resize(0);
      for (int i = 0; i < nLeft; ++i) idEnds.push_back( idPartons[i] );
    
      // Fragment off at random, but save nLeft/2 for final recombination.
      if (nNew > nLeft/2) {
        for (int i = 0; i < nNew - nLeft/2; ++i) {
          // When four quarks consider last one to be spectator.
          int iEnd = int( (nLeft - 1.) * Rndm::flat() );
          // Pick new flavour and form a new hadron.
          int idNew, idHad;
          do {
            idNew = flavSel.pick( idEnds[iEnd] );
            idHad = flavSel.combine( idEnds[iEnd], idNew);
          } while (idHad == 0);
          // Store new hadron and endpoint flavour.
          idProd.push_back( idHad);  
          idEnds[iEnd] = - idNew;
	}
      }
      
      // When only two quarks left, combine to form final hadron.
      if (nLeft == 2) {
        int idHad;
        if ( abs(idEnds[0]) > 8 && abs(idEnds[1]) > 8) 
          diquarkClash = true; 
	else { 
          do idHad = flavSel.combine( idEnds[0], idEnds[1]);
          while (idHad == 0);
          idProd.push_back( idHad); 
	} 

      // If four quarks, decide how to pair them up.
      } else {
        int iEnd1 = 0;
        int iEnd2 = 1;
        int iEnd3 = 2;
        int iEnd4 = 3;
        if ( Rndm::flat() < colRearrange) iEnd2 = 3;
        int relColSign = ( (idEnds[iEnd1] > 0 && idEnds[iEnd1] < 9)
          || idEnds[iEnd1] < -10 ) ? 1 : -1;
        if ( (idEnds[iEnd2] < 0 && idEnds[iEnd2] > -9)
          || idEnds[iEnd2] > 10 ) relColSign *= -1;
        if (relColSign == 1) iEnd2 = 2;
        if (iEnd2 == 2) iEnd3 = 1;
        if (iEnd2 == 3) iEnd4 = 1; 
        
        // Then combine to get final two hadrons.
        int idHad;
        if ( abs(idEnds[iEnd1]) > 8 && abs(idEnds[iEnd2]) > 8) 
          diquarkClash = true; 
	else { 
          do idHad = flavSel.combine( idEnds[iEnd1], idEnds[iEnd2]);
          while (idHad == 0);
          idProd.push_back( idHad);
	}  
        if ( abs(idEnds[iEnd3]) > 8 && abs(idEnds[iEnd4]) > 8) 
          diquarkClash = true; 
	else { 
          do idHad = flavSel.combine( idEnds[iEnd3], idEnds[iEnd4]);
          while (idHad == 0);
          idProd.push_back( idHad); 
	} 
      }

      // Find masses of the new hadrons.
      for (int i = nFilled; i < int(idProd.size()) ; ++i) {
        double mHad = ParticleDataTable::mass(idProd[i]);
        mProd.push_back( mHad);
        mDiff -= mHad;
      } 
    }

  // Keep on trying until enough phase space. 
  } while (mDiff < mSafety || diquarkClash);

  // Update particle multiplicity.
  mult = idProd.size() - 1;

  // Done.
  return true;

}

//*********

// Set colour flow and scale in a decay explicitly to partons.

bool ParticleDecays::setColours(Event& event) {

  // Decay to q qbar (or qbar q).
  if (mode == 32 && idProd[1] > 0 && idProd[1] < 9) {
    int newCol = event.nextColTag();
    cols[1] = newCol;
    acols[2] = newCol;
  } else if (mode == 32 && idProd[1] < 0 && idProd[1] > -9) {
    int newCol = event.nextColTag();
    cols[2] = newCol;
    acols[1] = newCol;

  // Decay to gg.
  } else if (mode == 32 && idProd[1] == 21) {
    int newCol1 = event.nextColTag();
    int newCol2 = event.nextColTag();
    cols[1] = newCol1;
    acols[1] = newCol2;
    cols[2] = newCol2;
    acols[2] = newCol1;

  // Decay to g g g.
  } else if (mode == 33 && idProd[1] == 21 && idProd[2] == 21 
    &&  idProd[3] == 21) { 
    int newCol1 = event.nextColTag();
    int newCol2 = event.nextColTag();
    int newCol3 = event.nextColTag();
    cols[1] = newCol1;
    acols[1] = newCol2;
    cols[2] = newCol2;
    acols[2] = newCol3;
    cols[3] = newCol3;
    acols[3] = newCol1;

  // Decay to g g gamma: locate which is gamma.
  } else if (mode == 33) {
    int iGlu1 = (idProd[1] == 21) ? 1 : 3;
    int iGlu2 = (idProd[2] == 21) ? 2 : 3;
    int newCol1 = event.nextColTag();
    int newCol2 = event.nextColTag();
    cols[iGlu1] = newCol1;
    acols[iGlu1] = newCol2;
    cols[iGlu2] = newCol2;
    acols[iGlu2] = newCol1;
   
  // Unknown decay mode means failure.
  } else return false;

  // Set maximum scale to be mass of decaying particle.
  double scale = mProd[0];
  for (int i = 1; i <= mult; ++i) event[iProd[i]].scale(scale);

  // Done.
  return true;
     
}

//**************************************************************************

} // end namespace Pythia8

