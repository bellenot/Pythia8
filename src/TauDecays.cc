// TauDecays.cc is a part of the PYTHIA event generator.
// Copyright (C) 2011 Philip Ilten, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the TauDecays class.

#include "TauDecays.h"

namespace Pythia8 {

//==========================================================================

// The TauDecays class.

//--------------------------------------------------------------------------

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.
  
  // Number of times to try a decay channel.
  const int    TauDecays::NTRYCHANNEL = 10;
  
  // Number of times to try a decay sampling.
  const int    TauDecays::NTRYDECAY   = 10000;

//--------------------------------------------------------------------------

  /*
    Initialize the TauDecays class with the necessary pointers: info,
    particle data, random number pointer, and Standard Model constants pointer.
    Additionally, the necessary matrix elements are initialized with the
    Standard Model constants, and particle data pointers.
   */
  void TauDecays::init(Info* INFO, Settings* SET, ParticleData* DATA,
		       Rndm* RNDM, Couplings *SM) {

    // Set the pointers.
    info = INFO; settings = SET; particleData = DATA; rndm = RNDM; smData = SM;

    // Initialize the hard matrix elements.
    hmeTwoFermions2W2TwoFermions     .initPointers(smData, particleData);
    hmeTwoFermions2Z2TwoFermions     .initPointers(smData, particleData);
    hmeTwoFermions2Gamma2TwoFermions .initPointers(smData, particleData);
    hmeTwoFermions2GammaZ2TwoFermions.initPointers(smData, particleData);
    hmeHiggsEven2TwoFermions         .initPointers(smData, particleData);
    hmeHiggsOdd2TwoFermions          .initPointers(smData, particleData);
    hmeHiggsCharged2TwoFermions      .initPointers(smData, particleData);
    hmeUnpolarized                   .initPointers(smData, particleData);

    // Initialize the tau decay matrix elements.
    hmeTau2Meson                   .initPointers(smData, particleData);
    hmeTau2TwoLeptons              .initPointers(smData, particleData);
    hmeTau2TwoMesonsViaVector      .initPointers(smData, particleData);
    hmeTau2TwoMesonsViaVectorScalar.initPointers(smData, particleData);
    hmeTau2ThreePions              .initPointers(smData, particleData);
    hmeTau2FourPions               .initPointers(smData, particleData);
    hmeTau2PhaseSpace              .initPointers(smData, particleData);
   }

//--------------------------------------------------------------------------

  /*
    Main method of the TauDecays class. Pass the index of the tau requested
    to be decayed along with the event record in which the tau exists. The
    tau is then decayed with proper spin correlations as well any partner. The
    children of the decays are written to the event record, and if the decays
    were succesful, a return value of true is supplied.
   */
  bool TauDecays::decay(int idxOut1, Event& event) {

    // Set the first outgoing particle of the hard process.
    out1 = HelicityParticle(event[idxOut1]); out1.idx = idxOut1;
    
    // Find the mediator of the hard process.
    int idxMediator = out1.mother1();
    int idxFirstOut1 = idxOut1;
    while(idxMediator > 0 && event[idxMediator].id() == out1.id()) {
      idxFirstOut1 = idxMediator;
      idxMediator = event[idxMediator].mother1();
    }
    mediator = HelicityParticle(event[idxMediator]);
    mediator.idx = idxMediator; mediator.direction = -1;
    
    // Find the incoming particles of the hard process.
    int idxIn1 = mediator.mother1();
    int idxIn2 = mediator.mother2();
    while(idxIn1 > 0 && event[idxIn1].id() == mediator.id()) {
      idxIn1 = event[idxIn1].mother1();
      idxIn2 = event[idxIn2].mother2();
    }
    in1 = HelicityParticle(event[idxIn1]); in1.idx = idxIn1; in1.direction = -1;
    in2 = HelicityParticle(event[idxIn2]); in2.idx = idxIn2; in2.direction = -1;
    
    // Find the second outgoing particle of the hard process.
    int idxOut2 = (mediator.daughter1() == idxFirstOut1)
      ? mediator.daughter2() : mediator.daughter1();
    while (idxOut2 > 0 && event[idxOut2].daughter1() != 0) {
      idxOut2 = event[idxOut2].daughter1();
    }
    out2 = HelicityParticle(event[idxOut2]); out2.idx = idxOut2;

    // Set the particles vector.
    particles.clear();
    particles.push_back(in1);
    particles.push_back(in2);
    particles.push_back(out1);
    particles.push_back(out2);
    
    // Set the hard matrix element.
    // Produced from a W.
    if (abs(mediator.id()) == 24 &&
	(abs(in1.id()) <= 18 || abs(in2.id()) <= 18)) {
      // S-channel production.
      if (abs(in1.id()) <= 18 && abs(in2.id()) <= 18)
	hardME = hmeTwoFermions2W2TwoFermions.initChannel(particles);
      // T-channel production. 
      else {
	bool fermion = abs(in1.id()) <= 18 ? 0 : 1;
	particles[!fermion] = event[particles[fermion].daughter1()].id() 
	  == mediator.id() ?
	  HelicityParticle(event[particles[fermion].daughter2()]) :
	  HelicityParticle(event[particles[fermion].daughter1()]);
	particles[!fermion].direction = 1;
	if (abs(particles[!fermion].id()) <= 18)
	  hardME = hmeTwoFermions2W2TwoFermions.initChannel(particles);
	else {					
	  info->errorMsg("Warning in TauDecays::decay: unknown "
			 "tau production, assuming unpolarized "
			 "and uncorrelated");
	  hardME = hmeUnpolarized.initChannel(particles);
	}
      }
      correlated = false;
    // Produced from a photon.
    } else if (abs(mediator.id()) == 22 && abs(in1.id()) <= 18) {
      particles.push_back(mediator);
      hardME = hmeTwoFermions2Gamma2TwoFermions.initChannel(particles);
      correlated = true;
    // Produced from a photon/Z.
    } else if (abs(mediator.id()) == 23 && abs(in1.id()) <= 18) {
      particles.push_back(mediator);
      if (settings->mode("WeakZ0:gmZmode") == 0)
	hardME = hmeTwoFermions2GammaZ2TwoFermions.initChannel(particles);
      else if (settings->mode("WeakZ0:gmZmode") == 1)
	hardME = hmeTwoFermions2Gamma2TwoFermions.initChannel(particles);
      else if (settings->mode("WeakZ0:gmZmode") == 2)
	hardME = hmeTwoFermions2Z2TwoFermions.initChannel(particles);
      correlated = true;
    // Produced from a CP even Higgs.
    } else if (abs(mediator.id()) == 25 || abs(mediator.id()) == 35) {
      hardME = hmeHiggsEven2TwoFermions.initChannel(particles);
      correlated = true;
    // Produced from a CP odd Higgs.
    } else if (abs(mediator.id()) == 36) {
      hardME = hmeHiggsOdd2TwoFermions.initChannel(particles);
      correlated = true;
    // Produced from a charged Higgs.
    } else if (abs(mediator.id()) == 37) {
      hardME = hmeHiggsCharged2TwoFermions.initChannel(particles);
      correlated = false;
    // Produced from a D or B meson decay with a single tau.
    } else if ((abs(mediator.id()) == 411 || abs(mediator.id()) == 431 ||
    		abs(mediator.id()) == 511 || abs(mediator.id()) == 521 ||
    		abs(mediator.id()) == 531 || abs(mediator.id()) == 541) &&
    	       abs(out2.id()) == 16) {
      particles[0] = HelicityParticle(mediator.id() > 0 ? -5 : 5, 0, 0, 0, 0,
				      0, 0, 0, 0., 0., 0., 0., 0., 0.,
				      particleData);
      particles[1] = HelicityParticle(mediator.id() > 0 ? 5 : -5, 0, 0, 0, 0,
				      0, 0, 0, 0., 0., 0., 0., 0., 0.,
				      particleData);
      particles[0].idx = 0; particles[1].idx = 1;
      // D or B meson decays into neutrino + tau + meson.
      if (mediator.daughter1()+2 == mediator.daughter2()) {
	particles[0].p(mediator.p());
      	particles[1].direction = 1; particles[1].id(-particles[1].id());
      	particles[1].p(particles[0].p() - particles[2].p() - particles[3].p());
      }
      // D or B meson decays into neutrino + tau.
      else {
    	particles[0].p(mediator.p()/2);
    	particles[1].p(mediator.p()/2);
      }
      hardME = hmeTwoFermions2W2TwoFermions.initChannel(particles);
      correlated = false;
    // Produced from a virtual photon with correlated taus.
    } else if (abs(out1.id()) == 15 && abs(out2.id()) == 15) {
      particles.push_back(mediator);
      particles[0] = HelicityParticle(-1, 0, 0, 0, 0, 0, 0,
				      0, mediator.p()/2, 0., 0., particleData);
      particles[1] = HelicityParticle(1, 0, 0, 0, 0, 0, 0,
				      0, mediator.p()/2, 0., 0., particleData);
      particles[0].direction = -1; particles[1].direction = -1;
      particles[0].idx = 0; particles[1].idx = 0;
      hardME = hmeTwoFermions2Gamma2TwoFermions.initChannel(particles);
      correlated = true;
    // Produced from an unknown process, assume unpolarized and uncorrelated.
    } else {
      info->errorMsg("Warning in TauDecays::decay: unknown "
		     "tau production, assuming unpolarized and uncorrelated");
      hardME = hmeUnpolarized.initChannel(particles);
      correlated = false;
    }

    // Pick the first tau to decay.
    HelicityParticle* tau;
    int idx;
    if (correlated) idx = rndm->flat() < 0.5 ? 2 : 3;
    else idx = abs(particles[2].id()) == 15 ? 2 : 3;
    tau = &particles[idx];

    // Calculate the density matrix and decay the tau.
    hardME->calculateRho(idx, particles);
    vector<HelicityParticle> children = createChildren(*tau);
    if (children.size() == 0) return false;

    // Decay the first tau.
    bool accepted = false;
    int  tries    = 0;
    while (!accepted) {
      isotropicDecay(children);
      double decayWeight    = decayME->decayWeight(children);
      double decayWeightMax = decayME->decayWeightMax(children);
      accepted = rndm->flat() < decayWeight / decayWeightMax;
      if (decayWeight > decayWeightMax)
	info->errorMsg("Warning in TauDecays::decay: maximum "
		       "decay weight exceeded in tau decay");
      if (tries > NTRYDECAY) {
	info->errorMsg("Warning in TauDecays::decay: maximum "
		       "number of decay attempts exceeded");
	break;
      }
      tries++;
    }
    writeDecay(event,children);

    // If a correlated second tau exists, decay that tau as well.
    if (correlated) {
      idx = idx == 2 ? 3 : 2;
      // Calculate the first tau decay matrix.
      decayME->calculateD(children);
      // Update the decay matrix for the tau.
      (*tau).D = children[0].D;
      // Switch the taus.
      tau = &particles[idx];
      // Calculate second tau's density matrix.
      hardME->calculateRho(idx, particles);
      // Decay the second tau.
      children.clear();
      children = createChildren(*tau);
      if (children.size() == 0) return false;
      accepted = false;
      tries    = 0;
      while (!accepted) {
    	isotropicDecay(children);
	double decayWeight    = decayME->decayWeight(children);
	double decayWeightMax = decayME->decayWeightMax(children);
	accepted = rndm->flat() < decayWeight / decayWeightMax;
	if (decayWeight > decayWeightMax)
	  info->errorMsg("Warning in TauDecays::decay: maximum "
			 "decay weight exceeded in correlated tau decay"); 
	if (tries > NTRYDECAY) {
	  info->errorMsg("Warning in TauDecays::decay: maximum "
			 "number of decay attempts exceeded");
	  break;
	}
	tries++;
      }
      writeDecay(event,children);
    }

    return true;

  }

//--------------------------------------------------------------------------

  /*
    Given a HelicityParticle parent, select the decay channel and return
    a vector of HelicityParticles containing the children, with the parent
    particle duplicated in the first entry of the vector.
   */
  vector<HelicityParticle> TauDecays::createChildren(HelicityParticle parent)
  {
    int meMode = 0;
    vector<HelicityParticle> children;

    // Set the parent as incoming.
    parent.direction = -1;

    // Setup decay data for the decaying particle.
    ParticleDataEntry decayData = parent.particleDataEntry();

    // Initialize the decay data.
    if (!decayData.preparePick(parent.id())) return children;

    // Try to pick a decay channel.
    bool decayed = false;
    int decayTries = 0;
    while (!decayed && decayTries < NTRYCHANNEL) {

      // Pick a decay channel.
      DecayChannel decayChannel = decayData.pickChannel();
      meMode = decayChannel.meMode();
      int decayMult = decayChannel.multiplicity();

      // Select children masses.
      bool allowed = false;
      int channelTries = 0;
      while (!allowed && channelTries < NTRYCHANNEL) {

	children.resize(0);
	children.push_back(parent);
      
        for (int i = 0; i < decayMult; ++i) {
	  // Grab child ID.
          int childId = decayChannel.product(i);
	  
	  // Flip sign for anti-particle decay.
          if (parent.id() < 0 && particleData->hasAnti(childId))
	    childId = -childId;
	  
	  // Grab child mass.
          double childMass = particleData->mass(childId);
	  
	  // Push back the child into the children vector.
	  children.push_back(HelicityParticle(childId, 91, parent.idx, 0, 0, 
					      0, 0, 0, 0., 0., 0., 0.,
					      childMass, 0., particleData));
	}
	
        // Check there is enough phase space for decay.
        if (decayMult > 1) {
          double massDiff = parent.m();
          for (int i = 0; i < decayMult; ++i) massDiff -= children[i].m();

	  // For now we just check kinematically available.
          if (massDiff > 0) { allowed = true; decayed = true; }
        }
	
	channelTries++;
      }

      decayTries++;
    }

    // Set the decay matrix element.
    
    // Two body decays.
    if (children.size() == 3) {
      if (meMode == 1521) {
	decayME = hmeTau2Meson.initChannel(children);
      }
      else decayME = hmeTau2PhaseSpace.initChannel(children);
    } 
    // Three body decays.
    else if (children.size() == 4) {
      // Leptonic decay.
      if (meMode == 1531)
	decayME = hmeTau2TwoLeptons.initChannel(children);
      // Two meson decay via vector meson.
      else if (meMode == 1532)
	decayME = hmeTau2TwoMesonsViaVector.initChannel(children);
      // Two meson decay via vector or scalar meson.
      else if (meMode == 1533)
	decayME = hmeTau2TwoMesonsViaVectorScalar.initChannel(children);
      // Flat phase space.
      else decayME = hmeTau2PhaseSpace.initChannel(children);
    }
    // Four body decays.
    else if (children.size() == 5) {
      // Three pion CLEO decay.
      if (meMode == 1541)
	decayME = hmeTau2ThreePions.initChannel(children);
      // Flat phase space.
      else decayME = hmeTau2PhaseSpace.initChannel(children);
    }
    // Five body decays.
    else if (children.size() == 6) {
      // Four pion Novosibirsk current.
      if (meMode == 1551)
	decayME = hmeTau2FourPions.initChannel(children);
      // Flat phase space.
      else decayME = hmeTau2PhaseSpace.initChannel(children);
    }
    // Flat phase space.
    else decayME = hmeTau2PhaseSpace.initChannel(children);
    return children;
  }

//--------------------------------------------------------------------------

  /*
    N-body decay using the M-generator algorithm described in
    "Monte Carlo Phase Space" by F. James in CERN 68-15, May 1968. Taken
    from ParticleDecays::mBody but modified to handle spin particles.

    Given a vector of HelicityParticles where the first particle is the mother,
    the remaining particles are decayed isotropically.
  */
  void TauDecays::isotropicDecay(vector<HelicityParticle>& p)
  {
    int decayMult = p.size() - 1;

    // Set the weight corrections.
    double WTCORRECTION[11] = { 1., 1., 1., 2., 5., 15., 60., 250., 1250.,
				7000., 50000. };

    // Mother and sum daughter masses.
    double m0      = p[0].m();
    double mSum    = p[1].m();
    for (int i = 2; i <= decayMult; ++i) mSum += p[i].m(); 
    double mDiff   = m0 - mSum;
      
    // Begin setup of intermediate invariant masses.
    vector<double> mInv;
    for (int i = 0; i <= decayMult; ++i) mInv.push_back( p[i].m());
      
    // Calculate the maximum weight in the decay.
    double wtPS;
    double wtPSmax = 1. / WTCORRECTION[decayMult];
    double mMax    = mDiff + p[decayMult].m();
    double mMin    = 0.; 
    for (int i = decayMult - 1; i > 0; --i) {
      mMax        += p[i].m();
      mMin        += p[i+1].m();
      double mNow  = p[i].m();
      wtPSmax *= 0.5 * sqrtpos( (mMax - mMin - mNow) * (mMax + mMin + mNow)
				* (mMax + mMin - mNow) * (mMax - mMin + mNow) 
				) / mMax;  
    }
      
    // Begin loop to find the set of intermediate invariant masses.
    vector<double> rndmOrd;
    do {
      wtPS  = 1.;
	  
      // Find and order random numbers in descending order.
      rndmOrd.clear();
      rndmOrd.push_back(1.);
      for (int i = 1; i < decayMult - 1; ++i) {
	double random = rndm->flat();
	rndmOrd.push_back(random);
	for (int j = i - 1; j > 0; --j) {
	  if (random > rndmOrd[j]) swap( rndmOrd[j], rndmOrd[j+1] );
	  else break;
	} 
      }
      rndmOrd.push_back(0.);
	  
      // Translate into intermediate masses and find weight.
      for (int i = decayMult - 1; i > 0; --i) {
	mInv[i] = mInv[i+1] + p[i].m() + (rndmOrd[i-1] - rndmOrd[i])
	  * mDiff; 
	wtPS   *= 0.5 * sqrtpos( (mInv[i] - mInv[i+1] - p[i].m()) 
				 * (mInv[i] + mInv[i+1] + p[i].m()) 
				 * (mInv[i] + mInv[i+1] - p[i].m()) 
				 * (mInv[i] - mInv[i+1] + p[i].m()) )
	  / mInv[i];  
      }
	  
      // If rejected, try again with new invariant masses.
    } while ( wtPS < rndm->flat() * wtPSmax ); 
	
    // Perform two-particle decays in the respective rest frame.
    vector<Vec4> pInv(decayMult + 1);
    for (int i = 1; i < decayMult; ++i) {
      double pAbs = 0.5 * sqrtpos( (mInv[i] - mInv[i+1] - p[i].m()) 
				   * (mInv[i] + mInv[i+1] + p[i].m())
				   * (mInv[i] + mInv[i+1] - p[i].m())
				   * (mInv[i] - mInv[i+1] + p[i].m()) )
	/ mInv[i]; 

      // Isotropic angles give three-momentum.
      double cosTheta = 2. * rndm->flat() - 1.;
      double sinTheta = sqrt(1. - cosTheta*cosTheta);
      double phi      = 2. * M_PI * rndm->flat();
      double pX       = pAbs * sinTheta * cos(phi);  
      double pY       = pAbs * sinTheta * sin(phi);  
      double pZ       = pAbs * cosTheta;  

      // Calculate energies, fill four-momenta.
      double eHad     = sqrt( p[i].m()*p[i].m() + pAbs*pAbs);
      double eInv     = sqrt( mInv[i+1]*mInv[i+1] + pAbs*pAbs);
      p[i].p( pX, pY, pZ, eHad);
      pInv[i+1].p( -pX, -pY, -pZ, eInv);
    }       
  
    // Boost decay products to the mother rest frame.
    p[decayMult].p( pInv[decayMult] );
    for (int iFrame = decayMult - 1; iFrame > 1; --iFrame) 
      for (int i = iFrame; i <= decayMult; ++i) 
	p[i].bst( pInv[iFrame], mInv[iFrame]);

    // Boost decay products to the current frame and return. 
    pInv[1].p( p[0].p() );
    for (int i = 1; i <= decayMult; ++i) p[i].bst( pInv[1], mInv[1] );
    return;
  }

//--------------------------------------------------------------------------

  /*
    Write the vector of HelicityParticles to the event record, excluding the
    first particle. Set the lifetime and production vertex of the particles
    and mark the first particle of the vector as decayed.
  */
  void TauDecays::writeDecay(Event& event, vector<HelicityParticle>& p) {
    int  decayMult   = p.size() - 1;
    Vec4 decayVertex = p[0].vDec();
    
    // Set additional information and append children to event.
    for (int i = 1; i <= decayMult; i++) {
      // Set child lifetime.
      p[i].tau(p[i].tau0() * rndm->exp());
      // Set child production vertex.
      p[i].vProd(decayVertex);
      // Append child to record.
      p[i].idx = event.append(p[i]);
    }
    
    // Mark the parent as decayed and set children.
    event[p[0].idx].statusNeg(); 
    event[p[0].idx].daughters(p[1].idx, p[decayMult].idx);
  }

//==========================================================================

} // end namespace Pythia8
