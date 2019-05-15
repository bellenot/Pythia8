// Function definitions (not found in the header) for the
// SpaceShower class.
// Copyright C 2006 Torbjorn Sjostrand 

#include "SpaceShower.h"

namespace Pythia8 {

//**************************************************************************

// The SpaceShower class.

//*********
 
// Definitions of static variables.
// (Values will be overwritten in initStatic call, so are purely dummy.)

bool SpaceShower::doQCDshower = true;
bool SpaceShower::doQEDshowerByQ = true;
bool SpaceShower::doQEDshowerByL = true;
int SpaceShower::pTmaxMatch = 0;
double SpaceShower::mc = 1.5;
double SpaceShower::mb = 4.8;
double SpaceShower::mc2 = 2.25;
double SpaceShower::mb2 = 23.04;
double SpaceShower::alphaSvalue = 0.127;
int SpaceShower::alphaSorder = 1;
bool SpaceShower::samePTasMI = true;
double SpaceShower::pT0Ref= 2.2; 
double SpaceShower::ecmRef = 1800.; 
double SpaceShower::ecmPow = 0.16; 
double SpaceShower::pTmin = 0.2; 
double SpaceShower::alphaEM = 0.00729735;
double SpaceShower::pTminChgQ = 0.5;
double SpaceShower::pTminChgL = 0.5e-3;
bool SpaceShower::doMEcorrections = true;
bool SpaceShower::doPhiPolAsym = true;
int SpaceShower::nQuark = 5;

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Switch to alternative (but equivalent) backwards evolution for
// g -> Q Qbar (Q = c or b) when below QTHRESHOLD * mQ2.
const double SpaceShower::CTHRESHOLD = 2.0; 
const double SpaceShower::BTHRESHOLD = 2.0; 

// Renew evaluation of PDF's when the pT2 step is bigger than this 
// (in addition to initial scale and c and b thresholds.)
const double SpaceShower::EVALPDFSTEP = 0.1;

// Lower limit on PDF value in order to avoid division by zero.
const double SpaceShower::TINYPDF = 1e-10;

// Lower limit on estimated evolution rate, below which stop.
const double SpaceShower::TINYKERNELPDF = 1e-5;

// Lower limit on pT2, below which branching is rejected. 
const double SpaceShower::TINYPT2 = 1e-6;

// No attempt to do backwards evolution of a heavy (c or b) quark 
// if evolution starts at a scale pT2 < HEAVYPT2EVOL * mQ2.
const double SpaceShower::HEAVYPT2EVOL = 1.1;

// No attempt to do backwards evolution of a heavy (c or b) quark 
// if evolution starts at a  x > HEAVYXEVOL * x_max, where 
// x_max is the largest possible x value for a g -> Q Qbar branching.
const double SpaceShower::HEAVYXEVOL = 0.9;
  
// When backwards evolution Q -> g + Q greates a heavy quark Q,
// an earlier branching g -> Q + Qbar will restrict kinematics
// to  M_{Q Qbar}^2 > EXTRASPACEQ * 4 m_Q^2. (Smarter to be found??) 
const double SpaceShower::EXTRASPACEQ = 2.0;

//*********

// Initialize static data members.

void SpaceShower::initStatic() {

  // Main flags to switch on and off branchings.
  doQCDshower = Settings::flag("SpaceShower:QCDshower");
  doQEDshowerByQ = Settings::flag("SpaceShower:QEDshowerByQ");
  doQEDshowerByL = Settings::flag("SpaceShower:QEDshowerByL");

  // Matching in pT of hard interaction to shower evolution.
  pTmaxMatch = Settings::mode("SpaceShower:pTmaxMatch"); 

  // Charm and bottom mass thresholds.
  mc = ParticleDataTable::m0(4); 
  mb = ParticleDataTable::m0(5); 
  mc2 = mc * mc;
  mb2 = mb * mb;

  // Parameters of alphaStrong generation.
  alphaSvalue = Settings::parameter("SpaceShower:alphaSvalue");
  alphaSorder = Settings::mode("SpaceShower:alphaSorder");
 
  // Regularization of QCD evolution for pT -> 0. Can be taken 
  // same as for multiple interactions, or be set separately.
  samePTasMI = Settings::flag("SpaceShower:samePTasMI"); 
  if (samePTasMI) {
    pT0Ref = Settings::parameter("MultipleInteractions:pT0Ref");
    ecmRef = Settings::parameter("MultipleInteractions:ecmRef");
    ecmPow = Settings::parameter("MultipleInteractions:ecmPow");
    pTmin = Settings::parameter("MultipleInteractions:pTmin");
  } else {
    pT0Ref = Settings::parameter("SpaceShower:pT0Ref");
    ecmRef = Settings::parameter("SpaceShower:ecmRef");
    ecmPow = Settings::parameter("SpaceShower:ecmPow");
    pTmin = Settings::parameter("SpaceShower:pTmin");
  }
 
  // Parameters of QED evolution.
  alphaEM = Settings::parameter("StandardModel:alphaEMfix");
  pTminChgQ = Settings::parameter("SpaceShower:pTminchgQ"); 
  pTminChgL = Settings::parameter("SpaceShower:pTminchgL"); 

  // Various other parameters. 
  doMEcorrections = Settings::flag("SpaceShower:MEcorrections");
  doPhiPolAsym = Settings::flag("SpaceShower:phiPolAsym");
  nQuark = Settings::mode("SpaceShower:nQuark");

}

//*********

// Initialize alphaStrong and pTmin parameters.

void SpaceShower::init( BeamParticle* beamAPtrIn, 
  BeamParticle* beamBPtrIn) {

  // Store input pointers for future use. 
  beamAPtr = beamAPtrIn;
  beamBPtr = beamBPtrIn;

  // Combinations of fixed couplings.
  alphaS2pi = 0.5 * alphaSvalue / M_PI;
  alphaEM2pi = 0.5 * alphaEM / M_PI;
  
  // Initialize alpha_strong generation.
  alphaS.init( alphaSvalue, alphaSorder); 
  
  // Lambda for 5, 4 and 3 flavours.
  Lambda3flav = alphaS.Lambda3(); 
  Lambda4flav = alphaS.Lambda4(); 
  Lambda5flav = alphaS.Lambda5(); 
  Lambda5flav2 = pow2(Lambda5flav);
  Lambda4flav2 = pow2(Lambda4flav);
  Lambda3flav2 = pow2(Lambda3flav);

  // Calculate invariant mass of system. Set current pT0 scale.
  sCM = m2( beamAPtr->p(), beamBPtr->p());
  eCM = sqrt(sCM);
  pT0 = pT0Ref * pow(eCM / ecmRef, ecmPow);

  // Derived parameters of QCD evolution.
  pT20 = pT0*pT0;
  pT2min = pTmin*pTmin;
  pT2minChgQ = pTminChgQ*pTminChgQ;
  pT2minChgL = pTminChgL*pTminChgL;

} 

//*********

// Find whether to limit maximum scale of emissions.

bool SpaceShower::limitPTmax( Event& event) {

  // User-set cases.
  if (pTmaxMatch == 1) return true;
  if (pTmaxMatch == 2) return false;
   
  // Look if any quark (u, d, s, c, b), gluon or photon in final state. 
  bool hasQGP = false;
  for (int i = 5; i < event.size(); ++i) {
    int idAbs=event[i].id();
    if (idAbs <= 5 || idAbs == 21 || idAbs == 22) hasQGP = true;
  }
  return (hasQGP) ? true : false;
 
}

//*********

// Prepare system for evolution; identify ME.
// Routine may be called after multiple interactions, 
// skipping part up to sizeOld, which was already processed.

void SpaceShower::prepare( Event& event, bool limitPTmax, int sizeOld) {

 // Reset list of (sub)systems first time around.
  if (sizeOld == 0) system.resize(0);

  // Find positions of incoming colliding partons.
  int in1 = 0;
  int in2 = 0;
  for (int i = sizeOld; i < int(event.size()); ++i) {  
    if (event[i].status() == -21 || event[i].status() == -31) {
      if (event[i].mother1() == 1) in1 = i;
      if (event[i].mother1() == 2) in2 = i;
    }
    if (in1 > 0 && in2 > 0) break;
  } 
  if (in1 <= 0 || in2 <= 0) {
    ErrorMessages::message("Error in SpaceShower::prepare: "
      "failed to locate the two incoming partons"); 
    return;
  }

  // Find list of all particles contained in incoming system.
  iParton.resize(0);
  iParton.push_back(in1);
  iParton.push_back(in2);
  for (int i = sizeOld; i < int(event.size()); ++i) {
    if (i != in1 && i != in2) {   
      int iMother = event[i].mother1();
      while ( iMother > 0 && iMother != in1 && iMother != in2 ) 
        iMother = event[iMother].mother1();  
      if (iMother == in1 || iMother == in2) iParton.push_back(i);
    }
  }

  // Store new system.
  system.push_back(iParton);
  sysNow = &system[system.size() - 1];  

  // Find further properties of the system.
  sysNow->id1 = event[in1].id();
  sysNow->id2 = event[in2].id();
  sysNow->m =  m( event[in1], event[in2] );
  sysNow->m2 = pow2(sysNow->m);
  sysNow->x1 = event[in1].pPlus() / event[0].e();  
  sysNow->x2 = event[in2].pMinus() / event[0].e();   

  // Find dipole ends for QCD radiation.
  double pTmax1 = (limitPTmax) ? event[in1].scale() : eCM;
  int colType1 = doQCDshower ? event[in1].colType() : 0;
  sysNow->dipEnd[0] = SpaceDipoleEnd( pTmax1, 1, colType1, 0, -1) ;
  double pTmax2 = (limitPTmax) ? event[in2].scale() : eCM;
  int colType2 = doQCDshower ? event[in2].colType() : 0;
  sysNow->dipEnd[1] = SpaceDipoleEnd( pTmax2, 2, colType2, 0, -1) ;

  // Find dipole ends for QED radiation. 
  int chgType1 = ( (event[in1].isQ() && doQEDshowerByQ)
    || (event[in1].isL() && doQEDshowerByL) )
    ? event[in1].icharge() : 0;
  sysNow->dipEnd[2] = SpaceDipoleEnd( pTmax1, 1, 0, chgType1, -1) ;
  int chgType2 = ( (event[in2].isQ() && doQEDshowerByQ)
    || (event[in2].isL() && doQEDshowerByL) )
    ? event[in2].icharge() : 0;
  sysNow->dipEnd[3] = SpaceDipoleEnd( pTmax2, 2, 0, chgType2, -1) ;

  // Now find matrix element corrections for system.
  findMEtype( event); 

}

//*********
 
// Select next pT in downwards evolution of the existing dipoles.

double SpaceShower::pTnext( double pTbegAll, double pTendAll) {
  
  // Starting values: no radiating dipole found.
  double pT2sel = pow2(pTendAll);
  sysSel = 0;
  dipEndSel = 0; 

  // Loop over all possible radiating systems.
  for (int iSystem = 0; iSystem < int(system.size()); ++iSystem) {
    iSysNow = iSystem;
    sysNow = &system[iSystem];  

    // Loop over four dipole ends for each radiating system.
    for (int iDipEnd = 0; iDipEnd < 4; ++iDipEnd) {
      dipEndNow = &sysNow->dipEnd[iDipEnd];        
      dipEndNow->pT2 = 0.;
   
      // Check whether dipole end should be allowed to shower. 
      double pT2begDip = pow2( min( pTbegAll, dipEndNow->pTmax ));
      if (pT2begDip > pT2sel 
        && ( dipEndNow->colType != 0 || dipEndNow->chgType != 0 ) ) {
        double pT2endDip = 0.;

        // Determine lower cut for evolution, for QCD or QED (q or l).      
        if (dipEndNow->colType != 0) pT2endDip = max( pT2sel, pT2min );   
        else if (abs(dipEndNow->chgType) != 3) pT2endDip 
          = max( pT2sel, pT2minChgQ );   
        else pT2endDip = max( pT2sel, pT2minChgL );  

        // Now do evolution in pT2, for QCD or QED 
        if (pT2begDip > pT2endDip) { 
          if (dipEndNow->colType != 0) pT2nextQCD( pT2begDip, pT2endDip);
          else pT2nextQED( pT2begDip, pT2endDip);
        }

        // Update if found larger pT than current maximum.
        if (dipEndNow->pT2 > pT2sel) {
          pT2sel = dipEndNow->pT2;
          iSysSel = iSysNow;
          sysSel = sysNow;
          dipEndSel = dipEndNow;
        }
      }

    // End loop over dipole ends and systems.
    }
  } 
  
  // Check if g -> Q + Qbar.
  if (dipEndSel != 0 && dipEndSel->idMother == 21 
    && ( abs(dipEndSel->idSister) == 4 || abs(dipEndSel->idSister) == 5) ) {
  } ; 

  // Return nonvanishing value if found pT is bigger than already found.
  return (sysSel == 0) ? 0. : sqrt(pT2sel); 
}

//*********

// Evolve a QCD dipole end. 

void SpaceShower::pT2nextQCD( double pT2begDip, double pT2endDip) { 

  // Read beam particle and incoming parton from appropriate side. 
  BeamParticle& beam = (dipEndNow->side == 1) ? *beamAPtr : *beamBPtr;
  int idDaughter = (dipEndNow->side == 1) ? sysNow->id1 : sysNow->id2;
  double xDaughter = (dipEndNow->side == 1) ? sysNow->x1 : sysNow->x2;
  bool isGluon = (idDaughter == 21) ? true : false;
  bool isValence = beam[iSysNow].isValence();
  int MEtype = dipEndNow->MEtype;

  // Some kinematical starting values.
  double pT2 = pT2begDip;
  double m2 = sysNow->m2;
  double xMaxAbs = beam.xMax(iSysNow);
  double zMinAbs = xDaughter / xMaxAbs;

  // Evolution below scale of massive quark or at large x is impossible.
  double idMassive = 0;
  if ( abs(idDaughter) == 4 ) idMassive = 4;
  if ( abs(idDaughter) == 5 ) idMassive = 5;
  bool isMassive = (idMassive > 0) ? true : false;
  double m2Massive = 0.;
  double mRatio = 0.;
  double zMaxMassive = 1.;
  double m2Threshold = pT2;
  if (isMassive) { 
    m2Massive = (idMassive == 4) ? mc2 : mb2;
    if (pT2 < HEAVYPT2EVOL * m2Massive) return;
    mRatio = sqrt( m2Massive / m2 );
    zMaxMassive = (1. -  mRatio) / ( 1. +  mRatio * (1. -  mRatio) ); 
    if (xDaughter > HEAVYXEVOL * zMaxMassive * xMaxAbs) return; 
  
    // Find threshold scale below which only g -> Q + Qbar will be allowed.
    m2Threshold = (idMassive == 4) ? min( pT2, CTHRESHOLD * mc2)
      : min( pT2, BTHRESHOLD * mb2); 
  }
  
  // Variables used inside evolution loop. (Mainly dummy start values.)
  int nFlavour = 3; 
  double b0 = 4.5;
  double Lambda2 = Lambda3flav2;
  double pT2minNow = pT2endDip; 
  int idMother = 0; 
  int idSister = 0;
  double z = 0.;
  double zMaxAbs = 0.; 
  double g2gInt = 0.; 
  double q2gInt = 0.; 
  double q2qInt = 0.;
  double g2qInt = 0.;
  double g2Qenhance = 0.;
  double xPDFdaughter = 0.;
  double xPDFmother[21] = {0.};
  double xPDFgMother = 0.;
  double xPDFmotherSum = 0.;
  double kernelPDF = 0.;
  double xMother = 0.;
  double wt = 0.;
  double Q2 = 0.;
  double mSister = 0.;
  double m2Sister = 0.;
  double pT2corr = 0.;
  double phi = 0.; 
  double pT2PDF = pT2;
  bool needNewPDF = true;

  // Begin evolution loop towards smaller pT values.
  do { 
    wt = 0.;

    // Initialize integrals of splitting kernels and evaluate parton 
    // densities at the beginning. Reinitialize after long evolution 
    // in pT2 or when crossing c and b flavour thresholds.
    if (needNewPDF || pT2 < EVALPDFSTEP * pT2PDF) {
      pT2PDF = pT2;

      // Determine overestimated z range; switch at c and b masses.
      if (pT2 > mb2) {
        nFlavour = 5;
        pT2minNow = mb2;
        b0 = 23./6.;
        Lambda2 = Lambda5flav2;
      } else if (pT2 > mc2) {
        nFlavour = 4;
        pT2minNow = mc2;
        b0 = 25./6.;
        Lambda2 = Lambda4flav2;
      } else { 
        nFlavour = 3;
        pT2minNow = pT2endDip;
        b0 = 27./6.;
        Lambda2 = Lambda3flav2;
      }
      zMaxAbs = 1. - 0.5 * (pT2minNow / m2) *
        ( sqrt( 1. + 4. * m2 / pT2minNow ) - 1. );
      if (isMassive) zMaxAbs = min( zMaxAbs, zMaxMassive);  

      // Go to another z range with lower mass scale if current is closed.
      if (zMinAbs > zMaxAbs) { 
        if (nFlavour == 3 || (idMassive == 4 && nFlavour == 4) 
          || idMassive == 5) return;
        pT2 = (nFlavour == 4) ? mc2 : mb2;
        continue;
      } 

      // Parton density of daughter at current scale. 
      xPDFdaughter = max( TINYPDF, 
        beam.xfISR(iSysNow, idDaughter, xDaughter, pT2) );

      // Integrals of splitting kernels for gluons: g -> g, q -> g.
      if (isGluon) {
        g2gInt = 6. * log(zMaxAbs * (1.-zMinAbs) 
          / (zMinAbs * (1.-zMaxAbs)));
        if (doMEcorrections) g2gInt *= calcMEmax(MEtype, 21);
        q2gInt = (16./3.) * (1./sqrt(zMinAbs) - 1./sqrt(zMaxAbs));
        if (doMEcorrections) q2gInt *= calcMEmax(MEtype, 1);

        // Parton density of potential quark mothers to a g.
        xPDFmotherSum = 0.;
        for (int i = -nQuark; i <= nQuark; ++i) {
          if (i == 0) {
            xPDFmother[10] = 0.;
          } else {
            xPDFmother[i+10] = beam.xfISR(iSysNow, i, xDaughter, pT2); 
            xPDFmotherSum += xPDFmother[i+10]; 
          }
        } 

        // Total QCD evolution coefficient for a gluon.
        kernelPDF = g2gInt + q2gInt * xPDFmotherSum / xPDFdaughter;

      // For valence quark only need consider q -> q g branchings.
      } else if (isValence) {
        q2qInt = (8./3.) * log( (1. - zMinAbs) / (1. - zMaxAbs) );
        if (doMEcorrections) q2qInt *= calcMEmax(MEtype, 1);
        g2qInt = 0.;
        kernelPDF = q2qInt; 

      // Integrals of splitting kernels for quarks: q -> q, g -> q.
      } else {
        q2qInt = (8./3.) * log( (1. - zMinAbs) / (1. - zMaxAbs) );
        if (doMEcorrections) q2qInt *= calcMEmax(MEtype, 1);
        g2qInt = 0.5 * (zMaxAbs - zMinAbs);
        if (doMEcorrections) g2qInt *= calcMEmax(MEtype, 21);

        // Increase estimated upper weight for g -> Q + Qbar.
        if (isMassive) {
          double m2log = log( m2Massive / Lambda2);
          g2Qenhance = log( log(pT2/Lambda2) / m2log ) 
            / log( log(m2Threshold/Lambda2) / m2log );
          g2qInt *= g2Qenhance;
	}

        // Parton density of a potential gluon mother to a q.
        xPDFgMother = beam.xfISR(iSysNow, 21, xDaughter, pT2);

        // Total QCD evolution coefficient for a quark.
        kernelPDF = q2qInt + g2qInt * xPDFgMother / xPDFdaughter;
      }

      // End evaluation of splitting kernels and parton densities.
      needNewPDF = false;
    }
    if (kernelPDF < TINYKERNELPDF) { pT2 = 0.; continue; }

    // Pick pT2 (in overestimated z range), for one of three different cases.
    // Assume form alphas(pT0^2 + pT^2) * dpT^2/(pT0^2 + pT^2).

    // Fixed alpha_strong.
    if (alphaSorder == 0) {
      pT2 = (pT2 + pT20) * pow( Rndm::flat(), 
        1. / (alphaS2pi * kernelPDF)) - pT20;

    // First-order alpha_strong.
    } else if (alphaSorder == 1) {
      pT2 = Lambda2 * pow( (pT2 + pT20) / Lambda2, 
        pow(Rndm::flat(), b0 / kernelPDF) ) - pT20;

    // For second order reject by second term in alpha_strong expression.
    } else {
      do pT2 = Lambda2 * pow( (pT2 + pT20) / Lambda2, 
        pow(Rndm::flat(), b0 / kernelPDF) ) - pT20;
      while (alphaS.alphaS2OrdCorr(pT2 + pT20) < Rndm::flat() 
        && pT2 > pT2minNow);
    } 

    // Check for pT2 values that prompt special action.

    // If fallen into b threshold region, force g -> b + bbar.
    if (idMassive == 5 && pT2 < m2Threshold) {
      pT2nearQCDthreshold( beam, m2Massive, m2Threshold, zMinAbs, 
        zMaxMassive );
      return;

    // If crossed b threshold, continue evolution from this threshold.
    } else if (nFlavour == 5 && pT2 < mb2) {  
      needNewPDF = true;
      pT2 = mb2;
      continue;

    // If fallen into c threshold region, force g -> c + cbar.
    } else if (idMassive == 4 && pT2 < m2Threshold) {
      pT2nearQCDthreshold( beam, m2Massive, m2Threshold, zMinAbs, 
        zMaxMassive );
      return; 

    // If crossed c threshold, continue evolution from this threshold.
    } else if (nFlavour == 4 && pT2 < mc2) { 
      needNewPDF = true;
      pT2 = mc2;
      continue;

    // Abort evolution if below cutoff scale, or below another branching.
    } else if (pT2 < pT2endDip) return; 

    // Select z value of branching to g, and corrective weight.
    if (isGluon) {
      // g -> g (+ g). 
      if (Rndm::flat() * kernelPDF < g2gInt) {
        idMother = 21;
        idSister = 21;
        z = 1. / ( 1. + ((1. - zMinAbs) / zMinAbs) * pow( (zMinAbs * 
          (1. - zMaxAbs)) / (zMaxAbs * (1. - zMinAbs)), Rndm::flat() ) );
        wt = pow2( 1. - z * (1. - z));
      } else {
      // q -> g (+ q): also select flavour. 
        double temp = xPDFmotherSum * Rndm::flat();
        idMother = -nQuark-1;
        do { temp -= xPDFmother[(++idMother) + 10]; } 
        while (temp > 0. && idMother < nQuark);  
        idSister = idMother;
        z = (zMinAbs * zMaxAbs) / pow2( sqrt(zMinAbs) + Rndm::flat() 
          * ( sqrt(zMaxAbs)- sqrt(zMinAbs) ));
        wt = 0.5 * (1. + pow2(1. - z)) * sqrt(z) 
          * xPDFdaughter / xPDFmother[idMother + 10];
      } 

    // Select z value of branching to q, and corrective weight.
    // Include massive kernel corrections for c and b quarks.
    } else {
      // q -> q (+ g). 
      if (isValence || Rndm::flat() * kernelPDF < q2qInt) {
        idMother = idDaughter;
        idSister = 21;
        // Valence more peaked at large z, like q -> g above ??
        z = 1. - (1. - zMinAbs) * pow( (1. - zMaxAbs) / (1. - zMinAbs),
          Rndm::flat() ); 
        if (!isMassive) { 
          wt = 0.5 * (1. + pow2(z));
	} else {
          wt = 0.5 * (1. + pow2(z) - z * pow2(1.-z) * m2Massive / pT2);
	}
      // g -> q (+ qbar). 
      } else {
        idMother = 21;
        idSister = - idDaughter; 
        z = zMinAbs + Rndm::flat() * (zMaxAbs - zMinAbs);
        if (!isMassive) { 
          wt = (pow2(z) + pow2(1.-z)) * xPDFdaughter / xPDFgMother ;
	} else {
          wt = (pow2(z) + pow2(1.-z) + 2. * z * (1.-z) * m2Massive / pT2) 
            * xPDFdaughter / (xPDFgMother * g2Qenhance) ;
	}
      }
    }

    // Derive Q2 and x of mother from pT2 and z. 
    Q2 = pT2 / (1.- z);
    xMother = xDaughter / z;
 
    // Forbidden emission if outside allowed z range for given pT2.
    mSister = ParticleDataTable::m0(idSister);
    m2Sister = pow2(mSister);
    pT2corr = Q2 - z * (m2 + Q2) * (Q2 + m2Sister) / m2;
    if(pT2corr < TINYPT2) { wt = 0.; continue; }

    // If creating heavy quark by Q -> g + Q then next need g -> Q + Qbar.
    // So minimum total mass2 is 4 * m2Sister, but use more to be safe.
    if ( isGluon && ( abs(idMother) == 4 || abs(idMother) == 5 )) {
      double m2QQsister =  EXTRASPACEQ * 4. * m2Sister;
      double pT2QQcorr = Q2 - z * (m2 + Q2) * (Q2 + m2QQsister) / m2;
      if(pT2QQcorr < TINYPT2) { wt = 0.; continue; }
    }  

    // Select phi angle of branching at random.
    phi = 2. * M_PI * Rndm::flat();

    // Evaluation of ME correction.
    if (doMEcorrections) wt *= calcMEcorr(MEtype, idMother, m2, z, Q2) 
      / calcMEmax(MEtype, idMother); 

    // Evaluation of new daughter and mother PDF's.
    double xPDFdaughterNew = max ( TINYPDF, 
      beam.xfISR(iSysNow, idDaughter, xDaughter, pT2) );
    double xPDFmotherNew = beam.xfISR(iSysNow, idMother, xMother, pT2);
    wt *= xPDFmotherNew / xPDFdaughterNew;

    // Check that valence step does not cause problem.
    if (wt > 1.) ErrorMessages::message("Warning in SpaceShower::"
      "pT2nextQCD: weight above unity"); 

  // Iterate until acceptable pT (or have fallen below pTmin).
  } while (wt < Rndm::flat()) ;

  // Save values for (so far) acceptable branching.
  dipEndNow->idMother = idMother;  
  dipEndNow->idSister = idSister;  
  dipEndNow->pT2 = pT2;  
  dipEndNow->z = z; 
  dipEndNow->Q2 = Q2;
  dipEndNow->mSister = mSister;  
  dipEndNow->m2Sister = m2Sister;  
  dipEndNow->pT2corr = pT2corr;  
  dipEndNow->phi = phi;  

}

//*********

// Evolve a QCD dipole end near threshold, with g -> Q + Qbar enforced.
// Note: No explicit Sudakov factor formalism here. Instead use that 
// df_Q(x, pT2) = (alpha_s/2pi) * (dT2/pT2) * ((gluon) * (splitting)).
// This implies that effects of Q -> Q + g are neglected in this range. 

void SpaceShower::pT2nearQCDthreshold( BeamParticle& beam, 
  double m2Massive, double m2Threshold, double zMinAbs, 
  double zMaxMassive) {

  // Read beam particle and incoming parton from appropriate side. 
  int idDaughter = (dipEndNow->side == 1) ? sysNow->id1 : sysNow->id2;
  double xDaughter = (dipEndNow->side == 1) ? sysNow->x1 : sysNow->x2;

  // Further initial values, to be used in kinematics and weighting.
  double m2 = sysNow->m2;
  double Lambda2 = (abs(idDaughter) == 4) ? Lambda4flav2 : Lambda5flav2;
  double logM2Lambda2 = log( m2Massive / Lambda2 );
  double xPDFmotherOld = beam.xfISR(iSysNow, 21, xDaughter, m2Threshold);

  // Variables used inside evolution loop. (Mainly dummy start values.)
  int loop = 0;
  double wt = 0.;
  double pT2 = 0.; 
  double z = 0.; 
  double Q2 = 0.; 
  double pT2corr = 0.;

  // Begin loop over tries to find acceptable g -> Q + Qbar branching. 
  do { 
    wt = 0.;

    // Check that not caught in infinite loop with impossible kinematics.
    if (++loop > 100) { 
      ErrorMessages::message("Error in SpaceShower::pT2nearQCDthreshold: "
        "stuck in loop"); 
      return; 
    }

    // Pick dpT2/pT2 in range [m2Massive,thresholdRatio * m2Massive]. 
    pT2 = m2Massive * pow( m2Threshold / m2Massive, Rndm::flat() ); 

    // Pick z flat in allowed range.
    z = zMinAbs + Rndm::flat() * (zMaxMassive - zMinAbs);

    // Check that kinematically possible choice.
    Q2 = pT2 / (1.-z) - m2Massive;
    pT2corr = Q2 - z * (m2 + Q2) * (Q2 + m2Massive) / m2;
    if(pT2corr < TINYPT2) continue;
    
    // Correction factor for running alpha_s.  ??
    wt = logM2Lambda2 / log( pT2 / Lambda2 ); 

    // Correction factor for splitting kernel.
    wt *= pow2(z) + pow2(1.-z) + 2. * z * (1.-z) * m2Massive / pT2;

    // Correction factor for gluon density.
    double xPDFmotherNew = beam.xfISR(iSysNow, 21, xDaughter/z, pT2);
    wt *= xPDFmotherNew / xPDFmotherOld;

  // Iterate until acceptable pT and z.
  } while (wt < Rndm::flat()) ;

  // Select phi angle of branching at random.
  double phi = 2. * M_PI * Rndm::flat();

  // Save values for (so far) acceptable branching.
  dipEndNow->idMother = 21;  
  dipEndNow->idSister = -idDaughter;  
  dipEndNow->pT2 = pT2;  
  dipEndNow->z = z; 
  dipEndNow->Q2 = Q2;
  dipEndNow->mSister = (abs(idDaughter) == 4) ? mc : mb;  
  dipEndNow->m2Sister = pow2(dipEndNow->mSister);  
  dipEndNow->pT2corr = pT2corr;  
  dipEndNow->phi = phi;  

}

//*********

// Evolve a QED dipole end. 

void SpaceShower::pT2nextQED( double pT2begDip, double pT2endDip) { 

  // Read beam particle and incoming parton from appropriate side. 
  BeamParticle& beam = (dipEndNow->side == 1) ? *beamAPtr : *beamBPtr;
  int idDaughter = (dipEndNow->side == 1) ? sysNow->id1 : sysNow->id2;
  double xDaughter = (dipEndNow->side == 1) ? sysNow->x1 : sysNow->x2;

  // Further starting values.
  int MEtype = dipEndNow->MEtype;
  bool isQuark = (abs(dipEndNow->chgType) < 3) ? true : false;
  bool isPhoton = (idDaughter == 22) ? true : false;
  double pT2 = pT2begDip;
  double m2 = sysNow->m2;
  double xMaxAbs = beam.xMax(iSysNow);
  double zMinAbs = xDaughter / xMaxAbs;

  // Determine overestimated z range
  double zMaxAbs = 1. - 0.5 * (pT2endDip / m2) *
      ( sqrt( 1. + 4. * m2 / pT2endDip ) - 1. );

  // Currently no f -> gamma branching implemented. ??
  if (isPhoton) return;

  // Integrals of splitting kernels for fermions: f -> f.
  double f2fInt = 0.;
  if (isQuark) f2fInt = (2. / 9.) * pow2(dipEndNow->chgType) 
    *  log( (1.-zMinAbs) /(1.-zMaxAbs) );
  if (doMEcorrections) f2fInt *= calcMEmax(MEtype, 1);
  double kernelPDF = alphaEM2pi * f2fInt;
  if (f2fInt < TINYKERNELPDF) return;
  
  // Variables used inside evolution loop. (Mainly dummy start values.)
  int idMother = 0;
  double z = 0.; 
  double xMother = 0.; 
  double wt = 0.; 
  double Q2 = 0.;
  double mSister = 0.; 
  double m2Sister = 0.;
  double pT2corr = 0.;
  double phi = 0.;
  
  // Begin evolution loop towards smaller pT values.
  do { 
    wt = 0.;

    // Pick pT2 (in overestimated z range).
    pT2 = pT2 * pow(Rndm::flat(), 1. / kernelPDF) ;

    // Abort evolution if below cutoff scale, or below another branching.
    if ( pT2 < pT2endDip) return; 

    // Select z value of branching q -> q + gamma, and corrective weight.
    idMother = idDaughter;
    z = 1. - (1. - zMinAbs) * pow( (1. - zMaxAbs) / (1. - zMinAbs),
      Rndm::flat() ); 
    wt = 0.5 * (1. + pow2(z));

    // Derive Q2 and x of mother from pT2 and z. 
    Q2 = pT2 / (1.- z);
    xMother = xDaughter / z;
 
    // Forbidden emission if outside allowed z range for given pT2.
    mSister = 0.;
    m2Sister = 0.;
    pT2corr = Q2 - z * (m2 + Q2) * (Q2 + m2Sister) / m2;
    if(pT2corr < TINYPT2) { wt = 0.; continue; }

    // Select phi angle of branching at random.
    phi = 2. * M_PI * Rndm::flat();

    // Evaluation of ME correction.
    if (doMEcorrections) wt *= calcMEcorr(MEtype, idMother, m2, z, Q2) 
      / calcMEmax(MEtype, idMother); 

    // Evaluation of new daughter and mother PDF's.
    double xPDFdaughterNew = max ( TINYPDF, 
      beam.xfISR(iSysNow, idDaughter, xDaughter, pT2) );
    double xPDFmotherNew = beam.xfISR(iSysNow, idMother, xMother, pT2);
    wt *= xPDFmotherNew / xPDFdaughterNew;

  // Iterate until acceptable pT (or have fallen below pTmin).
  } while (wt < Rndm::flat()) ;

  // Save values for (so far) acceptable branching.
  dipEndNow->idMother = idMother;  
  dipEndNow->idSister = 22;  
  dipEndNow->pT2 = pT2;  
  dipEndNow->z = z; 
  dipEndNow->Q2 = Q2;
  dipEndNow->mSister = mSister;  
  dipEndNow->m2Sister = m2Sister;  
  dipEndNow->pT2corr = pT2corr;  
  dipEndNow->phi = phi;  

}

//*********

// Kinematics of branching.
// Construct mother -> daughter + sister, with recoiler on other side. 

bool SpaceShower::branch( Event& event) {

  // Side on which branching occured.
  int side = dipEndSel->side;
  double sideSign = (side == 1) ? 1. : -1.;

  // Read in flavour and colour variables.
  int iDaughter = sysSel->iParton[side - 1];
  int iRecoiler = sysSel->iParton[2 - side];
  int idDaughter = (side == 1) ? sysSel->id1 : sysSel->id2;
  int idMother = dipEndSel->idMother;
  int idSister = dipEndSel->idSister;
  int colDaughter = event[iDaughter].col();
  int acolDaughter = event[iDaughter].acol();

  // Read in kinematical variables.
  double m = sysSel->m;
  double m2 = sysSel->m2;
  double x1 = sysSel->x1;
  double x2 = sysSel->x2;
  double pT2 = dipEndSel->pT2;
  double z = dipEndSel->z;
  double Q2 = dipEndSel->Q2; 
  double mSister = dipEndSel->mSister;
  double m2Sister = dipEndSel->m2Sister;
  double pT2corr = dipEndSel->pT2corr;
  double phi = dipEndSel->phi;

  // Take copy of existing system, to be given modified kinematics.
  int eventSizeOld = event.size();
  int systemSizeOld = sysSel->iParton.size();
  for ( int iCopy = 0; iCopy < systemSizeOld; ++iCopy) {
    int iOldCopy = sysSel->iParton[iCopy];
    int statusNew = (iOldCopy == iDaughter 
      || iOldCopy == iRecoiler) ? event[iOldCopy].status() : 44;
    event.copy(iOldCopy, statusNew);
  }
 
  // Define colour flow in branching.
  // Default corresponds to f -> f + gamma.
  int colMother = colDaughter;
  int acolMother = acolDaughter;
  int colSister = 0;
  int acolSister = 0; 
  if (idSister == 22) ; 
  // q -> q + g and 50% of g -> g + g; need new colour.
  else if (idSister == 21 && ( (idMother > 0 && idMother < 9)
  || (idMother == 21 && Rndm::flat() < 0.5) ) ) {  
    colMother = event.nextColTag();
    colSister = colMother;
    acolSister = colDaughter;
  // qbar -> qbar + g and other 50% of g -> g + g; need new colour.
  } else if (idSister == 21) {  
    acolMother = event.nextColTag();
    acolSister = acolMother;
    colSister = acolDaughter;
  // q -> g + q.
  } else if (idDaughter == 21 && idMother > 0) { 
    colMother = colDaughter;
    acolMother = 0;
    colSister = acolDaughter;
  // qbar -> g + qbar
  } else if (idDaughter == 21) {
    acolMother = acolDaughter;
    colMother = 0;
    acolSister = colDaughter;
  // g -> q + qbar.
  } else if (idDaughter > 0 && idDaughter < 9) {
    acolMother = event.nextColTag();
    acolSister = acolMother;
  // g -> qbar + q.
  } else if (idDaughter < 0 && idDaughter > -9) {
    colMother = event.nextColTag();
    colSister = colMother;
  // q -> gamma + q.
  } else if (idDaughter == 22 && idMother > 0) {
    colMother = event.nextColTag();
    colSister = colMother; 
   // qbar -> gamma + qbar.
  } else if (idDaughter == 22) {
    acolMother = event.nextColTag();
    acolSister = acolMother;
  }   

  // Construct kinematics of mother, sister and recoiler in old rest frame.
  double pTbranch = sqrt(pT2corr) * m2 / ( z * (m2 + Q2) );
  double pzMother = sideSign * 0.5 * m * ( (m2 - Q2) / ( z * (m2 + Q2) )
    + (Q2 + m2Sister) / m2 ); 
  double eMother = sqrt( pow2(pTbranch) + pow2(pzMother) );
  double pzSister = pzMother - sideSign * 0.5 * (m2 + Q2) / m;
  double eSister = sqrt( pow2(pTbranch) + pow2(pzSister) + m2Sister );
  double eNewRecoiler = 0.5 * (m2 + Q2) / m;
  Vec4 pMother( pTbranch, 0., pzMother, eMother );
  Vec4 pSister( pTbranch, 0., pzSister, eSister ); 
  Vec4 pNewRecoiler( 0., 0., -sideSign * eNewRecoiler, eNewRecoiler);

  // Indices of partons involved. Add new sister.
  int iMother = eventSizeOld + side - 1;
  int iNewRecoiler = eventSizeOld + 2 - side;
  int iSister = event.append( idSister, 43, iMother, 0, 0, 0,
     colSister, acolSister, pSister, mSister, sqrt(pT2) );

  // References to the partons involved.
  Particle& daughter = event[iDaughter];
  Particle& mother = event[iMother];
  Particle& newRecoiler = event[iNewRecoiler];
  Particle& sister = event.back();

  // Replace old by new mother; update new recoiler.
  mother.id( idMother );
  mother.status( -41);
  mother.cols( colMother, acolMother);
  mother.p( pMother);
  newRecoiler.status( -42);
  newRecoiler.p( pNewRecoiler);

  // Update mother and daughter pointers; also for beams.
  daughter.mothers( iMother, 0);
  mother.daughters( iSister, iDaughter); 
  event[1].daughter1( (side == 1) ? iMother : iNewRecoiler ); 
  event[2].daughter1( (side == 2) ? iMother : iNewRecoiler ); 

  // Find boost to old rest frame, and rotation -phi.
  RotBstMatrix Mtot;
  Mtot.bst(0., 0., (x2 - x1) / (x1 + x2) );
  Mtot.rot(0., -phi); 

  // Find boost from old rest frame to event cm frame.
  RotBstMatrix MfromRest;
  // The boost to the new rest frame.
  Vec4 sumNew = pMother + pNewRecoiler;
  double betaX = sumNew.px() / sumNew.e();
  double betaZ = sumNew.pz() / sumNew.e();
  MfromRest.bst( -betaX, 0., -betaZ);
  // Alignment of  radiator + recoiler to +- z axis, and rotation +phi.
  pMother.rotbst(MfromRest);  
  double theta = pMother.theta();
  if (side == 2) theta += M_PI;
  MfromRest.rot(-theta, phi); 
  // Longitudinal boost to radiator + recoiler in event cm frame.
  double x1New = (side == 1) ? x1 / z : x1;
  double x2New = (side == 2) ? x2 / z : x2;
  MfromRest.bst(0., 0., (x1New - x2New) / (x1New + x2New) );
  Mtot.rotbst(MfromRest);

  // Perform cumulative rotation/boost operation.
  // Mother, recoiler and sister from old rest frame to event cm frame.
  mother.rotbst(MfromRest);
  newRecoiler.rotbst(MfromRest);
  sister.rotbst(MfromRest);
  // The rest from (and to) event cm frame.
  for ( int i = eventSizeOld + 2; i < eventSizeOld + systemSizeOld; ++i) 
    event[i].rotbst(Mtot);  
 
  // Update list of partons in system; adding newly produced one.
  for ( int iCopy = 0; iCopy < systemSizeOld; ++iCopy) 
    sysSel->iParton[iCopy] = eventSizeOld + iCopy;
  sysSel->iParton.push_back(eventSizeOld + systemSizeOld);

  // Update info on system, including in beam remnants.
  if (side == 1) {
    sysSel->id1 = idMother;
    sysSel->x1 = x1New;
  } else {
    sysSel->id2 = idMother;
    sysSel->x2 = x2New;
  }
  sysSel->m2 = m2 / z;
  sysSel->m = sqrt(sysSel->m2);

  // Update info on dipole end.
  if (dipEndSel->colType != 0) dipEndSel->colType = mother.colType();
  if (dipEndSel->chgType != 0) dipEndSel->chgType = mother.icharge();
  dipEndSel->MEtype = 0;

  // Update info on beam remnants.
  BeamParticle& beamBranch = (side == 1) ? *beamAPtr : *beamBPtr;
  double xNew = (side == 1) ? x1New : x2New;
  beamBranch[iSysSel].lineidx( iMother, idMother, xNew);
  // Redo choice of companion kind whenever new flavour.
  if (idMother != idDaughter) {
    beamBranch.xfISR( iSysSel, idMother, xNew, pT2);
    beamBranch.pickValSeaComp();
  }
  BeamParticle& beamRecoil = (side == 1) ? *beamBPtr : *beamAPtr;
  beamRecoil[iSysSel].line( iNewRecoiler);

  // Done without any errors.
  return true;

}

//*********

// Find class of ME correction.

void SpaceShower::findMEtype( Event& event) {

  // Default values and no action.
  int MEtype = 0; 
  if (!doMEcorrections) ;

  // Identify systems producing a single resonance.
  else if (sysNow->iParton.size() == 3) {
    int idIn1 = event[sysNow->iParton[0]].id();
    int idIn2 = event[sysNow->iParton[1]].id();
    int idRes = event[sysNow->iParton[2]].id();

    // f + fbar -> vector boson. 
    if ( (idRes == 23 || abs(idRes) == 24 || idRes == 32 
      || idRes == 33 || abs(idRes) == 34 || abs(idRes) == 41)
      && abs(idIn1) < 20 && abs(idIn2) < 20 ) MEtype = 1;

    // g + g, gamma + gamma  -> Higgs boson.
    if ( (idRes == 25 || idRes == 35 || idRes == 36) 
       && ( ( idIn1 == 21 && idIn2 == 21 ) 
       || ( idIn1 == 22 && idIn2 == 22 ) ) ) MEtype = 2; 
  }

  // Set matrix element correction type. 
  for (int i = 0; i < 4; ++i) sysNow->dipEnd[i].MEtype = MEtype;

}

//*********

// Provide maximum of expected ME weight; for preweighting of evolution.

double SpaceShower::calcMEmax( int MEtype, int idMother) {

  // Currently only one non-unity case, so simplify.
  if (MEtype == 1 && idMother > 20 ) return 3.;
  return 1.;

}  

//*********

// Provide actual ME weight for current branching.

double SpaceShower::calcMEcorr(int MEtype, int idMother,
  double M2, double z, double Q2) {

  // Convert to Mandelstam variables.
  double sH = M2 / z;
  double tH = -Q2;
  double uH = Q2 - M2 * (1. - z) / z;

  // Corrections for f + fbar -> s-channel vector boson.
  if (MEtype == 1) {
    if (abs(idMother) < 20) {
      return (tH*tH + uH*uH + 2. * M2 * sH) / (sH*sH + M2*M2); 
    } else {
      return (sH*sH + uH*uH + 2. * M2 * tH) / (pow2(sH - M2) + M2*M2); 
    }

  // Corrections for g + g -> Higgs boson.
  } else if (MEtype == 2) {
    if (abs(idMother) < 20) {
      return (sH*sH + uH*uH) / (sH*sH + pow2(sH - M2)); 
    } else {
      return 0.5 * (pow4(sH) + pow4(tH) + pow4(uH) + pow4(M2)) 
        / pow2(sH*sH - M2 * (sH - M2)); 
    }    
  }

  return 1.;

}
 
//**************************************************************************

} // end namespace Pythia8
