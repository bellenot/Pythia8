// Function definitions (not found in the header) for the TimeShower class. 
// Copyright C 2006 Torbjorn Sjostrand

#include "TimeShower.h"

namespace Pythia8 {

//**************************************************************************

// The TimeShower class.

//*********
 
// Definitions of static variables.
// (Values will be overwritten in initStatic call, so are purely dummy.)

bool TimeShower::doQCDshower = true;
bool TimeShower::doQEDshowerByQ = true;
bool TimeShower::doQEDshowerByL = true;
bool TimeShower::doMEcorrections = true;
bool TimeShower::doPhiPolAsym = true;
int TimeShower::alphaSorder = 1;
int TimeShower::nQuark = 5;
double TimeShower::mc = 1.5;
double TimeShower::mb = 4.8;
double TimeShower::mc2 = 2.25;
double TimeShower::mb2 = 23.04;
double TimeShower::alphaSvalue = 0.1265;
double TimeShower::alphaS2pi = 0.02013;
double TimeShower::pTcolCutMin = 0.5;
double TimeShower::alphaEM = 0.00729735;
double TimeShower::alphaEM2pi = 1.62e-3;
double TimeShower::pTchgQCut = 0.5;
double TimeShower::pT2chgQCut = 0.25;
double TimeShower::pTchgLCut = 0.5e-3;
double TimeShower::pT2chgLCut = 0.25e-6;
double TimeShower::sin2thetaW = 0.232;
double TimeShower::mZ = 91.188;
double TimeShower::gammaZ = 2.478;

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// For small x approximate 1 - sqrt(1 - x) by x/2.
const double TimeShower::SIMPLIFYROOT = 1e-8;

// Do not allow x too close to 0 or 1 in matrix element expressions.
const double TimeShower::XMARGIN = 1e-10;

//*********

// Initialize static data members.

void TimeShower::initStatic() {

  // Main flags.
  doQCDshower = Settings::flag("TimeShower:QCDshower");
  doQEDshowerByQ = Settings::flag("TimeShower:QEDshowerByQ");
  doQEDshowerByL = Settings::flag("TimeShower:QEDshowerByL");
  doMEcorrections = Settings::flag("TimeShower:MEcorrections");
  doPhiPolAsym = Settings::flag("TimeShower:phiPolAsym"); 

  // Charm and bottom mass thresholds.
  mc = ParticleDataTable::m0(4); 
  mb = ParticleDataTable::m0(5); 
  mc2 = mc * mc;
  mb2 = mb * mb;

  // Parameters of alphaStrong generation .
  alphaSvalue = Settings::parameter("TimeShower:alphaSvalue");
  alphaSorder = Settings::mode("TimeShower:alphaSorder");
  alphaS2pi = 0.5 * alphaSvalue / M_PI;
 
  // Parameters of QCD evolution. 
  nQuark = Settings::mode("TimeShower:nQuark");
  pTcolCutMin = Settings::parameter("TimeShower:pTmin"); 
 
  // Parameters of QED evolution.
  alphaEM = Settings::parameter("StandardModel:alphaEMfix");
  alphaEM2pi = 0.5 * alphaEM / M_PI;
  pTchgQCut = Settings::parameter("TimeShower:pTminchgQ"); 
  pT2chgQCut = pow2(pTchgQCut);
  pTchgLCut = Settings::parameter("TimeShower:pTminchgL"); 
  pT2chgLCut = pow2(pTchgLCut);

  // Parameters needed for gamma/Z0 mixing.
  sin2thetaW = Settings::parameter("StandardModel:sin2thetaW");
  mZ = ParticleDataTable::m0(23);
  gammaZ = ParticleDataTable::width(23);

} 

//*********

// Initialize alphaStrong and related pTmin parameters.

void TimeShower::init() {

  // Initialize alphaStrong generation.
  alphaS.init( alphaSvalue, alphaSorder); 
  
  // Lambda for 5, 4 and 3 flavours.
  Lambda3flav = alphaS.Lambda3(); 
  Lambda4flav = alphaS.Lambda4(); 
  Lambda5flav = alphaS.Lambda5(); 
  Lambda5flav2 = pow2(Lambda5flav);
  Lambda4flav2 = pow2(Lambda4flav);
  Lambda3flav2 = pow2(Lambda3flav);
 
  // Parameters of QCD evolution. 
  pTcolCut = max( pTcolCutMin, 1.1 * Lambda3flav );
  pT2colCut = pow2(pTcolCut);  

}

//*********

// Top-level driver routine to do a single time-like shower. Used e.g. in 
// resonance decays, with showers decoupled from the rest of the event
// Input provided by range of involved partons.

void TimeShower::shower( Event& event, int iBeg, int iEnd, double pTmax) {

  // Prepare system for evolution.
  prepare( event, iBeg, iEnd);

  // Begin evolution down in pT from hard pT scale. 
  int loop =0;
  do {
    double pTtimes = pTnext( event, pTmax, 0.);
    ++loop;

    // Do a final-state emission (if allowed).
    if (pTtimes > 0) {
      branch( event); 
      pTmax = pTtimes;
    }
    
  // Keep on evolving until nothing is left to be done.
    else pTmax = 0.;
  } while (pTmax > 0.);   

}

//*********

// Prepare system for evolution; identify ME.
// Still to handle: resonance -> 3, strings to remnant, etc.

void TimeShower::prepare( Event& event, int iBegIn, int iEndIn) {

  // Reset list of radiating dipole ends.
  dipole.clear();

  // Range to be covered: provided or (default) all.
  int iBeg = (iBegIn > 0) ? iBegIn : 0;
  int iEnd = (iEndIn > 0) ? iEndIn + 1 : event.size();  

  // If only two decay products of resonance then relative recoilers.
  bool isRes2Two = (iEndIn - iBegIn == 1) ? true : false;  

  // Loop through event record to find possible dipole ends.
  for (int iRad = iBeg; iRad < iEnd; ++iRad) 
    if (event[iRad].remains()) {

    // Find dipole end formed by colour index.
    int colTag = event[iRad].col();     
    if (doQCDshower && colTag > 0) { 
      for (int iRec = iBeg; iRec < iEnd; ++iRec) {
        if (iRec != iRad && event[iRec].remains() 
        && (isRes2Two || event[iRec].acol() == colTag) ) {
          double pTmax = event[iRad].scale();
          int colType = 1;
          if (event[iRad].id() == 21) colType = 2;
          dipole.push_back( 
            TimeDipoleEnd(iRad, iRec, pTmax, colType, 0, -1) );
	}  
      }
    }

    // Find dipole end formed by anticolour index.
    int acolTag = event[iRad].acol();     
    if (doQCDshower && acolTag > 0) { 
      for (int iRec = iBeg; iRec < iEnd; ++iRec) {
        if (iRec != iRad && event[iRec].remains() 
        && (isRes2Two || event[iRec].col() == acolTag) ) {
          double pTmax = event[iRad].scale();
          int colType = -1;
          if (event[iRad].id() == 21) colType = -2;
          dipole.push_back( 
            TimeDipoleEnd(iRad, iRec, pTmax, colType, 0, -1) );
	}  
      }
    }

    // Find charge-dipole ends. Search procedure to be expanded??
    // Not correct e.g. for t -> b W!!
    int chgTag = event[iRad].icharge();     
    if (doQEDshowerByQ && chgTag != 0) { 
      for (int iRec = iBeg; iRec < iEnd; ++iRec) {
        if (iRec != iRad && event[iRec].remains() 
        && (isRes2Two || event[iRec].icharge() == -chgTag) ) {
          double pTmax = event[iRad].scale();
          dipole.push_back( 
            TimeDipoleEnd(iRad, iRec, pTmax, 0, chgTag, -1) );
	}  
      }
    }

  // End event record loop. Have now found the dipole ends.
  } 

  // Now loop through dipole ends to find matrix element corrections.
  for (int iDip = 0; iDip < int(dipole.size()); ++iDip) 
    findMEtype( event, dipole[iDip]); 

}

//*********

// Select next pT in downwards evolution of the existing dipoles.

double TimeShower::pTnext( Event& event, double pTbegAll, double pTendAll) {

  // Begin loop over all possible radiating dipole ends.
  double pT2sel = pTendAll*pTendAll;
  dipSel = 0;
  for (int iDip = 0; iDip < int(dipole.size()); ++iDip) {
    TimeDipoleEnd& dip = dipole[iDip]; 
   
    // Dipole properties. (Could partly be moved up to prepare??)
    dip.mRad = event[dip.iRadiator].m(); 
    dip.m2Rad = pow2(dip.mRad);
    dip.mRec = event[dip.iRecoiler].m(); 
    dip.m2Rec = pow2(dip.mRec);
    dip.mDip = m( event[dip.iRadiator], event[dip.iRecoiler] );
    dip.m2Dip = pow2(dip.mDip);

    // Find maximum evolution scale for dipole.
    dip.m2DipCorr = pow2(dip.mDip - dip.mRec) - dip.m2Rad; 
    double pTbegDip = min( pTbegAll, dip.pTmax ); 
    double pT2begDip = min( pow2(pTbegDip), 0.25 * dip.m2DipCorr);

    // Do QCD or QED (for q or l) evolution if it makes sense.
    dip.pT2 = 0.;
    if (pT2begDip > pT2sel) {
      if (dip.colType != 0) {
        double pT2endDip = max( pT2sel, pT2colCut );   
        if (pT2begDip > pT2endDip) pT2nextQCD(pT2begDip, pT2endDip, dip);
      } else if (dip.chgType != 0 && abs(dip.chgType) != 3) { 
        double pT2endDip = max( pT2sel, pT2chgQCut );   
        if (pT2begDip > pT2endDip) pT2nextQED(pT2begDip, pT2endDip, dip);
      } else if (abs(dip.chgType) == 3) { 
        double pT2endDip = max( pT2sel, pT2chgLCut );   
        if (pT2begDip > pT2endDip) pT2nextQED(pT2begDip, pT2endDip, dip);
      }

      // Update if found larger pT than current maximum found.
      if (dip.pT2 > pT2sel) {
        pT2sel = dip.pT2;
        dipSel = &dip;
      }
    } 
  } 

  // Return nonvanishing value if found pT bigger than already found.
  if (dipSel == 0) return 0.;
  else return sqrt(pT2sel); 
}

//*********

// Evolve a QCD dipole end. 

void TimeShower::pT2nextQCD(double pT2begDip, double pT2endDip, 
  TimeDipoleEnd& dip) { 

  // Upper estimate for matrix element weighting and colour factor.
  // Note that g -> g g and g -> q qbar are split on two sides.
  double wtPSglue = 2.;
  double colFac = (abs(dip.colType) == 1) ? 4./3. : 3./2.;
  double wtPSqqbar = (abs(dip.colType) == 2) ? 0.25 * nQuark : 0.;
  
  // Variables used inside evolution loop. (Mainly dummy start values.)
  dip.pT2 = pT2begDip;
  int nFlavour = 3;
  double zMinAbs = 0.5;
  double pT2min = pT2endDip;
  double b0 = 4.5;
  double Lambda2 = Lambda3flav2; 
  double emitCoefGlue = 0.;
  double emitCoefQqbar = 0.; 
  double emitCoefTot = 0.; 
  double wt = 0.; 
  bool mustFindRange = true;
  
  // Begin evolution loop towards smaller pT values.
  do { 

    // Initialize evolution coefficients at the beginning and
    // reinitialize when crossing c and b flavour thresholds.
    if (mustFindRange) {

      // Determine overestimated z range; switch at c and b masses.
      if (dip.pT2 > mb2) {
        nFlavour = 5;
        pT2min = mb2;
        b0 = 23./6.;
        Lambda2 = Lambda5flav2;
      } else if (dip.pT2 > mc2) {
        nFlavour = 4;
        pT2min = mc2;
        b0 = 25./6.;
        Lambda2 = Lambda4flav2;
      } else { 
        nFlavour = 3;
        pT2min = pT2endDip;
        b0 = 27./6.;
        Lambda2 = Lambda3flav2;
      }
      zMinAbs = 0.5 - sqrtpos( 0.25 - pT2min / dip.m2DipCorr );
      if (zMinAbs < SIMPLIFYROOT) zMinAbs = pT2min / dip.m2DipCorr;

      // Find emission coefficients for X -> X g and g -> q qbar.
      emitCoefGlue = wtPSglue * colFac * log(1. / zMinAbs - 1.);
      emitCoefTot = emitCoefGlue;
      if (abs(dip.colType) == 2) {
        emitCoefQqbar = wtPSqqbar * (1. - 2. * zMinAbs);
        emitCoefTot += emitCoefQqbar;
      }

      // Initialization done for current range.
      mustFindRange = false;
    } 

    // Pick pT2 (in overestimated z range) for fixed alpha_strong.
    if (alphaSorder == 0) {
      dip.pT2 = dip.pT2 * pow( Rndm::flat(), 
        1. / (alphaS2pi * emitCoefTot) );

    // Ditto for first-order alpha_strong.
    } else if (alphaSorder == 1) {
      dip.pT2 = Lambda2 * pow( dip.pT2 / Lambda2, 
        pow( Rndm::flat(), b0 / emitCoefTot) );

      // For second order reject by second term in alpha_strong expression.
    } else {
      do dip.pT2 = Lambda2 * pow( dip.pT2 / Lambda2, 
        pow( Rndm::flat(), b0 / emitCoefTot) );
      while (alphaS.alphaS2OrdCorr(dip.pT2) < Rndm::flat() 
        && dip.pT2 > pT2min);
    }
    wt = 0.;
  
    // If crossed c or b thresholds: continue evolution from threshold.
    if (nFlavour == 5 && dip.pT2 < mb2) {  
      mustFindRange = true;
      dip.pT2 = mb2;
    } else if ( nFlavour == 4 && dip.pT2 < mc2) { 
      mustFindRange = true;
      dip.pT2 = mc2;

    // Abort evolution if below cutoff scale, or below another branching.
    } else {
      if ( dip.pT2 < pT2endDip) { dip.pT2 = 0.; return; }

      // Pick kind of branching: X -> X g or g -> q qbar. 
      dip.flavour = 21;
      if (abs(dip.colType) == 2 && emitCoefQqbar > Rndm::flat() 
        * emitCoefTot) dip.flavour = 0; 

      // Pick z: either dz/(1-z) or flat dz.
      if (dip.flavour == 21) {
        dip.z = 1. - zMinAbs * pow( 1. / zMinAbs - 1., Rndm::flat() );
      } else { 
        dip.z = zMinAbs + (1. - 2. * zMinAbs) * Rndm::flat();   
      }
  
      // Do not accept branching if outside allowed z range.
      double zMin = 0.5 - sqrtpos( 0.25 - dip.pT2 / dip.m2DipCorr ); 
      if (zMin < SIMPLIFYROOT) zMin = dip.pT2 / dip.m2DipCorr;
      dip.m2 = dip.m2Rad + dip.pT2 / (dip.z * (1. - dip.z));
      if (dip.z > zMin && dip.z < 1. - zMin 
        && dip.m2 * dip.m2Dip < dip.z * (1. - dip.z) 
        * pow2(dip.m2Dip + dip.m2 - dip.m2Rec)) {

        // No z weight if to do matrix element corrections later on.
        if (dip.MEtype > 0) { 
          wt = 1.;

        // z weight for X -> Xg.
        } else if (dip.flavour == 21 && abs(dip.colType) == 1) {
          wt = (1. + pow2(dip.z)) / wtPSglue;
	} else if (dip.flavour == 21) {     
          wt = (1. + pow3(dip.z)) / wtPSglue;
           
        // Flavour choice and z weight for g -> q qbar.
        } else {
          dip.flavour = min(5, 1 + int(nQuark * Rndm::flat())); 
          dip.mFlavour = ParticleDataTable::m0(dip.flavour);
          double beta = sqrtpos( 1. - 4. * pow2(dip.mFlavour) / dip.m2 );
          wt = beta * ( pow2(dip.z) + pow2(1. - dip.z) );
        }
      }
    }

  // Iterate until acceptable pT (or have fallen below pTmin).
  } while (wt < Rndm::flat()) ;
}

//*********

// Evolve a QED dipole end. 

void TimeShower::pT2nextQED(double pT2begDip, double pT2endDip, 
  TimeDipoleEnd& dip) { 

  // Unique choice flavour.
  dip.flavour = 22;

  // Upper estimate for matrix element weighting and charge factor.
  double wtPSga = 2.;
  double chg = dip.chgType / 3.;

  // Determine overestimated z range. Find evolution coefficient.
  double zMinAbs = 0.5 - sqrtpos( 0.25 - pT2endDip / dip.m2DipCorr );
  if (zMinAbs < SIMPLIFYROOT) zMinAbs = pT2endDip / dip.m2DipCorr;
  double emitCoefTot = alphaEM2pi * pow2(chg) * wtPSga 
    * log(1. / zMinAbs - 1.);
  
  // Variables used inside evolution loop.
  dip.pT2 = pT2begDip;
  double wt; 
  
  // Begin evolution loop towards smaller pT values.
  do { 
 
    // Pick pT2 (in overestimated z range).
    dip.pT2 = dip.pT2 * pow(Rndm::flat(), 1. / emitCoefTot);
    wt = 0.;

    // Abort evolution if below cutoff scale, or below another branching.
    if ( dip.pT2 < pT2endDip) { dip.pT2 = 0.; return; }

    // Pick z according to dz/(1-z).
    dip.z = 1. - zMinAbs * pow( 1. / zMinAbs - 1., Rndm::flat() );
  
    // Do not accept branching if outside allowed z range.
    double zMin = 0.5 - sqrtpos( 0.25 - dip.pT2 / dip.m2DipCorr ); 
    if (zMin < SIMPLIFYROOT) zMin = dip.pT2 / dip.m2DipCorr;
    dip.m2 = dip.m2Rad + dip.pT2 / (dip.z * (1. - dip.z));
    if (dip.z > zMin && dip.z < 1. - zMin 
      && dip.m2 * dip.m2Dip < dip.z * (1. - dip.z) 
      * pow2(dip.m2Dip + dip.m2 - dip.m2Rec)) {
                      
      // No z weight if to do matrix element corrections later on.
      if (dip.MEtype > 0) { 
        wt = 1.;

      // z weight for X -> Xg.
      } else {
        wt = (1. + pow2(dip.z)) / wtPSga;
      }
    }

  // Iterate until acceptable pT (or have fallen below pTmin).
  } while (wt < Rndm::flat()) ;  
}

//*********

// ME corrections and kinematics that may give failure.
// Notation: radBef, recBef = radiator, recoiler before emission,
//           rad, rec, emt = radiator, recoiler, emitted efter emission.
//           (rad, emt distinguished by colour flow for g -> q qbar.) 

bool TimeShower::branch( Event& event) {

  // Find initial particles in dipole branching.
  int iRadBef = dipSel->iRadiator;
  int iRecBef = dipSel->iRecoiler;
  Particle& radBef = event[iRadBef]; 
  Particle& recBef = event[iRecBef];

  // Default flavours and colour tags for new particles in dipole branching. 
  int idRad = radBef.id();
  int idEmt = dipSel->flavour; 
  int colRad = radBef.col();
  int acolRad = radBef.acol();
  int colEmt = 0;
  int acolEmt = 0;

  // Default OK for photon emission.
  if (dipSel->flavour == 22) { 
  // New colour tag required for gluon emission.
  } else if (dipSel->flavour == 21 && dipSel->colType > 0) { 
    colEmt = colRad;  
    colRad = event.nextColTag();   
    acolEmt = colRad;
  } else if (dipSel->flavour == 21) { 
    acolEmt = acolRad;  
    acolRad = event.nextColTag();   
    colEmt = acolRad;
  // New flavours for g -> q qbar; split colours.
  } else if (dipSel->colType > 0) {
    idEmt = dipSel->flavour ;
    idRad = -idEmt;
    colEmt = colRad;
    colRad = 0; 
  } else {
    idEmt = -dipSel->flavour ;
    idRad = -idEmt;
    acolEmt = acolRad;
    acolRad = 0; 
  }

  // Construct kinematics in dipole rest frame: 
  // begin simple (like g -> g g).
  double eRadPlusEmt = 0.5 * (dipSel->m2Dip + dipSel->m2 - dipSel->m2Rec) 
    / dipSel->mDip;
  double e2RadPlusEmt = pow2(eRadPlusEmt);
  double pzRadPlusEmt = 0.5 * sqrtpos( pow2(dipSel->m2Dip - dipSel->m2 
    - dipSel->m2Rec) - 4. * dipSel->m2 * dipSel->m2Rec ) / dipSel->mDip;
  double pT2corr = dipSel->m2 * (e2RadPlusEmt * dipSel->z * (1. - dipSel->z)
    - 0.25 * dipSel->m2) / pow2(pzRadPlusEmt);
  double pTcorr = sqrtpos( pT2corr );
  double pzRad = (e2RadPlusEmt * dipSel->z - 0.5 * dipSel->m2) 
    / pzRadPlusEmt;
  double pzEmt = (e2RadPlusEmt * (1. - dipSel->z) - 0.5 * dipSel->m2) 
    / pzRadPlusEmt;
  double mRad = dipSel->mRad;
  double mEmt = 0.;

  // Kinematics reduction for q -> q g or q -> q gamma when m_q > 0. 
  if (abs(dipSel->colType) == 1 || dipSel->chgType != 0) { 
    pTcorr *= 1. - dipSel->m2Rad / dipSel->m2; 
    pzRad += pzEmt * dipSel->m2Rad / dipSel->m2;
    pzEmt *= 1. - dipSel->m2Rad / dipSel->m2;  
  // Kinematics reduction for g -> q qbar when m_q > 0;
  } else if (abs(dipSel->flavour) < 20) {
    mEmt = ParticleDataTable::m0(dipSel->flavour);
    mRad = mEmt;
    double beta = sqrtpos( 1. - 4. * pow2(mEmt) / dipSel->m2 );   
    pTcorr *= beta;
    pzRad = 0.5 * ( (1. + beta) * pzRad + (1. - beta) * pzEmt );
    pzEmt = pzRadPlusEmt - pzRad;
  } 

  // Find rest frame and angles of original dipole.
  RotBstMatrix M;
  M.fromCMframe(radBef.p(), recBef.p());

  // Evaluate coefficient of azimuthal asymmetry from gluon polarization.
  findAsymPol( event, dipSel);

  // Begin construction of new dipole kinematics: pick azimuthal angle.
  Vec4 pRad, pEmt, pRec;
  double wtPhi = 1.;
  do { 
    double phi = 2. * M_PI * Rndm::flat();

    // Define kinematics of branching in dipole rest frame.
    pRad = Vec4( pTcorr * cos(phi), pTcorr * sin(phi), pzRad, 
      sqrt( pow2(pTcorr) + pow2(pzRad) + pow2(mRad) ) );
    pEmt = Vec4( -pRad.px(), -pRad.py(), pzEmt,
      sqrt( pow2(pTcorr) + pow2(pzEmt) + pow2(mEmt) ) );
    pRec = Vec4( 0., 0., -pzRadPlusEmt, sqrt( pow2(pzRadPlusEmt) 
      + dipSel->m2Rec ) );

    // Rotate and boost dipole products to the event frame.
    pRad.rotbst(M);
    pEmt.rotbst(M);
    pRec.rotbst(M);

    // Azimuthal phi weighting: loop to new phi value if required.
    if (dipSel->asymPol != 0.) {
      Vec4 pRadBef = event[iRadBef].p();
      Vec4 pAunt = event[dipSel->iAunt].p();
      double cosPhi = cosphi( pRad, pAunt, pRadBef );
      wtPhi = ( 1. + dipSel->asymPol * (2. * pow2(cosPhi) - 1.) )
        / ( 1. + abs(dipSel->asymPol) );
    } 
  } while (wtPhi < Rndm::flat()) ;

  // Define new particles from dipole branching.
  double pTsel = sqrt(dipSel->pT2);
  Particle rad = Particle(idRad, 51, iRadBef, 0, 0, 0, 
    colRad, acolRad, pRad, mRad, pTsel); 
  Particle emt = Particle(idEmt, 51, iRadBef, 0, 0, 0,
    colEmt, acolEmt, pEmt, mEmt, pTsel);
  Particle rec = Particle(recBef.id(), 52, iRecBef, iRecBef, 0, 0, 
    recBef.col(), recBef.acol(), pRec, dipSel->mRec, pTsel); 

  // ME corrections can lead to branching being rejected.
  if (dipSel->MEtype > 0) {
    Particle& partner = (dipSel->iMEpartner == iRecBef) 
      ? rec : event[dipSel->iMEpartner];
    if ( findMEcorr( dipSel, rad, partner, emt) < Rndm::flat() ) 
      return false;
  }

  // Put new particles into the event record.
  int iRad = event.append(rad);
  int iEmt = event.append(emt);
  int iRec = event.append(rec);

  // Mark original dipole partons as branched and set daughters.
  event[dipSel->iRadiator].statusNeg();
  event[dipSel->iRadiator].daughters( iRad, iEmt); 
  event[dipSel->iRecoiler].statusNeg();
  event[dipSel->iRecoiler].daughters( iRec, iRec);
  
  // Photon emission: change info to reflect new dipole ends.
  if (dipSel->flavour == 22) { 
    dipSel->iRadiator = iRad;
    dipSel->iRecoiler = iRec;
    dipSel->pTmax = pTsel;
 
  // Gluon emission: update both dipole ends and add two new ones.
  } else if (dipSel->flavour == 21) { 
    dipSel->iRadiator = iRad;
    dipSel->iRecoiler = iEmt;
    dipSel->pTmax = pTsel;
    for (int i = 0; i < int(dipole.size()); ++i) {
      if (dipole[i].iRadiator == iRecBef && dipole[i].iRecoiler == iRadBef 
        && dipole[i].colType != 0) {
        dipole[i].iRadiator = iRec;
        dipole[i].iRecoiler = iEmt;
        dipole[i].pTmax = pTsel;
      }
    }
    int colType = (dipSel->colType > 0) ? 2 : -2 ;
    dipole.push_back( TimeDipoleEnd(iEmt, iRec, pTsel, colType, 0, 0));
    dipole.push_back( TimeDipoleEnd(iEmt, iRad, pTsel, -colType, 0, 0));

  // Gluon branching to q qbar: update current dipole and other of gluon.
  } else { 
    for (int i = 0; i < int(dipole.size()); ++i) {
      if (dipole[i].iRadiator == iRadBef && abs(dipole[i].colType) == 2) {
        dipole[i].colType /= 2;
        dipole[i].MEtype = 6;
        if ( &dipole[i] == dipSel ) dipole[i].iMEpartner = iRad;
        else dipole[i].iMEpartner = iEmt;
      }
    }
    dipSel->iRadiator = iEmt;
    dipSel->iRecoiler = iRec;
    dipSel->pTmax = pTsel;
  }

  // Now update other dipoles that also involved the radiator or recoiler.
  for (int i = 0; i < int(dipole.size()); ++i) {
    if (dipole[i].iRadiator == iRadBef) dipole[i].iRadiator = iRad;
    if (dipole[i].iRadiator == iRecBef) dipole[i].iRadiator = iRec;
    if (dipole[i].iRecoiler == iRadBef) dipole[i].iRecoiler = iRad;
    if (dipole[i].iRecoiler == iRecBef) dipole[i].iRecoiler = iRec;
    if (dipole[i].iMEpartner == iRadBef) dipole[i].iMEpartner = iRad;
    if (dipole[i].iMEpartner == iRecBef) dipole[i].iMEpartner = iRec;
  }

  // Done. 
  return true;

}

//*********

//void TimeShower::update( Event& event) {
// ;
//}

//*********

// Find class of QCD ME correction.
// MEtype classification follow codes in Norrbin article,
// additionally -1 = try to find type, 0 = no ME corrections.
// Warning: not yet tried out to do a correct assignment in 
// arbitrary multiparton configurations! ??

void TimeShower::findMEtype( Event& event, TimeDipoleEnd& dip) {

  // Return if no ME corrections to be applied.
  if (!doMEcorrections) { 
    dip.MEtype = 0;
    return;
  } 

  // Temporary solution: kill ME corrections for emissions off a gluon.??
  if (event[dip.iRadiator].id() == 21) {
    dip.MEtype = 0;
    return;
  } 

  // If no ME partner set, assume it is the recoiler.
  if (dip.iMEpartner < 0) dip.iMEpartner = dip.iRecoiler;

  // Now begin processing of colour dipole.
  if (dip.colType != 0) {

    // Find daughter types (may or may not be used later on).
    int dau1Type = findMEparticle(event[dip.iRadiator].id());
    int dau2Type = findMEparticle(event[dip.iMEpartner].id());
    int minDauType = min(dau1Type, dau2Type);
    int maxDauType = max(dau1Type, dau2Type);
    dip.MEorder = (dau2Type >= dau1Type) ? true : false;
    dip.MEsplit = (maxDauType <= 2 || maxDauType == 6) ? true : false; 
    dip.MEgluinoDau = (maxDauType == 6) ? true : false;
 
    // If type already set (or set not to have) then done.
    if (minDauType == 0 && dip.MEtype < 0) dip.MEtype = 0;
    if (dip.MEtype >= 0) return;
    dip.MEtype = 0;

    // Find mother type.
    int iMother = event[dip.iRadiator].mother1();
    int idMother = 0;
    if ( event[dip.iRecoiler].mother1() == iMother && iMother >= 0) 
      idMother = event[iMother].id();
    int motherType = (idMother != 0) ? findMEparticle(idMother) : 0;

    // Now start from default, which is no ME corrections, 
    // and try to find matching ME cases below.
    int MEkind = 0;
    int MEcombi = 4;
    dip.MEmix = 0.5;

    // Vector/axial vector -> q + qbar; q -> q + V.
    if (minDauType == 1 && maxDauType == 1 && 
      (motherType == 3 || motherType == 0) ) {
      MEkind = 2;
      if (idMother == 21 || idMother == 22) MEcombi = 1;
      else if (idMother == 23 || idMother == 0) {MEcombi = 3; 
        dip.MEmix = gammaZmix( event, iMother, dip.iRadiator, 
          dip.iRecoiler );}
      else if (idMother == 24) MEcombi = 4;
    }
    // For chi -> chi q qbar, use V/A -> q qbar as first approximation.??
    else if (minDauType == 1 && maxDauType == 1 && motherType == 5)
      MEkind =2;
    else if (minDauType == 1 && maxDauType == 3 && (motherType == 0
      || motherType == 1)) MEkind = 3;
 
    // Scalar/pseudoscalar -> q + qbar; q -> q + S.
    else if (minDauType == 1 && maxDauType == 1 && motherType == 4) {
      MEkind =4;
      if (idMother == 25 || idMother == 35 || idMother == 37) MEcombi = 1;
      else if (idMother == 36) MEcombi = 2;
    } 
    else if (minDauType == 1 && maxDauType == 4 && 
      (motherType == 0 || motherType == 1) ) MEkind = 5;
 
    // V -> ~q + ~qbar; ~q -> ~q + V; S -> ~q + ~qbar; ~q -> ~q + S.
    else if (minDauType == 2 && maxDauType == 2 && 
      (motherType == 0 || motherType == 3) ) MEkind = 6;
    else if (minDauType == 2 && maxDauType == 3 && 
      (motherType == 0 || motherType == 2) ) MEkind = 7;
    else if (minDauType == 2 && maxDauType == 2 && motherType == 4)
      MEkind = 8;
    else if (minDauType == 2 && maxDauType == 4 && 
      (motherType == 0 || motherType == 2) ) MEkind = 9;
 
    // chi -> q + ~qbar; ~q -> q + chi; q -> ~q + chi.
    else if (minDauType == 1 && maxDauType == 2 && 
      (motherType == 0 || motherType == 5) ) MEkind = 10;
    else if (minDauType == 1 && maxDauType == 5 && 
      (motherType == 0 || motherType == 2) ) MEkind = 11;
    else if (minDauType == 2 && maxDauType == 5 && 
      (motherType == 0 || motherType == 1) ) MEkind = 12;
 
    // ~g -> q + ~qbar; ~q -> q + ~g; q -> ~q + ~g.
    else if (minDauType == 1 && maxDauType == 2 && motherType == 6)
      MEkind = 13;
    else if (minDauType == 1 && maxDauType == 6 && 
      (motherType == 0 || motherType == 2) ) MEkind = 14;
    else if (minDauType == 2 && maxDauType == 6 && 
      (motherType == 0 || motherType == 1) ) MEkind = 15;

    // g -> ~g + ~g (eikonal approximation).
    else if (minDauType == 6 && maxDauType == 6 && motherType == 0)
      MEkind = 16;

    // Save ME type and gamma_5 admixture. 
    dip.MEtype = 5 * MEkind + MEcombi; 

  // Now begin processing of charge dipole - still primitive.
  } else if (dip.chgType != 0) {

    // Set defaults for QED case; then possibly done.
    dip.MEorder = true;
    dip.MEsplit = true; 
    dip.MEgluinoDau = false;
    if (dip.MEtype >= 0) return;

    // So far only ME corrections for q qbar or l lbar.
    int idDau1 = event[dip.iRadiator].id();
    int idDau2 = event[dip.iMEpartner].id();
    if (abs(idDau1) < 9 && abs(idDau2) < 9 && idDau1 * idDau2 < 0) ;
    else if (abs(idDau1) > 10 && abs(idDau1) < 19 && abs(idDau2) > 10
      && abs(idDau2) < 19 && idDau1 * idDau2 < 0) ;
    else { dip.MEtype = 0; return; }

    // Distinguish charge sum != 0 or = 0; in latter assume vector source.
    dip.MEtype = 101;
    if (idDau1 + idDau2 == 0) dip.MEtype = 102; 
    dip.MEmix = 1.;
  }

}

//*********
 
// Find type of particle for ME type: 0 = unknown, 1 = quark,
// 2 = squark, 3 = vector boson (also g), 4 = colourless scalar,
// 5 = colourless neutralino/chargino, 6 = gluino.

int TimeShower::findMEparticle( int id) {

  // find colour and spin of particle.
  int type = 0;
  int colType = ParticleDataTable::colType(id); 
  int spinType = ParticleDataTable::spinType(id);

  // Find particle type from colour and spin.
       if (colType == 1 && spinType == 2) type = 1;
  else if (colType == 1 && spinType == 1) type = 2;
  else if (colType == 0 && spinType == 3) type = 3;
  else if (colType == 2 && spinType == 3) type = 3;
  else if (colType == 0 && spinType == 1) type = 4;
  else if (colType == 0 && spinType == 2) type = 5;
  else if (colType == 2 && spinType == 2) type = 6;
  return type;
}  

//*********

// Find mixture of V and A in gamma/Z: energy- and flavour-dependent. 

double TimeShower::gammaZmix( Event& event, int iRes, int iDau1, int iDau2) {

  // Try to identify initial flavours; use e+e- as default.
  int idIn1 = -11;
  int idIn2 = 11;
  int iIn1 = (iRes >= 0) ? event[iRes].mother1() : -1;
  int iIn2 = (iRes >= 0) ? event[iRes].mother2() : -1;
  if (iIn1 >=0) idIn1 = event[iIn1].id();
  if (iIn2 >=0) idIn2 = event[iIn1].id();
         
  // In processes f + g/gamma -> f + Z only need find one fermion.
  if (idIn1 == 21 || idIn1 == 22) idIn1 = -idIn2;
  if (idIn2 == 21 || idIn2 == 22) idIn2 = -idIn1;
 
  // Initial flavours and couplings; return if don't make sense.
  if (idIn1 + idIn2 != 0 ) return 0.5;
  int idIn = abs(idIn1);
  if (idIn == 0 || idIn >18 ) return 0.5; 
  double eIn = ParticleDataTable::charge(idIn);
  double aIn = (eIn < -0.1) ? -1. : 1.;
  double vIn = aIn - 4. * sin2thetaW * eIn; 

  // Final flavours and couplings; return if don't make sense.
  if (event[iDau1].id() + event[iDau2].id() != 0) return 0.5;
  int idOut = abs(event[iDau1].id());
  if (idOut == 0 || idOut >18 ) return 0.5; 
  double eOut = ParticleDataTable::charge(idOut);
  double aOut = (eOut < -0.1) ? -1. : 1.;
  double vOut = aOut - 4. * sin2thetaW * eOut; 

  // Calculate gamma and Z0 propagators from kinematics.
  Vec4 psum = event[iDau1].p() + event[iDau2].p();
  double sHat = psum.m2Calc();
  double prop = 1. / ( pow2(sHat - mZ*mZ) + sHat * gammaZ*gammaZ); 
  double xwc = 1. / (16. * sin2thetaW * (1.-sin2thetaW));

  // Calculate vector and axial expressions and find mix.
  double vect = eIn*eIn * eOut*eOut 
    + 2. * eIn*vIn * eOut*vOut * xwc * sHat * (sHat - mZ*mZ) * prop 
    + (vIn*vIn + aIn*aIn) * vOut*vOut * xwc*xwc * sHat*sHat * prop;
  double axiv = (vIn*vIn + aIn*aIn) * aOut*aOut * xwc*xwc 
    * sHat*sHat * prop;
  return vect / (vect + axiv);
}

//*********

// Set up to calculate QCD ME correction with calcMEcorr.
// Normally for primary particles, but also from g/gamma -> f fbar.
  
double TimeShower::findMEcorr(TimeDipoleEnd* dip, Particle& rad, 
  Particle& partner, Particle& emt) {
  
  // Initial values and matrix element kind.
  double wtME = 1.;
  double wtPS = 1.; 
  int MEkind = dip->MEtype / 5;
  int MEcombi = dip->MEtype % 5;

  // Construct ME variables.
  Vec4 sum = rad.p() + partner.p() + emt.p();
  double eCM = sum.mCalc();
  double x1 = 2. * (sum * rad.p()) / pow2(eCM);
  double x2 = 2. * (sum * partner.p()) / pow2(eCM); 
  double r1 = rad.m() / eCM;
  double r2 = partner.m() / eCM; 

  // Derived ME variables, suitably protected.
  double x1minus = max(XMARGIN, 1. + r1*r1 - r2*r2 - x1);
  double x2minus = max(XMARGIN, 1. + r2*r2 - r1*r1 - x2) ;
  double x3 = max(XMARGIN, 2. - x1 - x2);

  // Begin processing of QCD dipoles.
  if (dip->colType !=0) {

    // Evaluate normal ME, for proper order of particles.
    if (dip->MEorder) 
         wtME = calcMEcorr(MEkind, MEcombi, dip->MEmix, x1, x2, r1, r2);
    else wtME = calcMEcorr(MEkind, MEcombi, dip->MEmix, x2, x1, r2, r1);

    // Split up total ME when two radiating particles.
    if (dip->MEsplit) wtME = wtME * x1minus / x3; 

    // Evaluate shower rate to be compared with.
    wtPS = 2. / ( x3 * x2minus );
    if (dip->MEgluinoDau) wtPS *= 9./4.;
  
  // For generic charge combination currently only massless expression.
  // (Masses included only to respect phase space boundaries.)
  } else if (dip->chgType !=0 && dip->MEtype == 101) {
    double chg1 = ParticleDataTable::charge(rad.id());
    double chg2 = ParticleDataTable::charge(partner.id());
    wtME = (x1*x1 + x2*x2) * pow2( chg1 * x1minus / x3 
      - chg2 * x2minus / x3 );
    wtPS = 2. * ( chg1*chg1 * x1minus / x3 + chg2*chg2 * x2minus / x3 ); 

  // For flavour neutral system assume vector source and include masses.
  } else if (dip->chgType !=0 && dip->MEtype == 102) {
    wtME = calcMEcorr(2, 1, dip->MEmix, x1, x2, r1, r2) * x1minus / x3;
    wtPS = 2. / ( x3 * x2minus );
  }
       
  // Return ratio of actual ME to assumed PS rate of emission.
  return wtME / wtPS; 
}

//*********

// Matrix elements for gluon (or photon) emission from
// a two-body state; to be used by the parton shower routine.
// Here x_i = 2 E_i/E_cm, r_i = m_i/E_cm and
// 1/sigma_0 d(sigma)/d(x_1)d(x_2) = (alpha-strong/2 pi) * C_F * (this),
// i.e. normalization is such that one recovers the familiar
// (x_1^2 + x_2^2)/((1-x_1)*(1-x_2)) for the massless case.
// Coupling structure:
// kind =  1 : eikonal soft-gluon expression (spin-independent)
//      =  2 : V -> q qbar (V = vector/axial vector colour singlet)
//      =  3 : q -> q V
//      =  4 : S -> q qbar (S = scalar/pseudoscalar colour singlet)
//      =  5 : q -> q S
//      =  6 : V -> ~q ~qbar (~q = squark)
//      =  7 : ~q -> ~q V
//      =  8 : S -> ~q ~qbar
//      =  9 : ~q -> ~q S
//      = 10 : chi -> q ~qbar (chi = neutralino/chargino)
//      = 11 : ~q -> q chi
//      = 12 : q -> ~q chi
//      = 13 : ~g -> q ~qbar
//      = 14 : ~q -> q ~g
//      = 15 : q -> ~q ~g
//      = 16 : (9/4)*(eikonal) for gg -> ~g ~g
// Note that the order of the decay products is important.
// combi = 1 : pure non-gamma5, i.e. vector/scalar/...
//       = 2 : pure gamma5, i.e. axial vector/pseudoscalar/....
//       = 3 : mixture mix*(combi=1) + (1-mix)*(combi=2)
//       = 4 : mixture (combi=1) +- (combi=2)

double TimeShower::calcMEcorr( int kind, int combiIn, double mixIn, 
  double x1, double x2, double r1, double r2) {

  // Frequent variable combinations.
  double x3 = 2. - x1 - x2;
  double x1s = x1 * x1;
  double x2s = x2 * x2;
  double x3s = x3 * x3;
  double x1c = x1 * x1s;
  double x2c = x2 * x2s;
  double x3c = x3 * x3s;
  double r1s = r1 * r1;
  double r2s = r2 * r2;
  double r1c = r1 * r1s;
  double r2c = r2 * r2s;
  double r1q = r1s * r1s;
  double r2q = r2s * r2s;
  double prop1 = 1. + r1s - r2s - x1; 
  double prop2 = 1. + r2s - r1s - x2;
  double prop1s = prop1 * prop1;
  double prop2s = prop2 * prop2;
  double prop12 = prop1 * prop2;
  double prop13 = prop1 * x3;
  double prop23 = prop2 * x3;

  // Check input values. Return zero outside allowed phase space.
  if (x1 - 2.*r1 < XMARGIN || prop1 < XMARGIN) return 0.;
  if (x2 - 2.*r2 < XMARGIN || prop2 < XMARGIN) return 0.;
  if (x1 + x2 - 1. - pow2(r1+r2) < XMARGIN) return 0.;
  if ((x1s - 4.*r1s) * (x2s - 4.*r2s) 
    - pow2( 2. * (1. - x1 - x2 + r1s + r2s) + x1*x2 ) < XMARGIN) return 0.;

  // Initial values; phase space.
  int combi = max(1, min(4, combiIn) ); 
  double mix = max(0., min(1., mixIn) );
  bool isSet1 = false;
  bool isSet2 = false;
  bool isSet4 = false;
  double ps = sqrtpos( pow2(1. - r1*r1 - r2*r2) - pow2(2. * r1 * r2) );
  double rLO = 0., rFO = 0., rLO1 = 0., rFO1 = 0., rLO2 = 0., 
    rFO2 = 0., rLO4 = 0., rFO4 = 0.;
  double offset = 0;
 
  // Select which kind of ME to use.
  switch (kind) {

    // case 1 is equal to default, i.e. eikonal expression.

    // V -> q qbar (V = gamma*/Z0/W+-/...).
    case 2:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(2.-r1s-r1q+6.*r1*r2-r2s+2.*r1s*r2s-r2q)/2.;
        rFO1 = -(3.+6.*r1s+r1q-6.*r1*r2+6.*r1c*r2-2.*r2s-6.*r1s*r2s
        +6.*r1*r2c+r2q-3.*x1+6.*r1*r2*x1+2.*r2s*x1+x1s-2.*r1s*x1s
        +3.*r1s*x3+6.*r1*r2*x3-r2s*x3-2.*x1*x3-5.*r1s*x1*x3
        +r2s*x1*x3+x1s*x3-3.*x3s-3.*r1s*x3s+r2s*x3s
        +2.*x1*x3s+x3c-x2)
        /prop2s
        -2.*(-3.+r1s-6.*r1*r2+6.*r1c*r2+3.*r2s-4.*r1s*r2s
        +6.*r1*r2c+2.*x1+3.*r1s*x1+r2s*x1-x1s-r1s*x1s
        -r2s*x1s+4.*x3+2.*r1s*x3+3.*r1*r2*x3-r2s*x3-3.*x1*x3
        -2.*r1s*x1*x3+x1s*x3-x3s-r1s*x3s+r1*r2*x3s+x1*x3s)
        /prop12
        -(-1.+2.*r1s+r1q+6.*r1*r2+6.*r1c*r2-2.*r2s-6.*r1s*r2s
        +6.*r1*r2c+r2q-x1-2.*r1s*x1-6.*r1*r2*x1+8.*r2s*x1+x1s
        -2.*r2s*x1s-r1s*x3+r2s*x3-r1s*x1*x3+r2s*x1*x3+x1s*x3+x2)
        /prop1s;
        rFO1 = rFO1/2.;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(2.-r1s-r1q-6.*r1*r2-r2s+2.*r1s*r2s-r2q)/2.;
        rFO2 = -(3.+6.*r1s+r1q+6.*r1*r2-6.*r1c*r2-2.*r2s-6.*r1s*r2s    
        -6.*r1*r2c+r2q-3.*x1-6.*r1*r2*x1+2.*r2s*x1+x1s-2.*r1s*x1s
        +3.*r1s*x3-6.*r1*r2*x3-r2s*x3-2.*x1*x3-5.*r1s*x1*x3
        +r2s*x1*x3+x1s*x3-3.*x3s-3.*r1s*x3s+r2s*x3s+2.*x1*x3s+x3c-x2)
        /prop2s
        -2.*(-3+r1s+6.*r1*r2-6.*r1c*r2+3.*r2s-4.*r1s*r2s-6.*r1*r2c
        +2.*x1+3.*r1s*x1+r2s*x1-x1s-r1s*x1s-r2s*x1s+4.*x3+2.*r1s*x3
        -3.*r1*r2*x3-r2s*x3-3.*x1*x3-2.*r1s*x1*x3+x1s*x3-x3s-r1s*x3s
        -r1*r2*x3s+x1*x3s)
        /prop12
        -(-1.+2.*r1s+r1q-6.*r1*r2-6.*r1c*r2-2.*r2s-6.*r1s*r2s
        -6.*r1*r2c+r2q-x1-2.*r1s*x1+6.*r1*r2*x1+8.*r2s*x1+x1s
        -2.*r2s*x1s-r1s*x3+r2s*x3-r1s*x1*x3+r2s*x1*x3+x1s*x3+x2)
        /prop1s;
        rFO2 = rFO2/2.;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(2.-r1s-r1q-r2s+2.*r1s*r2s-r2q)/2.;
        rFO4 = (1.-r1q+6.*r1s*r2s-r2q+x1+3.*r1s*x1-9.*r2s*x1-3.*x1s
        -r1s*x1s+3.*r2s*x1s+x1c-x2-r1s*x2+r2s*x2-r1s*x1*x2+r2s*x1*x2
        +x1s*x2)
        /prop1s 
        -2.*(1.+r1s+r2s-4.*r1s*r2s+r1s*x1+2.*r2s*x1-x1s-r2s*x1s
        +2.*r1s*x2+r2s*x2-3.*x1*x2+x1s*x2-x2s-r1s*x2s+x1*x2s)
        /prop12
        +(1.-r1q+6.*r1s*r2s-r2q-x1+r1s*x1-r2s*x1+x2-9.*r1s*x2
        +3.*r2s*x2+r1s*x1*x2-r2s*x1*x2-3.*x2s+3.*r1s*x2s-r2s*x2s
        +x1*x2s+x2c)
        /prop2s;
        rFO4 = rFO4/2.;
        isSet4 = true;
      }
      break; 
 
    // q -> q V.
    case 3:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-2.*r1s+r1q+r2s-6.*r1*r2s+r1s*r2s-2.*r2q);
        rFO1 = -2.*(-1.+r1-2.*r1s+2.*r1c-r1q+pow5(r1)-r2s+r1*r2s
        -5.*r1s*r2s+r1c*r2s-2.*r1*r2q+2.*x1-2.*r1*x1+2.*r1s*x1
        -2.*r1c*x1+2.*r2s*x1+5.*r1*r2s*x1+r1s*r2s*x1+2.*r2q*x1
        -x1s+r1*x1s-r2s*x1s+3.*x2+4.*r1s*x2+r1q*x2+2.*r2s*x2
        +2.*r1s*r2s*x2-4.*x1*x2-2.*r1s*x1*x2-r2s*x1*x2+x1s*x2
        -2.*x2s-2.*r1s*x2s+x1*x2s)
        /prop23
        +(2.*r2s+6.*r1*r2s-6.*r1s*r2s+6.*r1c*r2s+2.*r2q+6.*r1*r2q
        -r2s*x1+r1s*r2s*x1-r2q*x1+x2-r1q*x2-3.*r2s*x2-6.*r1*r2s*x2
        +9.*r1s*r2s*x2-2.*r2q*x2-x1*x2+r1s*x1*x2-x2s-3.*r1s*x2s
        +2.*r2s*x2s+x1*x2s)
        /prop2s
        +(-4.-8.*r1s-4.*r1q+4.*r2s-4.*r1s*r2s+8.*r2q+9.*x1+10.*r1s*x1
        +r1q*x1-3.*r2s*x1+6.*r1*r2s*x1+r1s*r2s*x1-2.*r2q*x1-6.*x1s-
        2.*r1s*x1s+x1c+7.*x2+8.*r1s*x2+r1q*x2-7.*r2s*x2+6.*r1*r2s*x2
        +r1s*r2s*x2-2.*r2q*x2-9.*x1*x2-3.*r1s*x1*x2+2.*r2s*x1*x2
        +2.*x1s*x2-3.*x2s-r1s*x2s+2.*r2s*x2s+x1*x2s)
	/x3s;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-2.*r1s+r1q+r2s+6.*r1*r2s+r1s*r2s-2.*r2q);
        rFO2 = 2*(1.+r1+2.*r1s+2.*r1c+r1q+pow5(r1)+r2s+r1*r2s
        +5.*r1s*r2s+r1c*r2s-2.*r1*r2q-2.*x1-2.*r1*x1-2.*r1s*x1
        -2.*r1c*x1-2.*r2s*x1+5.*r1*r2s*x1-r1s*r2s*x1-2.*r2q*x1+x1s
        +r1*x1s+r2s*x1s-3.*x2-4.*r1s*x2-r1q*x2-2.*r2s*x2
        -2.*r1s*r2s*x2+4.*x1*x2+2.*r1s*x1*x2+r2s*x1*x2-x1s*x2
        +2.*x2s+2.*r1s*x2s-x1*x2s)
        /prop23
        +(2.*r2s-6.*r1*r2s-6.*r1s*r2s-6.*r1c*r2s+2.*r2q-6.*r1*r2q
        -r2s*x1+r1s*r2s*x1-r2q*x1+x2-r1q*x2-3.*r2s*x2+6.*r1*r2s*x2
        +9.*r1s*r2s*x2-2.*r2q*x2-x1*x2+r1s*x1*x2-x2s-3.*r1s*x2s
        +2.*r2s*x2s+x1*x2s)
        /prop2s
        +(-4.-8.*r1s-4.*r1q+4.*r2s-4.*r1s*r2s+8.*r2q+9.*x1+10.*r1s*x1
        +r1q*x1-3.*r2s*x1-6.*r1*r2s*x1+r1s*r2s*x1-2.*r2q*x1-6.*x1s
        -2.*r1s*x1s+x1c+7.*x2+8.*r1s*x2+r1q*x2-7.*r2s*x2-6.*r1*r2s*x2
        +r1s*r2s*x2-2.*r2q*x2-9.*x1*x2-3.*r1s*x1*x2+2.*r2s*x1*x2
        +2.*x1s*x2-3.*x2s-r1s*x2s+2.*r2s*x2s+x1*x2s)
	/x3s;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-2.*r1s+r1q+r2s+r1s*r2s-2.*r2q);
        rFO4 = 2*(1.+2.*r1s+r1q+r2s+5.*r1s*r2s-2.*x1-2.*r1s*x1
        -2.*r2s*x1-r1s*r2s*x1-2.*r2q*x1+x1s+r2s*x1s-3.*x2-4.*r1s*x2
        -r1q*x2-2.*r2s*x2-2.*r1s*r2s*x2+4.*x1*x2+2.*r1s*x1*x2+r2s*x1*x2
        -x1s*x2+2.*x2s+2.*r1s*x2s-x1*x2s)
        /prop23
        +(2.*r2s-6.*r1s*r2s+2.*r2q-r2s*x1+r1s*r2s*x1-r2q*x1+x2-r1q*x2
        -3.*r2s*x2+9.*r1s*r2s*x2-2.*r2q*x2-x1*x2+r1s*x1*x2-x2s-3.*r1s*x2s
        +2.*r2s*x2s+x1*x2s)
        /prop2s
        +(-4.-8.*r1s-4.*r1q+4.*r2s-4.*r1s*r2s+8.*r2q+9.*x1+10.*r1s*x1
        +r1q*x1-3.*r2s*x1+r1s*r2s*x1-2.*r2q*x1-6.*x1s-2.*r1s*x1s+x1c
        +7.*x2+8.*r1s*x2+r1q*x2-7.*r2s*x2+r1s*r2s*x2-2.*r2q*x2-9.*x1*x2
        -3.*r1s*x1*x2+2.*r2s*x1*x2+2.*x1s*x2-3.*x2s-r1s*x2s+2.*r2s*x2s
        +x1*x2s)
        /x3s;
        isSet4 = true;
      }
      break; 
 
    // S -> q qbar    (S = h0/H0/A0/H+-/...).
    case 4:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-r1s-r2s-2.*r1*r2);
        rFO1 = -(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q+x1
        -r1s*x1+2.*r1*r2*x1+3.*r2s*x1+x2+r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -2.*(r1s+r1q-2.*r1c*r2+r2s-6.*r1s*r2s-2.*r1*r2c+r2q-r1s*x1
        +r1*r2*x1+2.*r2s*x1+2.*r1s*x2+r1*r2*x2-r2s*x2-x1*x2)
        /prop12
        -(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q+x1-r1s*x1
        +r2s*x1+x2+3.*r1s*x2+2.*r1*r2*x2-r2s*x2-x1*x2)
        /prop2s;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-r1s-r2s+2.*r1*r2);
        rFO2 = -(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1-2.*r1*r2*x1+3.*r2s*x1+x2+r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2-2.*r1*r2*x2-r2s*x2-x1*x2)
        /prop2s
        +2.*(-r1s-r1q-2.*r1c*r2-r2s+6.*r1s*r2s-2.*r1*r2c-r2q+r1s*x1
        +r1*r2*x1-2.*r2s*x1-2.*r1s*x2+r1*r2*x2+r2s*x2+x1*x2)
        /prop12;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s-r2s);
        rFO4 = -(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+3.*r2s*x1+x2
        +r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -2.*(r1s+r1q+r2s-6.*r1s*r2s+r2q-r1s*x1
        +2.*r2s*x1+2.*r1s*x2-r2s*x2-x1*x2)
        /prop12
        -(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+r2s*x1
        +x2+3.*r1s*x2-r2s*x2-x1*x2)
        /prop2s;
        isSet4 = true;
      }
      break; 
 
    // q -> q S.
    case 5:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.+r1s-r2s+2.*r1);
        rFO1 = (4.-4.*r1s+4.*r2s-3.*x1-2.*r1*x1+r1s*x1-r2s*x1-5.*x2
        -2.*r1*x2+r1s*x2-r2s*x2+x1*x2+x2s)
        /x3s
        -2.*(3.-r1-5.*r1s-r1c+3.*r2s+r1*r2s-2.*x1-r1*x1
        +r1s*x1-4.*x2+2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +(2.-2.*r1-6.*r1s-2.*r1c+2.*r2s-2.*r1*r2s-x1+r1s*x1
        -r2s*x1-3.*x2+2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.+r1s-r2s-2.*r1);
        rFO2 = (4.-4.*r1s+4.*r2s-3.*x1+2.*r1*x1+r1s*x1-r2s*x1-5.*x2
        +2.*r1*x2+r1s*x2-r2s*x2+x1*x2+x2s)
        /x3s
        -2.*(3.+r1-5.*r1s+r1c+3.*r2s-r1*r2s-2.*x1+r1*x1
        +r1s*x1-4.*x2+2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +(2.+2.*r1-6.*r1s+2.*r1c+2.*r2s+2.*r1*r2s-x1+r1s*x1
        -r2s*x1-3.*x2-2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.+r1s-r2s);
        rFO4 = (4.-4.*r1s+4.*r2s-3.*x1+r1s*x1-r2s*x1-5.*x2+r1s*x2
        -r2s*x2+x1*x2+x2s)
        /x3s
        -2.*(3.-5.*r1s+3.*r2s-2.*x1+r1s*x1-4.*x2+2.*r1s*x2
        -r2s*x2+x1*x2+x2s)
        /prop23
        +(2.-6.*r1s+2.*r2s-x1+r1s*x1-r2s*x1-3.*x2+3.*r1s*x2
        -r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet4 = true;
      }
      break; 
 
    // V -> ~q ~qbar  (~q = squark).
    case 6:
      rLO1 = ps*(1.-2.*r1s+r1q-2.*r2s-2.*r1s*r2s+r2q);
      rFO1 = 2.*3.+(1.+r1s+r2s-x1)*(4.*r1s-x1s)
      /prop1s
      +2.*(-1.-3.*r1s-r2s+x1+x1s*0.5+x2-x1*x2*0.5)
      /prop1
      +(1.+r1s+r2s-x2)*(4.*r2s-x2s)
      /prop2s
      +2.*(-1.-r1s-3.*r2s+x1+x2-x1*x2*0.5+x2s*0.5)
      /prop2
      -(-4.*r1s-4.*r1q-4.*r2s-8.*r1s*r2s-4.*r2q+2.*x1+6.*r1s*x1
      +6.*r2s*x1-2.*x1s+2.*x2+6.*r1s*x2+6.*r2s*x2-4.*x1*x2
      -2.*r1s*x1*x2-2.*r2s*x1*x2+x1s*x2-2.*x2s+x1*x2s)
      /prop12;
      isSet1 = true;
      break; 
 
    // ~q -> ~q V.
    case 7:
      rLO1 = ps*(1.-2.*r1s+r1q-2.*r2s-2.*r1s*r2s+r2q);
      rFO1 = 16.*r2s-8.*(4.*r2s+2.*r2s*x1+x2+r1s*x2+r2s*x2-x1*x2
      -2.*x2s)
      /(3.*prop2)
      +8.*(1.+r1s+r2s-x2)*(4.*r2s-x2s)
      /(3.*prop2s)
      +8.*(x1+x2)*(-1.-2.*r1s-r1q-2.*r2s+2.*r1s*r2s-r2q+2.*x1
      +2.*r1s*x1+2.*r2s*x1-x1s+2.*x2+2.*r1s*x2+2.*r2s*x2-2.*x1*x2-x2s)
      /(3.*x3s)
      +8.*(-1.-r1s+r2s-x1)*(2.*r2s*x1+x2+r1s*x2+r2s*x2-x1*x2-x2s)
      /(3.*prop2*x3)
      -8.*(1.+2.*r1s+r1q+2.*r2s-2.*r1s*r2s+r2q-2.*x1-2.*r1s*x1
      -4.*r2s*x1+x1s-3.*x2-3.*r1s*x2-3.*r2s*x2+3.*x1*x2+2.*x2s)
      /(3.*x3);
      rFO1 = 3.*rFO1/8.;
      isSet1 = true;
      break; 
 
    // S -> ~q ~qbar.
    case 8:
      rLO1 = ps;
      rFO1 = (-1.-2.*r1s-r1q-2.*r2s+2.*r1s*r2s-r2q+2.*x1+2.*r1s*x1
      +2.*r2s*x1-x1s-r2s*x1s+2.*x2+2.*r1s*x2+2.*r2s*x2-3.*x1*x2
      -r1s*x1*x2-r2s*x1*x2+x1s*x2-x2s-r1s*x2s+x1*x2s)
      /(prop1s*prop2s);
      rFO1 = 2.*rFO1;
      isSet1 = true;
      break; 
 
    // ~q -> ~q S.
    case 9:
      rLO1 = ps;
      rFO1 = (-1.-r1s-r2s+x2)
      /prop2s
      +(1.+r1s-r2s+x1)
      /prop23
      -(x1+x2)
      /x3s;
      isSet1 = true;
      break; 
 
    // chi -> q ~qbar   (chi = neutralino/chargino).
    case 10:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.+r1s-r2s+2.*r1);
        rFO1 = (2.*r1+x1)*(-1.-r1s-r2s+x1)
        /prop1s
        +2.*(-1.-r1s-2.*r1c-r2s-2.*r1*r2s+3.*x1*0.5+r1*x1
        -r1s*x1*0.5-r2s*x1*0.5+x2+r1*x2+r1s*x2-x1*x2*0.5)
        /prop12
        +(2.-2.*r1-6.*r1s-2.*r1c+2.*r2s-2.*r1*r2s-x1+r1s*x1
        -r2s*x1-3.*x2+2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-2.*r1+r1s-r2s);
        rFO2 = (2.*r1-x1)*(1.+r1s+r2s-x1)
        /prop1s
        +2.*(-1.-r1s+2.*r1c-r2s+2.*r1*r2s+3.*x1*0.5-r1*x1
        -r1s*x1*0.5-r2s*x1*0.5+x2-r1*x2+r1s*x2-x1*x2*0.5)
        /prop12
        +(2.+2.*r1-6.*r1s+2.*r1c+2.*r2s+2.*r1*r2s-x1+r1s*x1
        -r2s*x1-3.*x2-2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)/
        prop2s;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.+r1s-r2s);
        rFO4 = x1*(-1.-r1s-r2s+x1)
        /prop1s
        +2.*(-1.-r1s-r2s+3.*x1*0.5-r1s*x1*0.5-r2s*x1*0.5
        +x2+r1s*x2-x1*x2*0.5)
        /prop12
        +(2.-6.*r1s+2.*r2s-x1+r1s*x1-r2s*x1-3.*x2+3.*r1s*x2
        -r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet4 = true;
      }
      break; 
 
    // ~q -> q chi.
    case 11:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-pow2(r1+r2));
        rFO1 = (1.+r1s+2.*r1*r2+r2s-x1-x2)*(x1+x2)
        /x3s
        -(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2+2.*r1*r2*x2-r2s*x2-x1*x2)
        /prop2s
        +(-1.-2.*r1s-r1q-2.*r1*r2-2.*r1c*r2+2.*r1*r2c+r2q+x1+r1s*x1
        -2.*r1*r2*x1-3.*r2s*x1+2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /prop23;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-pow2(r1-r2));
        rFO2 = (1.+r1s-2.*r1*r2+r2s-x1-x2)*(x1+x2)
        /x3s
        -(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2-2.*r1*r2*x2-r2s*x2-x1*x2)
        /prop2s
        +(-1.-2.*r1s-r1q+2.*r1*r2+2.*r1c*r2-2.*r1*r2c+r2q+x1+r1s*x1
        +2.*r1*r2*x1-3.*r2s*x1+2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /prop23;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s-r2s);
        rFO4 = (1.+r1s+r2s-x1-x2)*(x1+x2)
        /x3s
        -(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+r2s*x1+x2
        +3.*r1s*x2-r2s*x2-x1*x2)
        /prop2s
        +(-1.-2.*r1s-r1q+r2q+x1+r1s*x1-3.*r2s*x1
        +2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /prop23;
        isSet4 = true;
      }
      break; 
 
    // q -> ~q chi.
    case 12:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-r1s+r2s+2.*r2);
        rFO1 = (2.*r2+x2)*(-1.-r1s-r2s+x2)
        /prop2s
        +(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1-2.*r2*x1+r2s*x1+x1s
        -3.*x2-r1s*x2-2.*r2*x2+r2s*x2+x1*x2)
        /x3s
        +2.*(-1.-r1s+r2+r1s*r2-r2s-r2c+x1+r2*x1+r2s*x1+2.*x2
        +r1s*x2-x1*x2*0.5-x2s*0.5)
        /prop23;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-r1s+r2s-2.*r2);
        rFO2 = (2.*r2-x2)*(1.+r1s+r2s-x2)
        /prop2s
        +(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1+2.*r2*x1+r2s*x1+x1s
        -3.*x2-r1s*x2+2.*r2*x2+r2s*x2+x1*x2)
        /x3s
        +2.*(-1.-r1s-r2-r1s*r2-r2s+r2c+x1-r2*x1+r2s*x1+2.*x2
        +r1s*x2-x1*x2*0.5-x2s*0.5)
        /prop23;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s+r2s);
        rFO4 = x2*(-1.-r1s-r2s+x2)
        /prop2s
        +(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1+r2s*x1+x1s
        -3.*x2-r1s*x2+r2s*x2+x1*x2)
        /x3s
        +2.*(-1.-r1s-r2s+x1+r2s*x1+2.*x2
        +r1s*x2-x1*x2*0.5-x2s*0.5)
        /prop23;
        isSet4 = true;
      }
      break; 
 
    // ~g -> q ~qbar.
    case 13:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.+r1s-r2s+2.*r1);
        rFO1 = 4.*(2.*r1+x1)*(-1.-r1s-r2s+x1)
        /(3.*prop1s)
        -(-1.-r1s-2.*r1c-r2s-2.*r1*r2s+3.*x1*0.5+r1*x1-r1s*x1*0.5
        -r2s*x1*0.5+x2+r1*x2+r1s*x2-x1*x2*0.5)
        /(3.*prop12)
        +3.*(-1.+r1-r1s-r1c-r2s+r1*r2s+2.*x1+r2s*x1-x1s*0.5+x2+r1*x2
        +r1s*x2-x1*x2*0.5)
        /prop13
        +3.*(4.-4.*r1s+4.*r2s-3.*x1-2.*r1*x1+r1s*x1-r2s*x1-5.*x2
        -2.*r1*x2+r1s*x2-r2s*x2+x1*x2+x2s)
        /x3s
        -3.*(3.-r1-5.*r1s-r1c+3.*r2s+r1*r2s-2.*x1-r1*x1+r1s*x1
        -4.*x2+2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +4.*(2.-2.*r1-6.*r1s-2.*r1c+2.*r2s-2.*r1*r2s-x1+r1s*x1-r2s*x1
        -3.*x2+2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /(3.*prop2s);
        rFO1 = 3.*rFO1/4.;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.+r1s-r2s-2.*r1);
        rFO2 = 4.*(2.*r1-x1)*(1.+r1s+r2s-x1)
        /(3.*prop1s)
        +3.*(-1.-r1-r1s+r1c-r2s-r1*r2s+2.*x1+r2s*x1-x1s*0.5
        +x2-r1*x2+r1s*x2-x1*x2*0.5)
        /prop13
        +(2.+2.*r1s-4.*r1c+2.*r2s-4.*r1*r2s-3.*x1+2.*r1*x1
        +r1s*x1+r2s*x1-2.*x2+2.*r1*x2-2.*r1s*x2+x1*x2)
        /(6.*prop12)
        +3.*(4.-4.*r1s+4.*r2s-3.*x1+2.*r1*x1+r1s*x1-r2s*x1-5.*x2
        +2.*r1*x2+r1s*x2-r2s*x2+x1*x2+x2s)
        /x3s
        -3.*(3.+r1-5.*r1s+r1c+3.*r2s-r1*r2s-2.*x1+r1*x1+r1s*x1-4.*x2
        +2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +4.*(2.+2.*r1-6.*r1s+2.*r1c+2.*r2s+2.*r1*r2s-x1+r1s*x1-r2s*x1
        -3.*x2-2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /(3.*prop2s);
        rFO2 = 3.*rFO2/4.;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.+r1s-r2s);
        rFO4 = 8.*x1*(-1.-r1s-r2s+x1)
        /(3.*prop1s)
        +6.*(-1-r1s-r2s+2.*x1+r2s*x1-x1s*0.5+x2+r1s*x2-x1*x2*0.5)
        /prop13
        +(2.+2.*r1s+2.*r2s-3.*x1+r1s*x1+r2s*x1-2.*x2-2.*r1s*x2+x1*x2)
        /(3.*prop12)
        +6.*(4.-4.*r1s+4.*r2s-3.*x1+r1s*x1-r2s*x1-5.*x2+r1s*x2-r2s*x2
        +x1*x2+x2s)
        /x3s
        -6.*(3.-5.*r1s+3.*r2s-2.*x1+r1s*x1-4.*x2+2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +8.*(2.-6.*r1s+2.*r2s-x1+r1s*x1-r2s*x1-3.*x2+3.*r1s*x2-r2s*x2
        +x1*x2+x2s)
        /(3.*prop2s);
        rFO4 = 3.*rFO4/8.;
        isSet4 = true;
      }
      break; 
 
    // ~q -> q ~g.
    case 14:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-r1s-r2s-2.*r1*r2);
        rFO1 = 64.*(1.+r1s+2.*r1*r2+r2s-x1-x2)*(x1+x2)
        /(9.*x3s)
        -16.*(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q
        +x1-r1s*x1+2.*r1*r2*x1+3.*r2s*x1+x2+r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -16.*(r1s+r1q-2.*r1c*r2+r2s-6.*r1s*r2s-2.*r1*r2c+r2q-r1s*x1
        +r1*r2*x1+2.*r2s*x1+2.*r1s*x2+r1*r2*x2-r2s*x2-x1*x2)
        /prop12
        -64.*(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2+2.*r1*r2*x2-r2s*x2-x1*x2)
        /(9.*prop2s)
        +8.*(-1.+r1q-2.*r1*r2+2.*r1c*r2-2.*r2s-2.*r1*r2c-r2q-2.*r1s*x1
        +2.*r2s*x1+x1s+x2-3.*r1s*x2-2.*r1*r2*x2+r2s*x2+x1*x2)
        /prop13
        -8.*(-1.-2.*r1s-r1q-2.*r1*r2-2.*r1c*r2+2.*r1*r2c+r2q+x1+r1s*x1
        -2.*r1*r2*x1-3.*r2s*x1+2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /(9.*prop23);
        rFO1 = 9.*rFO1/64.;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-r1s-r2s+2.*r1*r2);
        rFO2 = 64.*(1.+r1s-2.*r1*r2+r2s-x1-x2)*(x1+x2)
        /(9.*x3s)
        -16.*(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1-2.*r1*r2*x1+3.*r2s*x1+x2+r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -64.*(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2-2.*r1*r2*x2-r2s*x2-x1*x2)
        /(9.*prop2s)
        +16.*(-r1s-r1q-2.*r1c*r2-r2s+6.*r1s*r2s-2.*r1*r2c-r2q+r1s*x1
        +r1*r2*x1-2.*r2s*x1-2.*r1s*x2+r1*r2*x2+r2s*x2+x1*x2)
        /prop12
        +8.*(-1.+r1q+2.*r1*r2-2.*r1c*r2-2.*r2s+2.*r1*r2c-r2q-2.*r1s*x1
        +2.*r2s*x1+x1s+x2-3.*r1s*x2+2.*r1*r2*x2+r2s*x2+x1*x2)
        /prop13
        -8.*(-1.-2.*r1s-r1q+2.*r1*r2+2.*r1c*r2-2.*r1*r2c+r2q+x1+r1s*x1+
        2.*r1*r2*x1-3.*r2s*x1+2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /(9.*prop23);
        rFO2 = 9.*rFO2/64.;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s-r2s);
        rFO4 = 128.*(1.+r1s+r2s-x1-x2)*(x1+x2)
        /(9.*x3s)
        -32*(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+3.*r2s*x1+x2
        +r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -32.*(r1s+r1q+r2s-6.*r1s*r2s+r2q-r1s*x1+2.*r2s*x1+2.*r1s*x2
        -r2s*x2-x1*x2)
        /prop12
        -128.*(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+r2s*x1+x2+3.*r1s*x2
        -r2s*x2-x1*x2)
        /(9.*prop2s)
        +16.*(-1.+r1q-2.*r2s-r2q-2.*r1s*x1+2.*r2s*x1+x1s
        +x2-3.*r1s*x2+r2s*x2+x1*x2)
        /prop13
        -16.*(-1.-2.*r1s-r1q+r2q+x1+r1s*x1-3.*r2s*x1
        +2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /(9.*prop23);
        rFO4 = 9.*rFO4/128.;
        isSet4 = true;
      }
      break; 
 
    // q -> ~q ~g.
    case 15:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-r1s+r2s+2.*r2);
        rFO1 = 32*(2.*r2+x2)*(-1.-r1s-r2s+x2)
        /(9.*prop2s)
        +8.*(-1.-r1s-2.*r1s*r2-r2s-2.*r2c+x1+r2*x1+r2s*x1
        +3.*x2*0.5-r1s*x2*0.5+r2*x2-r2s*x2*0.5-x1*x2*0.5)
        /prop12
        +8.*(2.+2.*r1s-2.*r2-2.*r1s*r2-6.*r2s-2.*r2c-3.*x1-r1s*x1
        +2.*r2*x1+3.*r2s*x1+x1s-x2-r1s*x2+r2s*x2+x1*x2)
        /prop1s
        +32.*(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1-2.*r2*x1+r2s*x1+x1s
        -3.*x2-r1s*x2-2.*r2*x2+r2s*x2+x1*x2)
        /(9.*x3s)
        -8.*(3.+3.*r1s-r2+r1s*r2-5.*r2s-r2c-4.*x1-r1s*x1
        +2.*r2s*x1+x1s-2.*x2-r2*x2+r2s*x2+x1*x2)
        /prop13
        -8.*(-1.-r1s+r2+r1s*r2-r2s-r2c+x1+r2*x1+r2s*x1+2.*x2+r1s*x2
        -x1*x2*0.5-x2s*0.5)
        /(9.*prop23);
        rFO1 = 9.*rFO1/32.;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-r1s+r2s-2.*r2);
        rFO2 = 32*(2.*r2-x2)*(1.+r1s+r2s-x2)
        /(9.*prop2s)
        +8.*(-1.-r1s+2.*r1s*r2-r2s+2.*r2c+x1-r2*x1+r2s*x1
        +3.*x2*0.5-r1s*x2*0.5-r2*x2-r2s*x2*0.5-x1*x2*0.5)
        /prop12
        +8.*(2.+2.*r1s+2.*r2+2.*r1s*r2-6.*r2s+2.*r2c-3.*x1-r1s*x1
        -2.*r2*x1+3.*r2s*x1+x1s-x2-r1s*x2+r2s*x2+x1*x2)
        /prop1s
        -8.*(3.+3.*r1s+r2-r1s*r2-5.*r2s+r2c-4.*x1-r1s*x1+2.*r2s*x1+x1s
        -2.*x2+r2*x2+r2s*x2+x1*x2)
        /prop13
        +32*(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1+2.*r2*x1+r2s*x1
        +x1s-3.*x2-r1s*x2+2.*r2*x2+r2s*x2+x1*x2)
        /(9.*x3s)
        -8.*(-1.-r1s-r2-r1s*r2-r2s+r2c+x1-r2*x1+r2s*x1+2.*x2+r1s*x2
        -x1*x2*0.5-x2s*0.5)
        /(9.*prop23);
        rFO2 = 9.*rFO2/32.;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s+r2s);
        rFO4 = 64.*x2*(-1.-r1s-r2s+x2)
        /(9.*prop2s)
        +16.*(-1.-r1s-r2s+x1+r2s*x1+3.*x2*0.5-r1s*x2*0.5
        -r2s*x2*0.5-x1*x2*0.5)
        /prop12
        -16.*(3.+3.*r1s-5.*r2s-4.*x1-r1s*x1+2.*r2s*x1+x1s-2.*x2+r2s*x2
        +x1*x2)
        /prop13
        +64.*(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1+r2s*x1+x1s-3.*x2
        -r1s*x2+r2s*x2+x1*x2)
        /(9.*x3s)
        +16.*(2.+2.*r1s-6.*r2s-3.*x1-r1s*x1+3.*r2s*x1+x1s
        -x2-r1s*x2+r2s*x2+x1*x2)
        /prop1s
        -16.*(-1.-r1s-r2s+x1+r2s*x1+2.*x2+r1s*x2-x1*x2*0.5-x2s*0.5)
        /(9.*prop23);
        rFO4 = 9.*rFO4/64.;
        isSet4 = true;
      }
      break; 
 
    // g -> ~g ~g. Use (9/4)*eikonal. May be changed in the future.
    case 16:
      rLO = ps;
      if (combi == 2) offset = x3s;
      else if (combi == 3) offset = mix * x3s;
      else if (combi == 4) offset = 0.5 * x3s;
      rFO = ps * 4.5 * ( (x1+x2-1.+offset-r1s-r2s)/prop12 
      - r1s/prop2s - r2s/prop1s );
      break; 

    // Eikonal expression for kind == 1; also acts as default.
    default:
      rLO = ps;
      if (combi == 2) offset = x3s;
      else if (combi == 3) offset = mix * x3s;
      else if (combi == 4) offset = 0.5 * x3s;
      rFO = ps * 2. * ( (x1+x2-1.+offset-r1s-r2s)/prop12 
      - r1s/prop2s - r2s/prop1s );
      break;

  // End of ME cases. 
  }

  // Find relevant leading and first order expressions.
  if (combi == 1 && isSet1) {rLO = rLO1; rFO = rFO1;}     
  else if (combi == 2 && isSet2) {rLO = rLO2; rFO = rFO2;}     
  else if (combi == 3 && isSet1 && isSet2) {
    rLO = mix * rLO1 + (1.-mix) * rLO2; 
    rFO = mix * rFO1 + (1.-mix) * rFO2; }
  else if (isSet4) {rLO = rLO4; rFO = rFO4;}     
  else if (combi == 4 && isSet1 && isSet2) {
    rLO = 0.5 * (rLO1 + rLO2);
    rFO = 0.5 * (rFO1 + rFO2); }
  else if (isSet1) {rLO = rLO1; rFO = rFO1;} 

  // Return ratio of first to leading order cross section.     
  return rFO / rLO;
}  

//*********

// Find coefficient of azimuthal asymmetry from gluon polarization.

void TimeShower::findAsymPol( Event& event, TimeDipoleEnd* dip) {

  // Default is no asymmetry. Only gluons are studied.
  dip->asymPol = 0.;
  dip->iAunt = 0;
  int iRad = dip->iRadiator;
  if (!doPhiPolAsym || event[iRad].id() != 21) return;

  // Trace grandmother via possibly intermediate recoil copies.
  int iMother = event.iTopCopy(iRad);
  int iGrandM = event[iMother].mother1();

  // Check grandmother flavour and set aunt.
  if (!event[iGrandM].isQorG()) return;
  dip->iAunt = (event[iGrandM].daughter1() == iMother) 
    ? event[iGrandM].daughter2() : event[iGrandM].daughter1();

  // Coefficient from gluon production (approximate z by energy).
  double zProd = event[iRad].e() / (event[iRad].e() 
    + event[dip->iAunt].e());
  if (event[iGrandM].id() != 21) dip->asymPol = 2. * (1. - zProd) 
    / (1. + pow2(1. - zProd) );
  else dip->asymPol = pow2( (1. - zProd) / (1. - zProd * (1. - zProd) ) );

  // Coefficients from gluon decay.
  if (dip->flavour == 21) dip->asymPol *= pow2( (1. - dip->z) 
    / (1. - dip->z * (1. - dip->z) ) );
  else  dip->asymPol *= -2. * dip->z *( 1. - dip->z ) 
    / (1. - 2. * dip->z * (1. - dip->z) );

}

//*********

// Print the list of dipoles.

void TimeShower::list(ostream& os) {

  // Header.
  os << "\n --------  Dipole Listing  ---------------------------------"
     << "------------------- \n \n    i  rad  rec       pTmax  col  chg "  
     << "type  rec     mix  ord  spl  glu \n" << fixed << setprecision(3);
  
  // Loop over dipole list and print it.
  for (int i = 0; i < int(dipole.size()); ++i) {
    os << setw(5) << i << setw(5) << dipole[i].iRadiator << setw(5) 
       << dipole[i].iRecoiler << setw(12) << dipole[i].pTmax 
       << setw(5) << dipole[i].colType << setw(5) << dipole[i].chgType
       << setw(5) << dipole[i].MEtype << setw(5) << dipole[i].iMEpartner
       << setw(8) << dipole[i].MEmix << setw(5) << dipole[i].MEorder
       << setw(5) << dipole[i].MEsplit << setw(5) << dipole[i].MEgluinoDau
       << "\n";
  }
}

//**************************************************************************

} // end namespace Pythia8
