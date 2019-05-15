// Function definitions (not found in the header) for the 
// PhaseSpace and PhaseSpace2to2tauyz classes.
// Copyright C 2007 Torbjorn Sjostrand

#include "PhaseSpace.h"

namespace Pythia8 {

//**************************************************************************

// The PhaseSpace class.
// Base class for phase space generators.

//*********
 
// Definitions of static variables and functions.
// (Values will be overwritten in initStatic call, so are purely dummy.)

double PhaseSpace::mHatGlobalMin        = 4.;
double PhaseSpace::mHatGlobalMax        = -1.;
double PhaseSpace::pTHatMin             = 0.;
double PhaseSpace::pTHatMax             = -1.;
double PhaseSpace::pTHatMinDiverge      = 1.;
double PhaseSpace::pT2HatMin            = 0.;
double PhaseSpace::pT2HatMax            = 1.;
bool   PhaseSpace::useBreitWigners      = true;
double PhaseSpace::minWidthBreitWigners = 0.01;
bool   PhaseSpace::showSearch           = false;
bool   PhaseSpace::showViolation        = false;
int    PhaseSpace::gmZmodeGlob          = 0;

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Maximum cross section increase, just in case true maximum not found.
const double PhaseSpace::SAFETYMARGIN   = 1.05;

// Small number to avoid division by zero.
const double PhaseSpace::TINY           = 1e-20;

// Fraction of total weight that is shared evenly between all shapes.
const double PhaseSpace::EVENFRAC       = 0.4;

// Two cross sections with a small relative error are assumed same.
const double PhaseSpace::SAMESIGMA      = 1e-6;

// Do not include resonances peaked too far outside allowed mass region.
const double PhaseSpace::WIDTHMARGIN    = 20.;

// Special optimization treatment when two resonances at almost same mass.
const double PhaseSpace::SAMEMASS       = 0.01;

// Minimum phase space left when kinematics constraints are combined.
const double PhaseSpace::MASSMARGIN     = 0.01;

// When using Breit-Wigners in 2 -> 2 raise maximum weight estimate.
const double PhaseSpace::EXTRABWWTMAX   = 1.25;

// Size of Breit-Wigner threshold region, for mass selection biasing.
const double PhaseSpace::THRESHOLDSIZE  = 3.;

// Step size in optimal-mass search, for mass selection biasing.
const double PhaseSpace::THRESHOLDSTEP  = 0.2;

// Cutoff for f_e^e at x < 1 - 10^{-10} to be used in phase space selection.
// Note: the ...MIN quantities come from 1 - x_max or 1 - tau_max.
const double PhaseSpace::LEPTONXMIN     = 1e-10;
const double PhaseSpace::LEPTONXMAX     = 1. - 1e-10;
const double PhaseSpace::LEPTONXLOGMIN  = log(1e-10);
const double PhaseSpace::LEPTONXLOGMAX  = log(1. - 1e-10);
const double PhaseSpace::LEPTONTAUMIN   = 2e-10;

// Number of trial maxima around which maximum search is performed.
const int    PhaseSpace::NMAXTRY        = 2;

// Information on incoming beams.
BeamParticle* PhaseSpace::beamAPtr      = 0;
BeamParticle* PhaseSpace::beamBPtr      = 0;
int    PhaseSpace::idA                  = 0;
int    PhaseSpace::idB                  = 0;
double PhaseSpace::mA                   = 0.; 
double PhaseSpace::mB                   = 0.;
bool   PhaseSpace::hasLeptonBeams       = false;
bool   PhaseSpace::hasPointLeptons      = false;
  
// Pointer to the total/elastic/diffractive cross section object.
SigmaTotal* PhaseSpace::sigmaTotPtr     = 0;

// Pointer to userHooks object.
UserHooks* PhaseSpace::userHooksPtr     = 0;
bool   PhaseSpace::canModifySigma       = false;

//*********

// Initialize static data members.

void PhaseSpace::initStatic() {

  // Standard phase space cuts.
  mHatGlobalMin        = Settings::parm("PhaseSpace:mHatMin");
  mHatGlobalMax        = Settings::parm("PhaseSpace:mHatMax");
  pTHatMin             = Settings::parm("PhaseSpace:pTHatMin");
  pTHatMax             = Settings::parm("PhaseSpace:pTHatMax");
  pTHatMinDiverge      = Settings::parm("PhaseSpace:pTHatMinDiverge");

  // When to use Breit-Wigners.
  useBreitWigners      = Settings::flag("PhaseSpace:useBreitWigners");
  minWidthBreitWigners = Settings::parm("PhaseSpace:minWidthBreitWigners");

  // Quadratic limits often needed, so calculate them.
  pT2HatMin            = pow2(pTHatMin);
  pT2HatMax            = pow2(pTHatMax);

  // Print flag for maximization information.
  showSearch           = Settings::flag("PhaseSpace:showSearch");
  showViolation        = Settings::flag("PhaseSpace:showViolation");

  // Know whether a Z0 is pure Z0 or admixed with gamma*.
  gmZmodeGlob          = Settings::mode("SigmaProcess:gmZmode");  

}

//*********

// Store pointers to beams and SigmaTotal.
 
void PhaseSpace::setStaticPtrs( BeamParticle* beamAPtrIn, 
  BeamParticle* beamBPtrIn, SigmaTotal* sigmaTotPtrIn,
  UserHooks* userHooksPtrIn) {

  // Store input.
  beamAPtr        = beamAPtrIn;
  beamBPtr        = beamBPtrIn;
  sigmaTotPtr     = sigmaTotPtrIn;
  userHooksPtr    = userHooksPtrIn; 

  // Some commonly used beam information.
  idA             = beamAPtr->id(); 
  idB             = beamBPtr->id(); 
  mA              = beamAPtr->m(); 
  mB              = beamBPtr->m(); 

  // Flag if lepton beams, and if non-resolved ones.
  hasLeptonBeams  = ( beamAPtr->isLepton() || beamBPtr->isLepton() );
  hasPointLeptons = ( hasLeptonBeams 
    && (beamAPtr->isUnresolved() || beamBPtr->isUnresolved() ) );

  // Flag if user should be allow to reweight cross section.
  canModifySigma  = (userHooksPtr > 0) 
                  ? userHooksPtr->canModifySigma() : false; 

}

//*********

// Save pointers and values.

void PhaseSpace::initInfo(SigmaProcess* sigmaProcessPtrIn, double eCMIn) {

  // Store input pointers for future use. CM energy.
  sigmaProcessPtr = sigmaProcessPtrIn;
  eCM      = eCMIn;
  s        = eCM * eCM;

  // Default event-specific kinematics properties.
  x1H      = 1.;
  x2H      = 1.;
  m3       = 0.;
  m4       = 0.;
  s3       = m3 * m3;
  s4       = m4 * m4;
  mHat     = eCM;
  sH       = s;
  tH       = 0.;
  uH       = 0.;
  pTH      = 0.;
  theta    = 0.;
  phi      = 0.;
  runBW3H  = 1.;
  runBW4H  = 1.;

  // Default cross section information.
  sigmaNw  = 0.;
  sigmaMx  = 0.;
  sigmaNeg = 0.;

}

//*********

// Allow for nonisotropic decays when ME's available.

void PhaseSpace::decayKinematics( Event& process) {

  // Identify sets of sister partons. 
  int iResEnd = 5;
  for (int iResBeg = 5; iResBeg < process.size(); ++iResBeg) {
    if (iResBeg < iResEnd) continue;
    iResEnd = iResBeg + 1;
    while ( iResEnd < process.size() 
      && process[iResEnd].mother1() == process[iResBeg].mother1()
      && process[iResEnd].mother2() == process[iResBeg].mother2() )
      ++iResEnd;

    // Check that at least one of them is a resonance.
    bool hasRes = false;
    for (int iRes = iResBeg; iRes < iResEnd; ++iRes)
      if ( !process[iRes].isFinal() ) hasRes = true;
    if ( !hasRes ) continue; 

    // Evaluate matrix element and decide whether to keep kinematics.
    while ( sigmaProcessPtr->weightDecay( process, iResBeg, iResEnd) 
      < Rndm::flat() ) {

      // Find resonances for which to redo decay angles.
      for (int iRes = iResBeg; iRes < process.size(); ++iRes) {
        if ( process[iRes].isFinal() ) continue;
        int iResMother = iRes;
        while (iResMother >= iResEnd) 
          iResMother = process[iResMother].mother1();
        if (iResMother < iResBeg) continue;

        // Identify daughters. Find mother and daughter masses.
        int    i1 = process[iRes].daughter1();
        int    i2 = process[iRes].daughter2();
        double m0 = process[iRes].m();
        double m1 = process[i1].m();
        double m2 = process[i2].m();

        // Energies and absolute momentum in the rest frame.
        double e1   = 0.5 * (m0*m0 + m1*m1 - m2*m2) / m0;
        double e2   = 0.5 * (m0*m0 + m2*m2 - m1*m1) / m0;
        double pAbs = 0.5 * sqrtpos( (m0 - m1 - m2) * (m0 + m1 + m2)
          * (m0 + m1 - m2) * (m0 - m1 + m2) ) / m0;  

        // Pick isotropic angles to give three-momentum. 
        double cosTheta = 2. * Rndm::flat() - 1.;
        double sinTheta = sqrt(1. - cosTheta*cosTheta);
        double phi      = 2. * M_PI * Rndm::flat();
        double pX       = pAbs * sinTheta * cos(phi);  
        double pY       = pAbs * sinTheta * sin(phi);  
        double pZ       = pAbs * cosTheta;  

        // Fill four-momenta in mother rest frame and boost them. 
        Vec4 p1(  pX,  pY,  pZ, e1);
        Vec4 p2( -pX, -pY, -pZ, e2);
        p1.bst( process[iRes].p() );
        p2.bst( process[iRes].p() );
        process[i1].p( p1 );
        process[i2].p( p2 );

      // End loop over resonance decay chains.
      }

    // Ready to allow new test of matrix element.
    }

  // End loop over sets of sister resonances/partons. 
  }

}

//*********

// Determine how phase space should be sampled.

bool PhaseSpace::setupSampling1or2(bool is2, ostream& os) {

  // Should Z0 be treated as such or as gamma*/Z0?
  gmZmode         = gmZmodeGlob;
  int gmZmodeProc = sigmaProcessPtr->gmZmode();
  if (gmZmodeProc >= 0) gmZmode = gmZmodeProc;

  // Optional printout.
  if (showSearch) os <<  "\n Optimization printout for "  
    << sigmaProcessPtr->name() << "\n \n" << scientific << setprecision(3);

  // Check that open range in tau (+ set tauMin, tauMax).
  if (!limitTau(is2)) return false; 

  // Reset coefficients and matrices of equation system to solve.
  int binTau[8], binY[8], binZ[8];
  double vecTau[8], matTau[8][8], vecY[8], matY[8][8], vecZ[8], matZ[8][8];
  for (int i = 0; i < 8; ++i) {
    tauCoef[i] = 0.;
    yCoef[i]   = 0.;
    zCoef[i]   = 0.;
    binTau[i]  = 0;
    binY[i]    = 0;
    binZ[i]    = 0;
    vecTau[i]  = 0.;
    vecY[i]    = 0.;
    vecZ[i]    = 0.;
    for (int j = 0; j < 8; ++j) { 
      matTau[i][j] = 0.;
      matY[i][j]   = 0.;
      matZ[i][j]   = 0.;
    }  
  }
  sigmaMx  = 0.;
  sigmaNeg = 0.;
  
  // Number of used coefficients/points for each dimension: tau, y, c.
  nTau = (hasPointLeptons) ? 1 : 2;
  nY   = (hasPointLeptons) ? 1 : 3;
  nZ   = (is2) ? 5 : 1; 

  // Identify if any resonances contribute in s-channel.
  idResA = sigmaProcessPtr->resonanceA();
  if (idResA != 0) { 
     mResA = ParticleDataTable::m0(idResA);
     GammaResA = ParticleDataTable::mWidth(idResA);
     if (mHatMin > mResA + WIDTHMARGIN * GammaResA || (mHatMax > 0. 
       && mHatMax < mResA - WIDTHMARGIN * GammaResA) ) idResA = 0; 
  }
  idResB = sigmaProcessPtr->resonanceB();
  if (idResB != 0) { 
     mResB = ParticleDataTable::m0(idResB);
     GammaResB = ParticleDataTable::mWidth(idResB);
     if (mHatMin > mResB + WIDTHMARGIN * GammaResB || (mHatMax > 0.  
       && mHatMax < mResB - WIDTHMARGIN * GammaResB) ) idResB = 0; 
  }
  if (idResA == 0 && idResB != 0) {
    idResA = idResB;
    mResA = mResB;
    GammaResA = GammaResB;
    idResB = 0;
  }

  // More sampling in tau if resonances in s-channel.
  if (idResA !=0 && !hasPointLeptons) {
    nTau += 2;
    tauResA = mResA * mResA / s;
    widResA = mResA * GammaResA / s;
  }
  if (idResB != 0 && !hasPointLeptons) {
    nTau += 2;
    tauResB = mResB * mResB / s;
    widResB = mResB * GammaResB / s;
  }
  
  // More sampling in tau and y if incoming lepton beams.
  if (hasLeptonBeams && !hasPointLeptons) {
    ++nTau;
    nY += 2;
  }

  // Special case when both resonances have same mass.
  sameResMass = false;
  if (idResB != 0 && abs(mResA - mResB) < SAMEMASS * (GammaResA + GammaResB))
    sameResMass = true;
    
  // Default z value and weight required for 2 -> 1. Number of dimensions.
  z = 0.;
  wtZ = 1.;
  int nVar = (is2) ? 3 : 2;

  // Initial values, to be modified later.
  tauCoef[0] = 1.;
  yCoef[0]   = 0.5;
  yCoef[1]   = 0.5;
  zCoef[0]   = 1.; 

  // Step through grid in tau. Set limits on y and z generation. 
  for (int iTau = 0; iTau < nTau; ++iTau) {
    double posTau = 0.5;
    if (sameResMass && iTau > 1 && iTau < 6) posTau = (iTau < 4) ? 0.4 : 0.6;
    selectTau( iTau, posTau, is2);
    if (!limitY()) continue;
    if (is2 && !limitZ()) continue;

    // Step through grids in y and z. 
    for (int iY = 0; iY < nY; ++iY) {
      selectY( iY, 0.5);
      for (int iZ = 0; iZ < nZ; ++iZ) {
        if (is2) selectZ( iZ, 0.5);

        // Calculate cross section. Weight by phase-space volume
        // and Breit-Wigners for masses in 2 -> 2.
        if (is2) sigmaProcessPtr->set2Kin( x1H, x2H, sH, tH, m3, m4, 
                                           runBW3H, runBW4H);
        else     sigmaProcessPtr->set1Kin( x1H, x2H, sH);
        double sigmaNow = sigmaProcessPtr->sigmaPDF();
        sigmaNow *= wtTau * wtY * wtZ * wtBW; 

        // Allow possibility for user to modify cross section.
        if (canModifySigma) sigmaNow 
           *= userHooksPtr->multiplySigmaBy( sigmaProcessPtr, this, false);

        // Check if current maximum exceeded.
        if (sigmaNow > sigmaMx) sigmaMx = sigmaNow; 

        // Optional printout. Protect against negative cross sections.
        if (showSearch) os << " tau =" << setw(11) << tau << "  y =" 
	  << setw(11) << y << "  z =" << setw(11) << z
	  << "  sigma =" << setw(11) << sigmaNow << "\n";
        if (sigmaNow < 0.) sigmaNow = 0.; 

        // Sum up tau cross-section pieces in points used.
        if (!hasPointLeptons) {
          binTau[iTau] += 1;
          vecTau[iTau] += sigmaNow;
          matTau[iTau][0] += 1. / intTau0;
          matTau[iTau][1] += (1. / intTau1) / tau;
          if (idResA != 0) {
            matTau[iTau][2] += (1. / intTau2) / (tau + tauResA);
            matTau[iTau][3] += (1. / intTau3) 
              * tau / ( pow2(tau - tauResA) + pow2(widResA) );
          }
          if (idResB != 0) {
            matTau[iTau][4] += (1. / intTau4) / (tau + tauResB);
            matTau[iTau][5] += (1. / intTau5) 
              * tau / ( pow2(tau - tauResB) + pow2(widResB) );
          }
          if (hasLeptonBeams) matTau[iTau][nTau - 1] += (1. / intTau6) 
              * tau / max( LEPTONTAUMIN, 1. - tau);
        }

        // Sum up y cross-section pieces in points used.
        if (!hasPointLeptons) {
          binY[iY] += 1;
          vecY[iY] += sigmaNow;
          matY[iY][0] += (yMax / intY01) * (y + yMax);
          matY[iY][1] += (yMax / intY01) * (yMax - y);
          matY[iY][2] += (yMax / intY2) / cosh(y);
          if (hasLeptonBeams) {
            matY[iY][3] += (yMax / intY34) 
              / max( LEPTONXMIN, 1. - exp( y - yMax) );
            matY[iY][4] += (yMax / intY34) 
              / max( LEPTONXMIN, 1. - exp(-y - yMax) );
          }
	}

        // Integrals over z expressions at tauMax, to be used below.
        if (is2) {
          double p2AbsMax = 0.25 * (pow2(tauMax * s - s3 - s4) 
            - 4. * s3 * s4) / (tauMax * s);         
          double zMaxMax = sqrtpos( 1. - pT2HatMinNow / p2AbsMax );
          double zPosMaxMax = max(ratio34, unity34 + zMaxMax);
          double zNegMaxMax = max(ratio34, unity34 - zMaxMax);
          double intZ0Max = 2. * zMaxMax;
          double intZ12Max = log( zPosMaxMax / zNegMaxMax);
          double intZ34Max = 1. / zNegMaxMax - 1. / zPosMaxMax;  
  
          // Sum up z cross-section pieces in points used.
          binZ[iZ] += 1;
          vecZ[iZ] += sigmaNow;
          matZ[iZ][0] += 1.; 
          matZ[iZ][1] += (intZ0Max / intZ12Max) / zNeg;
          matZ[iZ][2] += (intZ0Max / intZ12Max) / zPos;
          matZ[iZ][3] += (intZ0Max / intZ34Max) / pow2(zNeg);
          matZ[iZ][4] += (intZ0Max / intZ34Max) / pow2(zPos);
	}

      // End of loops over phase space points. 
      }
    }
  }   

  // Fail if no non-vanishing cross sections.
  if (sigmaMx <= 0.) {
    sigmaMx = 0.;
    return false;
  }   

  // Solve respective equation system for better phase space coefficients.
  if (!hasPointLeptons) solveSys( nTau, binTau, vecTau, matTau, tauCoef);
  if (!hasPointLeptons) solveSys( nY, binY, vecY, matY, yCoef);
  if (is2)              solveSys( nZ, binZ, vecZ, matZ, zCoef);
  if (showSearch) os << "\n";

  // Provide cumulative sum of coefficients.
  tauCoefSum[0] = tauCoef[0];
    yCoefSum[0] =   yCoef[0];
    zCoefSum[0] =   zCoef[0];
  for (int i = 1; i < 8; ++ i) {
    tauCoefSum[i] = tauCoefSum[i - 1] + tauCoef[i]; 
      yCoefSum[i] =   yCoefSum[i - 1] +   yCoef[i]; 
      zCoefSum[i] =   zCoefSum[i - 1] +   zCoef[i]; 
  }
  // The last element should be > 1 to be on safe side in selection below.
  tauCoefSum[nTau - 1] = 2.;
    yCoefSum[nY   - 1] = 2.;
    zCoefSum[nZ   - 1] = 2.;
  
  
  // Begin find two most promising maxima among same points as before.
  int iMaxTau[NMAXTRY + 2], iMaxY[NMAXTRY + 2], iMaxZ[NMAXTRY + 2];
  double sigMax[NMAXTRY + 2];
  int nMax = 0;

  // Scan same grid as before in tau, y, z. 
  for (int iTau = 0; iTau < nTau; ++iTau) {
    double posTau = 0.5;
    if (sameResMass && iTau > 1 && iTau < 6) posTau = (iTau < 4) ? 0.4 : 0.6;
    selectTau( iTau, posTau, is2);
    if (!limitY()) continue;
    if (is2 && !limitZ()) continue;
    for (int iY = 0; iY < nY; ++iY) {
      selectY( iY, 0.5);
      for (int iZ = 0; iZ < nZ; ++iZ) {
        if (is2) selectZ( iZ, 0.5);

        // Calculate cross section. Weight by phase-space volume
        // and Breit-Wigners for masses in 2 -> 2.
        if (is2) sigmaProcessPtr->set2Kin( x1H, x2H, sH, tH, m3, m4,
                                           runBW3H, runBW4H);
        else     sigmaProcessPtr->set1Kin( x1H, x2H, sH);
        double sigmaNow = sigmaProcessPtr->sigmaPDF();
        sigmaNow *= wtTau * wtY * wtZ * wtBW;

        // Allow possibility for user to modify cross section.
        if (canModifySigma) sigmaNow 
          *= userHooksPtr->multiplySigmaBy( sigmaProcessPtr, this, false);

        // Optional printout. Protect against negative cross section.
        if (showSearch) os << " tau =" << setw(11) << tau << "  y =" 
	  << setw(11) << y << "  z =" << setw(11) << z
	  << "  sigma =" << setw(11) << sigmaNow << "\n";
        if (sigmaNow < 0.) sigmaNow = 0.; 

        // Check that point is not simply mirror of already found one.
        bool mirrorPoint = false;
        for (int iMove = 0; iMove < nMax; ++iMove)
          if (abs(sigmaNow - sigMax[iMove]) < SAMESIGMA 
	    * (sigmaNow + sigMax[iMove])) mirrorPoint = true; 

        // Add to or insert in maximum list. Only first two count.
	if (!mirrorPoint) {
          int iInsert = 0;
          for (int iMove = nMax - 1; iMove >= -1; --iMove) {
            iInsert = iMove + 1;
            if (iInsert == 0 || sigmaNow < sigMax[iMove]) break;
	    iMaxTau[iMove + 1] = iMaxTau[iMove];
	    iMaxY[iMove + 1] = iMaxY[iMove];
	    iMaxZ[iMove + 1] = iMaxZ[iMove];
	    sigMax[iMove + 1] = sigMax[iMove];
	  }
	  iMaxTau[iInsert] = iTau;
	  iMaxY[iInsert] = iY;
	  iMaxZ[iInsert] = iZ;
	  sigMax[iInsert] = sigmaNow;
          if (nMax < NMAXTRY) ++nMax;  
	}

      // Found two most promising maxima.
      }
    }
  }

  // Read out starting position for search.
  sigmaMx = sigMax[0]; 
  int beginVar = (hasPointLeptons) ? 2 : 0;
  for (int iMax = 0; iMax < nMax; ++iMax) {
    int iTau = iMaxTau[iMax];
    int iY = iMaxY[iMax];
    int iZ = iMaxZ[iMax];
    double tauVal = 0.5;
    double yVal = 0.5;
    double zVal = 0.5;
    int iGrid;
    double varVal, varNew, deltaVar, marginVar, sigGrid[3];

    // Starting point and step size in parameter space.
    for (int iRepeat = 0; iRepeat < 2; ++iRepeat) {
      // Run through (possibly a subset of) tau, y and z.
      for (int iVar = beginVar; iVar < nVar; ++iVar) {
        if (iVar == 0) varVal = tauVal;
        else if (iVar == 1) varVal = yVal;
        else varVal = zVal;
        deltaVar = (iRepeat == 0) ? 0.1 
          : max( 0.01, min( 0.05, min( varVal - 0.2, 0.98 - varVal) ) );
        marginVar = (iRepeat == 0) ? 0.02 : 0.002;
        int moveStart = (iRepeat == 0 && iVar == 0) ? 0 : 1;
        for (int move = moveStart; move < 9; ++move) {
 
          // Define new parameter-space point by step in one dimension.
          if (move == 0) {
            iGrid = 1;
            varNew = varVal;
          } else if (move == 1) {
            iGrid = 2;
            varNew = varVal + deltaVar;
          } else if (move == 2) {
            iGrid = 0;
            varNew = varVal - deltaVar;
          } else if (sigGrid[2] >= max( sigGrid[0], sigGrid[1]) 
            && varVal + 2. * deltaVar < 1. - marginVar) {
            varVal += deltaVar;
            sigGrid[0] = sigGrid[1];
            sigGrid[1] = sigGrid[2];
            iGrid = 2;
            varNew = varVal + deltaVar;
          } else if (sigGrid[0] >= max( sigGrid[1], sigGrid[2]) 
            && varVal - 2. * deltaVar > marginVar) {
            varVal -= deltaVar;
            sigGrid[2] = sigGrid[1];
            sigGrid[1] = sigGrid[0];
            iGrid = 0;
            varNew = varVal - deltaVar;
          } else if (sigGrid[2] >= sigGrid[0]) {
            deltaVar *= 0.5;
            varVal += deltaVar;
            sigGrid[0] = sigGrid[1];
            iGrid = 1;
            varNew = varVal;
          } else {
            deltaVar *= 0.5;
            varVal -= deltaVar;
            sigGrid[2] = sigGrid[1];
            iGrid = 1;
            varNew = varVal;
	  }
 
          // Convert to relevant variables and find derived new limits.
          bool insideLimits = true;
          if (iVar == 0) {
            tauVal = varNew;
            selectTau( iTau, tauVal, is2);
            if (!limitY()) insideLimits = false;
            if (is2 && !limitZ()) insideLimits = false; 
            if (insideLimits) {
              selectY( iY, yVal);
              if (is2) selectZ( iZ, zVal);
	    }
	  } else if (iVar == 1) {
            yVal = varNew;
            selectY( iY, yVal);
	  } else if (iVar == 2) {
            zVal = varNew;
            selectZ( iZ, zVal);
	  }

          // Evaluate cross-section. 
          double sigmaNow = 0.;
          if (insideLimits) {  

            // Calculate cross section. Weight by phase-space volume
            // and Breit-Wigners for masses in 2 -> 2.
            if (is2) sigmaProcessPtr->set2Kin( x1H, x2H, sH, tH, m3, m4,
                                               runBW3H, runBW4H);
            else     sigmaProcessPtr->set1Kin( x1H, x2H, sH);
            sigmaNow = sigmaProcessPtr->sigmaPDF();
            sigmaNow *= wtTau * wtY * wtZ * wtBW;

            // Allow possibility for user to modify cross section.
            if (canModifySigma) sigmaNow 
              *= userHooksPtr->multiplySigmaBy( sigmaProcessPtr, this, false);

            // Optional printout. Protect against negative cross section.
            if (showSearch) os << " tau =" << setw(11) << tau << "  y =" 
	      << setw(11) << y << "  z =" << setw(11) << z
	      << "  sigma =" << setw(11) << sigmaNow << "\n";
            if (sigmaNow < 0.) sigmaNow = 0.; 
          }

          // Save new maximum. Final maximum.
          sigGrid[iGrid] = sigmaNow;
          if (sigmaNow > sigmaMx) sigmaMx = sigmaNow;
	}
      }
    }
  }
  sigmaMx *= SAFETYMARGIN;

  // Optional printout.
  if (showSearch) os << "\n Final maximum = "  << setw(11) << sigmaMx << endl;

  // Done.
  return true;
}

//*********

// Select a trial kinematics phase space point.
// Note: by In is meant the integral over the quantity multiplying 
// coefficient cn. The sum of cn is normalized to unity.

bool PhaseSpace::trialKin1or2(bool is2, bool inEvent, ostream& os) {

  // Choose tau according to h1(tau)/tau, where
  // h1(tau) = c0/I0 + (c1/I1) * 1/tau 
  // + (c2/I2) / (tau + tauResA) 
  // + (c3/I3) * tau / ((tau - tauResA)^2 + widResA^2)
  // + (c4/I4) / (tau + tauResB) 
  // + (c5/I5) * tau / ((tau - tauResB)^2 + widResB^2)
  // + (c6/I6) * tau / (1 - tau).
  if (!limitTau(is2)) return false;
  int iTau = 0;
  if (!hasPointLeptons) {
    double rTau = Rndm::flat(); 
    while (rTau > tauCoefSum[iTau]) ++iTau; 
  }
  selectTau( iTau, Rndm::flat(), is2);

  // Choose y according to h2(y), where
  // h2(y) = (c0/I0) * (y-ymin) + (c1/I1) * (ymax-y) 
  // + (c2/I2) * 1/cosh(y) + (c3/I3) * 1 / (1 - exp(y-ymax)) 
  // + (c4/I4) * 1 / (1 - exp(ymin-y)).
  if (!limitY()) return false;
  int iY = 0;
  if (!hasPointLeptons) {
    double rY = Rndm::flat(); 
    while (rY > yCoefSum[iY]) ++iY; 
  }
  selectY( iY, Rndm::flat());

  // Choose z = cos(thetaHat) according to h3(z), where
  // h3(z) = c0/I0 + (c1/I1) * 1/(A - z) + (c2/I2) * 1/(A + z) 
  // + (c3/I3) * 1/(A - z)^2 + (c4/I4) * 1/(A + z)^2,
  // where A = 1 + 2*(m3*m4/sH)^2 (= 1 for massless products).
  if (is2) {
    if (!limitZ()) return false;
    int iZ = 0;
    double rZ = Rndm::flat(); 
    while (rZ > zCoefSum[iZ]) ++iZ; 
    selectZ( iZ, Rndm::flat());
  }
   
  // Calculate cross section. Weight by phase-space volume
  // and Breit-Wigners for masses in 2 -> 2.
  if (is2) sigmaProcessPtr->set2Kin( x1H, x2H, sH, tH, m3, m4, 
                                     runBW3H, runBW4H);
  else     sigmaProcessPtr->set1Kin( x1H, x2H, sH);
  sigmaNw = sigmaProcessPtr->sigmaPDF();
  sigmaNw *= wtTau * wtY * wtZ * wtBW;

  // Allow possibility for user to modify cross section.
  if (canModifySigma) sigmaNw 
    *= userHooksPtr->multiplySigmaBy( sigmaProcessPtr, this, inEvent);

  // Check if maximum violated.
  if (sigmaNw > sigmaMx) {
    ErrorMsg::message("Warning in PhaseSpace2to2tauyz::trialKin: "
      "maximum for cross section violated");
    double violFact = SAFETYMARGIN * sigmaNw / sigmaMx;
    sigmaMx = SAFETYMARGIN * sigmaNw; 
    // Optional printout of (all) violations.
    if (showViolation) { 
      if (violFact < 9.99) cout << fixed;
      else                 cout << scientific;
      os << " Maximum for " << sigmaProcessPtr->name() 
         << " increased by factor " << setprecision(3) << violFact 
         << " to " << scientific << sigmaMx << endl;
    }
  }

  // Check if negative cross section.
  if (sigmaNw < sigmaNeg) {
    ErrorMsg::message("Warning in PhaseSpace2to2tauyz::trialKin:"
      " negative cross section set 0", "for " +  sigmaProcessPtr->name() );
    sigmaNeg = sigmaNw;
    // Optional printout of (all) violations.
    if (showViolation) os << " Negative minimum for " 
      << sigmaProcessPtr->name() << " changed to " << scientific 
      << setprecision(3) << sigmaNeg << endl;
  }
  if (sigmaNw < 0.) sigmaNw = 0.;

  // Done.
  return true;
}

//*********

// Find range of allowed tau values.

bool PhaseSpace::limitTau(bool is2) {

  // Trivial reply for unresolved lepton beams.
  if (hasPointLeptons) {
    tauMin = 1.;
    tauMax = 1.;
    return true;
  }

  // Requirements from allowed mHat range.
  tauMin = sHatMin / s; 
  tauMax = (mHatMax < mHatMin) ? 1. : min( 1., sHatMax / s); 

  // Requirements from allowed pT range and masses.
  if (is2) {
    double mT3Min = sqrt(s3 + pT2HatMinNow);
    double mT4Min = sqrt(s4 + pT2HatMinNow);
    tauMin = max( tauMin, pow2(mT3Min + mT4Min) / s);
  }
  
  // Check that there is an open range.
  return (tauMax > tauMin);
}

//*********

// Find range of allowed y values.

bool PhaseSpace::limitY() {

  // Trivial reply for unresolved lepton beams.
  if (hasPointLeptons) {
    yMax = 1.;
    return true;
  }

  // Requirements from selected tau value.
  yMax = -0.5 * log(tau); 

  // For lepton beams requirements from cutoff for f_e^e.
  double yMaxMargin = (hasLeptonBeams) ? yMax + LEPTONXLOGMAX : yMax;

  // Check that there is an open range.
  return (yMaxMargin > 0.);
}

//*********

// Find range of allowed z = cos(theta) values.

bool PhaseSpace::limitZ() {

  // Default limits.
  zMin = 0.;
  zMax = 1.;

  // Requirements from pTHat limits.
  zMax = sqrtpos( 1. - pT2HatMinNow / p2Abs );
  if (pTHatMax > pTHatMin) zMin = sqrtpos( 1. - pT2HatMax / p2Abs );
 
  // Check that there is an open range.
  return (zMax > zMin);
}

//*********

// Select tau according to a choice of shapes.

void PhaseSpace::selectTau(int iTau, double tauVal, bool is2) {

  // Trivial reply for unresolved lepton beams.
  if (hasPointLeptons) {
    tau = 1.;
    wtTau = 1.;
    sH = s;
    mHat = sqrt(sH);
    if (is2) {
      p2Abs = 0.25 * (pow2(sH - s3 - s4) - 4. * s3 * s4) / sH; 
      pAbs = sqrtpos( p2Abs );
    }
    return;
  }

  // Contributions from s-channel resonances.
  double tRatA = 0.;
  double aLowA = 0.;
  double aUppA = 0.;
  if (idResA !=0) {
    tRatA = ((tauResA + tauMax) / (tauResA + tauMin)) * (tauMin / tauMax);
    aLowA = atan( (tauMin - tauResA) / widResA);
    aUppA = atan( (tauMax - tauResA) / widResA);
  }
  double tRatB = 0.;
  double aLowB = 0.;
  double aUppB = 0.;
  if (idResB != 0) {
    tRatB = ((tauResB + tauMax) / (tauResB + tauMin)) * (tauMin / tauMax);
    aLowB = atan( (tauMin - tauResB) / widResB);
    aUppB = atan( (tauMax - tauResB) / widResB);
  }
 
  // Contributions from 1 / (1 - tau)  for lepton beams.
  double aLowT = 0.;
  double aUppT = 0.;
  if (hasLeptonBeams) { 
    aLowT = log( max( LEPTONTAUMIN, 1. - tauMin) );
    aUppT = log( max( LEPTONTAUMIN, 1. - tauMax) ); 
    intTau6 = aLowT - aUppT;
  }  

  // Select according to 1/tau or 1/tau^2.
  if (iTau == 0) tau = tauMin * pow( tauMax / tauMin, tauVal);
  else if (iTau == 1) tau = tauMax * tauMin 
    / (tauMin + (tauMax - tauMin) * tauVal);  

  // Select according to 1 / (1 - tau) for lepton beams.
  else if (hasLeptonBeams && iTau == nTau - 1) 
    tau = 1. - exp( aUppT + intTau6 * tauVal );

  // Select according to 1 / (tau * (tau + tauRes)) or 
  // 1 / ((tau - tauRes)^2 + widRes^2) for resonances A and B.
  else if (iTau == 2) tau = tauResA * tauMin 
    / ((tauResA + tauMin) * pow( tRatA, tauVal) - tauMin);
  else if (iTau == 3) tau = tauResA + widResA 
    * tan( aLowA + (aUppA - aLowA) * tauVal);
  else if (iTau == 4) tau = tauResB * tauMin 
    / ((tauResB + tauMin) * pow( tRatB, tauVal) - tauMin);
  else if (iTau == 5) tau = tauResB + widResB 
    * tan( aLowB + (aUppB - aLowB) * tauVal);

  // Phase space weight in tau.
  intTau0 = log( tauMax / tauMin);
  intTau1 = (tauMax - tauMin) / (tauMax * tauMin);
  double invWtTau = (tauCoef[0] / intTau0) + (tauCoef[1] / intTau1) / tau;
  if (idResA != 0) {
    intTau2 = -log(tRatA) / tauResA;
    intTau3 = (aUppA - aLowA) / widResA; 
    invWtTau += (tauCoef[2] / intTau2) / (tau + tauResA) 
      + (tauCoef[3] / intTau3) * tau / ( pow2(tau - tauResA) + pow2(widResA) );
  }
  if (idResB != 0) {
    intTau4 = -log(tRatB) / tauResB;
    intTau5 = (aUppB - aLowB) / widResB; 
    invWtTau += (tauCoef[4] / intTau4) / (tau + tauResB) 
      + (tauCoef[5] / intTau5) * tau / ( pow2(tau - tauResB) + pow2(widResB) );
  }
  if (hasLeptonBeams) 
    invWtTau += (tauCoef[nTau - 1] / intTau6) 
      * tau / max( LEPTONTAUMIN, 1. - tau);
  wtTau = 1. / invWtTau;

  // Calculate sHat and absolute momentum of outgoing partons.
  sH = tau * s;
  mHat = sqrt(sH);
  if (is2) {
    p2Abs = 0.25 * (pow2(sH - s3 - s4) - 4. * s3 * s4) / sH; 
    pAbs = sqrtpos( p2Abs );
  }

}

//*********

// Select y according to a choice of shapes.

void PhaseSpace::selectY(int iY, double yVal) {

  // Trivial reply for unresolved lepton beams.
  if (hasPointLeptons) {
    y = 0.;
    wtY = 1.;
    x1H = 1.;
    x2H = 1.;
    return;
  }

  // Standard expressions used below.
  double atanMax = atan( exp(yMax) );  
  double atanMin = atan( exp(-yMax) );  
  double aUppY = (hasLeptonBeams) 
    ? log( max( LEPTONXMIN, LEPTONXMAX / tau - 1. ) ) : 0.;
  double aLowY = LEPTONXLOGMIN;

  // y - y_min or mirrored y_max - y.
  if (iY <= 1) y = yMax * (2. * sqrt(yVal) - 1.); 

  // 1 / cosh(y).
  else if (iY == 2) 
    y = log( tan( atanMin + (atanMax - atanMin) * yVal ) );

  // 1 / (1 - exp(y - y_max)) or mirrored 1 / (1 - exp(y_min - y)).
  else y = yMax - log( 1. + exp(aLowY + (aUppY - aLowY) * yVal) );

  // Mirror two cases. 
  if (iY == 1 || iY == 4) y = -y;

  // Phase space integral in y.
  intY01 = 0.5 * pow2(2. * yMax);
  intY2  = 2. * (atanMax - atanMin);
  intY34 = aUppY - aLowY;
  double invWtY = (yCoef[0] / intY01) * (y + yMax)
    + (yCoef[1] / intY01) * (yMax - y) + (yCoef[2] / intY2) / cosh(y);
  if (hasLeptonBeams) invWtY 
    += (yCoef[3] / intY34) / max( LEPTONXMIN, 1. - exp( y - yMax) )
    +  (yCoef[4] / intY34) / max( LEPTONXMIN, 1. - exp(-y - yMax) );  
  wtY = 1. / invWtY;

  // Calculate x1 and x2.
  x1H = sqrt(tau) * exp(y);
  x2H = sqrt(tau) * exp(-y);
}

//*********

// Select z = cos(theta) according to a choice of shapes.
// The selection is split in the positive- and negative-z regions,
// since a pTmax cut can remove the region around z = 0.

void PhaseSpace::selectZ(int iZ, double zVal) {

  // Mass-dependent dampening of pT -> 0 limit.
  ratio34 = max(TINY, 2. * s3 * s4 / pow2(sH));
  unity34 = 1. + ratio34;
  double ratiopT2 = 2. * pT2HatMinNow / sH;
  // ?? constant 0.0001; what if > this ?? 
  if (ratiopT2 < 0.0001) ratio34 = max( ratio34, ratiopT2);

  // Common expressions in z limits.
  double zPosMax = max(ratio34, unity34 + zMax);
  double zNegMax = max(ratio34, unity34 - zMax);
  double zPosMin = max(ratio34, unity34 + zMin);
  double zNegMin = max(ratio34, unity34 - zMin);

  // Flat in z.
  if (iZ == 0) {
    if (zVal < 0.5) z = -(zMax + (zMin - zMax) * 2. * zVal);
    else z = zMin + (zMax - zMin) * (2. * zVal - 1.);

  // 1 / (unity34 - z).
  } else if (iZ == 1) {
    double areaNeg = log(zPosMax / zPosMin);
    double areaPos = log(zNegMin / zNegMax); 
    double area = areaNeg + areaPos;
    if (zVal * area < areaNeg) {
      double zValMod = zVal * area / areaNeg;
      z = unity34 - zPosMax * pow(zPosMin / zPosMax, zValMod);
    } else {
      double zValMod = (zVal * area - areaNeg)/ areaPos;
      z = unity34 - zNegMin * pow(zNegMax / zNegMin, zValMod);
    }

  // 1 / (unity34 + z).
  } else if (iZ == 2) {
    double areaNeg = log(zNegMin / zNegMax);
    double areaPos = log(zPosMax / zPosMin);
    double area = areaNeg + areaPos;
    if (zVal * area < areaNeg) {
      double zValMod = zVal * area / areaNeg;
      z = zNegMax * pow(zNegMin / zNegMax, zValMod) - unity34;
    } else {
      double zValMod = (zVal * area - areaNeg)/ areaPos;
      z = zPosMin * pow(zPosMax / zPosMin, zValMod) - unity34;
    }

  // 1 / (unity34 - z)^2.
  } else if (iZ == 3) {
    double areaNeg = 1. / zPosMin - 1. / zPosMax;
    double areaPos = 1. / zNegMax - 1. / zNegMin; 
    double area = areaNeg + areaPos;
    if (zVal * area < areaNeg) {
      double zValMod = zVal * area / areaNeg;
      z = unity34 - 1. / (1./zPosMax + areaNeg * zValMod);
    } else {
      double zValMod = (zVal * area - areaNeg)/ areaPos;
      z = unity34 - 1. / (1./zNegMin + areaPos * zValMod);
    }

  // 1 / (unity34 + z)^2.
  } else if (iZ == 4) {
    double areaNeg = 1. / zNegMax - 1. / zNegMin;
    double areaPos = 1. / zPosMin - 1. / zPosMax; 
    double area = areaNeg + areaPos;
    if (zVal * area < areaNeg) {
      double zValMod = zVal * area / areaNeg;
      z = 1. / (1./zNegMax - areaNeg * zValMod) - unity34;
    } else {
      double zValMod = (zVal * area - areaNeg)/ areaPos;
      z = 1. / (1./zPosMin - areaPos * zValMod) - unity34;
    }
  }

  // Safety check for roundoff errors. Combinations with z.
  if (z < 0.) z = min(-zMin, max(-zMax, z));
  else z = min(zMax, max(zMin, z));
  zNeg = max(ratio34, unity34 - z);
  zPos = max(ratio34, unity34 + z);

  // Phase space integral in z.
  double intZ0 = 2. * (zMax - zMin);
  double intZ12 = log( (zPosMax * zNegMin) / (zPosMin * zNegMax) ); 
  double intZ34 = 1. / zPosMin - 1. / zPosMax + 1. / zNegMax 
    - 1. / zNegMin;
  wtZ = mHat * pAbs / ( (zCoef[0] / intZ0) + (zCoef[1] / intZ12) / zNeg
    + (zCoef[2] / intZ12) / zPos + (zCoef[3] / intZ34) / pow2(zNeg)
    + (zCoef[4] / intZ34) / pow2(zPos) );

  // Calculate tHat and uHat. Also gives pTHat.
  double sH34 = -0.5 * (sH - s3 - s4);
  tH  = sH34 + mHat * pAbs * z;
  uH  = sH34 - mHat * pAbs * z;
  pTH = sqrtpos( (tH * uH - s3 * s4) / sH); 

}

//*********

// Solve linear equation system for better phase space coefficients.
  
void PhaseSpace::solveSys( int n, int bin[8], double vec[8], 
  double mat[8][8], double coef[8], ostream& os) {

  // Optional printout.
  if (showSearch) {
    os << "\n Equation system: " << setw(5) << bin[0];  
    for (int j = 0; j < n; ++j) os << setw(12) << mat[0][j];
    os << setw(12) << vec[0] << "\n";
    for (int i = 1; i < n; ++i) {
      os << "                  " << setw(5) << bin[i];  
      for (int j = 0; j < n; ++j) os << setw(12) << mat[i][j];
      os << setw(12) << vec[i] << "\n";
    }
  }

  // Local variables.
  double vecNor[8], coefTmp[8];
  for (int i = 0; i < n; ++i) coefTmp[i] = 0.;

  // Check if equation system solvable.
  bool canSolve = true;  
  for (int i = 0; i < n; ++i) if (bin[i] == 0) canSolve = false;
  double vecSum = 0.;
  for (int i = 0; i < n; ++i) vecSum += vec[i];
  if (abs(vecSum) < TINY) canSolve = false;

  // Solve to find relative importance of cross-section pieces.  
  if (canSolve) {
    for (int i = 0; i < n; ++i) vecNor[i] = max( 0.1, vec[i] / vecSum);
    for (int k = 0; k < n - 1; ++k) {
      for (int i = k + 1; i < n; ++i) {
        if (abs(mat[k][k]) < TINY) {canSolve = false; break;}
        double ratio = mat[i][k] / mat[k][k];
        vec[i] -= ratio * vec[k];
        for (int j = k; j < n; ++j) mat[i][j] -= ratio * mat[k][j];
      }  
      if (!canSolve) break;
    }
    if (canSolve) {
      for (int k = n - 1; k >= 0; --k) {
        for (int j = k + 1; j < n; ++j) vec[k] -= mat[k][j] * coefTmp[j];
        coefTmp[k] = vec[k] / mat[k][k]; 
      }
    }
  }

  // Share evenly if failure.
  if (!canSolve) for (int i = 0; i < n; ++i) {
    coefTmp[i] = 1.;
    vecNor[i] = 0.1;
    if (vecSum > TINY) vecNor[i] = max(0.1, vec[i] / vecSum);
  }

  // Normalize coefficients, with piece shared democratically.
  double coefSum = 0.;
  vecSum = 0.;
  for (int i = 0; i < n; ++i) {  
    coefTmp[i] = max( 0., coefTmp[i]);
    coefSum += coefTmp[i];
    vecSum += vecNor[i];
  }
  if (coefSum > 0.) for (int i = 0; i < n; ++i) coef[i] = EVENFRAC / n 
    + (1. - EVENFRAC) * 0.5 * (coefTmp[i] / coefSum + vecNor[i] / vecSum); 
  else for (int i = 0; i < n; ++i) coef[i] = 1. / n;

  // Optional printout.
  if (showSearch) {
    os << " Solution:             ";  
    for (int i = 0; i < n; ++i) os << setw(12) << coef[i];
    os << "\n";
  }
}

//**************************************************************************

// PhaseSpace2to1tauy class.
// 2 -> 1 kinematics for normal subprocesses.

//*********

// Set limits for resonance mass selection.

bool PhaseSpace2to1tauy::setupMass() {

  // Mass limits for current resonance.
  int idRes = abs(sigmaProcessPtr->resonanceA());
  double mResMin = (idRes == 0) ? 0. : ParticleDataTable::mMin(idRes);
  double mResMax = (idRes == 0) ? 0. : ParticleDataTable::mMax(idRes);

  // Compare with global mass limits and pick tighter of them.
  mHatMin = max( mResMin, mHatGlobalMin);
  sHatMin = mHatMin*mHatMin;
  mHatMax = eCM;  
  if (mResMax > mResMin) mHatMax = min( mHatMax, mResMax);
  if (mHatGlobalMax > mHatGlobalMin) mHatMax = min( mHatMax, mHatGlobalMax);
  sHatMax = mHatMax*mHatMax;

  // Default Breit-Wigner weight.
  wtBW = 1.;

  // Fail if mass window (almost) closed.
  return (mHatMax > mHatMin + MASSMARGIN); 

}

//*********

// Construct the four-vector kinematics from the trial values. 

bool PhaseSpace2to1tauy::finalKin() {

  // Particle masses; incoming always on mass shell.
  mH[1] = 0.;
  mH[2] = 0.;
  mH[3] = mHat;

  // Incoming partons along beam axes. Outgoing has sum of momenta.
  pH[1] = Vec4( 0., 0., 0.5 * eCM * x1H, 0.5 * eCM * x1H); 
  pH[2] = Vec4( 0., 0., -0.5 * eCM * x2H, 0.5 * eCM * x2H); 
  pH[3] = pH[1] + pH[2];

  // Done.
  return true;
}

//**************************************************************************

// PhaseSpace2to2tauyz class.
// 2 -> 2 kinematics for normal subprocesses.

//*********

// Set up for fixed or Breit-Wigner mass selection.
  
bool PhaseSpace2to2tauyz::setupMasses() {

  // Set sHat limits - based on global limits only.
  mHatMin = mHatGlobalMin;
  sHatMin = mHatMin*mHatMin;
  mHatMax = eCM;  
  if (mHatGlobalMax > mHatGlobalMin) mHatMax = min( eCM, mHatGlobalMax);
  sHatMax = mHatMax*mHatMax;

  // Masses and widths of resonances. 
  // id3Mass, id4Mass = 0 also for light q so that is obtains no mass.
  // gmZmode == 1 means pure photon propagator; set at lower mass limit.
  id3Mass = abs(sigmaProcessPtr->id3Mass());
  m3Peak = (id3Mass == 0) ? 0. : ParticleDataTable::m0(id3Mass);
  m3Width = (id3Mass == 0) ? 0. : ParticleDataTable::mWidth(id3Mass);
  m3Min = (id3Mass == 0) ? 0. : ParticleDataTable::mMin(id3Mass);
  m3Max = (id3Mass == 0) ? 0. : ParticleDataTable::mMax(id3Mass);
  if(id3Mass == 23 && gmZmode == 1) m3Peak = m3Min;
  s3Peak = m3Peak*m3Peak;
  useBW3 = useBreitWigners && (m3Width > minWidthBreitWigners);
  if (!useBW3) m3Width = 0.;
  mw3 = m3Peak * m3Width;
  wm3Rat = m3Width / m3Peak;
  id4Mass = abs(sigmaProcessPtr->id4Mass());
  m4Peak = (id4Mass == 0) ? 0. : ParticleDataTable::m0(id4Mass);
  m4Width = (id4Mass == 0) ? 0. : ParticleDataTable::mWidth(id4Mass);
  m4Min = (id4Mass == 0) ? 0. : ParticleDataTable::mMin(id4Mass);
  m4Max = (id4Mass == 0) ? 0. : ParticleDataTable::mMax(id4Mass);
  if(id4Mass == 23 && gmZmode == 1) m4Peak = m4Min;
  s4Peak = m4Peak*m4Peak;
  useBW4 = useBreitWigners && (m4Width > minWidthBreitWigners);
  if (!useBW4) m4Width = 0.;
  mw4 = m4Peak * m4Width;
  wm4Rat = m4Width / m4Peak;

  // If either particle is massless then need extra pTHat cut.
  pTHatMinNow = pTHatMin;
  if (m3Peak < pTHatMinDiverge || m4Peak < pTHatMinDiverge)
    pTHatMinNow = max( pTHatMin, pTHatMinDiverge);
  pT2HatMinNow = pow2(pTHatMinNow);

  // Reduced mass range. Extra safety margin when Breit-Wigners.
  m34Max = (mHatMax < mHatMin) ? eCM : min( mHatMax, eCM);
  if (useBW3) {
    m3Lower = m3Min;
    m3Upper = m34Max;
    m3Upper -= (useBW4) ? m4Min : m4Peak;
    if (m3Max > m3Min) m3Upper = min( m3Upper, m3Max);
    s3Lower = m3Lower*m3Lower; 
    s3Upper = m3Upper*m3Upper;
  }
  if (useBW4) {
    m4Lower = m4Min;
    m4Upper = m34Max;
    m4Upper -= (useBW3) ? m3Min : m3Peak;
    if (m4Max > m4Min) m4Upper = min( m4Upper, m4Max); 
    s4Lower = m4Lower*m4Lower; 
    s4Upper = m4Upper*m4Upper;
  }

  // If closed phase space then unallowed process.
  bool physical = true;
  if (useBW3 && m3Upper < m3Lower + MASSMARGIN) physical = false;
  if (useBW4 && m4Upper < m4Lower + MASSMARGIN) physical = false;
  if (!useBW3 && !useBW4 && m34Max < m3Peak + m4Peak + MASSMARGIN)
    physical = false;  
  if (!physical) return false;

  // Prepare to select m3 by BW + flat + 1/s_3.
  // Determine relative coefficients by allowed mass range. 
  if (useBW3) {
    double distToThreshA = (m34Max - m3Peak - m4Peak) * m3Width
      / (pow2(m3Width) + pow2(m4Width)); 
    double distToThreshB = (m34Max - m3Peak - m4Min) / m3Width;
    double distToThresh = min( distToThreshA, distToThreshB);
    if (distToThresh > THRESHOLDSIZE) {
      frac3Flat = 0.1;
      frac3Inv  = 0.1;
    } else if (distToThresh > - THRESHOLDSIZE) {
      frac3Flat = 0.25 - 0.15 * distToThresh / THRESHOLDSIZE; 
      frac3Inv  = 0.15 - 0.05 * distToThresh / THRESHOLDSIZE; 
    } else {          
      frac3Flat = 0.4;
      frac3Inv  = 0.2;
    }
    // For gamma*/Z0: increase 1/s_i part and introduce 1/s_i^2 part.
    frac3Inv2 = 0.;
    if (id3Mass == 23 && gmZmode == 0) {
      frac3Flat *= 0.5;
      frac3Inv  = 0.5 * frac3Inv + 0.25;
      frac3Inv2 = 0.25;
    } else if (id3Mass == 23 && gmZmode == 1) {
      frac3Flat = 0.1;
      frac3Inv  = 0.4;
      frac3Inv2 = 0.4;
    }
    // Normalization integrals for the respective contribution.
    atan3Lower = atan( (s3Lower - s3Peak)/ mw3 ); 
    atan3Upper = atan( (s3Upper - s3Peak)/ mw3 ); 
    int3BW   = atan3Upper - atan3Lower;
    int3Flat = s3Upper - s3Lower;
    int3Inv  = log( s3Upper / s3Lower );
    int3Inv2 = 1./s3Lower - 1./s3Upper;
  }

  // Prepare to select m4 by BW + flat + 1/s_4.
  // Determine relative coefficients by allowed mass range. 
  if (useBW4) {
    double distToThreshA = (m34Max - m3Peak - m4Peak) * m4Width
      / (pow2(m3Width) + pow2(m4Width)); 
    double distToThreshB = (m34Max - m3Min - m4Peak) / m4Width;
    double distToThresh = min( distToThreshA, distToThreshB);
    if (distToThresh > THRESHOLDSIZE) {
      frac4Flat = 0.1;
      frac4Inv  = 0.1;
    } else if (distToThresh > - THRESHOLDSIZE) {
      frac4Flat = 0.25 - 0.15 * distToThresh / THRESHOLDSIZE; 
      frac4Inv  = 0.15 - 0.05 * distToThresh / THRESHOLDSIZE; 
    } else {          
      frac4Flat = 0.4;
      frac4Inv  = 0.2;
    }
    frac4Inv2 = 0.;
    // For gamma*/Z0: increase 1/s_i part and introduce 1/s_i^2 part.
    if (id4Mass == 23 && gmZmode == 0) {
      frac4Flat *= 0.5;
      frac4Inv  = 0.5 * frac4Inv + 0.25;
      frac4Inv2 = 0.25;

    } else if (id4Mass == 23 && gmZmode == 1) {
      frac4Flat = 0.1;
      frac4Inv  = 0.4;
      frac4Inv2 = 0.4;
    }
    // Normalization integrals for the respective contribution.
    atan4Lower = atan( (s4Lower - s4Peak)/ mw4 ); 
    atan4Upper = atan( (s4Upper - s4Peak)/ mw4 ); 
    int4BW   = atan4Upper - atan4Lower;
    int4Flat = s4Upper - s4Lower;
    int4Inv  = log( s4Upper / s4Lower );
    int4Inv2 = 1./s4Lower - 1./s4Upper;
  }

  // Initialization masses. Special cases when constrained phase space.
  m3 = (useBW3) ? min(m3Peak, m3Upper) : m3Peak;
  m4 = (useBW4) ? min(m4Peak, m4Upper) : m4Peak;
  if (m3 + m4 + THRESHOLDSIZE * (m3Width + m4Width) + MASSMARGIN > m34Max) {
    if (useBW3 && useBW4) physical = constrainedM3M4();
    else if (useBW3) physical = constrainedM3();
    else if (useBW4) physical = constrainedM4();
  }
  s3 = m3*m3;
  s4 = m4*m4;

  // Correct selected mass-spectrum to running-width Breit-Wigner.
  // Extra safety margin for maximum search.
  weightMasses();
  if (useBW3) wtBW *= EXTRABWWTMAX;
  if (useBW4) wtBW *= EXTRABWWTMAX;

  // Done.
  return physical;
  
}

//*********

// Special choice of m3 and m4 when m34Max push them off mass shell.
// Vary x in expression m3 + m4 = m34Max - x * (Gamma3 + Gamma4).
// For each x try to put either 3 or 4 as close to mass shell as possible.
// Maximize BW_3 * BW_4 * beta_34, where latter approximate phase space. 

bool PhaseSpace2to2tauyz::constrainedM3M4() {

  // Initial values.
  bool foundNonZero = false;
  double wtMassMax = 0.;
  double m3WtMax = 0.;
  double m4WtMax = 0.;
  double xMax = (m34Max - m3Lower - m4Lower) / (m3Width + m4Width);
  double xStep = THRESHOLDSTEP * min(1., xMax);
  double xNow = 0.;
  double wtMassXbin, wtMassMaxOld, m34, mT34Min, wtMassNow, 
    wtBW3Now, wtBW4Now, beta34Now;
 
  // Step through increasing x values.
  do {
    xNow += xStep;
    wtMassXbin = 0.;
    wtMassMaxOld = wtMassMax;
    m34 = m34Max - xNow * (m3Width + m4Width);

    // Study point where m3 as close as possible to on-shell.
    m3 = min( m3Upper, m34 - m4Lower);
    if (m3 > m3Peak) m3 = max( m3Lower, m3Peak);
    m4 = m34 - m3;
    if (m4 < m4Lower) {m4 = m4Lower; m3 = m34 - m4;} 

    // Check that inside phase space limit set by pTmin.
    mT34Min = sqrt(m3*m3 + pT2HatMinNow) + sqrt(m4*m4 + pT2HatMinNow);
    if (mT34Min < m34Max) {

      // Breit-Wigners and beta factor give total weight.
      wtMassNow = 0.;
      if (m3 > m3Lower && m3 < m3Upper && m4 > m4Lower && m4 < m4Upper) {
        wtBW3Now = mw3 / ( pow2(m3*m3 - s3Peak) + pow2(mw3) );
        wtBW4Now = mw4 / ( pow2(m4*m4 - s4Peak) + pow2(mw4) );
        beta34Now = sqrt( pow2(m34Max*m34Max - m3*m3 - m4*m4) 
          - pow2(2. * m3 * m4) ) / (m34Max*m34Max);
        wtMassNow = wtBW3Now * wtBW4Now * beta34Now;
      } 

      // Store new maximum, if any.
      if (wtMassNow > wtMassXbin) wtMassXbin = wtMassNow; 
      if (wtMassNow > wtMassMax) {
        foundNonZero = true;
        wtMassMax = wtMassNow;
        m3WtMax = m3;
        m4WtMax = m4;
      }
    }    

    // Study point where m4 as close as possible to on-shell.
    m4 = min( m4Upper, m34 - m3Lower);
    if (m4 > m4Peak) m4 = max( m4Lower, m4Peak);
    m3 = m34 - m4;
    if (m3 < m3Lower) {m3 = m3Lower; m4 = m34 - m3;}

    // Check that inside phase space limit set by pTmin.
    mT34Min = sqrt(m3*m3 + pT2HatMinNow) + sqrt(m4*m4 + pT2HatMinNow);
    if (mT34Min < m34Max) {

      // Breit-Wigners and beta factor give total weight.
      wtMassNow = 0.;
      if (m3 > m3Lower && m3 < m3Upper && m4 > m4Lower && m4 < m4Upper) {
        wtBW3Now = mw3 / ( pow2(m3*m3 - s3Peak) + pow2(mw3) );
        wtBW4Now = mw4 / ( pow2(m4*m4 - s4Peak) + pow2(mw4) );
        beta34Now = sqrt( pow2(m34Max*m34Max - m3*m3 - m4*m4) 
          - pow2(2. * m3 * m4) ) / (m34Max*m34Max);
        wtMassNow = wtBW3Now * wtBW4Now * beta34Now;
      } 

      // Store new maximum, if any.
      if (wtMassNow > wtMassXbin) wtMassXbin = wtMassNow; 
      if (wtMassNow > wtMassMax) {
        foundNonZero = true;
        wtMassMax = wtMassNow;
        m3WtMax = m3;
        m4WtMax = m4;
      }    
    } 

  // Continue stepping if increasing trend and more x range available.
  } while ( (!foundNonZero || wtMassXbin > wtMassMaxOld)
    && xNow < xMax - xStep); 

  // Restore best values for subsequent maximization. Return.
  m3 = m3WtMax;
  m4 = m4WtMax;
  return foundNonZero;

}

//*********

// Special choice of m3 when m34Max pushes it off mass shell.
// Vary x in expression m3 = m34Max - m4 - x * Gamma3.
// Maximize BW_3 * beta_34, where latter approximate phase space. 

bool PhaseSpace2to2tauyz::constrainedM3() {

  // Initial values.  
  bool foundNonZero = false;
  double wtMassMax = 0.;
  double m3WtMax = 0.;
  double mT4Min = sqrt(m4*m4 + pT2HatMinNow);
  double xMax = (m34Max - m3Lower - m4) / m3Width;
  double xStep = THRESHOLDSTEP * min(1., xMax);
  double xNow = 0.;
  double wtMassNow, mT34Min, wtBW3Now, beta34Now;
 
  // Step through increasing x values; gives m3 unambiguously.
  do {
    xNow += xStep;
    wtMassNow = 0.;
    m3 = m34Max - m4 - xNow * m3Width;

    // Check that inside phase space limit set by pTmin.
    mT34Min = sqrt(m3*m3 + pT2HatMinNow) + mT4Min;
    if (mT34Min < m34Max) {

      // Breit-Wigner and beta factor give total weight.
      wtBW3Now = mw3 / ( pow2(m3*m3 - s3Peak) + pow2(mw3) );
      beta34Now = sqrt( pow2(m34Max*m34Max - m3*m3 - m4*m4) 
        - pow2(2. * m3 * m4) ) / (m34Max*m34Max);
      wtMassNow = wtBW3Now * beta34Now;

      // Store new maximum, if any.
      if (wtMassNow > wtMassMax) {
        foundNonZero = true;
        wtMassMax = wtMassNow;
        m3WtMax = m3;
      }    
    }
     
  // Continue stepping if increasing trend and more x range available.
  } while ( (!foundNonZero || wtMassNow > wtMassMax) 
    && xNow < xMax - xStep); 

  // Restore best value for subsequent maximization. Return.
  m3 = m3WtMax;
  return foundNonZero;

}

//*********

// Special choice of m4 when m34Max pushes it off mass shell.
// Vary x in expression m4 = m34Max - m3 - x * Gamma4.
// Maximize BW_4 * beta_34, where latter approximate phase space. 

bool PhaseSpace2to2tauyz::constrainedM4() {

  // Initial values.  
  bool foundNonZero = false;
  double wtMassMax = 0.;
  double m4WtMax = 0.;
  double mT3Min = sqrt(m3*m3 + pT2HatMinNow);
  double xMax = (m34Max - m4Lower - m3) / m4Width;
  double xStep = THRESHOLDSTEP * min(1., xMax);
  double xNow = 0.;
  double wtMassNow, mT34Min, wtBW4Now, beta34Now;
 
  // Step through increasing x values; gives m4 unambiguously.
  do {
    xNow += xStep;
    wtMassNow = 0.;
    m4 = m34Max - m3 - xNow * m4Width;

    // Check that inside phase space limit set by pTmin.
    mT34Min = mT3Min + sqrt(m4*m4 + pT2HatMinNow);
    if (mT34Min < m34Max) {

      // Breit-Wigner and beta factor give total weight.
      wtBW4Now = mw4 / ( pow2(m4*m4 - s4Peak) + pow2(mw4) );
      beta34Now = sqrt( pow2(m34Max*m34Max - m3*m3 - m4*m4) 
        - pow2(2. * m3 * m4) ) / (m34Max*m34Max);
      wtMassNow = wtBW4Now * beta34Now;
 
      // Store new maximum, if any.
      if (wtMassNow > wtMassMax) {
        foundNonZero = true;
        wtMassMax = wtMassNow;
        m4WtMax = m4;
      }
    }    
 
  // Continue stepping if increasing trend and more x range available.
  } while ( (!foundNonZero || wtMassNow > wtMassMax) 
    && xNow < xMax - xStep); 

  // Restore best value for subsequent maximization.
  m4 = m4WtMax;
  return foundNonZero;

}

//*********

// Select Breit-Wigner-distributed or fixed masses.
  
bool PhaseSpace2to2tauyz::trialMasses() {

  // By default vanishing cross section.
  sigmaNw = 0.;

  // Distribution for m_3 is BW + flat + 1/s_3 + 1/s_3^2.
  if (useBW3) { 
    double pickForm = Rndm::flat();
    if (pickForm > frac3Flat + frac3Inv + frac3Inv2)
      s3 = s3Peak + mw3 * tan( atan3Lower + Rndm::flat() * int3BW ); 
    else if (pickForm > frac3Inv + frac3Inv2) 
      s3 = s3Lower + Rndm::flat() * (s3Upper - s3Lower);
    else if (pickForm > frac3Inv2) 
      s3 = s3Lower * pow( s3Upper / s3Lower, Rndm::flat() ); 
    else s3 = s3Lower * s3Upper 
      / (s3Lower + Rndm::flat() * (s3Upper - s3Lower)); 
    m3 = sqrt(s3);
  // Else m_3 is fixed at peak value.
  } else {
    m3 = m3Peak;
    s3 = s3Peak;
  }

  // Distribution for m_4 is BW + flat + 1/s_4 + 1/s_4^2.
  if (useBW4) { 
    double pickForm = Rndm::flat();
    if (pickForm > frac4Flat + frac4Inv + frac4Inv2)
      s4 = s4Peak + mw4 * tan( atan4Lower + Rndm::flat() * int4BW ); 
    else if (pickForm > frac4Inv + frac4Inv2) 
      s4 = s4Lower + Rndm::flat() * (s4Upper - s4Lower);
    else if (pickForm > frac4Inv2) 
      s4 = s4Lower * pow( s4Upper / s4Lower, Rndm::flat() );  
    else s4 = s4Lower * s4Upper 
      / (s4Lower + Rndm::flat() * (s4Upper - s4Lower)); 
    m4 = sqrt(s4);
  // Else m_4 is fixed at peak value.
  } else {
    m4 = m4Peak;
    s4 = s4Peak;
  }

  // If outside phase space then reject event.
  if (m3 + m4 + MASSMARGIN > m34Max) return false; 

  // Correct selected mass-spectrum to running-width Breit-Wigner.
  weightMasses();

  // Done.
  return true;
}

//*********

// Naively a fixed-width Breit-Wigner is used to pick the mass.
// Here come the correction factors for
// (i) preselection according to BW + flat in s_i + 1/s_i + 1/s_i^2,
// (ii) reduced allowed mass range,
// (iii) running width, i.e. m0*Gamma0 -> s*Gamma0/m0.
// In the end, the weighted distribution is a running-width BW.

void PhaseSpace2to2tauyz::weightMasses() {

  wtBW = 1.;
  if (useBW3) { 
    double genBW  = (1. - frac3Flat - frac3Inv - frac3Inv2) 
      * mw3 / ( (pow2(s3 - s3Peak) + pow2(mw3)) * int3BW)
      + frac3Flat / int3Flat + frac3Inv / (s3 * int3Inv)
      + frac3Inv2 / (s3*s3 * int3Inv2);
    double mw3Run = s3 * wm3Rat;
    runBW3H       = mw3Run / (pow2(s3 - s3Peak) + pow2(mw3Run)) / M_PI;
    wtBW         *= runBW3H / genBW ;
  }  
  if (useBW4) { 
    double genBW  = (1. - frac4Flat - frac4Inv - frac4Inv2) 
      * mw4 / ( (pow2(s4 - s4Peak) + pow2(mw4)) * int4BW)
      + frac4Flat / int4Flat + frac4Inv / (s4 * int4Inv)
      + frac4Inv2 / (s4*s4 * int4Inv2);
    double mw4Run = s4 * wm4Rat;
    runBW4H       = mw4Run / (pow2(s4 - s4Peak) + pow2(mw4Run)) / M_PI;
    wtBW         *= runBW4H / genBW;
  } 

}

//*********

// Construct the four-vector kinematics from the trial values. 

bool PhaseSpace2to2tauyz::finalKin() {

  // Assign masses to particles assumed massless in matrix elements.
  int id3 = sigmaProcessPtr->id(3);
  int id4 = sigmaProcessPtr->id(4);
  if (id3Mass == 0) { m3 = ParticleDataTable::m0(id3); s3 = m3*m3; }
  if (id4Mass == 0) { m4 = ParticleDataTable::m0(id4); s4 = m4*m4; }

  // Sometimes swap tHat <-> uHat to reflect chosen final-state order. 
  if (sigmaProcessPtr->swappedTU()) {
    swap(tH, uH);
    z = -z;
  }

  // Check that phase space still open after new mass assignment.
  if (m3 + m4 + MASSMARGIN > mHat) return false; 
  p2Abs = 0.25 * (pow2(sH - s3 - s4) - 4. * s3 * s4) / sH; 
  pAbs = sqrtpos( p2Abs );

  // Particle masses; incoming always on mass shell.
  mH[1] = 0.;
  mH[2] = 0.;
  mH[3] = m3;
  mH[4] = m4;

  // Incoming partons along beam axes.
  pH[1] = Vec4( 0., 0., 0.5 * eCM * x1H, 0.5 * eCM * x1H); 
  pH[2] = Vec4( 0., 0., -0.5 * eCM * x2H, 0.5 * eCM * x2H); 

  // Outgoing partons initially in collision CM frame along beam axes.
  pH[3] = Vec4( 0., 0.,  pAbs, 0.5 * (sH + s3 - s4) / mHat); 
  pH[4] = Vec4( 0., 0., -pAbs, 0.5 * (sH + s4 - s3) / mHat); 

  // Then rotate and boost them to overall CM frame
  theta = acos(z);
  phi = 2. * M_PI * Rndm::flat();
  betaZ = (x1H - x2H)/(x1H + x2H);   
  pH[3].rot( theta, phi);
  pH[4].rot( theta, phi);
  pH[3].bst( 0., 0., betaZ);
  pH[4].bst( 0., 0., betaZ);
  pTH = pAbs * sin(theta);

  // Done.
  return true;
}

//**************************************************************************

// PhaseSpace2to2eldiff class.
// 2 -> 2 kinematics set up for elastic and diffractive scattering.

//*********

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Number of tries to find acceptable (m^2, t) set.
const int PhaseSpace2to2eldiff::NTRY = 500;

// Maximum positive/negative argument for exponentiation.
const double PhaseSpace2to2eldiff::EXPMAX = 50.;

// Safety margin so sum of diffractive masses not too close to eCM.
const double PhaseSpace2to2eldiff::DIFFMASSMAX = 1e-8;

//*********

// Form of phase space sampling already fixed, so no optimization.
// However, need to read out relevant parameters from SigmaTotal.

bool PhaseSpace2to2eldiff::setupSampling() {

  // Find maximum = value of cross section.
  sigmaNw = sigmaProcessPtr->sigmaHat();
  sigmaMx = sigmaNw;

  // Masses of particles and minimal masses of diffractive states.
  m3ElDiff = (diffA) ? sigmaTotPtr->mMinXB()  : mA; 
  m4ElDiff = (diffB) ? sigmaTotPtr->mMinAX()  : mB; 
  s1 = mA * mA;
  s2 = mB * mB;
  s3 = pow2( m3ElDiff);
  s4 = pow2( m4ElDiff);

  // Parameters of low-mass-resonance diffractive enhancement.
  cRes = sigmaTotPtr->cRes();
  sResXB = pow2( sigmaTotPtr->mResXB());
  sResAX = pow2( sigmaTotPtr->mResAX());
  sProton = sigmaTotPtr->sProton();  

  // Elastic slope and lower limit diffractive slope.
  if (!diffA && !diffB) bMin = sigmaTotPtr->bSlopeEl();
  else if (!diffB) bMin = sigmaTotPtr->bMinSlopeXB();
  else if (!diffA) bMin = sigmaTotPtr->bMinSlopeAX();
  else bMin = sigmaTotPtr->bMinSlopeXX(); 
 
  // Determine maximum possible t range and coefficient of generation.
  lambda12 = sqrtpos( pow2( s - s1 - s2) - 4. * s1 * s2 );
  lambda34 = sqrtpos( pow2( s - s3 - s4) - 4. * s3 * s4 );
  double tempA = s - (s1 + s2 + s3 + s4) + (s1 - s2) * (s3 - s4) / s;
  double tempB = lambda12 *  lambda34 / s;
  double tempC = (s3 - s1) * (s4 - s2) + (s1 + s4 - s2 - s3)
    * (s1 * s4 - s2 * s3) / s;
  tLow = -0.5 * (tempA + tempB); 
  tUpp = tempC / tLow; 
  tAux = exp( max(-EXPMAX, bMin * (tLow - tUpp)) ) - 1.; 

  return true;

}

//*********

// Select a trial kinematics phase space point. Perform full
// Monte Carlo acceptance/rejection at this stage.

bool PhaseSpace2to2eldiff::trialKin( bool ) {

  // Loop over attempts to set up masses and t consistently.
  for (int loop = 0; ; ++loop) { 
    if (loop == NTRY) {
      ErrorMsg::message("Error in PhaseSpace2to2eldiff::trialKin: "
        " quit after repeated tries");
      return false;
    }
  
    // Select diffractive mass/masses according to dm^2/m^2.
    m3 = (diffA) ? m3ElDiff * pow( max(mA, eCM - m4ElDiff) / m3ElDiff,
      Rndm::flat()) : m3ElDiff;  
    m4 = (diffB) ? m4ElDiff * pow( max(mB, eCM - m3ElDiff) / m4ElDiff,
      Rndm::flat()) : m4ElDiff;
    s3 = m3 * m3;
    s4 = m4 * m4; 
 
    // Additional mass factors, including resonance enhancement.
    if (m3 + m4 >= eCM) continue;  
    if (diffA && !diffB) {
      double facXB = (1. - s3 / s)  
        * (1. + cRes * sResXB / (sResXB + s3));
      if (facXB < Rndm::flat() * (1. + cRes)) continue; 
    } else if (diffB && !diffA) {
      double facAX = (1. - s4 / s)  
        * (1. + cRes * sResAX / (sResAX + s4));
      if (facAX < Rndm::flat() * (1. + cRes)) continue; 
    } else if (diffA && diffB) {
      double facXX = (1. - pow2(m3 + m4) / s)  
        * (s * sProton / (s * sProton + s3 * s4))
        * (1. + cRes * sResXB / (sResXB + s3))
        * (1. + cRes * sResAX / (sResAX + s4));
      if (facXX < Rndm::flat() * pow2(1. + cRes)) continue; 
    }

    // Select t according to exp(bMin*t) and correct to right slope.
    tH = tUpp + log(1. + tAux * Rndm::flat()) / bMin;
    if (diffA || diffB) {
      double bDiff = 0.;
      if (diffA && !diffB) bDiff = sigmaTotPtr->bSlopeXB(s3) - bMin;
      else if (!diffA) bDiff = sigmaTotPtr->bSlopeAX(s4) - bMin;
      else bDiff = sigmaTotPtr->bSlopeXX(s3, s4) - bMin;
      bDiff = max(0., bDiff);
      if (exp( max(-50., bDiff * (tH - tUpp)) ) < Rndm::flat()) continue; 
    }
 
    // Check whether m^2 and t choices are consistent.
    lambda34 = sqrtpos( pow2( s - s3 - s4) - 4. * s3 * s4 );
    double tempA = s - (s1 + s2 + s3 + s4) + (s1 - s2) * (s3 - s4) / s;
    double tempB = lambda12 *  lambda34 / s;
    if (tempB < DIFFMASSMAX) continue;
    double tempC = (s3 - s1) * (s4 - s2) + (s1 + s4 - s2 - s3)
      * (s1 * s4 - s2 * s3) / s;
    double tLowNow = -0.5 * (tempA + tempB); 
    double tUppNow = tempC / tLowNow; 
    if (tH < tLowNow || tH > tUppNow) continue;

    // Careful reconstruction of scattering angle.
    double cosTheta = min(1., max(-1., (tempA + 2. * tH) / tempB));
    double sinTheta = 2. * sqrtpos( -(tempC + tempA * tH + tH * tH) )
      / tempB;
    theta = asin( min(1., sinTheta));
    if (cosTheta < 0.) theta = M_PI - theta;

    // Found acceptable kinematics, so no more looping.
    break;
  }

  return true;

}

//*********

// Construct the four-vector kinematics from the trial values. 

bool PhaseSpace2to2eldiff::finalKin() {

  // Particle masses; incoming always on mass shell.
  mH[1] = mA;
  mH[2] = mB;
  mH[3] = m3;
  mH[4] = m4;

  // Incoming particles along beam axes.
  pAbs = 0.5 * lambda12 / eCM;
  pH[1] = Vec4( 0., 0.,  pAbs, 0.5 * (s + s1 - s2) / eCM); 
  pH[2] = Vec4( 0., 0., -pAbs, 0.5 * (s + s2 - s1) / eCM); 

  // Outgoing particles initially along beam axes.
  pAbs = 0.5 * lambda34 / eCM;
  pH[3] = Vec4( 0., 0.,  pAbs, 0.5 * (s + s3 - s4) / eCM); 
  pH[4] = Vec4( 0., 0., -pAbs, 0.5 * (s + s4 - s3) / eCM); 

  // Then rotate them
  phi = 2. * M_PI * Rndm::flat();
  pH[3].rot( theta, phi);
  pH[4].rot( theta, phi);

  // Set some further info for completeness.
  x1H = 1.;
  x2H = 1.;
  sH = s;
  uH = s1 + s2 + s3 + s4 - sH - tH;
  mHat = eCM;
  p2Abs = pAbs * pAbs;
  betaZ = 0.;
  pTH = pAbs * sin(theta);

  // Done.
  return true;
}

//**************************************************************************

} // end namespace Pythia8
