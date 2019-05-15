// Function definitions (not found in the header) for the 
// PhaseSpace and PhaseSpace2to2tauyz classes.
// Copyright C 2006 Torbjorn Sjostrand

#include "PhaseSpace.h"

namespace Pythia8 {

//**************************************************************************

// The PhaseSpace class.
// Base class for phase space generators.

//*********
 
// Definitions of static variables and functions.
// (Values will be overwritten in initStatic call, so are purely dummy.)

double PhaseSpace::mHatMin = 4.;
double PhaseSpace::mHatMax = -1.;
double PhaseSpace::pTHatMin = 2.;
double PhaseSpace::pTHatMax = -1.;
double PhaseSpace::m3Min = 2.;
double PhaseSpace::m3Max = -1.;
double PhaseSpace::m4Min = 2.;
double PhaseSpace::m4Max = -1.;
double PhaseSpace::sHatMin = 16.;
double PhaseSpace::sHatMax = 1.;
double PhaseSpace::pT2HatMin = 4.;
double PhaseSpace::pT2HatMax = 1.;
double PhaseSpace::m3SMin = 4.;
double PhaseSpace::m3SMax = 1.;
double PhaseSpace::m4SMin = 4.;
double PhaseSpace::m4SMax = 1.;
bool PhaseSpace::showSearch = false;

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Maximum cross section increase, just in case true maximum not found.
const double PhaseSpace::SAFETYMARGIN = 1.05;

// Small number to avoid division by zero.
const double PhaseSpace::TINY = 1e-20;

// Fraction of total weight that is shared evenly between all shapes.
const double PhaseSpace::EVENFRAC = 0.4;

// Two cross sections with a small relative error are assumed same.
const double PhaseSpace::SAMESIGMA = 1e-6;

// Do not include resonances peaked too far outside allowed mass region.
const double PhaseSpace::WIDTHMARGIN = 20.;

// Special optimization treatment when two resonances at almost same mass.
const double PhaseSpace::SAMEMASS = 0.01;

//*********

// Initialize static data members.

void PhaseSpace::initStatic() {

  // Phase space cuts.
  mHatMin = Settings::parameter("PhaseSpace:mHatMin");
  mHatMax = Settings::parameter("PhaseSpace:mHatMax");
  pTHatMin = Settings::parameter("PhaseSpace:pTHatMin");
  pTHatMax = Settings::parameter("PhaseSpace:pTHatMax");
  m3Min = Settings::parameter("PhaseSpace:m3Min");
  m3Max = Settings::parameter("PhaseSpace:m3Max");
  m4Min = Settings::parameter("PhaseSpace:m4Min");
  m4Max = Settings::parameter("PhaseSpace:m4Max");

  // Quadratic limits often needed, so calculate them.
  sHatMin = pow2(mHatMin);
  sHatMax = pow2(mHatMax);
  pT2HatMin = pow2(pTHatMin);
  pT2HatMax = pow2(pTHatMax);
  m3SMin = pow2(m3Min);
  m3SMax = pow2(m3Max);
  m4SMin = pow2(m4Min);
  m4SMax = pow2(m4Max);

  // Print flag for maximization information.
  showSearch = Settings::flag("PhaseSpace:showSearch");
  
}

//*********

// Save pointers and values.

  void PhaseSpace::initInfo(SigmaProcess* sigmaProcessPtrIn, double eCMIn) {

  // Store input pointers for future use. CM energy.
  sigmaProcessPtr = sigmaProcessPtrIn;
  eCM = eCMIn;
  s = eCM * eCM;

  // Default event-specific kinematics properties.
  x1H = 1.;
  x2H = 1.;
  m3 = 0.;
  m4 = 0.;
  m3S = m3 * m3;
  m4S = m4 * m4;
  mHat = eCM;
  sH = s;
  tH = 0.;
  uH = 0.;
  theta = 0.;

  // Default cross section information.
  sigmaNw = 0.;
  sigmaMx = 0.;

}

//*********

// Determine how phase space should be sampled.

bool PhaseSpace::setupSampling1or2(bool is2) {

  // Optional printout.
  if (showSearch) cout <<  "\n Optimization printout for "  
    << sigmaProcessPtr->name() << "\n \n" << scientific << setprecision(3);

  // Set masses. Check that open range in tau (+ set tauMin, tauMax).
  if (is2) {
    m3 = sigmaProcessPtr->m(3);
    m3S = m3 * m3;
    m4 = sigmaProcessPtr->m(4);
    m4S = m4 * m4;
  }
  if (!limitTau(is2)) return false; 

  // Reset coefficients and matrices of equation system to solve.
  int binTau[8], binY[8], binZ[8];
  double vecTau[8], matTau[8][8], vecY[8], matY[8][8], vecZ[8], matZ[8][8];
  for (int i = 0; i < 8; ++i) {
    tauCoef[i] = 0.;
    yCoef[i] = 0.;
    zCoef[i] = 0.;
    binTau[i] = 0;
    binY[i] = 0;
    binZ[i] = 0;
    vecTau[i] = 0.;
    vecY[i] = 0.;
    vecZ[i] = 0.;
    for (int j = 0; j < 8; ++j) { 
      matTau[i][j] = 0.;
      matY[i][j] = 0.;
      matZ[i][j] = 0.;
    }  
  }
  sigmaMx = 0.;
  
  // Number of used coefficients/points for each dimension: tau, y, c.
  int nTau = 2;
  int nY = 3;
  int nZ = (is2)? 5 : 1; 

  // Identify if any resonances contribute in s-channel.
  idResA = sigmaProcessPtr->resonanceA();
  if (idResA != 0) { 
     mResA = ParticleDataTable::m0(idResA);
     GammaResA = ParticleDataTable::width(idResA);
     if (mHatMin > mResA + WIDTHMARGIN * GammaResA || (mHatMax > 0. 
       && mHatMax < mResA - WIDTHMARGIN * GammaResA) ) idResA = 0; 
  }
  idResB = sigmaProcessPtr->resonanceB();
  if (idResB != 0) { 
     mResB = ParticleDataTable::m0(idResB);
     GammaResB = ParticleDataTable::width(idResB);
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
  if (idResA !=0) {
    nTau += 2;
    tauResA = mResA * mResA / s;
    widResA = mResA * GammaResA / s;
    tRatA = ((tauResA + tauMax) / (tauResA + tauMin)) * (tauMin / tauMax);
    aUppA = atan( (tauMax - tauResA) / widResA);
    aLowA = atan( (tauMin - tauResA) / widResA);
  }
  if (idResB != 0) {
    nTau += 2;
    tauResB = mResB * mResB / s;
    widResB = mResB * GammaResB / s;
    tRatB = ((tauResB + tauMax) / (tauResB + tauMin)) * (tauMin / tauMax);
    aUppB = atan( (tauMax - tauResB) / widResB);
    aLowB = atan( (tauMin - tauResB) / widResB);
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
  yCoef[0] = 0.5;
  yCoef[1] = 0.5;
  zCoef[0] = 1.; 

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

        // Calculate cross section. Weight by phase-space volume.
        if (is2) sigmaProcessPtr->set2Kin( x1H, x2H, sH, tH);
        else     sigmaProcessPtr->set1Kin( x1H, x2H, sH);
        double sigmaNow = sigmaProcessPtr->sigmaPDF();
        sigmaNow *= wtTau * wtY * wtZ;
        if (sigmaNow > sigmaMx) sigmaMx = sigmaNow; 

        // Optional printout.
        if (showSearch) cout << " tau =" << setw(11) << tau << "  y =" 
	  << setw(11) << y << "  z =" << setw(11) << z
	  << "  sigma =" << setw(11) << sigmaNow << "\n";

        // Sum up tau cross-section pieces in points used.
        binTau[iTau] += 1;
        vecTau[iTau] += sigmaNow;
        matTau[iTau][0] += 1. / intTau0;
        matTau[iTau][1] += (1. /intTau1) / tau;
        if (idResA != 0) {
          matTau[iTau][2] += (1. /intTau2) / (tau + tauResA);
          matTau[iTau][3] += (1. /intTau3) 
            * tau / ( pow2(tau - tauResA) + pow2(widResA) );
        }
        if (idResB != 0) {
          matTau[iTau][4] += (1. /intTau4) / (tau + tauResB);
          matTau[iTau][5] += (1. /intTau5) 
            * tau / ( pow2(tau - tauResB) + pow2(widResB) );
        }
 
        // Sum up y cross-section pieces in points used.
        binY[iY] += 1;
        vecY[iY] += sigmaNow;
        matY[iY][0] += (1. / intY01) * (y + yMax);
        matY[iY][1] += (1. / intY01) * (yMax - y);
        matY[iY][2] += (1. / intY2) /cosh(y);

        // Integrals over z expressions at tauMax, to be used below.
        if (is2) {
          double p2AbsMax = 0.25 * (pow2(tauMax * s - m3S - m4S) 
            - 4. * m3S * m4S) / (tauMax * s);         
          double zMaxMax = sqrt( max( 0., 1. - pT2HatMin / p2AbsMax) );
          double zPosMaxMax = max(ratio34, unity34 + zMaxMax);
          double zNegMaxMax = max(ratio34, unity34 - zMaxMax);
          double intZ0Max = 2. * zMaxMax;
          double intZ12Max = log( zPosMaxMax / zNegMaxMax);
          double intZ34Max = 1. / zNegMaxMax - 1. / zPosMaxMax;  
  
          // Sum up z cross-section pieces in points used.
          binZ[iZ] += 1;
          vecZ[iZ] += sigmaNow;
          matZ[iZ][0] += 1. / intZ0Max; 
          matZ[iZ][1] += (1. / intZ12Max) / zNeg;
          matZ[iZ][2] += (1. / intZ12Max) / zPos;
          matZ[iZ][3] += (1. / intZ34Max) / pow2(zNeg);
          matZ[iZ][4] += (1. / intZ34Max) / pow2(zPos);
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
  solveSys( nTau, binTau, vecTau, matTau, tauCoef);
  solveSys( nY, binY, vecY, matY, yCoef);
  if (is2) solveSys( nZ, binZ, vecZ, matZ, zCoef);
  if (showSearch) cout << "\n";

  // Begin find two most promising maxima among same points as before.
  int iMaxTau[4], iMaxY[4], iMaxZ[4];
  double sigMax[4];
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

        // Calculate cross section. Weight by phase-space volume.
        if (is2) sigmaProcessPtr->set2Kin( x1H, x2H, sH, tH);
        else     sigmaProcessPtr->set1Kin( x1H, x2H, sH);
        double sigmaNow = sigmaProcessPtr->sigmaPDF();
        sigmaNow *= wtTau * wtY * wtZ;

        // Optional printout.
        if (showSearch) cout << " tau =" << setw(11) << tau << "  y =" 
	  << setw(11) << y << "  z =" << setw(11) << z
	  << "  sigma =" << setw(11) << sigmaNow << "\n";

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
          if (nMax < 2) ++nMax;  
	}

      // Found two most promising maxima.
      }
    }
  }

  // Read out starting position for search.
  sigmaMx = sigMax[0]; 
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
      for (int iVar = 0; iVar < nVar; ++iVar) {
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
            & varVal + 2. * deltaVar < 1. - marginVar) {
            varVal += deltaVar;
            sigGrid[0] = sigGrid[1];
            sigGrid[1] = sigGrid[2];
            iGrid = 2;
            varNew = varVal + deltaVar;
          } else if (sigGrid[0] >= max( sigGrid[1], sigGrid[2]) 
            & varVal - 2. * deltaVar > marginVar) {
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

            // Calculate cross section. Weight by phase-space volume.
            if (is2) sigmaProcessPtr->set2Kin( x1H, x2H, sH, tH);
            else     sigmaProcessPtr->set1Kin( x1H, x2H, sH);
            sigmaNow = sigmaProcessPtr->sigmaPDF();
            sigmaNow *= wtTau * wtY * wtZ;

            // Optional printout.
            if (showSearch) cout << " tau =" << setw(11) << tau << "  y =" 
	      << setw(11) << y << "  z =" << setw(11) << z
	      << "  sigma =" << setw(11) << sigmaNow << "\n";
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
  if (showSearch) cout << "\n Final maximum = "  << setw(11) << sigmaMx << endl;

  // Done.
  return true;
}

//*********

// Select a trial kinematics phase space point.
// Note: by In is meant the integral over the quantity multiplying 
// coefficient cn. The sum of cn is normalized to unity.

bool PhaseSpace::trialKin1or2(bool is2) {

  // Choose tau according to h1(tau)/tau, where
  // h1(tau) = c0/I0 + (c1/I1) * 1/tau 
  // + (c2/I2) / (tau + tauResA) + (c3/I3) * tau / ((tau - tauResA)^2 + widResA^2) 
  // + (c4/I4) / (tau + tauResB) + (c5/I5) * tau / ((tau - tauResB)^2 + widResB^2)
  if (!limitTau(is2)) return false;
  int iTau = 0;
  double rTau = Rndm::flat(); 
  if (rTau > tauCoef[0]) iTau = 1;
  if (rTau > tauCoef[0] +  tauCoef[1]) iTau = 2;
  if (rTau > tauCoef[0] +  tauCoef[1] + tauCoef[2]) iTau = 3;
  if (rTau > tauCoef[0] +  tauCoef[1] + tauCoef[2] + tauCoef[3]) iTau = 4;
  if (rTau > tauCoef[0] +  tauCoef[1] + tauCoef[2] + tauCoef[3] + tauCoef[4]) 
    iTau = 5;
  selectTau( iTau, Rndm::flat(), is2);

  // Choose y according to h2(y), where
  // h2(y) = (c0/I0) * (y-ymin) + (c1/I1) * (ymax-y) + (c2/I2) * 1/cosh(y).
  if (!limitY()) return false;
  int iY = 0;
  double rY = Rndm::flat(); 
  if (rY > yCoef[0]) iY = 1;
  if (rY > yCoef[0] + yCoef[1]) iY = 2;
  selectY( iY, Rndm::flat());

  // Choose z = cos(thetaHat) according to h3(z), where
  // h3(z) = c0/I0 + (c1/I1) * 1/(A - z) + (c2/I2) * 1/(A + z) 
  // + (c3/I3) * 1/(A - z)^2 + (c4/I4) * 1/(A + z)^2,
  // where A = 1 + 2*(m3*m4/sH)^2 (= 1 for massless products).
  if (is2) {
    if (!limitZ()) return false;
    int iZ = 0;
    double rZ = Rndm::flat(); 
    if (rZ > zCoef[0]) iZ = 1;
    if (rZ > zCoef[0] + zCoef[1]) iZ = 2;
    if (rZ > zCoef[0] + zCoef[1] + zCoef[2]) iZ = 3;
    if (rZ > zCoef[0] + zCoef[1] + zCoef[2] + zCoef[3]) iZ = 4;
    selectZ( iZ, Rndm::flat());
  }
   
  // Calculate cross section. Weight by phase-space volume.
  if (is2) sigmaProcessPtr->set2Kin( x1H, x2H, sH, tH);
  else     sigmaProcessPtr->set1Kin( x1H, x2H, sH);
  sigmaNw = sigmaProcessPtr->sigmaPDF();
  sigmaNw *= wtTau * wtY * wtZ;

  // Check if maximum violated.
  if (sigmaNw > sigmaMx) {
    ErrorMessages::message("Error in PhaseSpace2to2tauyz::trialKin: "
      "maximum for cross section violated");
    sigmaMx = SAFETYMARGIN * sigmaNw; 
  }

  // Done.
  return true;
}

//*********

// Find range of allowed tau values.

bool PhaseSpace::limitTau(bool is2) {

  // Requirements from allowed mHat range.
  tauMin = sHatMin / s; 
  tauMax = (mHatMax <= 0.) ? 1. : min( 1., sHatMax / s); 

  // Requirements from allowed pT range and masses.
  if (is2) {
    double mT3Min = sqrt(m3S + pT2HatMin);
    double mT4Min = sqrt(m4S + pT2HatMin);
    tauMin = max( tauMin, pow2(mT3Min + mT4Min) / s);
  }

  // Check that there is an open range.
  return (tauMax > tauMin) ? true : false;
}

//*********

// Find range of allowed y values.

bool PhaseSpace::limitY() {

  // Requirements from selected tau value.
  yMax = -0.5 * log(tau); 

  // Check that there is an open range.
  return (yMax > 0.) ? true : false;
}

//*********

// Find range of allowed z = cos(theta) values.

bool PhaseSpace::limitZ() {

  // Default limits.
  zMin = 0.;
  zMax = 1.;

  // Requirements from pTHat limits.
  zMax = sqrt( max( 0., 1. - pT2HatMin / p2Abs) );
  if (pTHatMax > 0.) zMin = sqrt( max( 0., 1. - pT2HatMax / p2Abs) );
 
  // Check that there is an open range.
  return (zMax > zMin) ? true : false;
}

//*********

// Select tau according to a choice of shapes.

void PhaseSpace::selectTau(int iTau, double tauVal, bool is2) {

  // Select according to 1/tau or 1/tau^2.
  if (iTau == 0) tau = tauMin * pow( tauMax / tauMin, tauVal);
  else if (iTau == 1) tau = tauMax * tauMin 
    / (tauMin + (tauMax - tauMin) * tauVal);  

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
  wtTau = 1. / invWtTau;

  // Calculate sHat and absolute momentum of outgoing partons.
  sH = tau * s;
  mHat = sqrt(sH);
  if (is2) {
    p2Abs = 0.25 * (pow2(sH - m3S - m4S) - 4. * m3S * m4S) / sH; 
    pAbs = sqrt( max(0., p2Abs) );
  }
}

//*********

// Select y according to a choice of shapes.

void PhaseSpace::selectY(int iY, double yVal) {

  // y - y_min, y_max - y.
  if (iY <= 1) {
    y = yMax * (2. * sqrt(yVal) - 1.);
    if (iY == 1) y = -y;

  // 1 / cosh(y).
  } else if (iY == 2) {
    double atanMax = atan( exp(yMax) );  
    double atanMin = atan( exp(-yMax) );  
    y = log( tan( atanMin + (atanMax - atanMin) * yVal ) );
  } 

  // Phase space integral in y.
  intY01 = 0.5 * pow2(2. * yMax);
  intY2 = 2. * (atan(exp(yMax)) - atan(exp(-yMax)));
  wtY = 1. / ( (yCoef[0] / intY01) * (y + yMax)
    + (yCoef[1] / intY01) * (yMax - y) + (yCoef[2] / intY2) / cosh(y) );

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
  ratio34 = max(TINY, 2. * m3S * m4S / pow2(sH));
  unity34 = 1. + ratio34;
  double ratiopT2 = 2. * pT2HatMin / sH;
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
    double areaNeg = log(zNegMax / zNegMin);
    double areaPos = log(zPosMin / zPosMax); 
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

  // Calculate tHat and uHat. 
  double sH34 = -0.5 * (sH - m3S - m4S);
  tH = sH34 + mHat * pAbs * z;
  uH = sH34 - mHat * pAbs * z;
}

//*********

// Solve linear equation system for better phase space coefficients.
  
void PhaseSpace::solveSys( int n, int bin[8], 
  double vec[8], double mat[8][8], double coef[8]) {

  // Optional printout.
  if (showSearch) {
    cout << "\n Equation system: " << setw(5) << bin[0];  
    for (int j = 0; j < n; ++j) cout << setw(12) << mat[0][j];
    cout << setw(12) << vec[0] << "\n";
    for (int i = 1; i < n; ++i) {
      cout << "                  " << setw(5) << bin[i];  
      for (int j = 0; j < n; ++j) cout << setw(12) << mat[i][j];
      cout << setw(12) << vec[i] << "\n";
    }
  }

  // Local variables.
  double vecNor[8], coefTmp[8];
  for (int i = 0; i < n; ++i) coefTmp[i] = 0.;

  // Check if equation sustem solvable.
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
    cout << " Solution:             ";  
    for (int i = 0; i < n; ++i) cout << setw(12) << coef[i];
    cout << "\n";
  }
}

//**************************************************************************

// PhaseSpace2to1tauy class.
// 2 -> 1 kinematics for normal subprocesses.

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

// Construct the four-vector kinematics from the trial values. 

bool PhaseSpace2to2tauyz::finalKin() {

  // Particle masses; incoming always on mass shell.
  mH[1] = 0.;
  mH[2] = 0.;
  mH[3] = m3;
  mH[4] = m4;

  // Incoming partons along beam axes.
  pH[1] = Vec4( 0., 0., 0.5 * eCM * x1H, 0.5 * eCM * x1H); 
  pH[2] = Vec4( 0., 0., -0.5 * eCM * x2H, 0.5 * eCM * x2H); 

  // Outgoing partons initially in collision CM frame along beam axes.
  pH[3] = Vec4( 0., 0.,  pAbs, 0.5 * (sH + m3S - m4S) / mHat); 
  pH[4] = Vec4( 0., 0., -pAbs, 0.5 * (sH + m4S - m3S) / mHat); 

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
  lambda12 = sqrt( max( 0., pow2( s - s1 - s2) - 4. * s1 * s2));
  lambda34 = sqrt( max( 0., pow2( s - s3 - s4) - 4. * s3 * s4));
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

bool PhaseSpace2to2eldiff::trialKin() {

  // Loop over attempts to set up masses and t consistently.
  for (int loop = 0; ; ++loop) { 
    if (loop == NTRY) {
      ErrorMessages::message("Error in PhaseSpace2to2eldiff::trialKin: "
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
    lambda34 = sqrt( max( 0., pow2( s - s3 - s4) - 4. * s3 * s4));
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
    double sinTheta = 2. * sqrt( max(0., -(tempC + tempA * tH + tH * tH)))
      / tempB;
    theta = asin( min(1., sinTheta));
    if (cosTheta < 0.) theta = M_PI - theta;

    // Found acceptable kinematics, so no more looping.
    break;
  }
  m3S = s3;
  m4S = s4;

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
