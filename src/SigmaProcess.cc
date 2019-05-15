// Function definitions (not found in the header) for the 
// SigmaProcess class, and classes derived from it.
// Copyright C 2006 Torbjorn Sjostrand

#include "SigmaProcess.h"

namespace Pythia8 {

//**************************************************************************

// The SigmaProcess class.
// Base class for cross sections.

//*********
 
// Definitions of static variables and functions.
// (Values will be overwritten in initStatic call, so are purely dummy.)

int SigmaProcess::alphaSorder = 1;
int SigmaProcess::nQuark = 3;
double SigmaProcess::alphaSvalue = 0.1265;
AlphaStrong SigmaProcess::alphaScalc;
AlphaEM SigmaProcess::alphaEMcalc;

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Conversion of GeV^{-2} to mb for cross section.
const double SigmaProcess::CONVERT2MB = 0.389380; 

// The sum of outgoing masses must not be too close to the cm energy.
const double SigmaProcess::MASSMARGIN = 0.1;

// Information on incoming beams.
int SigmaProcess::idA, SigmaProcess::idB;
double SigmaProcess::mA, SigmaProcess::mB; 
  
// Pointer to the total/elastic/diffractive cross section object.
SigmaTotal* SigmaProcess::sigmaTotPtr;

// Pointer to the SLHA object
SusyLesHouches* SigmaProcess::slha;

//*********

// Initialize static data members.

void SigmaProcess::initStatic() {

  // Parameters of alphaStrong generation .
  alphaSvalue = Settings::parm("SigmaProcess:alphaSvalue");
  alphaSorder = Settings::mode("SigmaProcess:alphaSorder");

  // Initialize alphaStrong generation.
  alphaScalc.init( alphaSvalue, alphaSorder); 

  // Maximum new quark flavour.
  nQuark = Settings::mode("SigmaProcess:nQuark");

}

//**************************************************************************

// The Sigma1Process class.
// Base class for resolved 2 -> 1 cross sections; derived from SigmaProcess.

//*********

// Input and complement kinematics for resolved 2 -> 1 process. 

bool Sigma1Process::set1Kin( double x1in, double x2in, double sHin) {

  // Default value only sensible for these processes.
  swapTU = false;

  // Incoming parton momentum fractions and sHat.
  x1 = x1in;
  x2 = x2in;
  sH = sHin;
  sH2 = sH * sH;

  // Use sHat as renormalization scale. Evaluate alpha_strong and alpha_EM.
  Q2RenH = sH;
  alpS = alphaScalc.alphaS(Q2RenH);  
  alpEM = alphaEMcalc.alphaEM(Q2RenH);  

  // Use sHat as factorization scale.
  Q2FacH = sH;

  // Done.
  return true;

}

//**************************************************************************

// The Sigma2Process class.
// Base class for resolved 2 -> 2 cross sections; derived from SigmaProcess.

//*********

// Input and complement kinematics for resolved 2 -> 2 process. 

bool Sigma2Process::set2Kin( double x1in, double x2in, double sHin, 
  double tHin, double m3in, double m4in) {

  // Default ordering of particles 3 and 4.
  swapTU = false;

  // Incoming parton momentum fractions.
  x1 = x1in;
  x2 = x2in;

  // Incoming masses and their squares.
  bool masslessKin = (id3Mass() == 0) && (id4Mass() == 0);
  if (masslessKin) {
    m3 = 0.;
    m4 = 0.;
  } else {
    m3 = m3in;
    m4 = m4in;
  }
  mH[3] = m3;
  mH[4] = m4;
  s3 = m3*m3;
  s4 = m4*m4;

  // Standard Mandelstam variables and their squares.
  sH = sHin;
  tH = tHin;
  uH = (masslessKin) ? -(sH + tH) : s3 + s4 - (sH + tH); 
  sH2 = sH * sH;
  tH2 = tH * tH;
  uH2 = uH * uH;

  // Calculate squared transverse momentum.
  pT2 = (masslessKin) ?  tH * uH / sH : (tH * uH - s3 * s4) / sH;

  // Use pT^2 as renormalization scale, generalized to m_T3 * m_T4
  // for massive case. Evaluate alpha_strong and alpha_EM.
  // Is this correct for processes with t-channel W/Z exchange??
  Q2RenH = (masslessKin) ? pT2 : sqrt((pT2 + s3) * (pT2 + s4));
  // For comparisons with Pythia 6.4 use scale similar to there.
  // Q2RenH = (masslessKin) ? pT2 : 0.5 * (s3 + s4) + pT2;
  alpS = alphaScalc.alphaS(Q2RenH);  
  alpEM = alphaEMcalc.alphaEM(Q2RenH);  

  // Use pT^2 as factorization scale, generalized to min( m_T3^2, m_T4^2)
  // for massive case. 
  // Is this correct for processes with t-channel W/Z exchange??
  Q2FacH = (masslessKin) ? pT2 : min (s3, s4) + pT2;
  // For comparisons with Pythia 6.4 use scale similar to there.
  // Q2FacH = (masslessKin) ? pT2 : 0.5 * (s3 + s4) + pT2;

  // Done.
  return true;

}

//*********

// As above, special kinematics for multiple interactions. 

bool Sigma2Process::set2KinMI( int id1in, int id2in, double x1in, double x2in,
  double sHin, double tHin, double uHin, double alpSin) {

  // Default ordering of particles 3 and 4.
  swapTU = false;
 
  // Incoming flavours and x values.
  id1 = id1in;
  id2 = id2in;
  x1 = x1in;
  x2 = x2in;

  // Standard Mandelstam variables and their squares.
  sH = sHin;
  tH = tHin;
  uH = uHin; 
  sH2 = sH * sH;
  tH2 = tH * tH;
  uH2 = uH * uH;

  // Strong coupling.
  alpS = alpSin;

  // Assume vanishing masses. (Will be modified in final kinematics.) 
  // Or are masses set elsewhere for MI??
  m3 = 0.;
  s3 = 0.;
  m4 = 0.;
  s4 = 0.; 

  //  Done.
  return true;

}

//*********

// Perform kinematics for a Multiple Interaction.

bool Sigma2Process::final2KinMI() {

  // Have to set flavours and colours.
  setIdColAcol();

  // Check that masses of outgoing particles not too big.
  m3 = ParticleDataTable::m0(idH[3]);
  m4 = ParticleDataTable::m0(idH[4]);
  double eCM = sqrt(sH);
  if (m3 + m4 + MASSMARGIN > eCM) return false;
  s3 = m3*m3;
  s4 = m4*m4;

  // Do kinematics of the decay.
  double eIn = 0.5 * eCM;
  double e3 = 0.5 * (sH + s3 - s4) / eCM;
  double e4 = 0.5 * (sH + s4 - s3) / eCM;
  double pAbs = sqrtpos( e3*e3 - s3 );
  double cosTheta = 1. + 2. * tH / sH;
  double sinTheta = 2. * sqrtpos( tH * uH ) / sH;
  thetaH = atan2( sinTheta, cosTheta);
  phiH = 2. * M_PI * Rndm::flat();
  double pZ = pAbs * cosTheta;
  double pX = pAbs * sinTheta * sin(phiH);
  double pY = pAbs * sinTheta * cos(phiH);
  double scale = eIn * sinTheta;

  // Fill particle info.
  parton[1] = Particle( idH[1], -31, 0, 0, 3, 4, colH[1], acolH[1],
    0., 0., eIn, eIn, 0., scale);
  parton[2] = Particle( idH[2], -31, 0, 0, 3, 4, colH[2], acolH[2],
    0., 0., -eIn, eIn, 0., scale);
  parton[3] = Particle( idH[3],  33, 1, 2, 0, 0, colH[3], acolH[3],
    pX, pY, pZ, e3, m3, scale);
  parton[4] = Particle( idH[4],  33, 1, 2, 0, 0, colH[4], acolH[4],
    -pX, -pY, -pZ, e4, m4, scale);

  // Boost particles from subprocess rest frame to event rest frame.
  double betaZ = (x1 - x2) / (x1 + x2);
  for (int i = 1; i <= 4; ++i) parton[i].bst(0., 0., betaZ);

  // Done.
  return true;

}  

//**************************************************************************

} // end namespace Pythia8
