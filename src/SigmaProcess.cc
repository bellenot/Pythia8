// Function definitions (not found in the header) for the 
// SigmaProcess class, and classes derived from it.
// Copyright C 2007 Torbjorn Sjostrand

#include "SigmaProcess.h"

namespace Pythia8 {

//**************************************************************************

// The SigmaProcess class.
// Base class for cross sections.

//*********
 
// Definitions of static variables and functions.
// (Values will be overwritten in initStatic call, so are purely dummy.)

int         SigmaProcess::alphaSorder   = 1;
int         SigmaProcess::alphaEMorder  = 1;
int         SigmaProcess::nQuark        = 3;
int         SigmaProcess::renormScale   = 0;
int         SigmaProcess::factorScale   = 0;
int         SigmaProcess::SMHiggsParity = 1;
double      SigmaProcess::alphaSvalue   = 0.1265;
double      SigmaProcess::renormMult    = 1.;
double      SigmaProcess::factorMult    = 1.;
double      SigmaProcess::SMHiggsEta    = 0.;
AlphaStrong SigmaProcess::alphaS;
AlphaEM     SigmaProcess::alphaEM;

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Conversion of GeV^{-2} to mb for cross section.
const double SigmaProcess::CONVERT2MB   = 0.389380; 

// The sum of outgoing masses must not be too close to the cm energy.
const double SigmaProcess::MASSMARGIN   = 0.1;

// Information on incoming beams.
int    SigmaProcess::idA, SigmaProcess::idB;
double SigmaProcess::mA,  SigmaProcess::mB; 
bool   SigmaProcess::hasLeptonBeams     = false;
  
// Pointer to the total/elastic/diffractive cross section object.
SigmaTotal* SigmaProcess::sigmaTotPtr;

// Pointer to the SLHA object.
SusyLesHouches* SigmaProcess::slha;

//*********

// Initialize static data members.

void SigmaProcess::initStatic() {

  // Parameters of alphaStrong generation .
  alphaSvalue   = Settings::parm("SigmaProcess:alphaSvalue");
  alphaSorder   = Settings::mode("SigmaProcess:alphaSorder");

  // Initialize alphaStrong generation.
  alphaS.init( alphaSvalue, alphaSorder); 

  // Parameters of alphaEM generation.
  alphaEMorder  = Settings::mode("SigmaProcess:alphaEMorder");

  // Initialize alphaEM generation.
  alphaEM.init( alphaEMorder); 

  // Maximum new quark flavour.
  nQuark        = Settings::mode("SigmaProcess:nQuark");

  // Renormalization scale choice.
  renormScale   = Settings::mode("SigmaProcess:renormScale"); 
  renormMult    = Settings::parm("SigmaProcess:renormMult"); 

  // Factorization scale choice.
  factorScale   = Settings::mode("SigmaProcess:factorScale"); 
  factorMult    = Settings::parm("SigmaProcess:factorMult"); 

  // Higgs parity assumption.
  SMHiggsParity = Settings::mode("SMHiggs:parity");
  SMHiggsEta    = Settings::parm("SMHiggs:etaParity");

}

//*********

// Evaluate weight for W decay distribution in t -> W b -> f fbar b.

double SigmaProcess::weightTopDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // If not pair W d/s/b and mother t then return unit weight.
  if (iResEnd - iResBeg != 2) return 1.;
  int iW = iResBeg;
  int iB = iResBeg + 1;
  int idW = process[iW].idAbs();
  int idB = process[iB].idAbs();
  if (idW != 24) {
    swap(iW, iB); 
    swap(idW, idB);
  } 
  if (idW != 24 || (idB != 1 && idB != 3 && idB != 5)) return 1.;
  int iT = process[iW].mother1(); 
  if (iT <= 0 || process[iT].idAbs() != 6) return 1.;

  // Find sign-matched order of W decay products. 
  int iF    = process[iW].daughter1(); 
  int iFbar = process[iW].daughter2();
  if (iFbar - iF != 1) return 1.; 
  if (process[iT].id() * process[iF].id() < 0) swap(iF, iFbar);

  // Weight and maximum weight.
  double wt    = (process[iT].p() * process[iFbar].p()) 
               * (process[iF].p() * process[iB].p());
  double wtMax = ( pow4(process[iT].m()) - pow4(process[iW].m()) ) / 8.;  

  // Done.
  return wt / wtMax;

}

//*********

// Evaluate weight for Z0/W+- decay distributions in 
// H -> Z0 Z0 or W+ W- -> f fbar f' fbar'.

double SigmaProcess::weightHiggsDecay( Event& process, int iResBeg, 
  int iResEnd) {

  // If not pair Z0 Z0 or W+ W- or not mother Higgs then return unit weight.
  if (iResEnd - iResBeg != 2) return 1.;
  int iZW1  = iResBeg;
  int iZW2  = iResBeg + 1;
  int idZW1 = process[iZW1].id();
  int idZW2 = process[iZW2].id();
  if (idZW1 < 0) {
    swap(iZW1, iZW2); 
    swap(idZW1, idZW2);
  } 
  if ( (idZW1 != 23 || idZW2 != 23) && (idZW1 != 24 || idZW2 != -24) )
    return 1.;
  int iH = process[iZW1].mother1(); 
  if (iH <= 0 || process[iH].id() != 25) return 1.;

  // Option with isotropic decays.
  if (SMHiggsParity == 0) return 1.;

  // Maximum and initial weight. 
  double wtMax = pow4(process[iH].m());
  double wt    = wtMax; 

  // Find sign-matched order of Z0/W+- decay products. 
  int i3 = process[iZW1].daughter1();
  int i4 = process[iZW1].daughter2();
  if (process[i3].id() < 0) swap( i3, i4); 
  int i5 = process[iZW2].daughter1();
  int i6 = process[iZW2].daughter2();
  if (process[i5].id() < 0) swap( i5, i6); 

  // Evaluate four-vector products and find masses..
  double p35  = 2. * process[i3].p() * process[i5].p(); 
  double p36  = 2. * process[i3].p() * process[i6].p(); 
  double p45  = 2. * process[i4].p() * process[i5].p(); 
  double p46  = 2. * process[i4].p() * process[i6].p(); 
  double p34  = 2. * process[i3].p() * process[i4].p(); 
  double p56  = 2. * process[i5].p() * process[i6].p(); 
  double mZW1 = process[iZW1].m();
  double mZW2 = process[iZW2].m();

  // For mixed CP states need epsilon product and gauge boson masses.
  double epsilonProd = 0.;
  if (SMHiggsParity == 3) {
    double p[4][4];
    for (int i = 0; i < 4; ++i) {
      int         ii = i3;
      if (i == 1) ii = i4;
      if (i == 2) ii = i5;
      if (i == 3) ii = i6;
      p[i][0] = process[ii].e();
      p[i][1] = process[ii].px();
      p[i][2] = process[ii].py();
      p[i][3] = process[ii].pz();
    }     
    epsilonProd 
      = p[0][0]*p[1][1]*p[2][2]*p[3][3] - p[0][0]*p[1][1]*p[2][3]*p[3][2] 
      - p[0][0]*p[1][2]*p[2][1]*p[3][3] + p[0][0]*p[1][2]*p[2][3]*p[3][1]
      + p[0][0]*p[1][3]*p[2][1]*p[3][2] - p[0][0]*p[1][3]*p[2][2]*p[3][1]
      - p[0][1]*p[1][0]*p[2][2]*p[3][3] + p[0][1]*p[1][0]*p[2][3]*p[3][2]
      + p[0][1]*p[1][2]*p[2][0]*p[3][3] - p[0][1]*p[1][2]*p[2][3]*p[3][0]
      - p[0][1]*p[1][3]*p[2][0]*p[3][2] + p[0][1]*p[1][3]*p[2][2]*p[3][0]
      + p[0][2]*p[1][0]*p[2][1]*p[3][3] - p[0][2]*p[1][0]*p[2][3]*p[3][1]
      - p[0][2]*p[1][1]*p[2][0]*p[3][3] + p[0][2]*p[1][1]*p[2][3]*p[3][0] 
      + p[0][2]*p[1][3]*p[2][0]*p[3][1] - p[0][2]*p[1][3]*p[2][1]*p[3][0]
      - p[0][3]*p[1][0]*p[2][1]*p[3][2] + p[0][3]*p[1][0]*p[2][2]*p[3][1] 
      + p[0][3]*p[1][1]*p[2][0]*p[3][2] - p[0][3]*p[1][1]*p[2][2]*p[3][0] 
      - p[0][3]*p[1][2]*p[2][0]*p[3][1] + p[0][3]*p[1][2]*p[2][1]*p[3][0];
  }

  // Z0 Z0 decay: vector and axial couplings of two fermion pairs.
  if (idZW1 == 23) {
    double vf1 = CoupEW::vf(process[i3].idAbs());
    double af1 = CoupEW::af(process[i3].idAbs());
    double vf2 = CoupEW::vf(process[i5].idAbs());
    double af2 = CoupEW::af(process[i5].idAbs());
    double va12asym = 4. * vf1 * af1 * vf2 * af2 
      / ( (vf1*vf1 + af1*af1) * (vf2*vf2 + af2*af2) );
    double etaMod = SMHiggsEta / pow2( ParticleDataTable::m0(23) );
    
    // Normal CP-even decay.
    if (SMHiggsParity == 1) wt = 8. * (1. + va12asym) * p35 * p46 
      + 8. * (1. - va12asym) * p36 * p45;

    // CP-odd decay.
    else if (SMHiggsParity == 2) wt = ( pow2(p35 + p46) 
      + pow2(p36 + p45) - 2. * p34 * p56
      - 2. * pow2(p35 * p46 - p36 * p45) / (p34 * p56) 
      + va12asym * (p35 + p36 - p45 - p46) * (p35 + p45 - p36 - p46) )
      / (1. +  va12asym);

    // Mixed CP states. 
    else wt = 32. * ( 0.25 * ( (1. + va12asym) * p35 * p46 
      + (1. - va12asym) * p36 * p45 ) - 0.5 * etaMod * epsilonProd
      * ( (1. + va12asym) * (p35 + p46) - (1. - va12asym) * (p36 + p45) )
      + 0.0625 * etaMod * etaMod * (-2. * pow2(p34 * p56) 
      - 2. * pow2(p35 * p46 - p36 * p45) 
      + p34 * p56 * (pow2(p35 + p46) + pow2(p36 + p45)) 
      + va12asym * p34 * p56 * (p35 + p36 - p45 - p46) 
      * (p35 + p45 - p36 - p46) ) ) / ( 1. * 2. * etaMod * mZW1 * mZW2 
      + 2. * pow2(etaMod * mZW1 * mZW2) * (1. + va12asym) );

  // W+ W- decay.
  } else if (idZW1 == 24) {
    double etaMod = SMHiggsEta / pow2( ParticleDataTable::m0(24) );
    
    // Normal CP-even decay.
    if (SMHiggsParity == 1) wt = 16. * p35 * p46; 

    // CP-odd decay.
    else if (SMHiggsParity == 2) wt = 0.5 * ( pow2(p35 + p46) 
      + pow2(p36 + p45) - 2. * p34 * p56  
      - 2. * pow2(p35 * p46 - p36 * p45) / (p34 * p56) 
      + (p35 + p36 - p45 - p46) * (p35 + p45 - p36 - p46) );

    // Mixed CP states. 
    else wt = 32. * ( 0.25 * 2. * p35 * p46 
      - 0.5 * etaMod * epsilonProd * 2. * (p35 + p46)
      + 0.0625 * etaMod * etaMod * (-2. * pow2(p34 * p56) 
      - 2. * pow2(p35 * p46 - p36 * p45) 
      + p34 * p56 * (pow2(p35 + p46) + pow2(p36 + p45)) 
      + p34 * p56 * (p35 + p36 - p45 - p46) * (p35 + p45 - p36 - p46) ) ) 
      / ( 1. * 2. * etaMod * mZW1 * mZW2 + 2. * pow2(etaMod * mZW1 * mZW2) );
  }

  // Done.
  return wt / wtMax;

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
  x1  = x1in;
  x2  = x2in;
  sH  = sHin;
  sH2 = sH * sH;

  // Use sHat as renormalization scale. Evaluate alpha_strong and alpha_EM.
  Q2RenH = renormMult * sH;
  alpS   = alphaS.alphaS(Q2RenH);  
  alpEM  = alphaEM.alphaEM(Q2RenH);  

  // Use sHat as factorization scale.
  Q2FacH = factorMult * sH;

  // Done.
  return true;

}

//**************************************************************************

// The Sigma2Process class.
// Base class for resolved 2 -> 2 cross sections; derived from SigmaProcess.

//*********

// Input and complement kinematics for resolved 2 -> 2 process. 

bool Sigma2Process::set2Kin( double x1in, double x2in, double sHin, 
  double tHin, double m3in, double m4in, double runBW3in, double runBW4in) {

  // Incoming flavours not known.
  id12IsSet = false;

  // Default ordering of particles 3 and 4.
  swapTU = false;

  // Incoming parton momentum fractions.
  x1     = x1in;
  x2     = x2in;

  // Incoming masses and their squares.
  bool masslessKin = (id3Mass() == 0) && (id4Mass() == 0);
  if (masslessKin) {
    m3   = 0.;
    m4   = 0.;
  } else {
    m3   = m3in;
    m4   = m4in;
  }
  mH[3]  = m3;
  mH[4]  = m4;
  s3     = m3 * m3;
  s4     = m4 * m4;

  // Standard Mandelstam variables and their squares.
  sH     = sHin;
  tH     = tHin;
  uH     = (masslessKin) ? -(sH + tH) : s3 + s4 - (sH + tH); 
  sH2    = sH * sH;
  tH2    = tH * tH;
  uH2    = uH * uH;

  // The nominal Breit-Wigner factors with running width.
  runBW3 = runBW3in;
  runBW4 = runBW4in; 

  // Calculate squared transverse momentum.
  pT2 = (masslessKin) ?  tH * uH / sH : (tH * uH - s3 * s4) / sH;

  // Use pT^2 as renormalization scale, generalized to m_T3 * m_T4
  // for massive case, alternatively PYTHIA6 scale.
  // Is this correct for processes with s/t-channel W/Z exchange??
  if (renormScale == 1) Q2RenH = (masslessKin) ? pT2 : 0.5 * (s3 + s4) + pT2;
  else Q2RenH = (masslessKin) ? pT2 : sqrt((pT2 + s3) * (pT2 + s4));
  Q2RenH *= renormMult;

  // Debug??
  Q2RenH = sH;

  // Evaluate alpha_strong and alpha_EM.
  alpS = alphaS.alphaS(Q2RenH);  
  alpEM = alphaEM.alphaEM(Q2RenH);  

  // Use pT^2 as factorization scale, generalized to min( m_T3^2, m_T4^2)
  // for massive case, alternatively PYTHIA6 scale. 
  // Is this correct for processes with s/t-channel W/Z exchange??
  if (factorScale == 1) Q2FacH = (masslessKin) ? pT2 : 0.5 * (s3 + s4) + pT2;
  else Q2FacH = (masslessKin) ? pT2 : min (s3, s4) + pT2;
  Q2FacH *= factorMult;

  // Debug??
  Q2FacH = sH;

  // Done.
  return true;

}

//*********

// As above, special kinematics for multiple interactions. 

bool Sigma2Process::set2KinMI( int id1in, int id2in, double x1in, double x2in,
  double sHin, double tHin, double uHin, double alpSin, double alpEMin,
  bool needMasses, double m3in, double m4in) {

  // Default ordering of particles 3 and 4.
  swapTU = false;
 
  // Incoming flavours and x values.
  id1       = id1in;
  id2       = id2in;
  id12IsSet = true;
  x1        = x1in;
  x2        = x2in;

  // Standard Mandelstam variables and their squares.
  sH        = sHin;
  tH        = tHin;
  uH        = uHin; 
  sH2       = sH * sH;
  tH2       = tH * tH;
  uH2       = uH * uH;

  // Strong and electroweak couplings.
  alpS      = alpSin;
  alpEM     = alpEMin;

  // Assume vanishing masses. (Will be modified in final kinematics.) 
  m3        = 0.;
  s3        = 0.;
  m4        = 0.;
  s4        = 0.;
  sHBeta    = sH; 

  // Scattering angle.
  cosTheta  = (tH - uH) / sH;
  sinTheta  = 2. * sqrtpos( tH * uH ) / sH;

  // In some cases must use masses and redefine meaning of tHat and uHat.
  if (needMasses) { 
    m3      = m3in;
    s3      = m3 * m3;
    m4      = m4in;
    s4      = m4 * m4;
    sHMass  = sH - s3 - s4;
    sHBeta  = sqrtpos(sHMass*sHMass - 4. * s3 * s4);   
    tH      = -0.5 * (sHMass - sHBeta * cosTheta); 
    uH      = -0.5 * (sHMass + sHBeta * cosTheta); 
    tH2     = tH * tH;
    uH2     = uH * uH;
  }

  // pT2 with masses (at this stage) included.
  pT2Mass   = 0.25 * sHBeta * pow2(sinTheta);

  //  Done.
  return true;

}

//*********

// Perform kinematics for a Multiple Interaction.

bool Sigma2Process::final2KinMI() {

  // Have to set flavours and colours.
  setIdColAcol();

  // Check that masses of outgoing particles not too big.
  m3           = ParticleDataTable::m0(idH[3]);
  m4           = ParticleDataTable::m0(idH[4]);
  double eCM   = sqrt(sH);
  if (m3 + m4 + MASSMARGIN > eCM) return false;
  s3           = m3 * m3;
  s4           = m4 * m4;

  // Do kinematics of the decay.
  double eIn   = 0.5 * eCM;
  double e3    = 0.5 * (sH + s3 - s4) / eCM;
  double e4    = 0.5 * (sH + s4 - s3) / eCM;
  double pAbs  = sqrtpos( e3*e3 - s3 );
  phi          = 2. * M_PI * Rndm::flat();
  double pZ    = pAbs * cosTheta;
  double pX    = pAbs * sinTheta * sin(phi);
  double pY    = pAbs * sinTheta * cos(phi);
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
