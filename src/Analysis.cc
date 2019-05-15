// Function definitions (not found in the header) for the 
// Sphericity and CellJet classes.
// Copyright C 2006 Torbjorn Sjostrand

#include "Analysis.h"

namespace Pythia8 {

//**************************************************************************

// Sphericity class.
// This class finds sphericity-related properties of an event.

//*********
 
// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Minimum number of particles to perform study.
const int Sphericity::NSTUDYMIN = 2;

// Assign mimimum squared momentum in weight to avoid division by zero. 
const double Sphericity::P2MIN = 1e-20;

// Second eigenvalue not too low or not possible to find eigenvectors.
const double Sphericity::EIGENVALUEMIN = 1e-10;

//*********
 
// Analyze event.

bool Sphericity::analyze(Event& event) {

  // Initial values, tensor and counters zero.
  eVal1 = eVal2 = eVal3 = 0.;
  eVec1 = eVec2 = eVec3 = 0.;
  double tt[4][4];
  for (int j = 1; j < 4; ++j) 
  for (int k = j; k < 4; ++k) tt[j][k] = 0.;
  int nStudy = 0;
  double denom = 0.;

  // Loop over desired particles in the event.
  for (int i = 0; i < event.size(); ++i) 
  if (event[i].remains()) {
    if (select > 2 && event[i].isNeutral() ) continue;
    if (select == 2 && event[i].isInvisible() ) continue;
    ++nStudy;

    // Calculate matrix to be diagonalized. Special cases for speed.
    double pNow[4];
    pNow[1] = event[i].px();
    pNow[2] = event[i].py();
    pNow[3] = event[i].pz();
    double p2Now = pNow[1]*pNow[1] + pNow[2]*pNow[2] + pNow[3]*pNow[3];
    double pWeight = 1.;
    if (powerInt == 1) pWeight = 1. / sqrt(max(P2MIN, p2Now));
    else if (powerInt == 0) pWeight = pow( max(P2MIN, p2Now), powerMod);
    for (int j = 1; j < 4; ++j)   
    for (int k = j; k < 4; ++k) tt[j][k] += pWeight * pNow[j] * pNow[k];
    denom += pWeight * p2Now;
  }

  // Very low multiplicities (0 or 1) not considered.
  if (nStudy < NSTUDYMIN) {
    ErrorMessages::message("Warning in Sphericity::analyze: "
    " too few particles"); 
    return false;
  }

  // Normalize tensor to trace = 1.
  for (int j = 1; j < 4; ++j) 
  for (int k = j; k < 4; ++k) tt[j][k] /= denom;
 
  // Find eigenvalues to matrix (third degree equation).
  double qCoef = ( tt[1][1] * tt[2][2] + tt[1][1] * tt[3][3] 
    + tt[2][2] * tt[3][3] - pow2(tt[1][2]) - pow2(tt[1][3]) 
    - pow2(tt[2][3]) ) / 3. - 1./9.;
  double qCoefRt = sqrt( -qCoef);
  double rCoef = -0.5 * ( qCoef + 1./9. + tt[1][1] * pow2(tt[2][3]) 
    + tt[2][2] * pow2(tt[1][3]) + tt[3][3] * pow2(tt[1][2]) 
    - tt[1][1] * tt[2][2] * tt[3][3] ) 
    + tt[1][2] * tt[1][3] * tt[2][3] + 1./27.; 
  double pTemp = max( min( rCoef / pow3(qCoefRt), 1.), -1.);
  double pCoef = cos( acos(pTemp) / 3.);
  double pCoefRt = sqrt( 3. * (1. - pow2(pCoef)) );
  eVal1 = 1./3. + qCoefRt * max( 2. * pCoef, pCoefRt - pCoef);
  eVal3 = 1./3. + qCoefRt * min( 2. * pCoef, -pCoefRt - pCoef);
  eVal2 = 1. - eVal1 - eVal3;

  // Begin find first and last eigenvector.
  for (int iVal = 0; iVal < 2; ++iVal) {
    double eVal = (iVal == 0) ? eVal1 : eVal3;

    // If all particles are back-to-back then only first axis meaningful.
    if (iVal > 1 && eVal2 < EIGENVALUEMIN) {
      ErrorMessages::message("Warning in Sphericity::analyze: "
      " particles too back-to-back"); 
      return false;
    }

    // Set up matrix to diagonalize.
    double dd[4][4];
    for (int j = 1; j < 4; ++j) {
      dd[j][j] = tt[j][j] - eVal;
      for (int k = j + 1; k < 4; ++k) {
        dd[j][k] = tt[j][k]; 
        dd[k][j] = tt[j][k]; 
      }
    }

    // Find largest = pivotal element in matrix.
    int jMax = 0;
    int kMax = 0;
    double ddMax = 0.;
    for (int j = 1; j < 4; ++j) 
    for (int k = 1; k < 4; ++k) 
    if (abs(dd[j][k]) > ddMax) {
      jMax = j;
      kMax = k;
      ddMax = abs(dd[j][k]);
    }

    // Subtract one row from the other two; find new largest element. 
    int jMax2 = 0;
    ddMax = 0.;
    for (int j = 1; j < 4; ++j) 
    if ( j != jMax) {
      double pivot = dd[j][kMax] / dd[jMax][kMax];
      for (int k = 1; k < 4; ++k) {
        dd[j][k] -= pivot * dd[jMax][k];
        if (abs(dd[j][k]) > ddMax) {
          jMax2 = j;
          ddMax = abs(dd[j][k]);
	}
      } 
    }

    // Construct eigenvector. Normalize to unit length. Random sign.
    int k1 = kMax + 1; if (k1 > 3) k1 -= 3;
    int k2 = kMax + 2; if (k2 > 3) k2 -= 3;
    double eVec[4];
    eVec[k1] = -dd[jMax2][k2];    
    eVec[k2] = dd[jMax2][k1];    
    eVec[kMax] = (dd[jMax][k1] * dd[jMax2][k2]
      - dd[jMax][k2] * dd[jMax2][k1]) / dd[jMax][kMax];
    double length = sqrt( pow2(eVec[1]) + pow2(eVec[2])
      + pow2(eVec[3]) );
    if (Rndm::flat() > 0.5) length = -length;

    // Store eigenvectors.
    if (iVal == 0) eVec1 = Vec4( eVec[1] / length,
      eVec[2] / length, eVec[3] / length, 0.);
    else eVec3 = Vec4( eVec[1] / length,
      eVec[2] / length, eVec[3] / length, 0.);
  }

  // Middle eigenvector is orthogonal to the other two.
  eVec2 = cross3( eVec1, eVec3);
  if (Rndm::flat() > 0.5) eVec2 = -eVec2;

  // Done.
  return true;
}

//*********

// Provide a listing of the info.
  
void Sphericity::list(ostream& os) {

  // Header.
  os << "\n --------  Pythia Sphericity Listing, power = " 
     << fixed << setprecision(3) << setw(6) << power 
     << "  ------ \n \n  no     lambda      e_x       e_y       e_z \n";

  // The three eigenvalues and eigenvectors.
  os << setprecision(5);
  os << "   1" << setw(11) << eVal1 << setw(11) << eVec1.px() 
     << setw(10) << eVec1.py() << setw(10) << eVec1.pz() << "\n";
  os << "   2" << setw(11) << eVal2 << setw(11) << eVec2.px() 
     << setw(10) << eVec2.py() << setw(10) << eVec2.pz() << "\n";
  os << "   3" << setw(11) << eVal3 << setw(11) << eVec3.px() 
     << setw(10) << eVec3.py() << setw(10) << eVec3.pz() << "\n";

  // Listing finished.
  os << "\n --------  End Pythia Sphericity Listing  ------------------"
     << endl;
}

//**************************************************************************

// CellJet class.
// This class performs a cone jet search in (eta, phi, E_T) space.

//*********
 
// Analyze event.

bool CellJet::analyze(Event& event) {

  // Initial values zero.
  jets.resize(0);
  vector<SingleCell> cells;

  // Loop over desired particles in the event.
  for (int i = 0; i < event.size(); ++i) 
  if (event[i].remains()) {
    if (select > 2 && event[i].isNeutral() ) continue;
    if (select == 2 && event[i].isInvisible() ) continue;

    // Find particle position in (eta, phi, pT) space.
    double etaNow = event[i].eta();
    if (abs(etaNow) > etaMax) continue;
    double phiNow = event[i].phi();
    double pTnow = event[i].pT();
    int iEtaNow = max(1, min( nEta, 1 + int(nEta * 0.5 
      * (1. + etaNow / etaMax) ) ) );
    int iPhiNow = max(1, min( nPhi, 1 + int(nPhi * 0.5
      * (1. + phiNow / M_PI) ) ) );
    int iCell = nPhi * iEtaNow + iPhiNow;

    // Add pT to cell already hit or book a new cell.
    bool found = false;
    for (int j = 0; j < int(cells.size()); ++j) {
      if (iCell == cells[j].iCell) { 
        found = true;
        ++cells[j].multiplicity;
        cells[j].eTcell += pTnow; 
        continue;
      }
    }
    if (!found) {
      double etaCell = (etaMax / nEta) * (2 * iEtaNow - 1 - nEta);
      double phiCell = (M_PI / nPhi) * (2 * iPhiNow - 1 - nPhi);
      cells.push_back( SingleCell( iCell, etaCell, phiCell, pTnow, 1) );
    }
  }

  // Smear true bin content by calorimeter resolution.
  if (smear > 0) 
  for (int j = 0; j < int(cells.size()); ++j) {
    double eTeConv = (smear < 2) ? 1. : cosh( cells[j].etaCell );
    double eBef = cells[j].eTcell * eTeConv; 
    double eAft = 0.;
    do eAft = eBef + resolution * sqrt(eBef) * Rndm::gauss();
    while (eAft < 0 || eAft > upperCut * eBef);
    cells[j].eTcell = eAft / eTeConv;
  }

  // Remove cells below threshold for seed or for use at all.
  for (int j = 0; j < int(cells.size()); ++j) { 
    if (cells[j].eTcell < eTseed) cells[j].canBeSeed = false;
    if (cells[j].eTcell < threshold) cells[j].isUsed = true;
  }

  // Find seed cell: the one with highest pT of not yet probed ones.
  for ( ; ; ) {
    int jMax = 0;
    double eTmax = 0.;
    for (int j = 0; j < int(cells.size()); ++j) 
    if (cells[j].canBeSeed && cells[j].eTcell > eTmax) {
      jMax = j;
      eTmax = cells[j].eTcell;
    }

    // If too small cell eT then done, else start new trial jet.  
    if (eTmax < eTseed) break;
    double etaCenter = cells[jMax].etaCell;
    double phiCenter = cells[jMax].phiCell;
    double eTjet = 0.;

    //  Sum up unused cells within required distance of seed.
    for (int j = 0; j < int(cells.size()); ++j) {
      if (cells[j].isUsed) continue;
      double dEta = abs( cells[j].etaCell - etaCenter );
      if (dEta > coneRadius) continue;
      double dPhi = abs( cells[j].phiCell - phiCenter );
      if (dPhi > M_PI) dPhi = 2. * M_PI - dPhi;
      if (dPhi > coneRadius) continue;
      if (pow2(dEta) + pow2(dPhi) > pow2(coneRadius)) continue;
      cells[j].isAssigned = true;
      eTjet += cells[j].eTcell;
    }

    // Reject cluster below minimum ET.
    if (eTjet < eTjetMin) {
      cells[jMax].canBeSeed = false; 
      for (int j = 0; j < int(cells.size()); ++j) 
        cells[j].isAssigned = false;

    // Else find new jet properties. 
    } else {
      double etaWeighted = 0.;
      double phiWeighted = 0.;
      int multiplicity = 0;
      Vec4 pMassive;
      for (int j = 0; j < int(cells.size()); ++j) 
      if (cells[j].isAssigned) {
        cells[j].canBeSeed = false; 
        cells[j].isUsed = true; 
        cells[j].isAssigned = false; 
        etaWeighted += cells[j].eTcell * cells[j].etaCell;
        double phiCell = cells[j].phiCell; 
        if (abs(phiCell - phiCenter) > M_PI) 
          phiCell += (phiCenter > 0.) ? 2. * M_PI : -2. * M_PI;
        phiWeighted += cells[j].eTcell * phiCell;
        multiplicity += cells[j].multiplicity;
        pMassive += cells[j].eTcell * Vec4( cos(cells[j].phiCell), 
          sin(cells[j].phiCell), sinh(cells[j].etaCell), 
          cosh(cells[j].etaCell) );
      } 
      etaWeighted /= eTjet;
      phiWeighted /= eTjet; 

      // Bookkeep new jet, in decreasing ET order.
      jets.push_back( SingleCellJet( eTjet, etaCenter, phiCenter,
        etaWeighted, phiWeighted, multiplicity, pMassive) ); 
      for (int i = int(jets.size()) - 1; i > 0; --i) {
        if (jets[i-1].eTjet > jets[i].eTjet) break;
        swap( jets[i-1], jets[i]);
      }
    }
  }

  // Done.
  return true;
}

//*********

// Provide a listing of the info.
  
void CellJet::list(ostream& os) {

  // Header.
  os << "\n --------  Pythia CellJet Listing, eTjetMin = " 
     << fixed << setprecision(3) << setw(8) << eTjetMin 
     << ", coneRadius = " << setw(5) << coneRadius 
     << "  ------------------------------ \n \n  no    "
     << " eTjet  etaCtr  phiCtr   etaWt   phiWt mult      p_x"
     << "        p_y        p_z         e          m \n";

  // The jets.
  for (int i = 0; i < int(jets.size()); ++i) {
    os << setw(4) << i << setw(10) << jets[i].eTjet << setw(8) 
       << jets[i].etaCenter << setw(8) << jets[i].phiCenter << setw(8) 
       << jets[i].etaWeighted << setw(8) << jets[i].phiWeighted 
       << setw(5) << jets[i].multiplicity << setw(11) 
       << jets[i].pMassive.px() << setw(11) << jets[i].pMassive.py()
       << setw(11) << jets[i].pMassive.pz() << setw(11) 
       << jets[i].pMassive.e() << setw(11)
       << jets[i].pMassive.mCalc() << "\n";  
  }

  // Listing finished.
  os << "\n --------  End Pythia CellJet Listing  ------------------"
     << "-------------------------------------------------"
     << endl;
}

//**************************************************************************

} // end namespace Pythia8
