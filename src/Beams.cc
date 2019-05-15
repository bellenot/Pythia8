// Function definitions (not found in the header) for the 
// BeamParticle and BeamRemnants classes.
// Copyright C 2006 Torbjorn Sjostrand

#include "Beams.h"

namespace Pythia8 {
 
//**************************************************************************

// The BeamParticle class.

//*********
 
// Definitions of static variables.
// (Values will be overwritten in initStatic call, so are purely dummy.)

int BeamParticle::maxValQuark = 3;
int BeamParticle::companionPower = 4;
bool BeamParticle::allowJunction = true;

//*********

// Initialize static data members.

void BeamParticle::initStatic() {

  // Maximum quark kind in allowed incoming beam hadrons.
  maxValQuark = Settings::mode("Beams:maxValQuark");

  // Assume g(x) ~ (1-x)^power/x to constrain companion to sea quark.
  companionPower = Settings::mode("Beams:companionPower");

  // Allow or not more than two valence quarks to be kicked out.
  allowJunction = Settings::flag("Beams:allowJunction");

}

//*********

// Initialize data on a beam particle.

void BeamParticle::init( int idIn, double pzIn, double eIn, double mIn, 
  PDF* pdfInPtr) {

  // Store info on the incoming beam.
  idBeam = idIn; initBeamKind(); 
  pBeam = Vec4( 0., 0., pzIn, eIn); 
  mBeam = mIn; 
  pdfBeamPtr = pdfInPtr; 

}

//*********

// Initialize kind and valence flavour content of incoming beam.
// For recognized hadrons one can generate multiple interactions.
// So far we do not handle diagonal mesons or K0S/K0L (or photons), 
// for which flavour content is only known after first interaction.

void BeamParticle::initBeamKind() {

  // Reset.
  idBeamAbs = abs(idBeam);
  isLeptonBeam = false;
  isHadronBeam = false;  
  isMesonBeam = false;  
  isBaryonBeam = false;  
  nValKinds = 0;

  // Check for leptons. Done if cannot be lowest-lying hadron state.
  if (idBeamAbs == 11 || idBeamAbs == 13 || idBeamAbs == 15) 
    isLeptonBeam = true; 
  if (idBeamAbs < 101 || idBeamAbs > 9999) return;
  
  // Resolve valence content for assumed meson. Flunk unallowed codes.
  if (idBeamAbs < 1000) {
    int id1 = idBeamAbs/100;    
    int id2 = (idBeamAbs/10)%10;
    if ( id1 < 1 || id1 > maxValQuark 
      || id2 < 1 || id2 > maxValQuark ) return;
    if (id2 == id1 || idBeamAbs == 130 || idBeamAbs == 310) return;
    isMesonBeam = true;
    
    // Store valence content of a confirmed meson.
    nValKinds = 2; nVal[0] = 1 ; nVal[1] = 1;
    if (id1%2 == 0) {idVal[0] = id1; idVal[1] = -id2;}
    else {idVal[0] = id2; idVal[1] = -id1;}      
  
  // Resolve valence content for assumed baryon. Flunk unallowed codes.
  } else { 
    int id1 = idBeamAbs/1000;
    int id2 = (idBeamAbs/100)%10;
    int id3 = (idBeamAbs/10)%10;
    if ( id1 < 1 || id1 > maxValQuark || id2 < 1 || id2 > maxValQuark 
      || id3 < 1 || id3 > maxValQuark) return;
    if (id2 > id1 || id3 > id1) return;
    isBaryonBeam = true;

    // Store valence content of a confirmed baryon.
    nValKinds = 1; idVal[0] = id1; nVal[0] = 1;
    if (id2 == id1) ++nVal[0];
    else {nValKinds = 2; idVal[1] = id2; nVal[1] = 1;}
    if (id3 == id1) ++nVal[0];
    else if (id3 == id2) ++nVal[1];
    else {idVal[nValKinds] = id3; nVal[nValKinds] = 1; ++nValKinds;}
  }
  
  // Flip flavours for antimeson or antibaryon, and then done.
  if (idBeam < 0) for (int i = 0; i < nValKinds; ++i) idVal[i] = -idVal[i];
  isHadronBeam = true;
  Q2ValFracSav = -1.;

}

//*********

double BeamParticle::xMax(int iSkip) {

  // Minimum requirement on remaining energy > nominal mass for hadron.
  double xLeft = 1.;
  if (isHadron()) xLeft -= m() / e();
  if (size() == 0) return xLeft;

  // Subtract what was carried away by initiators (to date).
  for (int i = 0; i < size(); ++i) 
    if (i != iSkip) xLeft -= resolved[i].x();
  return xLeft;

}

//*********

// Parton distributions, reshaped to take into account previous 
// multiple interactions. By picking a non-negative iSkip value,
// one particular interaction is skipped, as needed for ISR  

double BeamParticle::xfModified(int iSkip, int id, double x, double Q2) {

  // Initial values.
  idSave = id;
  iSkipSave = iSkip;
  xqVal = 0.;
  xqgSea = 0.;
  xqCompSum = 0.;

  // Fast procedure for first interaction. 
  if (size() == 0) {
    if (x >= 1.) return 0.;
    bool canBeVal = false;
    for (int i = 0; i < nValKinds; ++i) if (id == idVal[i]) canBeVal = true;
    if (canBeVal) { xqVal = xfVal( id, x, Q2); xqgSea = xfSea( id, x, Q2); }
    else xqgSea = xf( id, x, Q2);

  // More complicated procedure for non-first interaction.
  } else { 

    // Sum up the x already removed, and check that remaining x is enough.
    double xUsed = 0.;
    for (int i = 0; i < size(); ++i) 
      if (i != iSkip) xUsed += resolved[i].x();
    double xLeft = 1. - xUsed;
    if (x >= xLeft) return 0.;
    double xRescaled = x / xLeft;

    // Calculate total and remaining amount of x carried by valence quarks.
    double xValTot = 0.;
    double xValLeft = 0.;
    for (int i = 0; i < nValKinds; ++i) {
      nValLeft[i] = nVal[i];
      for (int j = 0; j < size(); ++j)  
      if (j != iSkip && resolved[j].isValence() 
        && resolved[j].id() == idVal[i]) --nValLeft[i];
      double xValNow =  xValFrac(i, Q2);
      xValTot += nVal[i] * xValNow;
      xValLeft += nValLeft[i] * xValNow;
    }

    // Calculate total amount of x carried by unmatched companion quarks.
    double xCompAdded = 0.;
    for (int i = 0; i < size(); ++i)  
    if (i != iSkip && resolved[i].isUnmatched()) xCompAdded 
      += xCompFrac( resolved[i].x() / (xLeft + resolved[i].x()) )
      // Typo warning: extrafactor missing in Skands&Sjostrand article;
      // <x> for companion refers to fraction of x left INCLUDING sea quark.
      // To be modified further??
      * (1. + resolved[i].x() / xLeft);
  
    // Calculate total rescaling factor and pdf for sea and gluon.
    double rescaleGS = max( 0., (1. - xValLeft - xCompAdded) 
      / (1. - xValTot) );
    xqgSea = rescaleGS * xfSea( id, xRescaled, Q2); 

    // Find valence part and rescale it to remaining number of quarks. 
    for (int i = 0; i < nValKinds; ++i) 
    if (id == idVal[i] && nValLeft[i] > 0) 
      xqVal = xfVal( id, xRescaled, Q2) 
      * double(nValLeft[i]) / double(nVal[i]); 
                                                                               
    // Find companion part, by adding all companion contributions.
    for (int i = 0; i < size(); ++i) 
    if (i != iSkip && resolved[i].id() == -id && resolved[i].isUnmatched()) {
      double xsRescaled = resolved[i].x() / (xLeft + resolved[i].x());
      double xcRescaled = x / (xLeft + resolved[i].x());  
      double xqCompNow = xCompDist( xcRescaled, xsRescaled); 
      resolved[i].xqCompanion( xqCompNow);
      xqCompSum += xqCompNow; 
    }
  }

  // Add total, but only return relevant part for ISR. More cases??
  // Watch out, e.g. g can come from either kind of quark.??
  xqgTot = xqVal + xqgSea + xqCompSum;
  if (iSkip >= 0) {
    if (resolved[iSkip].isValence()) return xqVal;
    if (resolved[iSkip].isUnmatched()) return xqgSea + xqCompSum; 
  }
  return xqgTot;   
  
}

//*********

// Decide whether a quark extracted from the beam is of valence, sea or
// companion kind; in the latter case also pick its companion.
// Assumes xfModified has already been called.

  void BeamParticle::pickValSeaComp() {

  // If parton already has a companion than reset code for this.
  int oldCompanion = resolved[iSkipSave].companion();
  if (oldCompanion >= 0) resolved[oldCompanion].companion(-2); 

  // Default assignment is sea.
  int vsc = -2;

  // For gluons or photons no sense of valence or sea.
  if (idSave == 21 || idSave == 22) vsc = -1;  

  // Decide if valence or sea quark.
  else {
    double xqRndm = xqgTot * Rndm::flat(); 
    if (xqRndm < xqVal) vsc = -3; 
    else if (xqRndm < xqVal + xqgSea) vsc = -2; 
 
    // If not either, loop over all possible companion quarks.
    else {
      xqRndm -= xqVal + xqgSea;
      for (int i = 0; i < size(); ++i) 
      if (i != iSkipSave && resolved[i].id() == -idSave 
        && resolved[i].isUnmatched()) {
        xqRndm -= resolved[i].xqCompanion();
        if (xqRndm < 0.) vsc = i; 
        break;
      }
    }
  }

  // Bookkeep assignment; for sea--companion pair both ways.  
  resolved[iSkipSave].companion(vsc);
  if (vsc >= 0) resolved[vsc].companion(iSkipSave);

} 

//*********

// Fraction of hadron momentum sitting in a valence quark distribution.
// Based on hardcoded parametrizations of CTEQ 5L numbers.

double BeamParticle::xValFrac(int j, double Q2) {

  // Only recalculate when required.
  if (Q2 != Q2ValFracSav) { 
    Q2ValFracSav = Q2;
     
    // Q2-dependence of log-log form; assume fixed Lambda = 0.2.
    double llQ2 = log( log( max( 1., Q2) / 0.04 ));

    // Fractions carried by u and d in proton.
    uValInt =  0.48 / (1. + 1.56 * llQ2);
    dValInt = 0.385 / (1. + 1.60 * llQ2);
  }

  // Baryon with three different quark kinds: (2 * u + d) / 3 of proton.  
  if (isBaryonBeam && nValKinds == 3) return (2. * uValInt + dValInt) / 3.;

  // Baryon with one or two identical: like d or u of proton.
  if (isBaryonBeam && nVal[j] == 1) return dValInt;
  if (isBaryonBeam && nVal[j] == 2) return uValInt;

  // Meson: (2 * u + d) / 2 of proton so same total valence quark fraction.
    return 0.5 * (2. * uValInt + dValInt);

}

//*********

// The momentum integral of a companion quark, with its partner at x_s, 
// using an approximate gluon density like (1 - x_g)^power / x_g.
// The value corresponds to an unrescaled range between 0 and 1 - x_s.

double BeamParticle::xCompFrac(double xs) {

  // Select case by power of gluon (1-x_g) shape.
  switch (companionPower) {

    case 0: 
       return xs * ( 5. + xs * (-9. - 2. * xs * (-3. + xs)) + 3. * log(xs) )
         / ( (-1. + xs) * (2. + xs * (-1. + 2. * xs)) );

    case 1:
       return -1. -3. * xs + ( 2. * pow2(-1. + xs) * (1. + xs + xs*xs))
         / ( 2. + xs*xs * (xs - 3.) + 3. * xs * log(xs) );

    case 2:
       return xs * ( (1. - xs) * (19. + xs * (43. + 4. * xs))
         + 6. * log(xs) * (1. + 6. * xs + 4.*xs*xs) ) /
        ( 4. * ( (xs - 1.) * (1. + xs * (4. + xs) )
        - 3. * xs * log(xs) * (1 + xs) ) );

    case 3:
      return 3. * xs * ( (xs - 1.) * (7. + xs * (28. + 13. * xs))
        - 2. * log(xs) * (1. + xs * (9. + 2. * xs * (6. + xs))) ) 
        / ( 4. + 27. * xs - 31. * pow3(xs) 
        + 6. * xs * log(xs) * (3. + 2. * xs * (3.+xs)) );

    default:
      return ( -9. * xs * (xs*xs - 1.) * (5. + xs * (24. + xs)) + 12. * xs 
        * log(xs) * (1. + 2. * xs) * (1. + 2. * xs * (5. + 2. * xs)) ) 
        / ( 8. * (1. + 2. * xs) * ((xs - 1.) * (1. + xs * (10. + xs))
        - 6. * xs * log(xs) * (1. + xs)) );

  }
}

//*********

// The x*f pdf of a companion quark at x_c, with its sea partner at x_s,
// using an approximate gluon density like (1 - x_g)^power / x_g. 
// The value corresponds to an unrescaled range between 0 and 1 - x_s.

double BeamParticle::xCompDist(double xc, double xs) {

  // Mother gluon momentum fraction. Check physical limit.
  double xg = xc + xs;
  if (xg > 1.) return 0.;

  // Common factor, including splitting kernel and part of gluon density
  // (and that it is x_c * f that is coded).
  double fac = 3. * xc * xs * (xc*xc + xs*xs) / pow4(xg);

  // Select case by power of gluon (1-x_g) shape.
  switch (companionPower) {

    case 0: 
      return fac / ( 2. - xs * (3. - xs * (3. - 2. * xs)) );

    case 1:
      return fac * (1. - xg) / ( 2. + xs*xs * (-3. + xs) + 3. * xs * log(xs) );

    case 2:
      return fac * pow2(1. - xg) / ( 2. * ((1. - xs) * (1. + xs * (4. + xs))
        + 3. * xs * (1. + xs) * log(xs)) );

    case 3:
      return fac * pow3(1. - xg) * 2. / ( 4. + 27. * xs - 31. * pow3(xs)
        + 6. * xs * log(xs) * (3. + 2. * xs * (3. + xs)) );

    default:
       return fac * pow4(1. - xg) / ( 2. * (1. + 2. * xs) * ((1. - xs) 
         * (1. + xs * (10. + xs)) + 6. * xs * log(xs) * (1. + xs)) );

  }
}

//*********

// Add required extra remnant flavour content. 

bool BeamParticle::remnantFlavours() {

  // A baryon will have a junction, unless a diquark is formed later.
  hasJunctionBeam = (isBaryon()) ? true : false;  

  // Store how many hard-scattering partons were removed from beam.
  nInit = size();

  // Find remaining valence quarks.
  for (int i = 0; i < nValKinds; ++i) {
    nValLeft[i] = nVal[i];
    for (int j = 0; j < nInit; ++j) if (resolved[j].isValence() 
      && resolved[j].id() == idVal[i]) --nValLeft[i];
    // Add remaining valence quarks to record. Partly temporary values.
    for (int k = 0; k < nValLeft[i]; ++k) append(0, idVal[i], 0., -3);
  }

  // If at least two valence quarks left in baryon remnant then form diquark.
  int nInitPlusVal = size(); 
  if (isBaryon() && nInitPlusVal - nInit >= 2) {

    // If three, pick two at random to form diquark, else trivial.  
    int iQ1 = nInit;
    int iQ2 = nInit + 1;
    if (nInitPlusVal - nInit == 3) {
      double pickDq = 3. * Rndm::flat();
      if (pickDq > 1.) iQ2 = nInit + 2;
      if (pickDq > 2.) iQ1 = nInit + 1;
    } 

    // Pick spin 0 or 1 according to SU(6) wave function factors.
    int idDq = StringFlav::makeDiquark( resolved[iQ1].id(),
      resolved[iQ2].id(), idBeam);

    // Overwrite with diquark flavour and remove one slot. No more junction.
    resolved[iQ1].id(idDq);
    if (nInitPlusVal - nInit == 3 && iQ2 == nInit + 1) 
      resolved[nInit + 1].id( resolved[nInit + 2].id() );      
    resolved.pop_back();
    hasJunctionBeam = false;
  } 

  // Find companion quarks to unmatched sea quarks.
  for (int i = 0; i < nInit; ++i)       
  if (resolved[i].isUnmatched()) {

    // Add companion quark to record; and bookkeep both ways.
    append(0, -resolved[i].id(), 0., i);
    resolved[i].companion(size() - 1);
  }

  // If no other remnants found, add a gluon to carry remaining momentum.
  if ( size() == nInit) append(0, 21, 1., -1);     

  // Set initiator and remnant masses.
  for (int i = 0; i < size(); ++i) { 
    if (i < nInit) resolved[i].m(0.);
    else resolved[i].m( ParticleDataTable::m0( resolved[i].id() ) );
  }

  // For debug purposes: reject beams with resolved junction topology.
  if (hasJunctionBeam && !allowJunction) return false; 

  // Done.
  return true;

}

//*********

// Correlate all initiators and remnants to make a colour singlet. 

bool BeamParticle::remnantColours(Event& event, vector<int>& colFrom,
  vector<int>& colTo) {

  // Copy initiator colour info from the event record to the beam.
  for (int i = 0; i < nInit; ++i) {
    int j =  resolved[i].line();
    resolved[i].cols( event[j].col(), event[j].acol()); 
  }

  // Pick initial colours for remnants.
  for (int i = nInit; i < size(); ++i) {
    int colType = ParticleDataTable::colType( resolved[i].id() );
    int col = (colType == 1 || colType == 2) ? event.nextColTag() : 0;
    int acol = (colType == -1 || colType == 2) ? event.nextColTag() : 0;
    resolved[i].cols( col, acol);
  }

  // Find number and position of valence quarks, of gluons, and
  // of sea-companion pairs (counted as gluons) in the beam remnants.
  vector<int> iVal;
  vector<int> iGlu;
  for (int i = 0; i < size(); ++i) {
    if ( resolved[i].isValence() ) iVal.push_back(i);
    if ( resolved[i].id() == 21 || (resolved[i].isCompanion() 
        && resolved[i].companion() > i) ) iGlu.push_back(i);
  }
      
  // Pick a valence quark to which gluons are attached. 
  // Do not resolve quarks in diquark. (More sophisticated??)
  int iValSel= iVal[0];
  if (iVal.size() == 2) {
    if ( abs(resolved[iValSel].id()) > 10 ) iValSel = iVal[1];
  } else {
    double rndmValSel = 3. * Rndm::flat();
    if (rndmValSel > 1.) iValSel= iVal[1]; 
    if (rndmValSel > 2.) iValSel= iVal[2]; 
  }

  // This valence quark defines initial (anti)colour.
  int iBeg = iValSel;
  bool hasCol = (resolved[iBeg].col() > 0) ? true : false; 
  int begCol = (hasCol) ? resolved[iBeg].col() : resolved[iBeg].acol();

  // Do random stepping through gluon/(sea+companion) list.
  vector<int> iGluRndm;
  for (int i = 0; i < int(iGlu.size()); ++i)
    iGluRndm.push_back( iGlu[i] );
  for (int iOrder = 0; iOrder < int(iGlu.size()); ++iOrder) {
    int iRndm = int( double(iGluRndm.size()) * Rndm::flat()); 
    int iGluSel = iGluRndm[iRndm];
    iGluRndm[iRndm] = iGluRndm[iGluRndm.size() - 1];
    iGluRndm.pop_back();

    // Find matching anticolour/colour to current colour/anticolour.
    int iEnd = iGluSel;
    int endCol = (hasCol) ? resolved[iEnd].acol() : resolved[iEnd].col();
    // Not gluon but sea+companion pair: go to other.
    if (endCol == 0) {
      iEnd = resolved[iEnd].companion();
      endCol = (hasCol) ? resolved[iEnd].acol() : resolved[iEnd].col();
    }

    // Collapse this colour-anticolour pair to the lowest one.
    if (begCol < endCol) {
      if (hasCol) resolved[iEnd].acol(begCol); 
      else resolved[iEnd].col(begCol);
      colFrom.push_back(endCol);
      colTo.push_back(begCol);
    } else {
      if (hasCol) resolved[iBeg].col(endCol); 
      else resolved[iBeg].acol(endCol);
      colFrom.push_back(begCol);
      colTo.push_back(endCol);
    }

    // Pick up the other colour of the recent gluon and repeat.
    iBeg = iEnd;  
    begCol = (hasCol) ? resolved[iBeg].col() : resolved[iBeg].acol();
    // Not gluon but sea+companion pair: go to other.
    if (begCol == 0) {
      iBeg = resolved[iBeg].companion();
      begCol = (hasCol) ? resolved[iBeg].col() : resolved[iBeg].acol();
    }

  // At end of gluon/(sea+companion) list.
  } 
 
  // Now begin checks, and also finding junction information.
  // Loop through remnant partons; isolate all colours and anticolours.  
  vector<int> colList;
  vector<int> acolList;
  for (int i = 0; i < size(); ++i) {  
    if (resolved[i].col() > 0) colList.push_back( resolved[i].col() ); 
    if (resolved[i].acol() > 0) acolList.push_back( resolved[i].acol() ); 
  }

  // Remove all matching colour-anticolour pairs.
  bool foundPair = true;
  while (foundPair && colList.size() > 0 && acolList.size() > 0) {
    foundPair = false;
    for (int iCol = 0; iCol < int(colList.size()); ++iCol) {
      for (int iAcol = 0; iAcol < int(acolList.size()); ++iAcol) {
	if (acolList[iAcol] == colList[iCol]) { 
          colList[iCol] = colList.back(); colList.pop_back();     
          acolList[iAcol] = acolList.back(); acolList.pop_back();     
          foundPair = true; break;
	}
      } if (foundPair) break;
    }
  } 

  // Usually one unmatched pair left to collapse.
  if (colList.size() == 1 && acolList.size() == 1) {
    int finalFrom = max( colList[0], acolList[0]);
    int finalTo = min( colList[0], acolList[0]);
    for (int i = 0; i < size(); ++i) {  
      if (resolved[i].col() == finalFrom) resolved[i].col(finalTo); 
      if (resolved[i].acol() == finalFrom) resolved[i].acol(finalTo); 
    }
    colFrom.push_back(finalFrom);
    colTo.push_back(finalTo);

  // Store an (anti)junction when three (anti)coloured daughters.
  } else if (hasJunctionBeam && colList.size() == 3 
    && acolList.size() == 0) {
    event.appendJunction( 1, colList[0], colList[1], colList[2]);
    junCol[0] = colList[0]; 
    junCol[1] = colList[1]; 
    junCol[2] = colList[2];
  } else if (hasJunctionBeam && acolList.size() == 3 
    && colList.size() == 0) {
    event.appendJunction( 2, acolList[0], acolList[1], acolList[2]);
    junCol[0] = acolList[0]; 
    junCol[1] = acolList[1]; 
    junCol[2] = acolList[2];

  // Any other nonvanishing values indicate failure.
  } else if (colList.size() > 0 || acolList.size() > 0) {
    ErrorMessages::message("Error in BeamParticle::remnantColours: "
      "leftover unmatched colours"); 
    return false;
  }

  // Done.
  return true;
}
   
//*********

// Print the list of resolved partons in a beam.

void BeamParticle::list(ostream& os) {

  // Header.
  os << "\n --------  Partons resolved in beam  -------------------" 
     << "-------------------------------------------------------\n "
     << "   i  line      id       x    comp   xqcomp      colours  " 
     << "    p_x        p_y        p_z         e          m \n";
  
  // Loop over list of removed partons and print it. 
  double xSum = 0.;
  Vec4 pSum;
  for (int i = 0; i < size(); ++i) {
    ResolvedParton res = resolved[i];
    os << fixed << setprecision(6) << setw(5) << i << setw(6) << res.line() 
       << setw(8) << res.id() << setw(10) << res.x() << setw(6) 
       << res.companion() << setw(10) << res.xqCompanion() 
       << setprecision(3) << setw(6) << res.col() << setw(6) << res.acol() 
       << setw(11) << res.px() << setw(11) << res.py() << setw(11) 
       << res.pz() << setw(11) << res.e() << setw(11) << res.m() << "\n";

    // Also find and print sum of x and p values.
    xSum += res.x();  
    pSum += res.p();
  }
  os << setprecision(6) << "             x sum:" << setw(10) << xSum 
     << setprecision(3) << "                      p sum:" << setw(11) 
     << pSum.px() << setw(11) << pSum.py() << setw(11) << pSum.pz() 
     << setw(11) << pSum.e() << "\n"; 

}

//**************************************************************************

// The BeamRemnants class.

//*********
 
// Definitions of static variables.
// (Values will be overwritten in initStatic call, so are purely dummy.)

double BeamRemnants::primordialKTwidth = 1.;
double BeamRemnants::valencePowerMeson = 0.8;
double BeamRemnants::valencePowerUinP = 3.5;
double BeamRemnants::valencePowerDinP = 2.0;
double BeamRemnants::valenceDiqEnhance = 2.0;
int BeamRemnants::companionPower = 4;

// Constants: could be changed here if desired, but normally should not.
// These are of technical nature, as described for each.

// Maximum number of tries to match colours and kinematics in the event.
const int BeamRemnants::NTRYCOLMATCH = 20; 
const int BeamRemnants::NTRYKINMATCH = 10; 

//*********

// Initialize static data members.

void BeamRemnants::initStatic() {

  // Width of primordial kT distribution.
  primordialKTwidth = Settings::parameter("Beams:primordialKTwidth");

  // Power of (1-x)^power/sqrt(x) for remnant valence quark distribution.
  valencePowerMeson = Settings::parameter("Beams:valencePowerMeson");
  valencePowerUinP = Settings::parameter("Beams:valencePowerUinP");
  valencePowerDinP = Settings::parameter("Beams:valencePowerDinP");

  // Enhancement factor of x of diquark.
  valenceDiqEnhance = Settings::parameter("Beams:valenceDiqEnhance");

  // Assume g(x) ~ (1-x)^power/x to constrain companion to sea quark.
  companionPower = Settings::mode("Beams:companionPower");

}

//*********

// Set up trial transverse and longitudinal kinematics for each beam 
// separately. Final decisions involve comparing the two beams.

bool BeamRemnants::add( BeamParticle& beamA, BeamParticle& beamB,
  Event& event) {

  // Temporary solution: do not treat remnants in e+e- annihilation
  // (and similar). Assumes no ISR so no photon remnants ??
  int idAabs = abs(beamA.id());
  int idBabs = abs(beamB.id());
  if (idAabs < 20 || idBabs < 20) return true;

  // Total and squared CM energy.
  double eCM = event[0].m();
  double s = eCM * eCM;

  // Add required extra remnant flavour content. 
  // Start all over if fails (in option where junctions not allowed).
  if (!beamA.remnantFlavours()) return false;
  if (!beamB.remnantFlavours()) return false;

  // Save current complete colour configuration for fast restoration.
  vector<int> colSave;
  vector<int> acolSave;
  // cout << " begin push_back " << endl;
  for (int i = 0; i < event.size(); ++i) {
    // cout << " push_back 1 i = " << i << endl;
    colSave.push_back( event[i].col() );
    // cout << " push_back 2 i = " << i << endl;
    acolSave.push_back( event[i].acol() );
  }
  //for (int i = 0; i < event.size(); ++i)
  //  colSave.push_back( event[i].col() );
  //for (int i = 0; i < event.size(); ++i)
  //  acolSave.push_back( event[i].acol() );
  // cout << " end push_back " << endl;
  event.saveJunctionSize();
  
  // Allow several tries to set colours of initiators and remnants.
  // Frequent "failures" since shortcutting colours separately on
  // the two event sides may give "colour singlet gluons" etc.
  bool physical;
  for (int iTry = 0; iTry < NTRYCOLMATCH; ++iTry) {
    physical = true;

    // Reset list of colour "collapses" (transformations).
    colFrom.resize(0);
    colTo.resize(0);      

    // First process each set of beam colours on its own.
    if (!beamA.remnantColours(event, colFrom, colTo)) physical = false;
    if (!beamB.remnantColours(event, colFrom, colTo)) physical = false;

    // Then check that colours and anticolours are matched in whole event.
    if (physical && !checkColours(beamA, beamB, event)) physical = false;     

    // If no problems then done, else restore and loop.
    if (physical) break;
    for (int i = 0; i < event.size(); ++i) {
      event[i].col( colSave[i] );
      event[i].acol( acolSave[i] );
    }
    event.restoreJunctionSize();

  // If no solution after several tries then failed.
  }
  if (!physical) {
    ErrorMessages::message("Error in BeamRemnants::add:"
      " colour tracing failed"); 
    return false;
  }

  // Allow ten tries to construct kinematics (but normally works first).
  double xSum[2], xInvM[2], w2Beam[2], wPosRem, wNegRem, w2Rem;
    for (int iTry = 0; iTry < NTRYKINMATCH; ++iTry) {
    physical = true;

    // Loop over the two beams. Sum px and py separately within each. 
    for (int iBeam = 0; iBeam < 2; ++iBeam) {
      BeamParticle& beam = (iBeam == 0) ? beamA : beamB; 
      double pxSum = 0.;
      double pySum = 0.;
 
      // Loop over the partons in a beam. Generate Gaussian pT.
      for (int i = 0; i < beam.size(); ++i) { 
        double px = primordialKTwidth * Rndm::gauss();
        double py = primordialKTwidth * Rndm::gauss();
        beam[i].px(px);
        beam[i].py(py);
        pxSum += px;
        pySum += py;
      }  

      // Share recoil evenly between all partons.
      double pxRec = -pxSum / beam.size();
      double pyRec = -pySum / beam.size();
      for (int i = 0; i < beam.size(); ++i) { 
        beam[i].px( beam[i].px() + pxRec );
        beam[i].py( beam[i].py() + pyRec );
      }

      // Pick unrescaled x values for remnants. Sum up (unscaled) p+ and p-.
      xSum[iBeam] = 0.;
      xInvM[iBeam] = 0.;
      for (int i = beam.sizeInit(); i < beam.size(); ++i) {
        double xPrel = xRemnant( i, beam); 
        beam[i].x(xPrel);
        xSum[iBeam] += xPrel;
        xInvM[iBeam] += beam[i].mT2()/xPrel;     
      }

      // Squared transverse mass for each beam, using lightcone x.
      w2Beam[iBeam] = xSum[iBeam] * xInvM[iBeam];
  
    // End separate treatment of the two beams. 
    } 
 
    // Recalculate kinematics of initiator systems with primordial kT.
    wPosRem = eCM;
    wNegRem = eCM;
    for (int i = 0; i < beamA.sizeInit(); ++i) { 
      double sHat = beamA[i].x() * beamB[i].x() * s;
      double sHatT = sHat 
        + pow2( beamA[i].px() + beamB[i].px()) 
        + pow2( beamA[i].py() + beamB[i].py()); 
      double rescale = sqrt( sHatT / sHat);
      wPosRem -= rescale * beamA[i].x() * eCM;
      wNegRem -= rescale * beamB[i].x() * eCM;

      // Kinematics forbids too large primordial kT.
      if (sqrt(sHatT) < beamA[i].pT() + beamB[i].pT()) physical = false;
    }

    // Check that remaining momentum is enough for remnants.
    if (wPosRem < 0. || wNegRem < 0.) physical = false;
    w2Rem = wPosRem * wNegRem;
    if (sqrt(w2Rem) < sqrt(w2Beam[0]) + sqrt(w2Beam[1])) physical = false;

    // End of loop over ten tries. Do not loop when solution found.  
    if (physical) break;
  }
  // If no solution after ten tries then failed.
  if (!physical) return false;

  // Construct energy and pz for each initiator pair.
  for (int i = 0; i < beamA.sizeInit(); ++i) { 
    double sHat = beamA[i].x() * beamB[i].x() * s;
    double sHatT = sHat + pow2( beamA[i].px() + beamB[i].px()) 
      + pow2( beamA[i].py() + beamB[i].py()); 
    double rescale = sqrt( sHatT / sHat);
    double wPos = rescale * beamA[i].x() * eCM;
    double wNeg = rescale * beamB[i].x() * eCM;
    double w2A = beamA[i].mT2();
    double w2B = beamB[i].mT2();
    double lambdaRoot = sqrtpos( pow2( sHatT - w2A - w2B) - 4. * w2A * w2B );
    double pPosA = 0.5 * (sHatT + w2A - w2B + lambdaRoot) / sHatT * wPos;
    beamA[i].e( 0.5 * (pPosA + w2A / pPosA) );
    beamA[i].pz( 0.5 * (pPosA - w2A / pPosA) );
    double pNegB = 0.5 * (sHatT + w2B - w2A + lambdaRoot) / sHatT * wNeg;
    beamB[i].e( 0.5 * (pNegB + w2B / pNegB) );
    beamB[i].pz( 0.5 * (w2B / pNegB - pNegB) );

    // Construct rotations and boosts caused by primordial kT.
    int iA = beamA[i].line();
    int iB = beamB[i].line();
    RotBstMatrix M;
    M.toCMframe( event[iA].p(), event[iB].p() );
    M.fromCMframe( beamA[i].p(), beamB[i].p() );
 
    // Copy initiators and their systems and boost them accordingly.
    int iAcopy = event.copy(iA, -61);
    event[iAcopy].rotbst(M);
    int iBcopy = event.copy(iB, -61);
    event[iBcopy].rotbst(M);
    for (int iAB = max(iA,iB) + 1; iAB < event.size(); ++iAB) {
      if (!event[iAB].remains()) break;
      int iABcopy = event.copy(iAB, 62);
      event[iABcopy].rotbst(M); 
    }

    // Update daughter info of mothers, i.e. of beams, for hardest interaction.
    if (i == 0) { 
      int mother = event[iAcopy].mother1();
      event[mother].daughter1(iAcopy);      
      mother = event[iBcopy].mother1();
      event[mother].daughter1(iBcopy);      
    }
  }

  // Construct x rescaling factors for the two remants.
  double lambdaRoot = sqrtpos( pow2(w2Rem - w2Beam[0] - w2Beam[1])
    - 4. * w2Beam[0] * w2Beam[1] );
  double rescaleA = (w2Rem + w2Beam[0] - w2Beam[1] + lambdaRoot)
    / (2. * w2Rem * xSum[0]) ;
  double rescaleB = (w2Rem + w2Beam[1] - w2Beam[0] + lambdaRoot)
    / (2. * w2Rem * xSum[1]) ;

  // Construct energy and pz for remnants in first beam.
  for (int i = beamA.sizeInit(); i < beamA.size(); ++i) {
    double pPos = rescaleA * beamA[i].x() * wPosRem;
    double pNeg = beamA[i].mT2() / pPos;
    beamA[i].e( 0.5 * (pPos + pNeg) );
    beamA[i].pz( 0.5 * (pPos - pNeg) );  

    // Add these partons to the normal event record.
    int iNew = event.append( beamA[i].id(), 63, 1, 0, 0, 0, 
      beamA[i].col(), beamA[i].acol(), beamA[i].p(), beamA[i].m() );  
    beamA[i].line( iNew);
  }

  // Construct energy and pz for remnants in second beam.
  for (int i = beamB.sizeInit(); i < beamB.size(); ++i) {
    double pNeg = rescaleB * beamB[i].x() * wNegRem;
    double pPos = beamB[i].mT2() / pNeg;
    beamB[i].e( 0.5 * (pPos + pNeg) );
    beamB[i].pz( 0.5 * (pPos - pNeg) );  

    // Add these partons to the normal event record.
    int iNew = event.append( beamB[i].id(), 63, 2, 0, 0, 0,  
      beamB[i].col(), beamB[i].acol(), beamB[i].p(), beamB[i].m() );  
    beamB[i].line( iNew);
  }

  // Done.
  return true;

}

//*********

// Pick unrescaled x values for beam remnant sharing.

double BeamRemnants::xRemnant( int i, BeamParticle& beam) {

  double x = 0.;

  // Calculation of x of valence quark or diquark, for latter as sum.
  if (beam[i].isValence()) {  

    // Resolve diquark into sum of two quarks.
    int id1 = beam[i].id();
    int id2 = 0;
    if (abs(id1) > 10) {
      id2 = (id1 > 0) ? (id1/100)%10 : -(((-id1)/100)%10);
      id1 = (id1 > 0) ? id1/1000 : -((-id1)/1000);
    }
 
    // Loop over (up to) two quarks; add their contributions.
    for (int iId = 0; iId < 2; ++iId) {
      int id = (iId == 0) ? id1 : id2;
      if (id == 0) break;
      double xPart = 0.; 

      // Assume form (1-x)^a / sqrt(x).
      double xPow = valencePowerMeson;
      if (beam.isBaryon()) {
        if (beam.nValenceKinds() == 3 || beam.nValenceKinds() == 1) 
          xPow = (3. * Rndm::flat() < 2.) 
            ? valencePowerUinP : valencePowerDinP ; 
        else if (beam.nValence(id) == 2) xPow = valencePowerUinP;
        else xPow = valencePowerDinP;
      }
      do xPart = pow2( Rndm::flat() );
      while ( pow(1. - xPart, xPow) < Rndm::flat() ); 

      // End loop over (up to) two quarks. Possibly enhancement for diquarks.
      x += xPart; 
    }
   if (id2 != 0) x *= valenceDiqEnhance;
      
  // Calculation of x of sea quark, based on companion association.
  } else if (beam[i].isCompanion()) {

    // Find rescaled x value of companion.
    double xLeft = 1.;
    for (int iInit = 0; iInit < beam.sizeInit(); ++iInit) 
      xLeft -= beam[iInit].x();
    double xCompanion = beam[ beam[i].companion() ].x();
    xCompanion /= (xLeft + xCompanion);  

    // Now use ansatz q(x; x_c) < N/(x +x_c) to pick x.
    do x = pow( xCompanion, Rndm::flat()) - xCompanion; 
    while ( pow( 1. - x - xCompanion, companionPower) * (pow2(x) 
      + pow2(xCompanion)) / pow2(x + xCompanion) < Rndm::flat() );

  // Else, rarely, a single gluon remnant, so value does not matter. 
  } else x = 1.;
  return x;

}

//*********

// Collapse colours and check that they are consistent.

bool BeamRemnants::checkColours( BeamParticle& beamA, BeamParticle& beamB,
  Event& event) {

  // Remove ambiguities when one colour collapses two ways.
  // Resolve chains where one colour is mapped to another.
  for (int iCol = 1; iCol < int(colFrom.size()); ++iCol) 
  for (int iColRef = 0; iColRef < iCol; ++iColRef) { 
    if (colFrom[iCol] == colFrom[iColRef]) {
      colFrom[iCol] = colTo[iCol];
      colTo[iCol] = colTo[iColRef]; 
    }
    if (colTo[iCol] == colFrom[iColRef]) colTo[iCol] = colTo[iColRef];
  }   
  
  // Transform event record colours from beam remnant colour collapses.
  for (int i = 0; i < event.size(); ++i) { 
    int col = event[i].col();
    int acol = event[i].acol(); 
    for (int iCol = 0; iCol < int(colFrom.size()); ++iCol) {
      if (col == colFrom[iCol]) {col = colTo[iCol]; event[i].col(col);} 
      if (acol == colFrom[iCol]) {acol = colTo[iCol]; event[i].acol(acol);} 
    }
  }

  // Transform junction colours from beam remnant colour collapses.
  for (int iJun = 0; iJun < event.sizeJunction(); ++iJun)
  for (int leg = 0; leg < 3; ++leg) {
    int col = event.colJunction(iJun, leg); 
    for (int iCol = 0; iCol < int(colFrom.size()); ++iCol) 
    if (col == colFrom[iCol]) {col = colTo[iCol]; 
      event.colJunction(iJun, leg, col);} 
  }

  // Transform beam remnant colours from beam remnant colour collapses.
  for (int iBeam = 0; iBeam < 2; ++iBeam) {
    BeamParticle& beam = (iBeam == 0) ? beamA : beamB; 
    for (int i = 0; i < beam.size(); ++i) { 
      int col = beam[i].col();
      int acol = beam[i].acol(); 
      for (int iCol = 0; iCol < int(colFrom.size()); ++iCol) {
        if (col == colFrom[iCol]) {col = colTo[iCol]; beam[i].col(col);} 
        if (acol == colFrom[iCol]) {acol = colTo[iCol]; beam[i].acol(acol);} 
      }
    }

    // Transform beam junction colours from beam remnant colour collapses.
    if (beam.hasJunction()) {
      for (int i = 0; i < 3; ++i) { 
        int col = beam.junctionCol(i);
        for (int iCol = 0; iCol < int(colFrom.size()); ++iCol)
        if (col == colFrom[iCol]) {col = colTo[iCol]; 
          beam.junctionCol(i, col);} 
      }
    }
  }

  // Arrays for current colours and anticolours.
  vector<int> colList;
  vector<int> acolList;

  // Find current colours and anticolours in the event record.
  for (int i = 0; i < event.size(); ++i) if (event[i].status() > 0) {
    int id = event[i].id();
    int col = event[i].col();
    int acol = event[i].acol(); 

    // Quarks must have colour set, antiquarks anticolour, gluons both.
    if ( (id > 0 && id < 9 && (col <= 0 || acol != 0) )
      || (id < 0 && id > -9 && (col != 0 || acol <= 0) )
      || (id == 21 && (col <= 0 || acol <= 0) ) ) {
      ErrorMessages::message("Error in BeamRemnants::checkColours: "
        "q/qbar/g has wrong colour slots set");
      return false;
    }

    // No gluon can have the same anticolour as colour.
    // Instead of failure: insert singlet gluon on nearest string piece, 
    // defined by (p_i*p_k)(p_k*p_j)/(p_i*p_j)??
    if (col > 0 && acol == col) return false;
    if (col > 0) colList.push_back( col );
    if (acol > 0) acolList.push_back( acol );
  }

  // Add colours from the two beam remnants, not yet in the event record. 
  for (int i = beamA.sizeInit(); i < beamA.size(); ++i) {
    int col = beamA[i].col();
    int acol = beamA[i].acol(); 
    // No gluon can have the same anticolour as colour.
    if (col > 0 && acol == col) return false;
    if (col > 0) colList.push_back( col );
    if (acol > 0) acolList.push_back( acol );
  }
  for (int i = beamB.sizeInit(); i < beamB.size(); ++i) {
    int col = beamB[i].col();
    int acol = beamB[i].acol(); 
    // No gluon can have the same anticolour as colour.
    if (col > 0 && acol == col) return false;
    if (col > 0) colList.push_back( col );
    if (acol > 0) acolList.push_back( acol );
  }

  // Check that not the same colour or anticolour appears twice.
  for (int iCol = 0; iCol < int(colList.size()) - 1; ++iCol) {
    int col = colList[iCol];
    for (int iCol2 = iCol + 1; iCol2 < int(colList.size()); ++iCol2) 
      if (colList[iCol2] == col) return false;
  }
  for (int iAcol = 0; iAcol < int(acolList.size()) - 1; ++iAcol) {
    int acol = acolList[iAcol];
    for (int iAcol2 = iAcol + 1; iAcol2 < int(acolList.size()); ++iAcol2) 
      if (acolList[iAcol2] == acol) return false;
  }

  // Remove all matching colour-anticolour pairs.
  bool foundPair = true;
  while (foundPair && colList.size() > 0 && acolList.size() > 0) {
    foundPair = false;
    for (int iCol = 0; iCol < int(colList.size()); ++iCol) {
      for (int iAcol = 0; iAcol < int(acolList.size()); ++iAcol) {
	if (acolList[iAcol] == colList[iCol]) { 
          colList[iCol] = colList.back(); colList.pop_back();     
          acolList[iAcol] = acolList.back(); acolList.pop_back();     
          foundPair = true; break;
	}
      } if (foundPair) break;
    }
  } 

  // Check that remaining (anti)colours are accounted for by junctions.
  for (int iBeam = 0; iBeam < 2; ++iBeam) {
    BeamParticle& beam = (iBeam == 0) ? beamA : beamB;
    if (beam.hasJunction() && beam.id() > 0) {
      for (int iColEnd = 0; iColEnd < 3; ++iColEnd) {
        int colEnd = beam.junctionCol(iColEnd);
        bool foundCol = false;
        for (int iCol = 0; iCol < int(colList.size()); ++iCol) 
	if (colList[iCol] == colEnd) { 
          colList[iCol] = colList.back(); colList.pop_back();     
          foundCol = true; break;
	}  
      } 
    } else if (beam.hasJunction()) {
      for (int iColEnd = 0; iColEnd < 3; ++iColEnd) {
        int colEnd = beam.junctionCol(iColEnd);
        bool foundCol = false;
        for (int iAcol = 0; iAcol < int(acolList.size()); ++iAcol) 
	if (acolList[iAcol] == colEnd) { 
          acolList[iAcol] = acolList.back(); acolList.pop_back();     
          foundCol = true; break;
	}
      }
    }
  }

  // Done.
  return (colList.size() == 0 && acolList.size() == 0) ? true : false;

}

//**************************************************************************

} // end namespace Pythia8
