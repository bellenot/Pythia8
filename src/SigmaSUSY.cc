// SigmaSUSY.cc is a part of the PYTHIA event generator.
// Copyright (C) 2009 Peter Skands, Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Function definitions (not found in the header) for the 
// supersymmetry simulation classes. 

#include "SigmaSUSY.h"

namespace Pythia8 {

//**************************************************************************

// Some explicit declarations of static member variables

bool            CoupSUSY::isInit                         = false;
SusyLesHouches* CoupSUSY::slhaPtr;
double          CoupSUSY::mW;
double          CoupSUSY::mZ;
double          CoupSUSY::mWpole;
double          CoupSUSY::wWpole;
double          CoupSUSY::mZpole;
double          CoupSUSY::wZpole;
double          CoupSUSY::sin2W;
double          CoupSUSY::sinW;
double          CoupSUSY::cosW;
double          CoupSUSY::tanb;
double          CoupSUSY::sinb;
double          CoupSUSY::cosb;

complex         CoupSUSY::Lsddg[7][4];
complex         CoupSUSY::Rsddg[7][4];
complex         CoupSUSY::Lsuug[7][4];
complex         CoupSUSY::Rsuug[7][4];

complex         CoupSUSY::OLpp[6][6]; 
complex         CoupSUSY::ORpp[6][6];
complex         CoupSUSY::OLp[3][3];
complex         CoupSUSY::ORp[3][3];
complex         CoupSUSY::OL[6][3];
complex         CoupSUSY::OR[6][3];

complex         CoupSUSY::LqqZ[7];
complex         CoupSUSY::RqqZ[7]; 
complex         CoupSUSY::LsdsdZ[7][7];
complex         CoupSUSY::RsdsdZ[7][7]; 
complex         CoupSUSY::LsusuZ[7][7]; 
complex         CoupSUSY::RsusuZ[7][7]; 
complex         CoupSUSY::LudW[4][4];
complex         CoupSUSY::RudW[4][4];
complex         CoupSUSY::LsusdW[7][7];
complex         CoupSUSY::RsusdW[7][7];
complex         CoupSUSY::LsddX[7][4][6];
complex         CoupSUSY::RsddX[7][4][6];
complex         CoupSUSY::LsuuX[7][4][6];
complex         CoupSUSY::RsuuX[7][4][6];
complex         CoupSUSY::LsduX[7][4][3]; 
complex         CoupSUSY::RsduX[7][4][3];
complex         CoupSUSY::LsudX[7][4][3];
complex         CoupSUSY::RsudX[7][4][3];

const bool      CoupSUSY::DEBUG                          = false;

//**************************************************************************

// CoupSUSY::initStatic 
// Initialize SM+SUSY couplings (only performed once)

void CoupSUSY::initStatic (SusyLesHouches* slhaPtrIn) {

  // Set pointer to SLHA
  slhaPtr = slhaPtrIn;
  
  // Is NMSSM switched on?
  bool nmssm = (slhaPtr->modsel(3) != 1 ? false : true);
  int nNeut = (nmssm ? 5 : 4);
  int nChar = 2;

  // Initialize pole masses 
  mZpole    = ParticleDataTable::m0(23);
  wZpole    = ParticleDataTable::mWidth(23);
  mWpole    = ParticleDataTable::m0(24);
  wWpole    = ParticleDataTable::mWidth(24);
  
  // Running masses and weak mixing angle 
  // (default to pole values if no running available)
  mW        = mWpole;
  mZ        = mZpole;
  sin2W     = 1.0 - pow(mW/mZ,2);  

  if (slhaPtr->gauge.exists(1) && slhaPtr->gauge.exists(2) 
      && slhaPtr->hmix.exists(3)) {
    double gp=slhaPtr->gauge(1);
    double g =slhaPtr->gauge(2);
    double v =slhaPtr->hmix(3);
    mW      = g * v / 2.0;
    mZ      = sqrt(pow(gp,2)+pow(g,2)) * v / 2.0;
    double tan2W   = pow2(gp)/pow2(g);
    if (DEBUG) cout << " tan2W = " << tan2W << endl;
    sin2W   = pow2(gp)/(pow2(g)+pow2(gp));  
  }
  sinW = sqrt(sin2W);
  cosW = sqrt(1.0-sin2W);

  // Tan(beta)
  // By default, use the running one in HMIX (if not found, use MINPAR)
  tanb = slhaPtr->hmix.exists(2) ? slhaPtr->hmix(2) : slhaPtr->minpar(3);
  cosb = sqrt( 1.0 / (1.0 + tanb*tanb) );
  sinb = sqrt(max(0.0,1.0-cosb*cosb));
  
  // tmp : verbose output
  if (DEBUG) {
    cout << " sin2W(Q) = " << sin2W << "  mW(Q) = " << mW 
         << "  mZ(Q) = " << mZ << endl;
    cout << " vev(Q) = " << slhaPtr->hmix(3) << " tanb(Q) = " << tanb
         << endl;
    for (int i=1;i<=3;i++) {
      for (int j=1;j<=3;j++) {
	cout << " VCKM  [" << i << "][" << j << "] = " 
             << scientific << setw(10) << VCKM::Vgen(i,j) << endl;
      }
    }
  }  
  
  // Shorthand for squark mixing matrices 
  SusyLesHouches::matrixblock<6> Ru(slhaPtr->usqmix);
  SusyLesHouches::matrixblock<6> Rd(slhaPtr->dsqmix);
  SusyLesHouches::matrixblock<6> imRu(slhaPtr->imusqmix);
  SusyLesHouches::matrixblock<6> imRd(slhaPtr->imusqmix);
  
  // Construct ~g couplings
  for (int i=1 ; i<=6 ; i++) {
    for (int j=1 ; j<=3 ; j++) {
      Lsddg[i][j] = complex( Rd(i,j)  ,  imRd(i,j));
      Rsddg[i][j] = complex(-Rd(i,j+3), -imRd(i,j+3));
      Lsuug[i][j] = complex( Ru(i,j)  ,  imRu(i,j));
      Rsuug[i][j] = complex(-Ru(i,j+3), -imRu(i,j+3));
    }
  }
  
  // Construct qqZ couplings
  for (int i=1 ; i<=6 ; i++) {
    
    // q[i] q[i] Z (def with extra factor 2 compared to [Okun])
    LqqZ[i] = CoupEW::af(i) - 2.0*CoupEW::ef(i)*sin2W ;
    RqqZ[i] =               - 2.0*CoupEW::ef(i)*sin2W ;

    // tmp: verbose output
    if (DEBUG) {
      cout << " LqqZ  [" << i << "][" << i << "] = " 
           << scientific << setw(10) << LqqZ[i] 
           << " RqqZ  [" << i << "][" << i  << "] = " 
           << scientific << setw(10) << RqqZ[i] << endl;
    }
  }

  // Construct ~q~qZ couplings
  for (int i=1 ; i<=6 ; i++) {

    // Squarks can have off-diagonal couplings as well
    for (int j=1 ; j<=6 ; j++) {
      
      // ~d[i] ~d[j] Z
      LsdsdZ[i][j] = 0.0;
      RsdsdZ[i][j] = 0.0;
      for (int k=1;k<=3;k++) {
	complex Rdik  = complex(Rd(i,k),  imRd(i,k)  );
	complex Rdjk  = complex(Rd(j,k),  imRd(j,k)  );
	complex Rdik3 = complex(Rd(i,k+3),imRd(i,k+3));
	complex Rdjk3 = complex(Rd(j,k+3),imRd(j,k+3));
	LsdsdZ[i][j] += LqqZ[1] * (Rdik*conj(Rdjk)); 
	RsdsdZ[i][j] += RqqZ[1] * (Rdik3*conj(Rdjk3)); 
      }
      
      // ~u[i] ~u[j] Z
      LsusuZ[i][j] = 0.0;
      RsusuZ[i][j] = 0.0; 
      for (int k=1;k<=3;k++) {
	complex Ruik  = complex(Ru(i,k)  ,imRu(i,k)  );
	complex Rujk  = complex(Ru(j,k)  ,imRu(j,k)  );
	complex Ruik3 = complex(Ru(i,k+3),imRu(i,k+3));
	complex Rujk3 = complex(Ru(j,k+3),imRu(j,k+3));
	LsusuZ[i][j] += LqqZ[2] * (Ruik*conj(Rujk)); 
	RsusuZ[i][j] += RqqZ[2] * (Ruik3*conj(Rujk3)); 
      }
      
    // tmp: verbose output
      if (DEBUG) {
	if (max(abs(LsdsdZ[i][j]),abs(RsdsdZ[i][j])) > 1e-6) {
	  cout << " LsdsdZ[" << i << "][" << j << "] = " 
               << scientific << setw(10) << LsdsdZ[i][j]
	       << " RsdsdZ[" << i << "][" << j << "] = " 
               << scientific << setw(10) << RsdsdZ[i][j] << endl;
	}
	if (max(abs(LsusuZ[i][j]),abs(RsusuZ[i][j]))> 1e-6) {
	  cout << " LsusuZ[" << i << "][" << j << "] = " 
               << scientific << setw(10) << LsusuZ[i][j]
	       << " RsusuZ[" << i << "][" << j << "] = " 
               << scientific << setw(10) << RsusuZ[i][j] << endl;
	}
      }
	
    }
    
  }
  
  // Construct udW couplings
  
  // Loop over up [i] and down [j] quark generation
  for (int i=1;i<=3;i++) {
    for (int j=1;j<=3;j++) {
      
      // CKM matrix (use Pythia one if no SLHA)
      // (NB: could also try input one if no running one found, but
      // would then need to compute from Wolfenstein)
      complex Vij=VCKM::Vgen(i,j);
      if (slhaPtr->vckm.exists()) {
	Vij=complex(slhaPtr->vckm(i,j),slhaPtr->imvckm(i,j));
      }
      
      // u[i] d[j] W  
      LudW[i][j] = sqrt(2.0) * cosW * Vij;
      RudW[i][j] = 0.0;
            
      // tmp: verbose output
      if (DEBUG) {
	cout << " LudW  [" << i << "][" << j << "] = " 
             << scientific << setw(10) << LudW[i][j]
	     << " RudW  [" << i << "][" << j << "] = " 
             << scientific << setw(10) << RudW[i][j] << endl;
      }
    }
  }

  // Construct ~u~dW couplings

  // Loop over ~u[k] and ~d[l] flavours
  for (int k=1;k<=6;k++) {
    for (int l=1;l<=6;l++) {
	  
      LsusdW[k][l]=0.0; 
      RsusdW[k][l]=0.0;

      // Loop over u[i] and d[j] flavours
      for (int i=1;i<=3;i++) { 
	for (int j=1;j<=3;j++) {
	  
	  // CKM matrix (use Pythia one if no SLHA)
	  // (NB: could also try input one if no running one found, but
	  // would then need to compute from Wolfenstein)
	  complex Vij=VCKM::Vgen(i,j);
	  if (slhaPtr->vckm.exists()) {
	    Vij=complex(slhaPtr->vckm(i,j),slhaPtr->imvckm(i,j));
	  }
      
	  // ~u[k] ~d[l] W (add one term for each quark flavour i,j)
	  complex Ruki = complex(Ru(k,i),imRu(k,i));
	  complex Rdlj = complex(Rd(l,j),imRd(l,j));
	  LsusdW[k][l] += sqrt(2.0) * cosW * Vij * Ruki * conj(Rdlj);
	  RsusdW[k][l] += 0.0;
	  
	}
      }

      // tmp: verbose output
      if (DEBUG) {
	if (max(abs(LsusdW[k][l]),abs(RsusdW[k][l]))> 1e-6) {
	  cout << " LsusdW[" << k << "][" << l << "] = " 
               << scientific << setw(10) << LsusdW[k][l]
	       << " RsusdW[" << k << "][" << l << "] = " 
               << scientific << setw(10) << RsusdW[k][l] << endl;
	}
      }

    }
  }
  
  // Now we come to the ones with really many indices
  
  // Construct ~chi0 couplings (allow for 5 neutralinos in NMSSM)
  for (int i=1;i<=nNeut;i++) {
    
    // Ni1, Ni2, Ni3, Ni4, Ni5
    complex ni1,ni2,ni3,ni4,ni5;
    if (not nmssm) {	
      ni1=complex( slhaPtr->nmix(i,1), slhaPtr->imnmix(i,1) );
      ni2=complex( slhaPtr->nmix(i,2), slhaPtr->imnmix(i,2) );
      ni3=complex( slhaPtr->nmix(i,3), slhaPtr->imnmix(i,3) );
      ni4=complex( slhaPtr->nmix(i,4), slhaPtr->imnmix(i,4) );
      ni5=complex( 0.0, 0.0);
    } else {
      ni1=complex( slhaPtr->nmnmix(i,1), slhaPtr->imnmnmix(i,1) );
      ni2=complex( slhaPtr->nmnmix(i,2), slhaPtr->imnmnmix(i,2) );
      ni3=complex( slhaPtr->nmnmix(i,3), slhaPtr->imnmnmix(i,3) );
      ni4=complex( slhaPtr->nmnmix(i,4), slhaPtr->imnmnmix(i,4) );
      ni5=complex( slhaPtr->nmnmix(i,5), slhaPtr->imnmnmix(i,5) );
    }
    
    // Change to positive mass convention
    complex iRot( 0., 1.);
    if (slhaPtr->mass(idNeut(i)) < 0.) {
      ni1 *= iRot;
      ni2 *= iRot;
      ni3 *= iRot;
      ni4 *= iRot;
      ni5 *= iRot;
    }
    
    // ~chi0 [i] ~chi0 [j] Z : loop over [j]
    for (int j=1; j<=nNeut; j++) {
      
      // neutralino [j] higgsino components
      complex nj3, nj4;
      if (not nmssm) {
	nj3=complex( slhaPtr->nmix(j,3), slhaPtr->imnmix(j,3) );
	nj4=complex( slhaPtr->nmix(j,4), slhaPtr->imnmix(j,4) );
      } else {
	nj3=complex( slhaPtr->nmnmix(j,3), slhaPtr->imnmnmix(j,3) );
	nj4=complex( slhaPtr->nmnmix(j,4), slhaPtr->imnmnmix(j,4) );
      }
      // Change to positive mass convention
      if (slhaPtr->mass(idNeut(j)) < 0.) {
	nj3 *= iRot;
	nj4 *= iRot;
      }
      
      // ~chi0 [i] ~chi0 [j] Z : couplings
      OLpp[i][j] = -0.5 * ni3 * conj(nj3) + 0.5 * ni4 * conj(nj4);
      ORpp[i][j] =  0.5 * conj(ni3) * nj3 - 0.5 * conj(ni4) * nj4;
      
    // tmp: verbose output
      if (DEBUG) {
	cout << " OL''  [" << i << "][" << j << "] = " 
             << scientific << setw(10) << OLpp[i][j]
	     << " OR''  [" << i << "][" << j << "] = " 
             << scientific << setw(10) << ORpp[i][j] << endl;
      }
	
    }
    
    // ~chi0 [i] ~chi+ [j] W : loop over [j]
    for (int j=1; j<=nChar; j++) {
      
      // Chargino mixing
      complex uj1, uj2, vj1, vj2;      
      uj1=complex( slhaPtr->umix(j,1), slhaPtr->imumix(j,1) );
      uj2=complex( slhaPtr->umix(j,2), slhaPtr->imumix(j,2) );
      vj1=complex( slhaPtr->vmix(j,1), slhaPtr->imvmix(j,1) );
      vj2=complex( slhaPtr->vmix(j,2), slhaPtr->imvmix(j,2) );
      
      // ~chi0 [i] ~chi+ [j] W : couplings
      OL[i][j] = -1.0/sqrt(2.0)*ni4*conj(vj2)+ni2*conj(vj1);
      OR[i][j] = 1.0/sqrt(2.0)*conj(ni3)*uj2+conj(ni2)*uj1;
      
    // tmp: verbose output
      if (DEBUG) {
	cout << " OL    [" << i << "][" << j << "] = " 
             << scientific << setw(10) << OL[i][j]
	     << " OR    [" << i << "][" << j << "] = " 
             << scientific << setw(10) << OR[i][j] << endl;
      }
    }
    
    // Charges
    double ed  = -1.0/3.0;
    double T3d = -0.5;
    double eu  =  2.0/3.0;
    double T3u =  0.5;
    
    // Loop over quark [k] generation
    for (int k=1;k<=3;k++) {
      
      // Set quark masses
      // Initial guess 0,0,0,mc,mb,mt with the latter from the PDT
      double mu = ParticleDataTable::m0(2*k);
      double md = ParticleDataTable::m0(2*k-1);
      if (k == 1) { mu=0.0 ; md=0.0; }
      if (k == 2) { md=0.0 ; mu=0.0; } 
      
      // Compute running mass from Yukawas and vevs if possible.
      if (slhaPtr->yd.exists() && slhaPtr->hmix.exists(3)) {
	double ykk=slhaPtr->yd(k,k);
	double v1=slhaPtr->hmix(3)/sqrt(1+pow(tanb,2));
	if (ykk > 0.0) md = ykk * v1 / sqrt(2.0) ;
      }
      if (slhaPtr->yu.exists() && slhaPtr->hmix.exists(3)) {
	double ykk=slhaPtr->yu(k,k);
	double v2=slhaPtr->hmix(3)/sqrt(1.0+1.0/pow(tanb,2));
	if (ykk > 0.0) mu = ykk * v2 / sqrt(2.0) ;
      }
      
      // tmp: verbose output
      if (DEBUG) {
	cout  <<  " Gen = " << k << " mu = " << mu << " md = " << md 
              << " yUU,DD = " << slhaPtr->yu(k,k) << "," 
              << slhaPtr->yd(k,k) << endl;
      }
      
      // Loop over squark [j] flavour
      for (int j=1;j<=6;j++) {
	
	// Squark mixing
	complex Rdjk  = complex(Rd(j,k),  imRd(j,k)  );
	complex Rdjk3 = complex(Rd(j,k+3),imRd(j,k+3));
	complex Rujk  = complex(Ru(j,k),  imRu(j,k)  );
	complex Rujk3 = complex(Ru(j,k+3),imRu(j,k+3));
	
	// ~d[j] d[k] ~chi0[i]
	double rt2 = sqrt(2.0);
	LsddX[j][k][i] = ((ed-T3d)*sinW/cosW*ni1 + T3d*ni2)*conj(Rdjk)/rt2
	  + md*ni3*conj(Rdjk3)/2.0/rt2/mW/cosb; 
	RsddX[j][k][i] = -ed*sinW/cosW*conj(ni1)*conj(Rdjk3)/rt2 
	  + md*conj(ni3)*conj(Rdjk)/2.0/rt2/mW/cosb;

	// ~u[j] u[k] ~chi0[i]
	LsuuX[j][k][i] = ((eu-T3u)*sinW/cosW*ni1 + T3u*ni2)*conj(Rujk)/rt2
	  + mu*ni4*conj(Rujk3)/2.0/rt2/mW/sinb;
	RsuuX[j][k][i] = -eu*sinW/cosW*conj(ni1)*conj(Rujk3)/rt2
	  + mu*conj(ni4)*conj(Rujk)/2.0/rt2/mW/sinb;

	if (DEBUG) {
	  if (abs(LsddX[j][k][i]) > 1e-6) {
	    // tmp: verbose output
	    cout << " LsddX[" << j << "][" << k << "][" << i << "] = "
		 << scientific << setw(10) << LsddX[j][k][i] << endl;
	  }
	  if (abs(RsddX[j][k][i]) > 1e-6) {
	    // tmp: verbose output
	    cout << " RsddX[" << j << "][" << k << "][" << i << "] = "
		 << scientific << setw(10) << RsddX[j][k][i] << endl;
	  }
	  if (abs(LsuuX[j][k][i]) > 1e-6) {
	    // tmp: verbose output
	    cout << " LsuuX[" << j << "][" << k << "][" << i << "] = "
		 << scientific << setw(10) << LsuuX[j][k][i] << endl;
	  }
	  if (abs(RsuuX[j][k][i]) > 1e-6) {
	    // tmp: verbose output
	    cout << " RsuuX[" << j << "][" << k << "][" << i << "] = "
		 << scientific << setw(10) << RsuuX[j][k][i] << endl;
	  }
	}
      }
      
    }
    
  }
  
  // Construct ~chi+ couplings 
  for (int i=1;i<=nChar;i++) {
    
    // Ui1, Ui2, Vi1, Vi2
    complex ui1,ui2,vi1,vi2;    
    ui1=complex( slhaPtr->umix(i,1), slhaPtr->imumix(i,1) );
    ui2=complex( slhaPtr->umix(i,2), slhaPtr->imumix(i,2) );
    vi1=complex( slhaPtr->vmix(i,1), slhaPtr->imvmix(i,1) );
    vi2=complex( slhaPtr->vmix(i,2), slhaPtr->imvmix(i,2) );
    
    // ~chi+ [i] ~chi- [j] Z : loop over [j]
    for (int j=1; j<=nChar; j++) {
      
      // Chargino mixing
      complex uj1, uj2, vj1, vj2;      
      uj1=complex( slhaPtr->umix(j,1), slhaPtr->imumix(j,1) );
      uj2=complex( slhaPtr->umix(j,2), slhaPtr->imumix(j,2) );
      vj1=complex( slhaPtr->vmix(j,1), slhaPtr->imvmix(j,1) );
      vj2=complex( slhaPtr->vmix(j,2), slhaPtr->imvmix(j,2) );
      
      // ~chi+ [i] ~chi- [j] Z : couplings
      OLp[i][j] = -vi1*conj(vj1) - 0.5*vi2*conj(vj2) 
	+ ( (i == j) ? sin2W : 0.0);
      ORp[i][j] = -conj(ui1)*uj1 - 0.5*conj(ui2)*uj2 
	+ ( (i == j) ? sin2W : 0.0);
      
      if (DEBUG) {
	// tmp: verbose output
	cout << " OL'   [" << i << "][" << j << "] = " 
             << scientific << setw(10) << OLp[i][j]
	     << " OR'   [" << i << "][" << j << "] = " 
             << scientific << setw(10) << ORp[i][j] << endl;
      }
    }
    
    // Loop over quark [l] flavour
    for (int l=1;l<=3;l++) {
      
      // Set quark [l] masses 
      // Initial guess 0,0,0,mc,mb,mt with the latter from the PDT
      double mul = ParticleDataTable::m0(2*l);
      double mdl = ParticleDataTable::m0(2*l-1);
      if (l == 1) { mul=0.0 ; mdl=0.0; }
      if (l == 2) { mdl=0.0 ; mul=0.0; } 
      
      // Compute running mass from Yukawas and vevs if possible.
      if (slhaPtr->yd.exists() && slhaPtr->hmix.exists(3)) {
	double yll=slhaPtr->yd(l,l);
	double v1=slhaPtr->hmix(3)/sqrt(1+pow(tanb,2));
	if (yll > 0.0) mdl = yll * v1 / sqrt(2.0) ;
      }
      if (slhaPtr->yu.exists() && slhaPtr->hmix.exists(3)) {
	double yll=slhaPtr->yu(l,l);
	double v2=slhaPtr->hmix(3)/sqrt(1.0+1.0/pow(tanb,2));
	if (yll > 0.0) mul = yll * v2 / sqrt(2.0) ;
      }
      
      // Loop over squark [j] flavour
      for (int j=1;j<=6;j++) {

	// Loop over off-diagonal quark [k] generation
	for (int k=1;k<=3;k++) {

	  // Set quark [k] masses 
	  // Initial guess 0,0,0,0,mb,mt with the latter from the PDT
	  double muk = ParticleDataTable::m0(2*k);
	  double mdk = ParticleDataTable::m0(2*k-1);
	  if (k == 1) { muk=0.0 ; mdk=0.0; }
	  if (k == 2) { mdk=0.0 ; muk=0.0; } 
      
	  // Compute running mass from Yukawas and vevs if possible.
	  if (slhaPtr->yd.exists() && slhaPtr->hmix.exists(3)) {
	    double ykk=slhaPtr->yd(k,k);
	    double v1=slhaPtr->hmix(3)/sqrt(1+pow(tanb,2));
	    if (ykk > 0.0) mdk = ykk * v1 / sqrt(2.0) ;
	  }
	  if (slhaPtr->yu.exists() && slhaPtr->hmix.exists(3)) {
	    double ykk=slhaPtr->yu(k,k);
	    double v2=slhaPtr->hmix(3)/sqrt(1.0+1.0/pow(tanb,2));
	    if (ykk > 0.0) muk = ykk * v2 / sqrt(2.0) ;
	  }	  

	  // CKM matrix (use Pythia one if no SLHA)
	  // (NB: could also try input one if no running one found, but
	  // would then need to compute from Wolfenstein)
	  complex Vlk=VCKM::Vgen(l,k);
	  complex Vkl=VCKM::Vgen(k,l);
	  if (slhaPtr->vckm.exists()) {
	    Vlk=complex(slhaPtr->vckm(l,k),slhaPtr->imvckm(l,k));
	    Vkl=complex(slhaPtr->vckm(k,l),slhaPtr->imvckm(k,l));
	  }
      
	  // Squark mixing
	  complex Rdjk  = complex(Rd(j,k),  imRd(j,k)  );
	  complex Rdjk3 = complex(Rd(j,k+3),imRd(j,k+3));
	  complex Rujk  = complex(Ru(j,k),  imRu(j,k)  );
	  complex Rujk3 = complex(Ru(j,k+3),imRu(j,k+3));

	  // sqrt(2)
	  double rt2 = sqrt(2.0);
	  
	  // ~d[j] u[l] ~chi+[i]
	  LsduX[j][l][i] += 0.5*(ui1*conj(Rdjk) 
				 - mdk*ui2*conj(Rdjk3)/rt2/mW/cosb)*Vlk;
	  RsduX[j][l][i] -= 0.5*mul*conj(vi2)*Vlk*Rdjk/rt2/mW/sinb; 

	  // ~u[j] d[l] ~chi+[i]
	  LsudX[j][l][i] += 0.5*(vi1*conj(Rujk)
				 - muk*vi2*conj(Rujk3)/rt2/mW/sinb)*conj(Vkl);
	  RsudX[j][l][i] -= 0.5*mdl*conj(ui2)*conj(Vkl)*conj(Rujk)/rt2/mW/cosb;

	}

	if (DEBUG) {
	  if (max(abs(LsduX[j][l][i]),abs(RsduX[j][l][i])) > 1e-6) {
	    // tmp: verbose output
	    cout << " LsduX[" << j << "][" << l << "][" << i << "] = "
		 << scientific << setw(10) << LsduX[j][l][i];
	    cout << " RsduX[" << j << "][" << l << "][" << i << "] = "
		 << scientific << setw(10) << RsduX[j][l][i] << endl;
	  }
	  if (max(abs(LsudX[j][l][i]),abs(RsudX[j][l][i])) > 1e-6) {
	    // tmp: verbose output
	    cout << " LsudX[" << j << "][" << l << "][" << i << "] = "
		 << scientific << setw(10) << LsudX[j][l][i];
	    cout << " RsudX[" << j << "][" << l << "][" << i << "] = "
		 << scientific << setw(10) << RsudX[j][l][i] << endl;
	  }
	}
	
      }
      
      
    }
    
  }
  
  // Let everyone know we are ready
  isInit = true;
  
}
  
//**************************************************************************

// Sigma2qqbar2chi0chi0 
// Cross section for gaugino pair production: neutralino pair

//*********

// Initialize process. 
  
void Sigma2qqbar2chi0chi0::initProc() {

  // First make sure CoupSUSY is initialized
  if (not CoupSUSY::isInit) CoupSUSY::initStatic(slhaPtr);

  // Construct name of process. 
  nameSave = "q q'bar -> " + ParticleDataTable::name(id3) + " " 
    + ParticleDataTable::name(id4) + " (q,q'=d,u,s,c,b)";

}

//*********

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour. 

void Sigma2qqbar2chi0chi0::sigmaKin() {

  // Common flavour-independent factor.
  sigma0 = 4 * (M_PI / sH2 / pow2(CoupSUSY::sin2W)) * pow2(alpEM) ; 

  // Factor 1/2 for identical final particles.
  if (id3 == id4) sigma0 *= 0.5;

  // Auxiliary factors for use below
  ui       = uH - s3;
  uj       = uH - s4;
  ti       = tH - s3;
  tj       = tH - s4;
  double sV= sH - pow2(CoupSUSY::mZpole);
  double d = pow2(sV) + pow2(CoupSUSY::mZpole * CoupSUSY::wZpole);
  propZ    = complex( sV / d, CoupSUSY::mZpole * CoupSUSY::wZpole / d);

}

//*********

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence. 

double Sigma2qqbar2chi0chi0::sigmaHat() {

  // Only allow quark-antiquark incoming states
  if (id1*id2 >= 0) {
    return 0.0;    
  }
  
  // Only allow incoming states with sum(charge) = 0
  if ((id1+id2) % 2 != 0) {
    return 0.0;    
  }

  // Shorthands
  int idAbs1    = abs(id1);  
  int idAbs2    = abs(id2);  

  // Flavour-dependent kinematics-dependent couplings.
  complex QuLL(0.0),QtLL(0.0),QuRR(0.0),QtRR(0.0);
  complex QuLR(0.0),QtLR(0.0),QuRL(0.0),QtRL(0.0);
  if (idAbs1 == idAbs2) {
    QuLL = CoupSUSY::LqqZ[idAbs1] * CoupSUSY::OLpp[id3chi][id4chi];
    QtLL = CoupSUSY::LqqZ[idAbs1] * CoupSUSY::ORpp[id3chi][id4chi];
    QuRR = CoupSUSY::RqqZ[idAbs1] * CoupSUSY::ORpp[id3chi][id4chi];
    QtRR = CoupSUSY::RqqZ[idAbs1] * CoupSUSY::OLpp[id3chi][id4chi];
    QuLL *= propZ / 4.0 / (1.0-CoupSUSY::sin2W);
    QtLL *= propZ / 4.0 / (1.0-CoupSUSY::sin2W);
    QuRR *= propZ / 4.0 / (1.0-CoupSUSY::sin2W);
    QtRR *= propZ / 4.0 / (1.0-CoupSUSY::sin2W);  
  }

  // Add t-channel squark flavour sums to QmXY couplings
  for (int ksq=1; ksq<=6; ksq++) {    

    // Flavour indices
    int ifl1 = (idAbs1+1) / 2;
    int ifl2 = (idAbs2+1) / 2;

    // squark id and squark-subtracted u and t
    int idsq=((ksq+2)/3)*1000000 + 2*((ksq-1) % 3) + (idAbs1+1) % 2 + 1;
    double msq2=pow(ParticleDataTable::m0(idsq),2);
    double usq     = uH - msq2;
    double tsq     = tH - msq2;

    // Couplings
    complex Lsqq1X3 = CoupSUSY::LsuuX[ksq][ifl1][id3chi];
    complex Lsqq1X4 = CoupSUSY::LsuuX[ksq][ifl1][id4chi];
    complex Lsqq2X3 = CoupSUSY::LsuuX[ksq][ifl2][id3chi];
    complex Lsqq2X4 = CoupSUSY::LsuuX[ksq][ifl2][id4chi];
    complex Rsqq1X3 = CoupSUSY::RsuuX[ksq][ifl1][id3chi];
    complex Rsqq1X4 = CoupSUSY::RsuuX[ksq][ifl1][id4chi];
    complex Rsqq2X3 = CoupSUSY::RsuuX[ksq][ifl2][id3chi];
    complex Rsqq2X4 = CoupSUSY::RsuuX[ksq][ifl2][id4chi];
    if (idAbs1 % 2 != 0) {
      Lsqq1X3 = CoupSUSY::LsddX[ksq][ifl1][id3chi];
      Lsqq1X4 = CoupSUSY::LsddX[ksq][ifl1][id4chi];
      Lsqq2X3 = CoupSUSY::LsddX[ksq][ifl2][id3chi];
      Lsqq2X4 = CoupSUSY::LsddX[ksq][ifl2][id4chi];
      Rsqq1X3 = CoupSUSY::RsddX[ksq][ifl1][id3chi];
      Rsqq1X4 = CoupSUSY::RsddX[ksq][ifl1][id4chi];
      Rsqq2X3 = CoupSUSY::RsddX[ksq][ifl2][id3chi];
      Rsqq2X4 = CoupSUSY::RsddX[ksq][ifl2][id4chi];      
    }

    // QuXY
    QuLL += conj(Lsqq1X4)*Lsqq2X3/usq;
    QuLR += conj(Lsqq1X4)*Rsqq2X3/usq;
    QuRL += conj(Rsqq1X4)*Lsqq2X3/usq;
    QuRR += conj(Rsqq1X4)*Rsqq2X3/usq;

    // QtXY
    QtLL -= conj(Lsqq1X3)*Lsqq2X4/tsq;
    QtLR -= conj(Lsqq1X3)*Rsqq2X4/tsq;
    QtRL -= conj(Rsqq1X3)*Lsqq2X4/tsq;
    QtRR -= conj(Rsqq1X3)*Rsqq2X4/tsq;

  }

  // Compute matrix element weight
  double weight = 0;
  // Average over separate helicity contributions
  // LL (ha = -1, hb = +1) (divided by 4 for average)            
  weight += norm(QuLL) * ui * uj + norm(QtLL) * ti * tj
    + 2 * real(conj(QuLL) * QtLL) * m3 * m4 * sH;
  // RR (ha =  1, hb = -1) (divided by 4 for average)        
  weight += norm(QtRR) * ti * tj + norm(QuRR) * ui * uj  
    + 2 * real(conj(QuRR) * QtRR) * m3 * m4 * sH;
  // RL (ha =  1, hb =  1) (divided by 4 for average)        
  weight += norm(QuRL) * ui * uj + norm(QtRL) * ti * tj
    - real(conj(QuRL) * QtRL) * (uH * tH - s3 * s4);
  // LR (ha = -1, hb = -1) (divided by 4 for average)        
  weight += norm(QuLR) * ui * uj + norm(QtLR) * ti * tj
    - real(conj(QuLR) * QtLR) * (uH * tH - s3 * s4);

  // Cross section, including colour factor.
  double sigma = sigma0 * weight;
  if (idAbs1 < 9) sigma /= 3.;

  // Answer.
  return sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2qqbar2chi0chi0::setIdColAcol() {

  // Set flavours.
  setId( id1, id2, id3, id4);

  // Colour flow topologies. Swap when antiquarks.
  if (abs(id1) < 9) setColAcol( 1, 0, 0, 1, 0, 0, 0, 0);
  else              setColAcol( 0, 0, 0, 0, 0, 0, 0, 0);
  if (id1 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2qqbar2charchi0
// Cross section for gaugino pair production: neutralino-chargino

//*********

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour. 

void Sigma2qqbar2charchi0::sigmaKin() {

  // Common flavour-independent factor.
  sigma0 = 4 * (M_PI / sH2 / pow2(CoupSUSY::sin2W) ) * pow2(alpEM) ; 

  // Auxiliary factors for use below
  ui        = uH - s3;
  uj        = uH - s4;
  ti        = tH - s3;
  tj        = tH - s4;
  double sW = sH - pow2(CoupSUSY::mWpole);
  double d  = pow2(sW) + pow2(CoupSUSY::mWpole * CoupSUSY::wWpole);
  propW     = complex( sW / d, CoupSUSY::mWpole * CoupSUSY::wWpole / d);

}

//*********

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence. 

double Sigma2qqbar2charchi0::sigmaHat() {

  // Only allow particle-antiparticle incoming states
  if (id1*id2 >= 0) {
    return 0.0;    
  }
  
  // Only allow incoming states with sum(charge) = final state
  if (abs(id1) % 2 == abs(id2) % 2) return 0.0;
  int isPos  = (id3chi > 0 ? 1 : 0);
  if (id1 < 0 && id1 > -10 && abs(id1) % 2 == 1-isPos ) return 0.0;
  else if (id1 > 0 && id1 < 10 && abs(id1) % 2 == isPos ) return 0.0;

  // Flavour-dependent kinematics-dependent couplings.
  int idAbs1    = abs(id1);  
  complex QuLL(0.0),QtLL(0.0),QuRR(0.0),QtRR(0.0);
  complex QuLR(0.0),QtLR(0.0),QuRL(0.0),QtRL(0.0);
  
  // Calculate everything from udbar -> ~chi+ ~chi0 template process
  int    id1tmp = id1;
  int    id2tmp = id2;
  double uHtmp  = uH;
  double tHtmp  = tH;
  double uitmp  = ui;
  double ujtmp  = uj;
  double titmp  = ti;
  double tjtmp  = tj;

  // u dbar , ubar d : do nothing
  if (idAbs1 % 2 == 0) {
  }

  // dbar u , d ubar : swap 1<->2 and t<->u
  else if (idAbs1 % 2 == 1) {
    id1tmp = id2;
    id2tmp = id1;
    uHtmp  = tH;
    tHtmp  = uH;
    uitmp  = ti;
    ujtmp  = tj;
    titmp  = ui;
    tjtmp  = uj;
  }

  // Generation indices
  int iGu = abs(id1tmp)/2;
  int iGd = (abs(id2tmp)+1)/2;

  // s-channel W contribution
  QuLL = conj(CoupSUSY::LudW[iGu][iGd]) 
    * conj(CoupSUSY::OL[id4chi][abs(id3chi)])
    * propW / 4.0 / CoupSUSY::cosW;
  QtLL = conj(CoupSUSY::LudW[iGu][iGd]) 
    * conj(CoupSUSY::OR[id4chi][abs(id3chi)])
    * propW / 4.0 / CoupSUSY::cosW;

  // Add t-channel squark flavour sums to QmXY couplings
  for (int jsq=1; jsq<=6; jsq++) {    

    int idsu=((jsq+2)/3)*1000000 + 2*((jsq-1) % 3) + 2;
    int idsd=((jsq+2)/3)*1000000 + 2*((jsq-1) % 3) + 1;
    double msd2 = pow(ParticleDataTable::m0(idsd),2);
    double msu2 = pow(ParticleDataTable::m0(idsu),2);
    double tsq  = tHtmp - msd2;
    double usq  = uHtmp - msu2;

    QuLL += conj(CoupSUSY::LsuuX[jsq][iGu][id4chi])
      *CoupSUSY::LsudX[jsq][iGd][abs(id3chi)]/usq;
    QuRR += conj(CoupSUSY::RsuuX[jsq][iGu][id4chi])
      *CoupSUSY::RsudX[jsq][iGd][abs(id3chi)]/usq;
    QuLR += conj(CoupSUSY::LsuuX[jsq][iGu][id4chi])
      *CoupSUSY::RsudX[jsq][iGd][abs(id3chi)]/usq;
    QuRL += conj(CoupSUSY::RsuuX[jsq][iGu][id4chi])
      *CoupSUSY::LsudX[jsq][iGd][abs(id3chi)]/usq;
    QtLL -= conj(CoupSUSY::LsduX[jsq][iGu][abs(id3chi)])
      *CoupSUSY::LsddX[jsq][iGd][id4chi]/tsq;
    QtRR -= conj(CoupSUSY::RsduX[jsq][iGu][abs(id3chi)])
      *CoupSUSY::RsddX[jsq][iGd][id4chi]/tsq;
    QtLR -= conj(CoupSUSY::LsduX[jsq][iGu][abs(id3chi)])
      *CoupSUSY::RsddX[jsq][iGd][id4chi]/tsq;
    QtRL -= conj(CoupSUSY::RsduX[jsq][iGu][abs(id3chi)])
      *CoupSUSY::LsddX[jsq][iGd][id4chi]/tsq;
    
  }

  // Compute matrix element weight
  double weight = 0;

  // Average over separate helicity contributions
  // (if swapped, swap ha, hb if computing polarized cross sections)
  // LL (ha = -1, hb = +1) (divided by 4 for average)            
  weight += norm(QuLL) * uitmp * ujtmp + norm(QtLL) * titmp * tjtmp
    + 2 * real(conj(QuLL) * QtLL) * m3 * m4 * sH;
  // RR (ha =  1, hb = -1) (divided by 4 for average)        
  weight += norm(QtRR) * titmp * tjtmp + norm(QuRR) * uitmp * ujtmp  
    + 2 * real(conj(QuRR) * QtRR) * m3 * m4 * sH;
  // RL (ha =  1, hb =  1) (divided by 4 for average)        
  weight += norm(QuRL) * uitmp * ujtmp + norm(QtRL) * titmp * tjtmp
    - real(conj(QuRL) * QtRL) * (uHtmp * tHtmp - s3 * s4);
  // LR (ha = -1, hb = -1) (divided by 4 for average)        
  weight += norm(QuLR) * uitmp * ujtmp + norm(QtLR) * titmp * tjtmp
    - real(conj(QuLR) * QtLR) * (uHtmp * tHtmp - s3 * s4);

  // Cross section, including colour factor.
  double sigma = sigma0 * weight;
  if (idAbs1 < 9) sigma /= 3.; 

  // Answer.
  return sigma;    

}

//**************************************************************************

// Sigma2qqbar2charchar
// Cross section for gaugino pair production: chargino-chargino

//*********

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour. 

void Sigma2qqbar2charchar::sigmaKin() {

  // Common flavour-independent factor.
  sigma0 = 4 * (M_PI / sH2 / pow2(CoupSUSY::sin2W)) * pow2(alpEM) ; 

  // Auxiliary factors for use below
  ui       = uH - s3;
  uj       = uH - s4;
  ti       = tH - s3;
  tj       = tH - s4;
  double sV= sH - pow2(CoupSUSY::mZpole);
  double d = pow2(sV) + pow2(CoupSUSY::mZpole * CoupSUSY::wZpole);
  propZ    = complex( sV / d, CoupSUSY::mZpole * CoupSUSY::wZpole / d);

}

//*********

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence. 

double Sigma2qqbar2charchar::sigmaHat() { 

  // Only allow quark-antiquark incoming states
  if (id1*id2 >= 0) {
    return 0.0;
  }
  
  // Only allow incoming states with sum(charge) = 0
  if ((id1+id2) % 2 != 0) {
    return 0.0;    
  }
  
  //if (id1 > 0 || id1==-1 || id1==-3 || id1==-5) return 0.0;
  //if (id1 < 0 || id1==1 || id1==3 || id1==5) return 0.0;
  
  // Flavour-dependent kinematics-dependent couplings.
  int idAbs1    = abs(id1);  
  int idAbs2    = abs(id2);  
  int i3        = abs(id3chi);
  int i4        = abs(id4chi);
  
  // Flavour-dependent kinematics-dependent couplings.
  complex QuLL(0.0),QtLL(0.0),QuRR(0.0),QtRR(0.0);
  complex QuLR(0.0),QtLR(0.0),QuRL(0.0),QtRL(0.0);

  // Add Z/gamma* for same-flavour in-quarks
  if (idAbs1 == idAbs2) {

    // s-channel Z ( L<->R for qbar q incoming )
    QuLL = CoupSUSY::LqqZ[idAbs1] * (id1 > 0 ? CoupSUSY::OLp[i3][i4]
				     : CoupSUSY::ORp[i3][i4]);
    QtLL = CoupSUSY::LqqZ[idAbs1] * (id1 > 0 ? CoupSUSY::ORp[i3][i4]
				     : CoupSUSY::OLp[i3][i4]);    
    QuRR = CoupSUSY::RqqZ[idAbs1] * (id1 > 0 ? CoupSUSY::ORp[i3][i4]
				     : CoupSUSY::OLp[i3][i4]);
    QtRR = CoupSUSY::RqqZ[idAbs1] * (id1 > 0 ? CoupSUSY::OLp[i3][i4]
				     : CoupSUSY::ORp[i3][i4]);
    QuLL *= propZ / 4.0 / (1.0-CoupSUSY::sin2W);
    QtLL *= propZ / 4.0 / (1.0-CoupSUSY::sin2W);
    QuRR *= propZ / 4.0 / (1.0-CoupSUSY::sin2W);
    QtRR *= propZ / 4.0 / (1.0-CoupSUSY::sin2W);  

    // s-channel gamma* (only for same-type charginos)
    if (i3 == i4) {

      // Charge of in-particles
      double q = 2.0/3.0;
      if (idAbs1 % 2 == 1) q = -1.0/3.0;      
      QuLL -= q * CoupSUSY::sin2W / 2.0 / sH;
      QuRR -= q * CoupSUSY::sin2W / 2.0 / sH;
      QtLL -= q * CoupSUSY::sin2W / 2.0 / sH;
      QtRR -= q * CoupSUSY::sin2W / 2.0 / sH;

    }	
  }

  // Add t- or u-channel squark flavour sums to QmXY couplings
  for (int ksq=1; ksq<=6; ksq++) {    

    // Positive-sign on side 1: t-channel diagrams only
    if (ParticleDataTable::chargeType(id1) > 0) {

      // u ubar -> chi+ chi- : add t-channel ~d
      if (id1 % 2 == 0) {
	int iG1    = (abs(id1)+1)/2;
	int iG2    = (abs(id2)+1)/2;
	int idsd   = ((ksq+2)/3)*1000000 + 2*((ksq-1) % 3) + 1;
	double msq = ParticleDataTable::m0(idsd);
	double tsq = tH - pow2(msq);
	QtLL -= conj(CoupSUSY::LsduX[ksq][iG1][i3])
	  *CoupSUSY::LsduX[ksq][iG2][i4]/tsq;
	QtLR -= conj(CoupSUSY::LsduX[ksq][iG1][i3])
	  *CoupSUSY::RsduX[ksq][iG2][i4]/tsq;
	QtRL -= conj(CoupSUSY::RsduX[ksq][iG1][i3])
	  *CoupSUSY::LsduX[ksq][iG2][i4]/tsq;
	QtRR -= conj(CoupSUSY::RsduX[ksq][iG1][i3])
	  *CoupSUSY::RsduX[ksq][iG2][i4]/tsq;
      }

      // dbar d -> chi+ chi- : add t-channel ~u
      // But with QtXY <-> -conj(QtYX) since qbar incoming
     else {
	int iG1    = (abs(id1)+1)/2;
	int iG2    = (abs(id2)+1)/2;
	int idsu   = ((ksq+2)/3)*1000000 + 2*((ksq-1) % 3) + 2;
	double msq = ParticleDataTable::m0(idsu);
	double tsq = tH - pow2(msq);
	QtLL += conj(CoupSUSY::LsudX[ksq][iG2][i4])
	  *CoupSUSY::LsudX[ksq][iG1][i3]/tsq;
	QtLR += conj(CoupSUSY::LsudX[ksq][iG2][i4])
	  *CoupSUSY::RsudX[ksq][iG1][i3]/tsq;
	QtRL += conj(CoupSUSY::RsudX[ksq][iG2][i4])
	  *CoupSUSY::LsudX[ksq][iG1][i3]/tsq;
	QtRR += conj(CoupSUSY::RsudX[ksq][iG2][i4])
	  *CoupSUSY::RsudX[ksq][iG1][i3]/tsq;
      }

    } 

    // Negative-sign on side 1: u-channel diagrams only
    else { 
    
      // ubar u -> chi+ chi- : add u-channel ~d
      // But with QuXY <-> -conj(QuYX) since qbar incoming
      if (id1 % 2 == 0) {
	int iG1    = abs(id1)/2;
	int iG2    = abs(id2)/2;
	int idsd   = ((ksq+2)/3)*1000000 + 2*((ksq-1) % 3) + 1;
	double msq = ParticleDataTable::m0(idsd);
	double usq = uH - pow2(msq);
	QuLL -= conj(CoupSUSY::LsduX[ksq][iG2][i3])
	  *CoupSUSY::LsduX[ksq][iG1][i4]/usq;
	QuLR -= conj(CoupSUSY::LsduX[ksq][iG2][i3])
	  *CoupSUSY::RsduX[ksq][iG1][i4]/usq;
	QuRL -= conj(CoupSUSY::RsduX[ksq][iG2][i3])
	  *CoupSUSY::LsduX[ksq][iG1][i4]/usq;
	QuRR -= conj(CoupSUSY::RsduX[ksq][iG2][i3])
	  *CoupSUSY::RsduX[ksq][iG1][i4]/usq;
      }

      // d dbar -> chi+ chi- : add u-channel ~u
      else {
	int iG1    = (abs(id1)+1)/2;
	int iG2    = (abs(id2)+1)/2;
	int idsu   = ((ksq+2)/3)*1000000 + 2*((ksq-1) % 3) + 2;
	double msq = ParticleDataTable::m0(idsu);
	double usq = uH - pow2(msq);
	QuLL += conj(CoupSUSY::LsudX[ksq][iG1][i4])
	  *CoupSUSY::LsudX[ksq][iG2][i3]/usq;
	QuLR += conj(CoupSUSY::LsudX[ksq][iG1][i4])
	  *CoupSUSY::RsudX[ksq][iG2][i3]/usq;
	QuRL += conj(CoupSUSY::RsudX[ksq][iG1][i4])
	  *CoupSUSY::LsudX[ksq][iG2][i3]/usq;
	QuRR += conj(CoupSUSY::RsudX[ksq][iG1][i4])
	  *CoupSUSY::RsudX[ksq][iG2][i3]/usq;
      }

    }
  }

  // Compute matrix element weight
  double weight = 0;

  // Average over separate helicity contributions
  // LL (ha = -1, hb = +1) (divided by 4 for average)            
  weight += norm(QuLL) * ui * uj + norm(QtLL) * ti * tj
    + 2 * real(conj(QuLL) * QtLL) * m3 * m4 * sH;
  // RR (ha =  1, hb = -1) (divided by 4 for average)        
  weight += norm(QtRR) * ti * tj + norm(QuRR) * ui * uj  
    + 2 * real(conj(QuRR) * QtRR) * m3 * m4 * sH;
  // RL (ha =  1, hb =  1) (divided by 4 for average)        
  weight += norm(QuRL) * ui * uj + norm(QtRL) * ti * tj
    - real(conj(QuRL) * QtRL) * (uH * tH - s3 * s4);
  // LR (ha = -1, hb = -1) (divided by 4 for average)        
  weight += norm(QuLR) * ui * uj + norm(QtLR) * ti * tj
    - real(conj(QuLR) * QtLR) * (uH * tH - s3 * s4);

  // Cross section, including colour factor.
  double sigma = sigma0 * weight;
  if (idAbs1 < 9) sigma /= 3.; 

  // Answer.
  return sigma;    

}

//**************************************************************************

// Sigma2qgchi0squark 
// Cross section for gaugino-squark production: neutralino-squark

//*********

// Initialize process. 
  
void Sigma2qg2chi0squark::initProc() {

  // First make sure CoupSUSY is initialized
  if (not CoupSUSY::isInit) CoupSUSY::initStatic(slhaPtr);

  // Construct name of process. 
  if (id4 % 2 == 0) {
    nameSave = "q g -> " + ParticleDataTable::name(id3) + " " 
      + ParticleDataTable::name(id4) + " + c.c. (q=u,c)";
  } 
  else {
    nameSave = "q g -> " + ParticleDataTable::name(id3) + " " 
      + ParticleDataTable::name(id4) + " + c.c. (q=d,s,b)";
  }
}

//*********

// Evaluate d(sigmaHat)/d(tHat), part independent of incoming flavour. 

void Sigma2qg2chi0squark::sigmaKin() {

  // Common flavour-independent factor.
  // tmp: alphaS = 0.1 for counter-checks
  sigma0 = M_PI / sH2 / CoupSUSY::sin2W * alpEM * 0.1;//alpS;

  // Auxiliary factors for use below
  ui       = uH - s3;
  uj       = uH - s4;
  ti       = tH - s3;
  tj       = tH - s4;

}

//*********

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence. 

double Sigma2qg2chi0squark::sigmaHat() {

  // Antiquark -> antisquark
  int idq = id1;
  if (id1 == 21 || id1 == 22) idq = id2;
  if (idq < 0) {
    id4 = -abs(id4);
  } else {
    id4 = abs(id4);
  }

  // tmp: only allow incoming quarks on side 1
  //  if (id1 < 0 || id1 == 21) return 0.0;

  // Generation index
  int iGq = (abs(idq)+1)/2;

  // Only accept u(bar) -> ~u(bar) and d(bar) -> ~d(bar)
  if (ParticleDataTable::chargeType(idq) != ParticleDataTable::chargeType(id4))
    return 0.0;
  
  // Couplings
  complex LsqqX, RsqqX;
  if (idq % 2 == 0) {
    LsqqX = CoupSUSY::LsuuX[id4sq][iGq][id3chi];
    RsqqX = CoupSUSY::RsuuX[id4sq][iGq][id3chi];
  }
  else { 
    LsqqX = CoupSUSY::LsddX[id4sq][iGq][id3chi];
    RsqqX = CoupSUSY::RsddX[id4sq][iGq][id3chi];
  }  

  // Prefactors : swap u and t if gq instead of qg
  double fac1, fac2;
  if (idq == id1) {
    fac1 = -ui/sH + 2.0 * ( uH*tH - s4*s3 )/sH/tj;
    fac2 = ti/tj * ( (tH + s4)/tj + (ti - uj)/sH );
  } else {
    fac1 = -ti/sH + 2.0 * ( uH*tH - s4*s3 )/sH/uj;
    fac2 = ui/uj * ( (uH + s4)/uj + (ui - tj)/sH );
  }

  // Compute matrix element weight
  double weight = 0.0;

  // Average over separate helicity contributions
  // (for qbar g : ha -> -ha )
  // LL (ha = -1, hb = +1) (divided by 4 for average)            
  weight += fac2 * norm(LsqqX) / 2.0;
  // RR (ha =  1, hb = -1) (divided by 4 for average)        
  weight += fac2 * norm(RsqqX) / 2.0;
  // RL (ha =  1, hb =  1) (divided by 4 for average)        
  weight += fac2 * norm(RsqqX) / 2.0 + fac1 * norm(RsqqX);
  // LR (ha = -1, hb = -1) (divided by 4 for average)        
  weight += fac2 * norm(LsqqX) / 2.0 + fac1 * norm(LsqqX);

  double sigma = sigma0 * weight;
  if (abs(idq) < 9) sigma /= 3.;

  // Answer.
  return sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2qg2chi0squark::setIdColAcol() {

  // Set flavours.
  setId( id1, id2, id3, (id1*id2 > 0 ? abs(id4) : -abs(id4)));

  // Colour flow topology. Swap if first is gluon, or when antiquark.
  if (id1 != 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else setColAcol( 1, 2, 2, 0, 0, 0, 1, 0);
  if (id1*id2 < 0) swapColAcol();

}

//**************************************************************************

// Sigma2qg2charsquark
// Cross section for gaugino-squark production: chargino-squark

//*********

// Initialize process. 
  
void Sigma2qg2charsquark::initProc() {

  // First make sure CoupSUSY is initialized
  if (not CoupSUSY::isInit) CoupSUSY::initStatic(slhaPtr);

  // Construct name of process. 
  if (id4 % 2 == 0) {
    nameSave = "q g -> " + ParticleDataTable::name(id3) + " " 
      + ParticleDataTable::name(id4) + " + c.c. (q=d,s,b)";
  } 
  else {
    nameSave = "q g -> " + ParticleDataTable::name(id3) + " " 
      + ParticleDataTable::name(id4) + " + c.c. (q=u,c)";
  }
}

//*********

// Evaluate d(sigmaHat)/d(tHat), including incoming flavour dependence. 

double Sigma2qg2charsquark::sigmaHat() {

  // Antiquark -> antisquark
  int idq = id1;
  if (id1 == 21) idq = id2;
  if (idq > 0) {
    id3 = id3Sav;
    id4 = id4Sav;
  } else {
    id3 = -id3Sav;
    id4 = -id4Sav;
  }

  // Only accept u(bar) -> ~d(bar) and d(bar) -> ~u(bar)
  if (ParticleDataTable::chargeType(idq) == ParticleDataTable::chargeType(id4))
    return 0.0;
  
  // Generation index
  int iGq = (abs(idq)+1)/2;

  // Couplings
  complex LsqqX, RsqqX;
  if (idq % 2 == 0) {
    LsqqX = CoupSUSY::LsduX[id4sq][iGq][id3chi];
    RsqqX = CoupSUSY::RsduX[id4sq][iGq][id3chi];
  }
  else { 
    LsqqX = CoupSUSY::LsduX[id4sq][iGq][id3chi];
    RsqqX = CoupSUSY::RsduX[id4sq][iGq][id3chi];
  }  

  // Prefactors : swap u and t if gq instead of qg
  double fac1, fac2;
  if (idq == id1) {
    fac1 = -ui/sH + 2.0 * ( uH*tH - s4*s3 )/sH/tj;
    fac2 = ti/tj * ( (tH + s4)/tj + (ti - uj)/sH );
  } else {
    fac1 = -ti/sH + 2.0 * ( uH*tH - s4*s3 )/sH/uj;
    fac2 = ui/uj * ( (uH + s4)/uj + (ui - tj)/sH );
  }

  // Compute matrix element weight
  double weight = 0.0;

  // Average over separate helicity contributions
  // (a, b refers to qg configuration)
  // LL (ha = -1, hb = +1) (divided by 4 for average)            
  weight += fac2 * norm(LsqqX) / 2.0;
  // RR (ha =  1, hb = -1) (divided by 4 for average)        
  weight += fac2 * norm(RsqqX) / 2.0;
  // RL (ha =  1, hb =  1) (divided by 4 for average)        
  weight += fac2 * norm(RsqqX) / 2.0 + fac1 * norm(RsqqX);
  // LR (ha = -1, hb = -1) (divided by 4 for average)        
  weight += fac2 * norm(LsqqX) / 2.0 + fac1 * norm(LsqqX);

  double sigma = sigma0 * weight;
  if (abs(idq) < 9) sigma /= 3.;

  // Answer.
  return sigma;    

}

//*********

// Select identity, colour and anticolour.

void Sigma2qg2charsquark::setIdColAcol() {

  // Set flavours.
  if (id1 > 0 && id2 > 0) {
    setId( id1, id2, id3Sav, id4Sav);
  } else {
    setId( id1, id2,-id3Sav,-id4Sav);
  }

  // Colour flow topology. Swap if first is gluon, or when antiquark.
  if (id1 != 21) setColAcol( 1, 0, 2, 1, 0, 0, 2, 0);
  else setColAcol( 1, 2, 2, 0, 0, 0, 1, 0);
  if (id1 < 0 || id2 < 0) swapColAcol();

}

//**************************************************************************

} // end namespace Pythia8

