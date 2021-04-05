//==========================================================================
// This file has been automatically generated for Pythia 8 by
// MadGraph5_aMC@NLO v. 2.7.2, 2020-03-17
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <iostream> 
#include "Parameters_heft.h"
#include "Pythia8/PythiaStdlib.h"

using namespace Pythia8; 

// Initialize static instance
Parameters_heft * Parameters_heft::instance = 0; 

// Function to get static instance - only one instance per program
Parameters_heft * Parameters_heft::getInstance()
{
  if (instance == 0)
    instance = new Parameters_heft(); 

  return instance; 
}

void Parameters_heft::setIndependentParameters(ParticleData * & pd, CoupSM * &
    csm, SusyLesHouches * & slhaPtr)
{
  mdl_WH1 = pd->mWidth(9000006); 
  mdl_WH = pd->mWidth(25); 
  mdl_WW = pd->mWidth(24); 
  mdl_WZ = pd->mWidth(23); 
  mdl_WT = pd->mWidth(6); 
  mdl_MP = pd->m0(9000006); 
  mdl_MTA = pd->m0(15); 
  mdl_MH = pd->m0(25); 
  mdl_MZ = pd->m0(23); 
  mdl_MB = pd->m0(5); 
  mdl_MT = pd->m0(6); 
  mdl_ymtau = pd->mRun(15, pd->m0(24)); 
  mdl_ymt = pd->mRun(6, pd->m0(24)); 
  mdl_ymb = pd->mRun(5, pd->m0(24)); 
  mdl_Gf = M_PI * csm->alphaEM(((pd->m0(23)) * (pd->m0(23)))) * ((pd->m0(23)) *
      (pd->m0(23)))/(sqrt(2.) * ((pd->m0(24)) * (pd->m0(24))) * (((pd->m0(23))
      * (pd->m0(23))) - ((pd->m0(24)) * (pd->m0(24)))));
  aEWM1 = 1./csm->alphaEM(((pd->m0(23)) * (pd->m0(23)))); 
  mdl_conjg__CKM3x3 = 1.; 
  ZERO = 0.; 
  mdl_complexi = std::complex<double> (0., 1.); 
  mdl_MZ__exp__2 = ((mdl_MZ) * (mdl_MZ)); 
  mdl_MZ__exp__4 = ((mdl_MZ) * (mdl_MZ) * (mdl_MZ) * (mdl_MZ)); 
  mdl_sqrt__2 = sqrt(2.); 
  mdl_MH__exp__4 = ((mdl_MH) * (mdl_MH) * (mdl_MH) * (mdl_MH)); 
  mdl_MT__exp__4 = ((mdl_MT) * (mdl_MT) * (mdl_MT) * (mdl_MT)); 
  mdl_MH__exp__2 = ((mdl_MH) * (mdl_MH)); 
  mdl_MT__exp__2 = ((mdl_MT) * (mdl_MT)); 
  mdl_MH__exp__12 = pow(mdl_MH, 12.); 
  mdl_MH__exp__10 = pow(mdl_MH, 10.); 
  mdl_MH__exp__8 = pow(mdl_MH, 8.); 
  mdl_MH__exp__6 = pow(mdl_MH, 6.); 
  mdl_MT__exp__6 = pow(mdl_MT, 6.); 
  mdl_aEW = 1./aEWM1; 
  mdl_MW = sqrt(mdl_MZ__exp__2/2. + sqrt(mdl_MZ__exp__4/4. - (mdl_aEW * M_PI *
      mdl_MZ__exp__2)/(mdl_Gf * mdl_sqrt__2)));
  mdl_sqrt__aEW = sqrt(mdl_aEW); 
  mdl_ee = 2. * mdl_sqrt__aEW * sqrt(M_PI); 
  mdl_MW__exp__2 = ((mdl_MW) * (mdl_MW)); 
  mdl_sw2 = 1. - mdl_MW__exp__2/mdl_MZ__exp__2; 
  mdl_cw = sqrt(1. - mdl_sw2); 
  mdl_sqrt__sw2 = sqrt(mdl_sw2); 
  mdl_sw = mdl_sqrt__sw2; 
  mdl_g1 = mdl_ee/mdl_cw; 
  mdl_gw = mdl_ee/mdl_sw; 
  mdl_v = (2. * mdl_MW * mdl_sw)/mdl_ee; 
  mdl_ee__exp__2 = ((mdl_ee) * (mdl_ee)); 
  mdl_MW__exp__12 = pow(mdl_MW, 12.); 
  mdl_MW__exp__10 = pow(mdl_MW, 10.); 
  mdl_MW__exp__8 = pow(mdl_MW, 8.); 
  mdl_MW__exp__6 = pow(mdl_MW, 6.); 
  mdl_MW__exp__4 = ((mdl_MW) * (mdl_MW) * (mdl_MW) * (mdl_MW)); 
  mdl_AH = (47. * mdl_ee__exp__2 * (1. - (2. * mdl_MH__exp__4)/(987. *
      mdl_MT__exp__4) - (14. * mdl_MH__exp__2)/(705. * mdl_MT__exp__2) + (213.
      * mdl_MH__exp__12)/(2.634632e7 * mdl_MW__exp__12) + (5. *
      mdl_MH__exp__10)/(119756. * mdl_MW__exp__10) + (41. *
      mdl_MH__exp__8)/(180950. * mdl_MW__exp__8) + (87. *
      mdl_MH__exp__6)/(65800. * mdl_MW__exp__6) + (57. * mdl_MH__exp__4)/(6580.
      * mdl_MW__exp__4) + (33. * mdl_MH__exp__2)/(470. * mdl_MW__exp__2)))/(72.
      * ((M_PI) * (M_PI)) * mdl_v);
  mdl_v__exp__2 = ((mdl_v) * (mdl_v)); 
  mdl_lam = mdl_MH__exp__2/(2. * mdl_v__exp__2); 
  mdl_yb = (mdl_ymb * mdl_sqrt__2)/mdl_v; 
  mdl_yt = (mdl_ymt * mdl_sqrt__2)/mdl_v; 
  mdl_ytau = (mdl_ymtau * mdl_sqrt__2)/mdl_v; 
  mdl_muH = sqrt(mdl_lam * mdl_v__exp__2); 
  mdl_gw__exp__2 = ((mdl_gw) * (mdl_gw)); 
  mdl_cw__exp__2 = ((mdl_cw) * (mdl_cw)); 
  mdl_sw__exp__2 = ((mdl_sw) * (mdl_sw)); 
}
void Parameters_heft::setIndependentCouplings()
{
  GC_1 = -(mdl_AH * mdl_complexi); 
  GC_2 = -(mdl_ee * mdl_complexi)/3.; 
  GC_3 = (2. * mdl_ee * mdl_complexi)/3.; 
  GC_5 = mdl_ee * mdl_complexi; 
  GC_6 = 2. * mdl_ee__exp__2 * mdl_complexi; 
  GC_7 = -mdl_ee__exp__2/(2. * mdl_cw); 
  GC_8 = (mdl_ee__exp__2 * mdl_complexi)/(2. * mdl_cw); 
  GC_18 = mdl_cw * mdl_complexi * mdl_gw; 
  GC_19 = -(mdl_complexi * mdl_gw__exp__2); 
  GC_20 = mdl_cw__exp__2 * mdl_complexi * mdl_gw__exp__2; 
  GC_21 = -2. * mdl_complexi * mdl_lam; 
  GC_22 = -4. * mdl_complexi * mdl_lam; 
  GC_23 = -6. * mdl_complexi * mdl_lam; 
  GC_24 = -(mdl_ee * mdl_MW); 
  GC_26 = (mdl_ee__exp__2 * mdl_complexi)/(2. * mdl_sw__exp__2); 
  GC_28 = (mdl_ee * mdl_complexi)/(2. * mdl_sw); 
  GC_29 = mdl_ee/(2. * mdl_sw); 
  GC_40 = -(mdl_cw * mdl_ee * mdl_complexi)/(2. * mdl_sw); 
  GC_42 = -((mdl_cw * mdl_ee * mdl_complexi)/mdl_sw); 
  GC_44 = -mdl_ee__exp__2/(2. * mdl_sw); 
  GC_45 = -(mdl_ee__exp__2 * mdl_complexi)/(2. * mdl_sw); 
  GC_48 = -(mdl_ee * mdl_complexi * mdl_MW)/(2. * mdl_sw); 
  GC_49 = (mdl_ee * mdl_MW)/(2. * mdl_sw); 
  GC_51 = (mdl_ee * mdl_MZ)/(2. * mdl_sw); 
  GC_52 = -(mdl_ee * mdl_complexi * mdl_MZ)/(2. * mdl_cw * mdl_sw); 
  GC_53 = -(mdl_ee * mdl_complexi * mdl_sw)/(6. * mdl_cw); 
  GC_54 = (mdl_ee * mdl_complexi * mdl_sw)/(2. * mdl_cw); 
  GC_55 = mdl_complexi * mdl_gw * mdl_sw; 
  GC_56 = -2. * mdl_cw * mdl_complexi * mdl_gw__exp__2 * mdl_sw; 
  GC_57 = mdl_complexi * mdl_gw__exp__2 * mdl_sw__exp__2; 
  GC_58 = -(mdl_cw * mdl_ee * mdl_complexi)/(2. * mdl_sw) + (mdl_ee *
      mdl_complexi * mdl_sw)/(2. * mdl_cw);
  GC_59 = (mdl_cw * mdl_ee * mdl_complexi)/(2. * mdl_sw) + (mdl_ee *
      mdl_complexi * mdl_sw)/(2. * mdl_cw);
  GC_60 = (mdl_cw * mdl_ee)/(2. * mdl_sw) + (mdl_ee * mdl_sw)/(2. * mdl_cw); 
  GC_61 = (mdl_cw * mdl_ee__exp__2 * mdl_complexi)/mdl_sw - (mdl_ee__exp__2 *
      mdl_complexi * mdl_sw)/mdl_cw;
  GC_62 = (mdl_cw * mdl_ee * mdl_MW)/(2. * mdl_sw) - (mdl_ee * mdl_MW *
      mdl_sw)/(2. * mdl_cw);
  GC_64 = -(mdl_ee__exp__2 * mdl_complexi) + (mdl_cw__exp__2 * mdl_ee__exp__2 *
      mdl_complexi)/(2. * mdl_sw__exp__2) + (mdl_ee__exp__2 * mdl_complexi *
      mdl_sw__exp__2)/(2. * mdl_cw__exp__2);
  GC_65 = mdl_ee__exp__2 * mdl_complexi + (mdl_cw__exp__2 * mdl_ee__exp__2 *
      mdl_complexi)/(2. * mdl_sw__exp__2) + (mdl_ee__exp__2 * mdl_complexi *
      mdl_sw__exp__2)/(2. * mdl_cw__exp__2);
  GC_66 = -(mdl_ee__exp__2 * mdl_v)/(2. * mdl_cw); 
  GC_68 = -2. * mdl_complexi * mdl_lam * mdl_v; 
  GC_69 = -6. * mdl_complexi * mdl_lam * mdl_v; 
  GC_70 = (mdl_ee__exp__2 * mdl_complexi * mdl_v)/(2. * mdl_sw__exp__2); 
  GC_73 = mdl_ee__exp__2 * mdl_complexi * mdl_v + (mdl_cw__exp__2 *
      mdl_ee__exp__2 * mdl_complexi * mdl_v)/(2. * mdl_sw__exp__2) +
      (mdl_ee__exp__2 * mdl_complexi * mdl_sw__exp__2 * mdl_v)/(2. *
      mdl_cw__exp__2);
  GC_74 = -((mdl_complexi * mdl_yb)/mdl_sqrt__2); 
  GC_75 = mdl_yb/mdl_sqrt__2; 
  GC_79 = -(mdl_yt/mdl_sqrt__2); 
  GC_80 = -((mdl_complexi * mdl_yt)/mdl_sqrt__2); 
  GC_84 = -mdl_ytau; 
  GC_86 = -((mdl_complexi * mdl_ytau)/mdl_sqrt__2); 
  GC_87 = mdl_ytau/mdl_sqrt__2; 
  GC_100 = (mdl_ee * mdl_complexi * mdl_conjg__CKM3x3)/(mdl_sw * mdl_sqrt__2); 
  GC_101 = mdl_yb * mdl_conjg__CKM3x3; 
  GC_102 = -(mdl_yt * mdl_conjg__CKM3x3); 
}

void Parameters_heft::setDependentParameters(ParticleData * & pd, CoupSM * &
    csm, SusyLesHouches * & slhaPtr, double alpS)
{
  aS = alpS; 
  mdl_sqrt__aS = sqrt(aS); 
  G = 2. * mdl_sqrt__aS * sqrt(M_PI); 
  mdl_G__exp__2 = ((G) * (G)); 
  mdl_GH = -(mdl_G__exp__2 * (1. + (13. * mdl_MH__exp__6)/(16800. *
      mdl_MT__exp__6) + mdl_MH__exp__4/(168. * mdl_MT__exp__4) + (7. *
      mdl_MH__exp__2)/(120. * mdl_MT__exp__2)))/(12. * ((M_PI) * (M_PI)) *
      mdl_v);
  mdl_Gphi = -(mdl_G__exp__2 * (1. + mdl_MH__exp__6/(560. * mdl_MT__exp__6) +
      mdl_MH__exp__4/(90. * mdl_MT__exp__4) + mdl_MH__exp__2/(12. *
      mdl_MT__exp__2)))/(8. * ((M_PI) * (M_PI)) * mdl_v);
}


void Parameters_heft::setDependentCouplings()
{
  GC_17 = -(G * mdl_Gphi); 
  GC_16 = (mdl_complexi * mdl_Gphi)/8.; 
  GC_15 = mdl_complexi * mdl_G__exp__2 * mdl_GH; 
  GC_14 = -(G * mdl_GH); 
  GC_13 = -(mdl_complexi * mdl_GH); 
  GC_12 = mdl_complexi * mdl_G__exp__2; 
  GC_11 = mdl_complexi * G; 
  GC_10 = -G; 
}

// Routines for printing out parameters
void Parameters_heft::printIndependentParameters()
{
  cout <<  "heft model parameters independent of event kinematics:" << endl; 
  cout << setw(20) <<  "mdl_WH1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WH1 << endl;
  cout << setw(20) <<  "mdl_WH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WH << endl;
  cout << setw(20) <<  "mdl_WW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WW << endl;
  cout << setw(20) <<  "mdl_WZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WZ << endl;
  cout << setw(20) <<  "mdl_WT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WT << endl;
  cout << setw(20) <<  "mdl_MP " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MP << endl;
  cout << setw(20) <<  "mdl_MTA " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MTA << endl;
  cout << setw(20) <<  "mdl_MH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MH << endl;
  cout << setw(20) <<  "mdl_MZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MZ << endl;
  cout << setw(20) <<  "mdl_MB " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MB << endl;
  cout << setw(20) <<  "mdl_MT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MT << endl;
  cout << setw(20) <<  "mdl_ymtau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymtau << endl;
  cout << setw(20) <<  "mdl_ymt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymt << endl;
  cout << setw(20) <<  "mdl_ymb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymb << endl;
  cout << setw(20) <<  "mdl_Gf " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_Gf << endl;
  cout << setw(20) <<  "aEWM1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aEWM1 << endl;
  cout << setw(20) <<  "mdl_conjg__CKM3x3 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM3x3 << endl;
  cout << setw(20) <<  "ZERO " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << ZERO << endl;
  cout << setw(20) <<  "mdl_complexi " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_complexi << endl;
  cout << setw(20) <<  "mdl_MZ__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MZ__exp__2 << endl;
  cout << setw(20) <<  "mdl_MZ__exp__4 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MZ__exp__4 << endl;
  cout << setw(20) <<  "mdl_sqrt__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_sqrt__2 << endl;
  cout << setw(20) <<  "mdl_MH__exp__4 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MH__exp__4 << endl;
  cout << setw(20) <<  "mdl_MT__exp__4 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MT__exp__4 << endl;
  cout << setw(20) <<  "mdl_MH__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MH__exp__2 << endl;
  cout << setw(20) <<  "mdl_MT__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MT__exp__2 << endl;
  cout << setw(20) <<  "mdl_MH__exp__12 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MH__exp__12 << endl;
  cout << setw(20) <<  "mdl_MH__exp__10 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MH__exp__10 << endl;
  cout << setw(20) <<  "mdl_MH__exp__8 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MH__exp__8 << endl;
  cout << setw(20) <<  "mdl_MH__exp__6 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MH__exp__6 << endl;
  cout << setw(20) <<  "mdl_MT__exp__6 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MT__exp__6 << endl;
  cout << setw(20) <<  "mdl_aEW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_aEW << endl;
  cout << setw(20) <<  "mdl_MW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MW << endl;
  cout << setw(20) <<  "mdl_sqrt__aEW " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sqrt__aEW << endl;
  cout << setw(20) <<  "mdl_ee " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ee << endl;
  cout << setw(20) <<  "mdl_MW__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MW__exp__2 << endl;
  cout << setw(20) <<  "mdl_sw2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_sw2 << endl;
  cout << setw(20) <<  "mdl_cw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_cw << endl;
  cout << setw(20) <<  "mdl_sqrt__sw2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sqrt__sw2 << endl;
  cout << setw(20) <<  "mdl_sw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_sw << endl;
  cout << setw(20) <<  "mdl_g1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_g1 << endl;
  cout << setw(20) <<  "mdl_gw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_gw << endl;
  cout << setw(20) <<  "mdl_v " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_v << endl;
  cout << setw(20) <<  "mdl_ee__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_ee__exp__2 << endl;
  cout << setw(20) <<  "mdl_MW__exp__12 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MW__exp__12 << endl;
  cout << setw(20) <<  "mdl_MW__exp__10 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MW__exp__10 << endl;
  cout << setw(20) <<  "mdl_MW__exp__8 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MW__exp__8 << endl;
  cout << setw(20) <<  "mdl_MW__exp__6 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MW__exp__6 << endl;
  cout << setw(20) <<  "mdl_MW__exp__4 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MW__exp__4 << endl;
  cout << setw(20) <<  "mdl_AH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_AH << endl;
  cout << setw(20) <<  "mdl_v__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_v__exp__2 << endl;
  cout << setw(20) <<  "mdl_lam " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_lam << endl;
  cout << setw(20) <<  "mdl_yb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yb << endl;
  cout << setw(20) <<  "mdl_yt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yt << endl;
  cout << setw(20) <<  "mdl_ytau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ytau << endl;
  cout << setw(20) <<  "mdl_muH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_muH << endl;
  cout << setw(20) <<  "mdl_gw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_gw__exp__2 << endl;
  cout << setw(20) <<  "mdl_cw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_cw__exp__2 << endl;
  cout << setw(20) <<  "mdl_sw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sw__exp__2 << endl;
}
void Parameters_heft::printIndependentCouplings()
{
  cout <<  "heft model couplings independent of event kinematics:" << endl; 
  cout << setw(20) <<  "GC_1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_1 << endl;
  cout << setw(20) <<  "GC_2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_2 << endl;
  cout << setw(20) <<  "GC_3 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_3 << endl;
  cout << setw(20) <<  "GC_5 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_5 << endl;
  cout << setw(20) <<  "GC_6 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_6 << endl;
  cout << setw(20) <<  "GC_7 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_7 << endl;
  cout << setw(20) <<  "GC_8 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_8 << endl;
  cout << setw(20) <<  "GC_18 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_18 << endl;
  cout << setw(20) <<  "GC_19 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_19 << endl;
  cout << setw(20) <<  "GC_20 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_20 << endl;
  cout << setw(20) <<  "GC_21 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_21 << endl;
  cout << setw(20) <<  "GC_22 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_22 << endl;
  cout << setw(20) <<  "GC_23 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_23 << endl;
  cout << setw(20) <<  "GC_24 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_24 << endl;
  cout << setw(20) <<  "GC_26 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_26 << endl;
  cout << setw(20) <<  "GC_28 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_28 << endl;
  cout << setw(20) <<  "GC_29 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_29 << endl;
  cout << setw(20) <<  "GC_40 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_40 << endl;
  cout << setw(20) <<  "GC_42 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_42 << endl;
  cout << setw(20) <<  "GC_44 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_44 << endl;
  cout << setw(20) <<  "GC_45 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_45 << endl;
  cout << setw(20) <<  "GC_48 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_48 << endl;
  cout << setw(20) <<  "GC_49 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_49 << endl;
  cout << setw(20) <<  "GC_51 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_51 << endl;
  cout << setw(20) <<  "GC_52 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_52 << endl;
  cout << setw(20) <<  "GC_53 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_53 << endl;
  cout << setw(20) <<  "GC_54 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_54 << endl;
  cout << setw(20) <<  "GC_55 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_55 << endl;
  cout << setw(20) <<  "GC_56 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_56 << endl;
  cout << setw(20) <<  "GC_57 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_57 << endl;
  cout << setw(20) <<  "GC_58 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_58 << endl;
  cout << setw(20) <<  "GC_59 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_59 << endl;
  cout << setw(20) <<  "GC_60 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_60 << endl;
  cout << setw(20) <<  "GC_61 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_61 << endl;
  cout << setw(20) <<  "GC_62 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_62 << endl;
  cout << setw(20) <<  "GC_64 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_64 << endl;
  cout << setw(20) <<  "GC_65 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_65 << endl;
  cout << setw(20) <<  "GC_66 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_66 << endl;
  cout << setw(20) <<  "GC_68 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_68 << endl;
  cout << setw(20) <<  "GC_69 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_69 << endl;
  cout << setw(20) <<  "GC_70 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_70 << endl;
  cout << setw(20) <<  "GC_73 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_73 << endl;
  cout << setw(20) <<  "GC_74 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_74 << endl;
  cout << setw(20) <<  "GC_75 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_75 << endl;
  cout << setw(20) <<  "GC_79 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_79 << endl;
  cout << setw(20) <<  "GC_80 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_80 << endl;
  cout << setw(20) <<  "GC_84 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_84 << endl;
  cout << setw(20) <<  "GC_86 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_86 << endl;
  cout << setw(20) <<  "GC_87 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_87 << endl;
  cout << setw(20) <<  "GC_100 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_100 << endl;
  cout << setw(20) <<  "GC_101 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_101 << endl;
  cout << setw(20) <<  "GC_102 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_102 << endl;
}
void Parameters_heft::printDependentParameters()
{
  cout <<  "heft model parameters dependent on event kinematics:" << endl; 
  cout << setw(20) <<  "aS " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aS << endl;
  cout << setw(20) <<  "mdl_sqrt__aS " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_sqrt__aS << endl;
  cout << setw(20) <<  "G " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << G << endl;
  cout << setw(20) <<  "mdl_G__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_G__exp__2 << endl;
  cout << setw(20) <<  "mdl_GH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_GH << endl;
  cout << setw(20) <<  "mdl_Gphi " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_Gphi << endl;
}
void Parameters_heft::printDependentCouplings()
{
  cout <<  "heft model couplings dependent on event kinematics:" << endl; 
  cout << setw(20) <<  "GC_17 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_17 << endl;
  cout << setw(20) <<  "GC_16 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_16 << endl;
  cout << setw(20) <<  "GC_15 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_15 << endl;
  cout << setw(20) <<  "GC_14 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_14 << endl;
  cout << setw(20) <<  "GC_13 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_13 << endl;
  cout << setw(20) <<  "GC_12 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_12 << endl;
  cout << setw(20) <<  "GC_11 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_11 << endl;
  cout << setw(20) <<  "GC_10 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_10 << endl;
}


