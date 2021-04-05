//==========================================================================
// This file has been automatically generated for Pythia 8
// MadGraph5_aMC@NLO v. 2.7.2, 2020-03-17
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Pythia8_parameters_heft_H
#define Pythia8_parameters_heft_H

#include <complex> 

#include "Pythia8/ParticleData.h"
#include "Pythia8/StandardModel.h"
#include "Pythia8/SusyLesHouches.h"

using namespace std; 

using namespace Pythia8; 

class Parameters_heft
{
  public:

    static Parameters_heft * getInstance(); 

    // Model parameters independent of aS
    double mdl_WH1, mdl_WH, mdl_WW, mdl_WZ, mdl_WT, mdl_MP, mdl_MTA, mdl_MH,
        mdl_MZ, mdl_MB, mdl_MT, mdl_ymtau, mdl_ymt, mdl_ymb, mdl_Gf, aEWM1,
        mdl_conjg__CKM3x3, ZERO, mdl_MZ__exp__2, mdl_MZ__exp__4, mdl_sqrt__2,
        mdl_MH__exp__4, mdl_MT__exp__4, mdl_MH__exp__2, mdl_MT__exp__2,
        mdl_MH__exp__12, mdl_MH__exp__10, mdl_MH__exp__8, mdl_MH__exp__6,
        mdl_MT__exp__6, mdl_aEW, mdl_MW, mdl_sqrt__aEW, mdl_ee, mdl_MW__exp__2,
        mdl_sw2, mdl_cw, mdl_sqrt__sw2, mdl_sw, mdl_g1, mdl_gw, mdl_v,
        mdl_ee__exp__2, mdl_MW__exp__12, mdl_MW__exp__10, mdl_MW__exp__8,
        mdl_MW__exp__6, mdl_MW__exp__4, mdl_AH, mdl_v__exp__2, mdl_lam, mdl_yb,
        mdl_yt, mdl_ytau, mdl_muH, mdl_gw__exp__2, mdl_cw__exp__2,
        mdl_sw__exp__2;
    std::complex<double> mdl_complexi; 
    // Model parameters dependent on aS
    double aS, mdl_sqrt__aS, G, mdl_G__exp__2, mdl_GH, mdl_Gphi; 
    // Model couplings independent of aS
    std::complex<double> GC_1, GC_2, GC_3, GC_5, GC_6, GC_7, GC_8, GC_18,
        GC_19, GC_20, GC_21, GC_22, GC_23, GC_24, GC_26, GC_28, GC_29, GC_40,
        GC_42, GC_44, GC_45, GC_48, GC_49, GC_51, GC_52, GC_53, GC_54, GC_55,
        GC_56, GC_57, GC_58, GC_59, GC_60, GC_61, GC_62, GC_64, GC_65, GC_66,
        GC_68, GC_69, GC_70, GC_73, GC_74, GC_75, GC_79, GC_80, GC_84, GC_86,
        GC_87, GC_100, GC_101, GC_102;
    // Model couplings dependent on aS
    std::complex<double> GC_17, GC_16, GC_15, GC_14, GC_13, GC_12, GC_11,
        GC_10;

    // Set parameters that are unchanged during the run
    void setIndependentParameters(ParticleData * & pd, CoupSM * & csm,
        SusyLesHouches * & slhaPtr);
    // Set couplings that are unchanged during the run
    void setIndependentCouplings(); 
    // Set parameters that are changed event by event
    void setDependentParameters(ParticleData * & pd, CoupSM * & csm,
        SusyLesHouches * & slhaPtr, double alpS);
    // TMP: hardcoded bogus implementation with no arguments since this
    // is being called from within the matrix elements.
    void setDependentParameters() {}; 

    // Set couplings that are changed event by event
    void setDependentCouplings(); 

    // Print parameters that are unchanged during the run
    void printIndependentParameters(); 
    // Print couplings that are unchanged during the run
    void printIndependentCouplings(); 
    // Print parameters that are changed event by event
    void printDependentParameters(); 
    // Print couplings that are changed event by event
    void printDependentCouplings(); 

  private:
    static Parameters_heft * instance; 
}; 

#endif  // Pythia8_parameters_heft_H

