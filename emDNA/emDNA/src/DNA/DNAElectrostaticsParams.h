// DNAElectrostaticsParams header
// Juan Wei, Nicolas Clauvelin


// this header file defines all the constants and parameters related to the
// treatment of the electrostatics of DNA


#ifndef emDNA_DNAElectrostaticsParams_h
#define emDNA_DNAElectrostaticsParams_h


#include <emDNA_Includes.h>


namespace DNAElec {

    // physical constants related to electrostatics treatment
    // this assumes 2 charges per base pair with a 76% concentration (Manning
    // theory)
    // for a 10mL concentration the Debye length is 30.4 Angstroms
    const Real WaterPermittivity = Real(77.4); // for T=300 K
    const Real VacuumPermittivity = Real(8.854188e-22);
    const Real ChargeDNA = Real(2.0*0.48*1.6021773e-19);
    const Real SolutionConcentration = Real(0.10);
    const Real kappaDebye = Real(0.329)*std::sqrt(SolutionConcentration);
    const Real kbT = Real(1.380658e-23*300.0);
    const Real dhConstant =
    ChargeDNA*ChargeDNA/(4.0*M_PI*
                         WaterPermittivity*
                         VacuumPermittivity*
                         kbT);

    // neighbor exclusion range
    const Size ExclusionRange = 2;

};


#endif  // emDNA_DNAElectrostaticsParams_h
