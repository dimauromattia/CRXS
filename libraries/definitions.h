#ifndef DEFINITIONS_H
#define DEFINITIONS_H


#define dMass_proton         0.9382720813;     // PDG Review 2016
#define dMass_neutron        0.939565379;      // PDG Review 2016
#define dMass_deuteron       1.875612928;      // PDG Review 2016

// Integration type:
#define     TRAPEZSUM   0
#define     UPPERSUM    1
#define     LOWERSUM    2


#define     NO_OUT      0
#define     WARN_OUT    1
#define     ALL_OUT     2

#define     warnout                         "\033[9;31mWarning: "
#define     warnend                         "\033[0m"

#define     warning_momentum_integration            "\033[9;31mWarning: In the momoentum integration the error is more than 1%. \033[0m"
#define     warning_eta_integration                 "\033[9;31mWarning: In the eta integration the range might be too small. \033[0m"
#define     warning_parametrization_not_available   "\033[9;31mWarning: Your requested parametrisation is not available. Return -1.\033[0m"
#define     warning_sourceterm_integration          "\033[9;31mWarning: In the source term integration the error is more than 1%.\033[0m"

#ifdef DEBUG
#define varOut(x) std::cout << "  "; std::cout.width(70); std::cout << std::left << #x <<"    " << x << std::endl;
#define funOut(x) std::cout << "*******************************************************************************************" << std::endl; std::cout << "**  " << #x << std::endl; std::cout << "*******************************************************************************************" << std::endl;
#else
#define varOut(x)
#define funOut(x)
#endif

#define out(x) std::cout << "  "; std::cout.width(70); std::cout << std::left << #x <<"    " << x << std::endl;


namespace CRACS {
    static double fMass_proton         = 0.9382720813;     // PDG Review 2016
    static double fMass_neutron        = fMass_proton;     // 0.939565379;      // PDG Review 2016
    static double fMass_deuteron       = 2*fMass_proton;   // 1.875612928;      // PDG Review 2016
    static double fMass_helium3        = 3*fMass_proton;   // FIXME
    static double fMass_helium4        = 4.002602*0.9314940954;   // Wikipedia (4He)
    static double fMass_helium         = 4.002602*0.9314940954;
    
    
    
    enum spectrumType{
        RIGIDITY    = 1,
        EKINPERN    = 2,
        RATIO       = 3,
        RATIO_EKIN  = 4
    };
    
}


#endif
