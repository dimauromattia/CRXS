#include "math.h"

#include "clike.h"

double fMass_proton         = 0.9382720813;     // PDG Review 2016

#define C_array_to_double(NAM) double C##NAM = C_array[NAM];


//! Parametrization of the invariant antiproton production (pbar) crosssection from pp in CMF
/*!
 *  All cross sections are given in mbarn
 *  All enegies, momenta, and masses have unit GeV.
 *
 *  Taken from:     Winkler, M. W.; 2016;
 *                  Cosmic Ray Antiprotons at High Energies;
 *                  arXiv:1701.04866
 *
 *  \f[  E_{\bar{p}} \frac{ d^3 \sigma_{pp}^{(\bar{p})} }{d^3 p_{\bar{p}} }  (s, E_{\bar{p}}, p_{T,\bar{p}})  \f]
 *
 *  \param double  s         CM energy.
 *  \param doulbe  E_pbar    Energy of the produced antiproton in CMF
 *  \param doulbe  pT_pbar   Transverse momentum of the produced antiproton in CMF.
 *  \param doulbe* C_array   Ci=(1,...,16), parameters. C1=C_array[1], C2=C_array[2], ... (C_array[0] is not used)
 * */
double inv_pp_pbar_CM__Winkler(double s, double E_pbar_d, double pT_pbar, double* C_array, int len_C_array){
    
    C_array_to_double( 5);
    C_array_to_double( 6);
    C_array_to_double( 7);
    C_array_to_double( 8);
    C_array_to_double( 9);
    C_array_to_double(10);
    C_array_to_double(11);
    C_array_to_double(12);
    C_array_to_double(13);
    
    double E_pbar = fabs(E_pbar_d);
    if (s<16*fMass_proton*fMass_proton){
        return 0;
    }
    if ( pow(pT_pbar, 2.) > pow(E_pbar, 2.) - pow(fMass_proton, 2.) ){
        return 0;
    }
    double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
    double x_R          =   E_pbar/E_pbar_Max;
    if ( x_R > 1. )
    return  0.;
    
    double m_T = sqrt(  pT_pbar*pT_pbar  +  fMass_proton*fMass_proton  );
    
    double R = 1.;
    if (sqrt(s)<10) {
        R = (1 +C9*pow(10-sqrt(s),5))  *  exp(C10*(10-sqrt(s))*pow(x_R-fMass_proton/E_pbar_Max, 2));
    }
    double sigma_in = C11 +C12*log(sqrt(s)) + C13*pow(log(sqrt(s)), 2);
    double X = C8 * pow(log(sqrt(s)/4./fMass_proton),2);
    
    //double f0_p = R * sigma_in * C5 * pow(1-x_R, C6) * exp(-m_T/C7); // 2014 paper
    double f0_p = R * sigma_in * C5 * pow(1-x_R, C6) * pow( 1+X*(m_T-fMass_proton), -1./X/C7 );
    double invCsCM = f0_p;
    return invCsCM;
}


double deltaHyperon(double s, double* C_array, int len_C_array){
    
    C_array_to_double( 1);
    C_array_to_double( 2);
    C_array_to_double( 3);
    C_array_to_double( 4);
    
    double factor   = 0.81;
    
    double hyperon = C1 + C2/(1+pow(C3/s,C4));
    hyperon       *= factor;  // branching ratio to Lambda and Sigma  ###### UNCERTAINTY??
    return hyperon;
}

double deltaIsospin(double s,  double* C_array, int len_C_array){
    
    C_array_to_double(14);
    C_array_to_double(15);
    C_array_to_double(16);
    
    return C14/(1+pow(s/C15, C16));
    
}



//! Parametrization of the invariant antiproton production (pbar) crosssection from pp in CMF
/*!
 *  Taken from:     di Mauro, et al.; 2014;
 *                  A new evaluation of the antiproton production cross section for cosmic ray studies;
 *                  DOI: 10.1103/PhysRevD.90.085017
 *
 *  \f[  E_{\bar{p}} \frac{ d^3 \sigma_{pp}^{(\bar{p})} }{d^3 p_{\bar{p}} }  (s, E_{\bar{p}}, p_{T,\bar{p}})  \f]
 *
 *  \param double s         CM energy.
 *  \param doulbe E_pbar    Energy of the produced antiproton in CMF
 *  \param doulbe pT_pbar   Transverse momentum of the produced antiproton in CMF.
 *  \param doulbe* C_array   Ci=(1,...,16), parameters. C1=C_array[1], C2=C_array[2], ... (C_array[0] is not used)
 * */
double inv_pp_pbar_CM__diMauro( double s, double E_pbar, double pT_pbar, double* C_array, int len_C_array ){
    
    C_array_to_double( 1);
    C_array_to_double( 2);
    C_array_to_double( 3);
    C_array_to_double( 4);
    C_array_to_double( 5);
    C_array_to_double( 6);
    C_array_to_double( 7);
    C_array_to_double( 8);
    C_array_to_double( 9);
    C_array_to_double(10);
    C_array_to_double(11);
    
    if (s<16*fMass_proton*fMass_proton){
    return 0;
    }
    if ( pow(pT_pbar*0.9, 2.) > pow(E_pbar, 2.) - pow(fMass_proton, 2.) ){
        return 0;
    }
    double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
    double x_R          =   E_pbar/E_pbar_Max;
    if ( x_R > 1. ){
    return  0.;
    }
    double sigma_in = tot_pp__diMauro(s) - el_pp__diMauro(s);
    double invCsCM      = sigma_in *
    pow(1 - x_R, C1) *
    exp(-C2 * x_R)*
    fabs(C3 * pow( s, C4 /2. ) * exp( -C5 *pT_pbar                 ) +
         C6 * pow( s, C7 /2. ) * exp( -C8 *pT_pbar*pT_pbar         ) +
         C9 * pow( s, C10/2. ) * exp( -C11*pT_pbar*pT_pbar*pT_pbar )
         );
    return invCsCM;
}

//! Parametrization of the total pp cross section.
/*!
 *  Taken from:     di Mauro, et al.; 2014;
 *                  A new evaluation of the antiproton production cross section for cosmic ray studies;
 *                  DOI: 10.1103/PhysRevD.90.085017
 *
 *  \param double s         CM energy.
 *
 * */
double tot_pp__diMauro(double s){
    if (s<0) return 0;
    double Zpp  = 33.44;
    double Y1pp = 13.53;
    double Y2pp = 6.38 ;
    double n1   = 0.324;
    double n2   = 0.324;
    double hbar2= 0.38937966; // GeV^2 mbarn (PDG)
    double M    = 2.06;
    double Bpp  = M_PI * hbar2/M/M;
    double sM   = pow(2*fMass_proton+M, 2);
    double sigmaPP = Zpp + Bpp*pow( log(s/sM), 2) + Y1pp * pow(sM/s, n1) - Y2pp * pow(sM/s, n2);
    
    return sigmaPP;
}



//! Parametrization of the elastic pp cross section.
/*!
 *  Taken from:     di Mauro, et al.; 2014;
 *                  A new evaluation of the antiproton production cross section for cosmic ray studies;
 *                  DOI: 10.1103/PhysRevD.90.085017
 *
 *  \param double s         CM energy.
 *
 * */
double el_pp__diMauro(double s){
    if (s<0) return 0;
    double Zpp  = 144.98;
    double Y1pp = 2.64;
    double Y2pp = 137.27 ;
    double n1   = 1.57;
    double n2   = -4.65e-3;
    double hbar2= 0.38937966; // GeV^2 mbarn (PDG)
    double M    = 3.06;
    double Bpp  = M_PI * hbar2/M/M;
    double sM   = pow(2*fMass_proton+M, 2);
    double sigmaPP = Zpp + Bpp*pow( log(s/sM), 2) + Y1pp * pow(sM/s, n1) - Y2pp * pow(sM/s, n2);
    
    return sigmaPP;
}


//! Parametrization of the target and projectile overlap function in pbar production
/*!
 *  Taken from:     NA49;
 *                  Inclusive production of protons, anti-protons, neutrons, deuterons and tritons in p+C collisions at 158 GeV/c beam momentum;
 *                  arXiv:1207.6520v3
 *
 *        Fig. 69, Tab. 14
 *
 *  \param   double x_F         Feynman parameter p/p_max.
 *  \return  double             Overlap function for projectile
 **/

double pbar_overlap_function_projectile(double x_F){
    
    int n = 25;
    
    if( x_F < -0.2499 ) return 0;
    if( x_F >  0.25   ) return 1;
    
    double xF[] = {-0.25, -0.225, -0.2, -0.175, -0.15, -0.125, -0.1, -0.075, -0.05, -0.025, 0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25};
    double F[]  = {0., 0.0003, 0.0008, 0.0027, 0.010, 0.035, 0.110, 0.197, 0.295, 0.4, 0.5, 0.6, 0.705, 0.803, 0.890, 0.965, 0.990, 0.9973, 0.9992, 0.9997, 1.0};
    
    double xl = 0;
    double xu = 1;
    double Fl = -1;
    double Fu = -1;
    
    for (int i =0  ; i<n-1; i++) {
        xl = xF[i];
        xu = xF[i+1];
        Fl = F[i];
        Fu = F[i+1];
        if (xu>x_F) break;
    }
    
    double ret = Fl + (Fu-Fl)*(x_F-xl)/(xu-xl);
    
    
    return ret;
    
}

//! Parametrization of the target and target overlap function in pbar production
/*!
 *  Taken from:     NA49;
 *                  Inclusive production of protons, anti-protons, neutrons, deuterons and tritons in p+C collisions at 158 GeV/c beam momentum;
 *                  arXiv:1207.6520v3
 *
 *        Fig. 69, Tab. 14
 *
 *  \param   double x_F         Feynman parameter p/p_max.
 *  \return  double             Overlap function for target
 **/
double pbar_overlap_function_target(double x_F){
    return 1.-pbar_overlap_function_projectile(x_F);
}

