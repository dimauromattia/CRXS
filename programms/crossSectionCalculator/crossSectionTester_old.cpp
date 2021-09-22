#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <unistd.h>

#include <time.h>

#include <configHandler.h>
#include <fileTools.h>

#include "QString"


#define warnout std::cout << "\033[9;31m"
#define warnend "\033[0m" << std::endl


std::string fProgramDescription = "Test Antiproton production cross sections";

//
// Definition of Constants
//
//  All enegies, momenta, and masses are given in   GeV
//  All cross sections are given in                 mbarn
//


double fMass_proton         = 0.9382720813;     // PDG Review 2016
double fMass_neutron        = 0.939565379;      // PDG Review 2016
double fMass_deuteron       = 1.875612928;      // PDG Review 2016

enum fIntegrationType{
    TRAPEZSUM=0,
    UPPERSUM=1,
    LOWERSUM=2
};

enum fErrorOutput{
    NO_OUT=0,
    WARN_OUT=1,
    ALL_OUT=2
};


//! Parametrization of the invariant antiproton production crosssection from pp in CMF
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
 *  \param doulbe Ci        i=(1,...,11), parameters (Default values taken according to the best fit from Eq. 13 of upper reference)
 * */
double InvariantAntiprotonProductionCrossSection_CM_diMauro(double s, double E_pbar, double pT_pbar,
                                                            double C1,
                                                            double C2,
                                                            double C3,
                                                            double C4,
                                                            double C5,
                                                            double C6,
                                                            double C7,
                                                            double C8,
                                                            double C9,
                                                            double C10,
                                                            double C11
                                                            ){
    if (s<16*fMass_proton)
        return 0;
    if ( pow(pT_pbar, 2.) > pow(E_pbar, 2.) - pow(fMass_proton, 2.) )
        return 0;
    double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
    double x_R          =   E_pbar/E_pbar_Max;
    if ( x_R > 1. )
        return  0.;
    double E_inc        =   s/( 2.*fMass_proton ) - 2.*fMass_proton;
    double sigma_0      =   44.4;   //Cross section normalization in mbarn
    double sigma_in     =   sigma_0 * (  1. - 0.62 * exp( - E_inc/0.2  ) * sin( 10.9/pow( E_inc, 0.28 ) )  );
    double invCsCM      = sigma_in *
            pow(1 - x_R, C1) *
            exp(-C2 * x_R)*
            fabs(   C3 * pow( s, C4 /2. ) * exp( -C5 *pT_pbar                 ) +
                    C6 * pow( s, C7 /2. ) * exp( -C8 *pT_pbar*pT_pbar         ) +
                    C9 * pow( s, C10/2. ) * exp( -C11*pT_pbar*pT_pbar*pT_pbar )
                 );
    //    std::cout << "E_pbar_Max "  << E_pbar_Max       << std::endl;
    //    std::cout << "x_R        "  << x_R              << std::endl;
    //    std::cout << "E_inc      "  << E_inc            << std::endl;
    //    std::cout << "sigma_in   "  << sigma_in         << std::endl;
    //    std::cout << "invCsCM    "  << invCsCM          << std::endl;
    return invCsCM;
}
double InvariantAntiprotonProductionCrossSection_CM_diMauro(double s, double E_pbar, double pT_pbar){
    double C1 = 4.448;
    double C2 = 3.75;
    double C3 = 0.00502;
    double C4 = 0.708;
    double C5 = 3.527;
    double C6 = 0.236;
    double C7 = -0.729;
    double C8 = 2.517;
    double C9 = -1.822e-11;
    double C10 = 3.527;
    double C11 = 0.384;
    return InvariantAntiprotonProductionCrossSection_CM_diMauro(s, E_pbar, pT_pbar, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11);
}

double InvariantAntineutronProductionCrossSection_CM_diMauro(double s, double E_pbar, double pT_pbar){
    return 1.3 * InvariantAntiprotonProductionCrossSection_CM_diMauro(s, E_pbar, pT_pbar);
}

double TotalProtonProtonCrossSection_diMauro(double s){
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

double InvariantAntideuteronProductionCrossSection_CM( double s, double E_dbar, double pT_dbar, double p_coalescence,
                                                              double (*totalProtonProtonCrossSection)(double),
                                                              double (*invariantAntiprotonProductionCrossSection)(double, double, double),
                                                              double (*invariantAntineutronProductionCrossSection)(double, double, double) ){
    
    double sReduced = pow( sqrt(s) - E_dbar, 2 );
    double CS         = fMass_deuteron/fMass_proton/fMass_neutron *
                        (4./3. * M_PI * pow(p_coalescence,3)) *
                        1./2. / totalProtonProtonCrossSection(s)*
                        (
                         invariantAntiprotonProductionCrossSection  (s, E_dbar/2, pT_dbar/2)*
                         invariantAntineutronProductionCrossSection (sReduced, E_dbar/2, pT_dbar/2) +
                         invariantAntineutronProductionCrossSection (s, E_dbar/2, pT_dbar/2)*
                         invariantAntiprotonProductionCrossSection  (sReduced, E_dbar/2, pT_dbar/2)
                        );
    return CS;
    
}

double InvariantAntideuteronProductionCrossSection_CM_diMauro(double s, double E_dbar, double pT_dbar){
    
    return InvariantAntideuteronProductionCrossSection_CM( s, E_dbar, pT_dbar, 0.100, TotalProtonProtonCrossSection_diMauro, InvariantAntiprotonProductionCrossSection_CM_diMauro, InvariantAntineutronProductionCrossSection_CM_diMauro);
    
}



//!  Production cross section in Lab Frame from pp collison
/*!
 *  Lab frame:     One proton is at rest and the other one has incident energy E_p
 *
 *  \f[  \frac{ d^3 \sigma_{pp}^{(X), LAB} }{d^3 p_{X} }  (E_p, E_{X}, p_{T,X})  \f]
 *
 *  \param double E_p           Incident proton energy in LAB frame.
 *  \param doulbe E_product     Energy of the produced particle in LAB frame.
 *  \param doulbe pT_product    Transverse momentum of the produced particle in LAB frame.
 *  \param doulbe m_product     Mass of the produced particle.
 *  \param double (*invariantProductionCrossSection)(double, doulbe, double) Corresponding invariant particle production CS in CMF.
 * */

double d3p_ProductionCrossSection_LAB(double E_p_LAB, double E_product_LAB, double cos_theta_product_LAB, double m_product, double (*invariantProductionCrossSection)(double, double, double) )
{
    double s = 2*fMass_proton*fMass_proton + 2 * E_p_LAB * fMass_proton;
    
    //Lorentz Transformation for pbar in CM frame
    double beta          = (E_p_LAB - fMass_proton)/sqrt(E_p_LAB*E_p_LAB - fMass_proton*fMass_proton);
    double gamma         = 1./sqrt(1 - beta*beta);
    double gammabeta     = gamma * beta;
    double p_product_LAB = sqrt(  pow( E_product_LAB, 2 ) - pow( m_product, 2 )  );
    
    double E_product   = gamma * E_product_LAB - gammabeta * p_product_LAB*cos_theta_product_LAB;
    double pT_product  = p_product_LAB*sqrt(1-cos_theta_product_LAB*cos_theta_product_LAB);
    
    //    std::cout << "s        "  << s          << std::endl;
    //    std::cout << "beta     "  << beta       << std::endl;
    //    std::cout << "gamma    "  << gamma      << std::endl;
    //    std::cout << "E_pbar   "  << E_pbar     << std::endl;
    
    return 1.0/E_product_LAB * invariantProductionCrossSection(s, E_product, pT_product);
}


//! Antiproton production cross section in Lab Frame from pp collison
/*!
 *  Lab frame:     One proton is at rest and the other one has incident energy E_p
 *
 *  \f[  \frac{ d^3 \sigma_{pp}^{(\bar{p}, LAB)} }{d^3 p_{\bar{p}} }  (E_p, E_{\bar{p}}, p_{T,\bar{p}})  \f]
 *
 *  \param double E_p       Incident proton energy in LAB frame.
 *  \param doulbe E_pbar    Energy of the produced antiproton in LAB frame.
 *  \param doulbe pT_pbar   Transverse momentum of the produced antiproton in LAB frame.
 *  \param double (*invariantAntiprotonProductionCrossSection)(double, doulbe, double) 
 *                          Corresponding invariant pabr production CS in CMF, Default is di Mauro.
 * */
double d3p_AntiprotonProductionCrossSection_LAB(double E_p_LAB, double E_pbar_LAB, double cos_theta_pbar_LAB, double (*invariantAntiprotonProductionCrossSection)(double, double, double)=InvariantAntiprotonProductionCrossSection_CM_diMauro ){
    return d3p_ProductionCrossSection_LAB(E_p_LAB, E_pbar_LAB, cos_theta_pbar_LAB, fMass_proton, invariantAntiprotonProductionCrossSection);
}


//!  Energy differential production cross section in Lab Frame from pp collison
/*!
 *  Lab frame:     One proton is at rest and the other one has incident energy E_p
 *
 *  \f[  \frac{ d \sigma_{pp}^{(X), LAB} }{d E_{X} }  (E_p, E_{X})  \f]
 *
 *  \param double E_p           Incident proton energy in LAB frame.
 *  \param doulbe E_product     Energy of the produced particle in LAB frame.
 *  \param doulbe m_product     Mass of the produced particle.
 *  \param double (*invariantProductionCrossSection)(double, doulbe, double) 
 *                              Corresponding invariant particle production CS in CMF.
 *  \param int    precision     Number of steps for the integration over cos(theta).
 * */
double dE_ProductionCrossSection_LAB(double E_p_LAB, double E_product_LAB, double m_product, double (*invariantProductionCrossSection)(double, double, double), double m_remainder, int precision=100000, int output=WARN_OUT ){
    
    double p_product_Lab = sqrt(E_product_LAB*E_product_LAB-m_product*m_product);
    
    double phi_integration              = 2*M_PI;
    double Jacobian_and_conversion      = p_product_Lab*E_product_LAB;
    //     Jacobian_determinant         = p_product_Lab*p_product_Lab;
    //     dp_to_dE_conversion          = E_product_LAB/p_product_Lab;
    
    //Calculate maximal cos_theta
    
    double s = 2*fMass_proton*fMass_proton + 2 * E_p_LAB * fMass_proton;
    
    double beta          = (E_p_LAB - fMass_proton)/sqrt(E_p_LAB*E_p_LAB - fMass_proton*fMass_proton);
    double gamma         = 1./sqrt(1 - beta*beta);
    double gammabeta     = gamma * beta;
    double p_product_LAB = sqrt(  pow( E_product_LAB, 2 ) - pow( m_product, 2 )  );
    
    double cos_theta_max = (  gamma * E_product_LAB -  (s + m_product*m_product - m_remainder*m_remainder)/2/sqrt(s) ) / (gammabeta * p_product_LAB);
    
    //    std::cout << "s             "  << s               << std::endl;
    //    std::cout << "beta          "  << beta            << std::endl;
    //    std::cout << "gamma         "  << gamma           << std::endl;
    //    std::cout << "E_pbar        "  << E_pbar          << std::endl;
    //    std::cout << "cos_theta_m   "  << cos_theta_max   << std::endl;
    
    
    //Numerical integration for theta integration, cos(theta) in 'precision' steps
    double integral_Trapez  = 0;
    double integral_Upper   = 0;
    double integral_Lower   = 0;
    double delta_cos_theta  = (1.-cos_theta_max)/precision;
    double cos_theta        = cos_theta_max;
    
    double lastCS = d3p_ProductionCrossSection_LAB(E_p_LAB, E_product_LAB, cos_theta, m_product, invariantProductionCrossSection);
    for ( ; cos_theta<=1; cos_theta+=delta_cos_theta) {
        
        double cos_theta_product_LAB = cos_theta;
        double E_product   = gamma * E_product_LAB - gammabeta * p_product_LAB*cos_theta_product_LAB;
        double pT_product  = p_product_LAB*sqrt(1-cos_theta_product_LAB*cos_theta_product_LAB);
        double d3p_ProductionCrossSection_LAB = 1.0/E_product_LAB * invariantProductionCrossSection(s, E_product, pT_product);
        
        double CS = d3p_ProductionCrossSection_LAB;
        
        integral_Trapez  += (lastCS+CS)/2;
        integral_Upper   += std::max(lastCS, CS);
        integral_Lower   += std::min(lastCS, CS);
        
        lastCS = CS;
    }
    
    integral_Trapez    *=  delta_cos_theta*phi_integration*Jacobian_and_conversion;
    integral_Upper     *=  delta_cos_theta*phi_integration*Jacobian_and_conversion;
    integral_Lower     *=  delta_cos_theta*phi_integration*Jacobian_and_conversion;
    
    double error = std::max( integral_Upper-integral_Trapez, integral_Trapez-integral_Lower );
    
    if(integral_Trapez>0 && error/integral_Trapez>0.01 && output>NO_OUT){
        warnout << "Warning: Error is more than 1%  " << warnend;
    }
    if (output==ALL_OUT) {
        std::cout << "Trapez Sum:      " << integral_Trapez << std::endl;
        std::cout << "Upper Sum:       " << integral_Upper  << std::endl;
        std::cout << "Lower Sum:       " << integral_Lower  << std::endl;
        std::cout << "Relative Error:  " << error/integral_Trapez  << std::endl;
    }
    
    return integral_Trapez;
}



double dE_AntiprotonProductionCrossSection_LAB(double E_p_LAB, double E_pbar_LAB, int precision=100000, int output=WARN_OUT, double (*invariantAntiprotonProductionCrossSection)(double, double, double)=InvariantAntiprotonProductionCrossSection_CM_diMauro ){
    return dE_ProductionCrossSection_LAB(E_p_LAB, E_pbar_LAB, fMass_proton, invariantAntiprotonProductionCrossSection, 3*fMass_proton, precision, output);
}



double dE_AntideuteronProductionCrossSection_LAB(double E_p_LAB, double E_dbar_LAB, int precision=100000, int output=WARN_OUT,
                                                 double (* invariantAntideuteronProductionCrossSection)(double, double, double)=InvariantAntideuteronProductionCrossSection_CM_diMauro ){
    return dE_ProductionCrossSection_LAB(E_p_LAB, E_dbar_LAB, fMass_deuteron, invariantAntideuteronProductionCrossSection, 4*fMass_proton, precision, output);
}



int main(int argc, char *argv[])
{
    // Config handling
    CRACS::ConfigHandler* config = CRACS::ConfigHandler::GetInstance();
    config->SetInput(argc, argv, fProgramDescription);
    
    int step = 1;
    config->AddOptionInt     (  "step",        step,       "Step    Default: 1"                     );
    double s;
    double E;
    double Ep;
    double pT;
    int n = 100000;
    int o = 1;
    config->AddOptionDouble  (  "s",  s  );
    config->AddOptionDouble  (  "Ep", Ep );
    config->AddOptionDouble  (  "E",  E  );
    config->AddOptionDouble  (  "pT", pT );
    config->AddOptionInt     (  "n",  n  );
    config->AddOptionInt     (  "o",  o  );
    config->CheckOption();
    
    
    if (step==1) {
        std::cout << InvariantAntiprotonProductionCrossSection_CM_diMauro(s, E, pT) << std::endl;
    }
    
    if (step==2) {
        std::cout << d3p_AntiprotonProductionCrossSection_LAB(Ep, E, pT) << std::endl;
    }
    
    if (step==3) {
        std::cout << dE_AntiprotonProductionCrossSection_LAB(Ep+fMass_proton, E+fMass_proton, n, o) << std::endl;
    }
    
    if (step==4) {
        for (double dE=0; dE<=8; dE+=0.1) {
            double Ekin_p = pow(10, dE);
            for (double dx=-6; dx<=0; dx+=0.1) {
                double Ekin_pbar = Ekin_p*pow(10, dx);
                std::cout << Ekin_p << "    " << Ekin_pbar << "    " << dE_AntiprotonProductionCrossSection_LAB(Ekin_p+fMass_proton, Ekin_pbar+fMass_proton) << std::endl;
            }
        }
    }
    
    if (step==5) {
        double Ekin_p = 300.081;
        for (double dx=-6; dx<=0; dx+=0.05) {
            double Ekin_pbar = Ekin_p*pow(10, dx);
            std::cout << Ekin_p << "    " << Ekin_pbar << "    " << dE_AntiprotonProductionCrossSection_LAB(Ekin_p+fMass_proton, Ekin_pbar+fMass_proton) << std::endl;
        }
        
    }
    
    if (step==10) {
        CRACS::FileTool f("sigma_bestfitalldef_eq13.dat");
        f.ExtractNumberTable(3, " ");
        for (int i=0; i<f.GetNumberOfLines()-1; i++) {
            double Ekin_p       = f.NumberTable(i, 0);
            double Ekin_pbar    = f.NumberTable(i, 1);
            double cs_serpico   = f.NumberTable(i, 2)*1e31;
            double cs_my        = dE_AntiprotonProductionCrossSection_LAB(Ekin_p+fMass_proton, Ekin_pbar+fMass_proton, n, o);
            std:: cout << Ekin_p << "    " << Ekin_pbar << "    " << cs_serpico << "    " << cs_my << std::endl;
        }
    }
    
    if (step==20) {
        double Ekin_p = Ep;
        for (double dx=-3; dx<=0; dx+=0.05) {
            double Ekin_dbar = Ekin_p*pow(10, dx);
            std::cout << Ekin_p << "    " << Ekin_dbar/2 << "    " << dE_AntideuteronProductionCrossSection_LAB(Ekin_p+fMass_proton, Ekin_dbar+2*fMass_proton, n, o) << std::endl;
        }
    }
    
    if (step==30) {
        for (double ds=-2; ds<=4; ds+=0.1) {
            double s = pow(10, ds);
            std::cout << s << "    " << TotalProtonProtonCrossSection_diMauro(s) << std::endl;
        }
    }
    
    
    return 0;
    
    
}
