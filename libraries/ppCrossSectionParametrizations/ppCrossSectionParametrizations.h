#ifndef PPCROSSSECTIONPARAMETRIZATIONS_H
#define PPCROSSSECTIONPARAMETRIZATIONS_H

#include <iostream>
#include <math.h>
#include <definitions.h>




namespace CRACS {
    
    

    
    /*! \brief Parametrization of various \f p + p \rightarrow f + X \f cross sections.
     *
     *  Differential cross sections (CS) are given in the center of mass (CM) frame.
     *  They shall be given in the invariant form
     *
     *  \f[  E \frac{ d^3 \sigma_{pp}^{(f)} }{d^3 p_{f} }  (s, E_{f}, p_{f})  \f]
     *
     *  All cross sections are given in mbarn
     *  All enegies, momenta, and masses have unit GeV.
     *
     */

    class ppCSParametrizations{
        
    public:
        
       
        
        //! Parametrization of the total pp cross section.
        /*!
         *  Taken from:     di Mauro, et al.; 2014;
         *                  A new evaluation of the antiproton production cross section for cosmic ray studies;
         *                  DOI: 10.1103/PhysRevD.90.085017
         *
         *  \param double s         CM energy.
         *
         * */
        static double tot_pp__diMauro(double s){
            if (s<0)
                return 0;
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
        static double el_pp__diMauro(double s){
            if (s<0)
                return 0;
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
        
        //! inv_pA_pbar_CM()
        /*!
         *
         *  Transform a pp cross section to pA accroding to the scaling introduced by
         *
         *  Ref:            R.P. Duperray, et al.; 2003;
         *                  Parameterization of the antiproton inclusive production cross section on nuclei;
         *                  DOI: 10.1103/PhysRevD.68.094017
         *                  Eq. 4
         *
         *  \param double s         CM energy.
         *  \param doulbe E_pbar    Energy of the produced antiproton in CMF, encode backward scattering by negative E_pbar
         *  \param doulbe pT_pbar   Transverse momentum of the produced antiproton in CMF.
         *  \param double inv_pA_pbar_CM*  Invariant ross section in CM fram as function onf s, E_pbar, pT_pbar.
         *  \param int    A         Mass Number
         * */
        static double inv_pA_pbar_CM(  double s, double E_pbar_d, double pT_pbar, double (*inv_pp_CM)(double, double, double), int A){
            
            double E_pbar = fabs(E_pbar_d);
            
            double CS = inv_pp_CM(s, E_pbar, pT_pbar);
            
            if (CS<1e-90)
                return 0;
            
            
            double D1  = 0.7;
            double D2  = 0.16990;
            double D3  = 10.28;
            
            double sigma_0      =   (  1+0.016*sin( 5.3-2.63*log(A) )  );
            double sigma_0_1    =   (  1+0.016*sin( 5.3-2.63*log(1) )  );
            double factor       =   pow( A, D1 + D2*log( sqrt(s)/D3)*pT_pbar  )*sigma_0/sigma_0_1;

            double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
            double x_R_d        =   E_pbar_d/E_pbar_Max;
            
            double D4 = 0.25;
            double factor2 = A*D4*x_R_d;
            
            
                       
            //out(factor)
            return factor * factor2 * CS;
        }
        
        static double inv_pA_pbar_CM__fromDuperray(  double s, double E_pbar_d, double pT_pbar, double (*inv_pp_CM)(double, double, double), int A){
            
            double E_pbar = fabs(E_pbar_d);
            
            double CS = inv_pp_CM(s, E_pbar, pT_pbar);
            
            if (CS<1e-90)
            return 0;
            
            
            double D1  = 0.7;
            double D2  = 0.16990;
            double D3  = 10.28;
            
            double sigma_0      =   (  1+0.016*sin( 5.3-2.63*log(A) )  );
            double sigma_0_1    =   (  1+0.016*sin( 5.3-2.63*log(1) )  );
            double factor       =   pow( A, D1 + D2*log( sqrt(s)/D3)*pT_pbar  )*sigma_0/sigma_0_1;
            
            //out(factor)
            
            if(factor!=factor) return 0;
            
            return factor * CS;
        }
        
        
////////////////////////
//  Tan and Ng
////////////////////////
        
        
        //! Parametrization of the invariant antiproton production (pbar) crosssection from pp in CMF
        /*!
         *  Taken from:     L C Tan and L K Ng; 1983;
         *                  Calculation of the equilibrium antiproton spectrum;
         *                  Journal of Physics G: Nuclear Physics
         *
         *  \f[  E_{\bar{p}} \frac{ d^3 \sigma_{pp}^{(\bar{p})} }{d^3 p_{\bar{p}} }  (s, E_{\bar{p}}, p_{T,\bar{p}})  \f]
         *
         *  \param double s         CM energy.
         *  \param doulbe E_pbar    Energy of the produced antiproton in CMF
         *  \param doulbe pT_pbar   Transverse momentum of the produced antiproton in CMF.
         * */
        static double inv_pp_pbar_CM__Tan_Ng(double s, double E_pbar, double pT_pbar){
            if (s<16*fMass_proton)
                return 0;
            if ( pow(pT_pbar, 2.) > pow(E_pbar, 2.) - pow(fMass_proton, 2.) )
                return 0;
            double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
            double x_R          =   E_pbar/E_pbar_Max;
            if ( x_R > 1. )
                return  0.;
            
            
            double theta_0_5__minus__x_R = 1;
            if (0.5-x_R<0){
                theta_0_5__minus__x_R = 0;
            }
            
            double f = 3.34 * exp( -17.6*x_R )*theta_0_5__minus__x_R + 2.10* pow(1-x_R,7.80);
            double A = 3.95 * exp( -2.76*x_R );
            double B = 40.5 * exp( -3.21*x_R )* pow( x_R, 2.13 );
            
            double invCsCM      = f * exp( -A*pT_pbar - B*pT_pbar*pT_pbar );
            
            return invCsCM;
        }
        
        
        static double inv_pHe_pbar_CM__Tan_Ng(double s, double E_pbar, double pT_pbar){
            return inv_pA_pbar_CM__fromDuperray(s, E_pbar, pT_pbar, inv_pp_pbar_CM__Tan_Ng, 4);
        }

        
////////////////////////
//  di Mauro
////////////////////////
        
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
         *  \param doulbe Ci        i=(1,...,11), parameters
         * */
        static double inv_pp_pbar_CM__diMauro(double s, double E_pbar, double pT_pbar,
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
            E_pbar = fabs(E_pbar);
//            out(s)
//            out(pow(pT_pbar, 2.))
//            out(pow(E_pbar, 2.) - pow(fMass_proton, 2.))
            if (s<16*fMass_proton*fMass_proton)
                return 0;
            if ( pow(pT_pbar*0.9, 2.) > pow(E_pbar, 2.) - pow(fMass_proton, 2.) ){
                varOut(pT_pbar)
                varOut(sqrt(pow(E_pbar, 2.) - pow(fMass_proton, 2.)))
                return 0;
            }
            double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
            double x_R          =   E_pbar/E_pbar_Max;
            if ( x_R > 1. )
                return  0.;
            //double E_inc        =   s/( 2.*fMass_proton ) - 2.*fMass_proton;
            //double sigma_0      =   44.4;   //Cross section normalization in mbarn
            //double sigma_in     =   sigma_0 * (  1. - 0.62 * exp( - E_inc/0.2  ) * sin( 10.9/pow( E_inc, 0.28 ) )  );
            double sigma_in = tot_pp__diMauro(s) - el_pp__diMauro(s);
            double invCsCM      = sigma_in *
            pow(1 - x_R, C1) *
            exp(-C2 * x_R)*
            fabs(C3 * pow( s, C4 /2. ) * exp( -C5 *pT_pbar                 ) +
                 C6 * pow( s, C7 /2. ) * exp( -C8 *pT_pbar*pT_pbar         ) +
                 C9 * pow( s, C10/2. ) * exp( -C11*pT_pbar*pT_pbar*pT_pbar )
                 );
//            out(invCsCM)
            return invCsCM;
        }
        
        
        //! It is the same as inv_pp_pbar_CM__diMauro() with default values.
        /*!
         *  Default values for Ci are taken according to the best fit from Eq. 13 of upper reference
         *
         *  \param double s         CM energy.
         *  \param doulbe E_pbar    Energy of the produced antiproton in CMF
         *  \param doulbe pT_pbar   Transverse momentum of the produced antiproton in CMF.
         * */
        static double inv_pp_pbar_CM__diMauro(double s, double E_pbar, double pT_pbar){
            double C1 = 4.448;
            double C2 = 3.735;
            double C3 = 0.00502;
            double C4 = 0.708;
            double C5 = 3.527;
            double C6 = 0.236;
            double C7 = -0.729;
            double C8 = 2.517;
            double C9 = -1.822e-11;
            double C10 = 3.527;
            double C11 = 0.384;
            return inv_pp_pbar_CM__diMauro(s, E_pbar, pT_pbar, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11);
        }
        
        static double inv_pp_nbar_CM__diMauro(double s, double E_pbar, double pT_pbar){
            return 1.3*inv_pp_pbar_CM__diMauro(s, E_pbar, pT_pbar);
        }
        
        static double inv_pHe_pbar_CM__diMauro(double s, double E_pbar, double pT_pbar){
            return inv_pA_pbar_CM__fromDuperray(s, E_pbar, pT_pbar, inv_pp_pbar_CM__diMauro, 4);
        }
        
        static double inv_HeHe_pbar_CM__diMauro(double s, double E_pbar, double pT_pbar){
            return inv_pA_pbar_CM__fromDuperray(s, E_pbar, pT_pbar, inv_pHe_pbar_CM__diMauro, 4);
        }
        

        
////////////////////////
//  di Mauro 12
////////////////////////
        
        
        //! It is the same as inv_pp_pbar_CM__diMauro() with default values.
        /*!
         *  Default values for Ci are taken according to the best fit from Eq. 12 of upper reference
         *
         *  \param double s         CM energy.
         *  \param doulbe E_pbar    Energy of the produced antiproton in CMF
         *  \param doulbe pT_pbar   Transverse momentum of the produced antiproton in CMF.
         * */
        static double inv_pp_pbar_CM__diMauro12(double s, double E_pbar, double pT_pbar){
            double C1 = 4.499;
            double C2 = 3.41;
            double C3 = 0.00942;
            double C4 = 0.445;
            double C5 = 3.502;
            double C6 = 0.0622;
            double C7 = -0.247;
            double C8 = 2.576;
            double C9 = 0;
            double C10 = 0;
            double C11 = 0;
            return inv_pp_pbar_CM__diMauro(s, E_pbar, pT_pbar, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11);
        }
        
        static double inv_pp_nbar_CM__diMauro12(double s, double E_pbar, double pT_pbar){
            return 1.3*inv_pp_pbar_CM__diMauro12(s, E_pbar, pT_pbar);
        }
        
        static double inv_pHe_pbar_CM__diMauro12(double s, double E_pbar, double pT_pbar){
            return inv_pA_pbar_CM__fromDuperray(s, E_pbar, pT_pbar, inv_pp_pbar_CM__diMauro12, 4);
        }
        
        static double inv_HeHe_pbar_CM__diMauro12(double s, double E_pbar, double pT_pbar){
            return inv_pA_pbar_CM__fromDuperray(s, E_pbar, pT_pbar, inv_pHe_pbar_CM__diMauro12, 4);
        }
        
        
        
////////////////////////
//  Korsmeier
////////////////////////
        
        
        //! It is the same as inv_pp_pbar_CM__diMauro() with default values, refitted Korsmeier.
        /*!
         *  Default values for Ci are taken from the best from "crossSectionFitter.cpp"
         *
         *  \param double s         CM energy.
         *  \param doulbe E_pbar    Energy of the produced antiproton in CMF
         *  \param doulbe pT_pbar   Transverse momentum of the produced antiproton in CMF.
         * */
        static double inv_pp_pbar_CM__Korsmeier(double s, double E_pbar, double pT_pbar){
            double C1 = 4.43390e+00;
            double C2 = 4.18682e+00;
            double C3 = 5.82113e-03;
            double C4 = 5.88957e-01;
            double C5 = 3.38266e+00;
            double C6 = 1.47670e-01;
            double C7 =-5.37479e-01;
            double C8 = 2.38055e+00;
            double C9 = 0;
            double C10= 0;
            double C11= 0;
            return inv_pp_pbar_CM__diMauro(s, E_pbar, pT_pbar, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11);
        }
        
        static double inv_pp_nbar_CM__Korsmeier(double s, double E_pbar, double pT_pbar){
            return 1.3*inv_pp_pbar_CM__Korsmeier(s, E_pbar, pT_pbar);
        }
        
        static double inv_pHe_pbar_CM__Korsmeier(double s, double E_pbar, double pT_pbar){
            return inv_pA_pbar_CM__fromDuperray(s, E_pbar, pT_pbar, inv_pp_pbar_CM__Korsmeier, 4);
        }
        
        
////////////////////////
//  Winkler (Kappl)
////////////////////////
        
        //! Parametrization of the invariant antiproton production (pbar) crosssection from pp in CMF
        /*!
         *  Taken from:     Winkler, M. W.; 2016;
         *                  Cosmic Ray Antiprotons at High Energies;
         *                  arXiv:1701.04866
         *
         *  \f[  E_{\bar{p}} \frac{ d^3 \sigma_{pp}^{(\bar{p})} }{d^3 p_{\bar{p}} }  (s, E_{\bar{p}}, p_{T,\bar{p}})  \f]
         *
         *  \param double s         CM energy.
         *  \param doulbe E_pbar    Energy of the produced antiproton in CMF
         *  \param doulbe pT_pbar   Transverse momentum of the produced antiproton in CMF.
         *  \param doulbe Ci        i=(1,...,16), parameters
         * */
        static double inv_pp_pbar_CM__Winkler(double s, double E_pbar_d, double pT_pbar, bool with_antineutrons, bool with_hyperons,
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
                                              double C11,
                                              double C12,
                                              double C13,
                                              double C14,
                                              double C15,
                                              double C16
                                              ){
            double E_pbar = fabs(E_pbar_d);
            if (s<16*fMass_proton*fMass_proton)
                return 0;
            if ( pow(pT_pbar*0.9, 2.) > pow(E_pbar, 2.) - pow(fMass_proton, 2.) ){
                varOut(pT_pbar)
                varOut(sqrt(pow(E_pbar, 2.) - pow(fMass_proton, 2.)))
                return 0;
            }
            double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
            double x_R          =   E_pbar/E_pbar_Max;
            if ( x_R > 1. )
                return  0.;
             
            double m_T = sqrt(  pT_pbar*pT_pbar  +  fMass_proton*fMass_proton  );
            
            double hyperon = C1 + C2/(1+pow(C3/s,C4));
            hyperon       *= 0.81;  // branching ratio to Lambda and Sigma
            
            double isospin = C14/(1+pow(s/C15, C16));
            
            double R = 1.;
            if (sqrt(s)<10) {
                R = (1 +C9*pow(10-sqrt(s),5))  *  exp(C10*(10-sqrt(s))*pow(x_R-fMass_proton/E_pbar_Max, 2));
            }
            
            double sigma_in = C11 +C12*log(sqrt(s)) + C13*pow(log(sqrt(s)), 2);
            
            
            double X = C8 * pow(log(sqrt(s)/4./fMass_proton),2);
            
            
            //double f0_p = R * sigma_in * C5 * pow(1-x_R, C6) * exp(-m_T/C7); // 2014 paper
            double f0_p = R * sigma_in * C5 * pow(1-x_R, C6) * pow( 1+X*(m_T-fMass_proton), -1./X/C7 );//*exp(-fMass_proton/C7);
            
            
            double invCsCM = f0_p;
            
            if (with_antineutrons){
                invCsCM *= (2+isospin);
            }
            
            if (with_hyperons)
                invCsCM *= (1+hyperon);
            
            
            
            return invCsCM;
        }
        
        static double inv_pp_pbar_CM__Winkler(double s, double E_pbar, double pT_pbar, bool with_antineutron, bool with_hyperon){
            double C1 = 0.31;
            double C2 = 0.30;
            double C3 = 146.*146.;
            double C4 = 0.9;
            
            double C5 = 0.047;
            double C6 = 7.76;
            double C7 = 0.168;
            
            double C8 = 0.038;
            
            double C9 = 1.0e-3;
            double C10= 0.7;
            
            double C11 = 30.9;
            double C12 = -1.74;
            double C13 = 0.71;
            
            double C14 = 0.114;
            double C15 = 144*144;
            double C16 = 0.51;
            return inv_pp_pbar_CM__Winkler(s, E_pbar, pT_pbar, with_antineutron, with_hyperon, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15, C16);
        }
        
        static double inv_pp_pbar_CM__Winkler(double s, double E_pbar, double pT_pbar){
            double res = inv_pp_pbar_CM__Winkler(s, E_pbar, pT_pbar, false, false);
            if (res!=res) return 0;
            return res;
        }
        
        static double inv_pp_pbar_CM__WinklerWithHypWithNbar(double s, double E_pbar, double pT_pbar){
            double res = inv_pp_pbar_CM__Winkler(s, E_pbar, pT_pbar, true, true);
            if (res!=res) return 0;
            return res;
        }
        
        static double IS_winkler(double s, double C14 = 0.114, double C15 = 20736, double C16 = 0.51){
            return  C14/(1+pow(s/C15, C16));
        }
        
        static double inv_pHe_pbar_CM__WinklerWithHypWithNbar(double s, double E_pbar_d, double pT_pbar){
            
            double E_pbar = fabs(E_pbar_d);
            double p  = sqrt( E_pbar*E_pbar-fMass_proton*fMass_proton );
            double pL = sqrt( p*p-pT_pbar*pT_pbar );
            
            if (E_pbar_d<0) pL*=-1;
            
            double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
            double p_pbar_Max   =   sqrt( E_pbar_Max*E_pbar_Max-fMass_proton*fMass_proton );
            
            double xF = pL*2./sqrt(s);
            
            double nu_He = 1.25;
            double factor = 4./nu_He * ( nu_He*(1.+IS_winkler(s)*0.5)*pbar_overlap_function_target(xF) +  pbar_overlap_function_projectile(xF) );
            
            if (factor!=factor) return 0;
            
            return inv_pp_pbar_CM__Winkler(s, E_pbar_d, pT_pbar, true, true)*factor;
        }
        
        
        static double inv_pHe_pbar_CM__Winkler(double s, double E_pbar_d, double pT_pbar){
            
            double E_pbar = fabs(E_pbar_d);
            double p  = sqrt( E_pbar*E_pbar-fMass_proton*fMass_proton );
            double pL = sqrt( p*p-pT_pbar*pT_pbar );
            if (pL!=pL) return 0;
            
            if (E_pbar_d<0) pL*=-1;
            
            double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
            double p_pbar_Max   =   sqrt( E_pbar_Max*E_pbar_Max-fMass_proton*fMass_proton );
            
            double xF = pL*2./sqrt(s);
            
            double nu_He = 1.25;
            double factor = 4./nu_He * ( nu_He*(1+IS_winkler(s)*0.5)*pbar_overlap_function_target(xF) +  pbar_overlap_function_projectile(xF) );
            
            
            return inv_pp_pbar_CM__Winkler(s, E_pbar_d, pT_pbar, false, false)*factor;
        }
        

        static double inv_Hep_pbar_CM__WinklerWithHypWithNbar(double s, double E_pbar_d, double pT_pbar){
            
            double E_pbar = fabs(E_pbar_d);
            double p  = sqrt( E_pbar*E_pbar-fMass_proton*fMass_proton );
            double pL = sqrt( p*p-pT_pbar*pT_pbar );
            
            if (E_pbar_d<0) pL*=-1;
            
            double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
            double p_pbar_Max   =   sqrt( E_pbar_Max*E_pbar_Max-fMass_proton*fMass_proton );
            
            double xF = pL*2./sqrt(s);
            
            double nu_He = 1.25;
            double factor = 4./nu_He * ( pbar_overlap_function_target(xF) +  nu_He*(1+IS_winkler(s)*0.5)*pbar_overlap_function_projectile(xF) );
            
            if (factor!=factor) return 0;
            
            return inv_pp_pbar_CM__Winkler(s, E_pbar_d, pT_pbar, true, true)*factor;
        }
        
        static double inv_Hep_pbar_CM__Winkler(double s, double E_pbar_d, double pT_pbar){
            
            double E_pbar = fabs(E_pbar_d);
            double p  = sqrt( E_pbar*E_pbar-fMass_proton*fMass_proton );
            double pL = sqrt( p*p-pT_pbar*pT_pbar );
            
            if (E_pbar_d<0) pL*=-1;
            
            double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
            double p_pbar_Max   =   sqrt( E_pbar_Max*E_pbar_Max-fMass_proton*fMass_proton );
            
            double xF = pL*2./sqrt(s);
            
            double nu_He = 1.25;
            double factor = 4./nu_He * ( pbar_overlap_function_target(xF) +  nu_He*(1+IS_winkler(s)*0.5)*pbar_overlap_function_projectile(xF) );
            
            if (factor!=factor) return 0;
            
            return inv_pp_pbar_CM__Winkler(s, E_pbar_d, pT_pbar, false, false)*factor;
        }
        

        static double inv_HeHe_pbar_CM__WinklerWithHypWithNbar(double s, double E_pbar_d, double pT_pbar){
            
            double E_pbar = fabs(E_pbar_d);
            double p  = sqrt( E_pbar*E_pbar-fMass_proton*fMass_proton );
            double pL = sqrt( p*p-pT_pbar*pT_pbar );
            
            
            if (E_pbar_d<0) pL*=-1;
            
            double xF = pL*2./sqrt(s);
            
            double nu_He = 1.25;
            double factor = pow(4./nu_He,2) *(1+IS_winkler(s)*0.5)* ( nu_He*pbar_overlap_function_target(xF) +  nu_He*pbar_overlap_function_projectile(xF) );
            
            if (factor!=factor) return 0;
            
            return inv_pp_pbar_CM__Winkler(s, E_pbar_d, pT_pbar, true, true)*factor;
        }
        
        static double inv_HeHe_pbar_CM__Winkler(double s, double E_pbar_d, double pT_pbar){
            
            double E_pbar = fabs(E_pbar_d);
            double p  = sqrt( E_pbar*E_pbar-fMass_proton*fMass_proton );
            double pL = sqrt( p*p-pT_pbar*pT_pbar );
            
            if (E_pbar_d<0) pL*=-1;
            
            double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
            
            double xF = pL*2./sqrt(s);
            
            double nu_He = 1.25;
            double factor = pow(4./nu_He,2) *(1+IS_winkler(s)*0.5)* ( nu_He*pbar_overlap_function_target(xF) +  nu_He*pbar_overlap_function_projectile(xF) );
            
            if (factor!=factor) return 0;
            
            return inv_pp_pbar_CM__Winkler(s, E_pbar_d, pT_pbar, false, false)*factor;
        }
        
        
        //! Parametrization of the target and projectile overlap function in pbar production
        /*!
         *  Taken from:     NA49;
         *                  Inclusive production of protons, anti-protons, neutrons, deuterons and tritons in p+C collisions at 158 GeV/c beam momentum;
         *                  arXiv:1207.6520v3
         *
         *        Fig. 69, Tab. 14
         *
         *  \param double x_F         Feynman parameter p/p_max.
         **/
        static double pbar_overlap_function_projectile(double x_F){
            
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
        
        static double pbar_overlap_function_target(double x_F){
            
            return 1.-pbar_overlap_function_projectile(x_F);
            
        }
        
////////////////////////
//  Duperray
////////////////////////
        
        
        //! Parametrization of the invariant antiproton production (pbar) crosssection from p p in CMF
        /*!
         *  Taken from:     R.P. Duperray, et al.; 2003;
         *                  Parameterization of the antiproton inclusive production cross section on nuclei;
         *                  DOI: 10.1103/PhysRevD.68.094017
         *                  Eq. 6
         *
         *  \f[  E_{\bar{p}} \frac{ d^3 \sigma_{pA}^{(\bar{p})} }{d^3 p_{\bar{p}} }  (s, E_{\bar{p}}, p_{T,\bar{p}})  \f]
         *
         *  The CMF referrs to the nucleon nucleon system!
         *
         *  \param double s         CM energy.
         *  \param doulbe E_pbar    Energy of the produced antiproton in CMF
         *  \param doulbe pT_pbar   Transverse momentum of the produced antiproton in CMF.
         *  \param doulbe Ci        i=(1,...,10), parameters
         * */
        static double inv_pp_pbar_CM__Duperray(double s, double E_pbar, double pT_pbar,
                                               double D1,
                                               double D2,
                                               double D3,
                                               double D4,
                                               double D5,
                                               double D6,
                                               double D7
                                               ){
            E_pbar = fabs(E_pbar);
            if (s<16*fMass_proton)
                return 0;
            if ( pow(pT_pbar, 2.) > pow(E_pbar, 2.) - pow(fMass_proton, 2.) )
                return 0;
            double E_pbar_Max   =   ( s-8*fMass_proton*fMass_proton )/2./sqrt( s );
            double x_R          =   E_pbar/E_pbar_Max;
            if ( x_R > 1. )
                return  0.;
            double E_inc        =   (  s + 2*fMass_proton*fMass_proton )/( 2.*fMass_proton );
            double sigma_0      =   45*(  1+0.016*sin( 5.3 )  );
            double sigma_in     =   sigma_0 * (  1. - 0.62 * exp( - E_inc*1000./200.  ) * sin( 10.9/pow( E_inc*1000, 0.28 ) )  );
            double invCsCM      =   sigma_in *
            pow( 1 - x_R, D1 ) *
            exp(-D2 * x_R)*
            (   D3 * pow( s, D4 /2. ) * exp( -D5 *pT_pbar          ) +
             D6 *                    exp( -D7 *pT_pbar*pT_pbar  )
             );
            if (invCsCM!=invCsCM) return 0;
            return invCsCM;
        }
        
        //! Parametrization of the invariant antiproton production (pbar) crosssection from p p in CMF
        /*!
         *  Taken from:     R.P. Duperray, et al.; 2003;
         *                  Parameterization of the antiproton inclusive production cross section on nuclei;
         *                  DOI: 10.1103/PhysRevD.68.094017
         *                  Eq. 6, with the given best fit parameters Di (D1 and D2 have to be interchanged, see di Mauro et al. 2014)
         *
         *  \f[  E_{\bar{p}} \frac{ d^3 \sigma_{pA}^{(\bar{p})} }{d^3 p_{\bar{p}} }  (s, E_{\bar{p}}, p_{T,\bar{p}})  \f]
         *
         *  The CMF referrs to the nucleon nucleon system!
         *
         *  \param double s         CM energy.
         *  \param doulbe E_pbar    Energy of the produced antiproton in CMF
         *  \param doulbe pT_pbar   Transverse momentum of the produced antiproton in CMF.
         * */
        static double inv_pp_pbar_CM__Duperray(double s, double E_pbar, double pT_pbar){
            double D1   = 4.34;
            double D2   = 3.461;
            double D3   = 0.007855;
            double D4   = 0.5121;
            double D5   = 3.6620;
            double D6   = 0.02307;
            double D7   = 3.2540;
            return inv_pp_pbar_CM__Duperray( s, E_pbar, pT_pbar, D1, D2, D3, D4, D5, D6, D7 );
        }


        //! Parametrization of the invariant antiproton production (pbar) crosssection from p A in CMF
        /*!
         *  Taken from:     R.P. Duperray, et al.; 2003;
         *                  Parameterization of the antiproton inclusive production cross section on nuclei;
         *                  DOI: 10.1103/PhysRevD.68.094017
         *                  Eq. 4
         *
         *  \f[  E_{\bar{p}} \frac{ d^3 \sigma_{pA}^{(\bar{p})} }{d^3 p_{\bar{p}} }  (s, E_{\bar{p}}, p_{T,\bar{p}})  \f]
         *
         *  The CMF referrs to the nucleon nucleon system!
         *
         *  \param double s         CM energy.
         *  \param doulbe E_pbar    Energy of the produced antiproton in CMF
         *  \param doulbe pT_pbar   Transverse momentum of the produced antiproton in CMF.
         *  \param int    A         Mass Number
         *  \param doulbe m_A       Mass of particle A
         *  \param doulbe Ci        i=(1,...,10), parameters
         * */
        static double inv_pA_pbar_CM__Duperray(double s, double E_pbar, double pT_pbar, int A, double m_A,
                                               double C1,
                                               double C2,
                                               double C3,
                                               double C4,
                                               double C5,
                                               double C6,
                                               double C7,
                                               double C8,
                                               double C9,
                                               double C10
                                               ){
            E_pbar = fabs(E_pbar);
            if (s<16*fMass_proton)
                return 0;
            if ( pow(pT_pbar, 2.) > pow(E_pbar, 2.) - pow(fMass_proton, 2.) )
                return 0;
            double E_pbar_Max   =   ( s+fMass_proton*fMass_proton-(2*fMass_proton+m_A/A)*(2*fMass_proton+m_A/A) )/2./sqrt( s );
            double x_R          =   E_pbar/E_pbar_Max;
            if ( x_R > 1. )
                return  0.;
            double E_inc        =   (  s + m_A*m_A/A/A + fMass_proton*fMass_proton )/( 2.*m_A/A );
            double sigma_0      =   45*pow(A, 0.7)*(  1+0.016*sin( 5.3-2.63*log(A) )  );   //Cross section normalization in mbarn
            double sigma_in     =   sigma_0 * (  1. - 0.62 * exp( - E_inc*1000./200.  ) * sin( 10.9/pow( E_inc*1000, 0.28 ) )  );
            double invCsCM      =   sigma_in *
            pow( A      , C1*log( sqrt(s)/C2)*pT_pbar  ) *
            pow( 1 - x_R, C3*log( sqrt(s)   )          ) *
            exp(-C4 * x_R)*
            (   C5 * pow( s, C6 /2. ) * exp( -C7 *pT_pbar                 ) +
                C8 * pow( s, C9 /2. ) * exp( -C10*pT_pbar*pT_pbar         )
                 );
            if (invCsCM!=invCsCM) return 0;
            return invCsCM;
        }
        
        
        //! Parametrization of the invariant antiproton production (pbar) crosssection from p A in CMF
        /*!
         *  Taken from:     R.P. Duperray, et al.; 2003;
         *                  Parameterization of the antiproton inclusive production cross section on nuclei;
         *                  DOI: 10.1103/PhysRevD.68.094017
         *                  Eq. 4 (with best Fit values)
         *
         *  \f[  E_{\bar{p}} \frac{ d^3 \sigma_{pA}^{(\bar{p})} }{d^3 p_{\bar{p}} }  (s, E_{\bar{p}}, p_{T,\bar{p}})  \f]
         *
         *  The CMF referrs to the nucleon nucleon system!
         *
         *  \param double s         CM energy.
         *  \param doulbe E_pbar    Energy of the produced antiproton in CMF
         *  \param doulbe pT_pbar   Transverse momentum of the produced antiproton in CMF.
         *  \param int    A         Mass Number
         *  \param doulbe m_A       Mass of particle A
         * */
        static double inv_pA_pbar_CM__Duperray(double s, double E_pbar, double pT_pbar, int A=1, double m_A=fMass_proton){
            double C1  = 0.16990;
            double C2  = 10.28;
            double C3  = 2.269;
            double C4  = 3.707;
            double C5  = 0.009205;
            double C6  = 0.4812;
            double C7  = 3.3600;
            double C8  = 0.06394;
            double C9  = -0.1824;
            double C10 = 2.4850;
            //return C1;
            return inv_pA_pbar_CM__Duperray(s, E_pbar, pT_pbar, A, m_A, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10);
        }
        
        
        static double inv_pHe_pbar_CM__Duperray(double s, double E_pbar, double pT_pbar){
            return inv_pA_pbar_CM__Duperray(s, E_pbar, pT_pbar, 4);
        }
        static double inv_HeHe_pbar_CM__Duperray(double s, double E_pbar, double pT_pbar){
            return inv_pA_pbar_CM__fromDuperray(s, E_pbar, pT_pbar, inv_pHe_pbar_CM__Duperray, 4);
        }
        
////////////////////////
//  Anderson
////////////////////////

        
        //! Parametrization of the invariant proton scattering crosssection from pp in CMF
        /*!
         *  Taken from:     Anderson, et al.; 1967;
         *                  PROTON AND PION SPECTRA FROM PROTON-PROTON INTERACTIONS AT 10, 20, AND 30 BeV/c*;
         *                  DOI: https://doi.org/10.1103/PhysRevLett.19.198
         *
         *
         *  \param double s         CM energy.
         *  \param doulbe E_p       Energy of the scattered proton in CMF
         *  \param doulbe pT_p      Transverse momentum of the scattered proton in CMF.
         * */
        static double inv_pp_p_CM__Anderson(double s, double E_p, double pT_p){
            E_p = fabs(E_p);
            if (s<4*E_p*E_p)
                return 0;
            if ( pow(pT_p, 2.) > pow(E_p, 2.) - pow(fMass_proton, 2.) )
                return 0;
            return E_p/2./M_PI*pT_p*610*exp(-pT_p/0.166);
        }

////////////////////////
//  Dbar Cross Section
////////////////////////

        //! Parametrization of the invariant antideuteron production (Dbar) crosssection from pp in CMF
        /*!
         *  Calculation accoriding to:
         *                  Chardonnet, et al.; 1997;
         *                  The production of anti-matter in our galaxy;
         *                  DOI: 10.1016/S0370-2693(97)00870-8
         *
         *
         *  \param double s                 CM energy.
         *  \param doulbe E_dbar            Energy of the produced antineutron in CMF
         *  \param doulbe pT_dbar           Transverse momentum of the produced antineutron in CMF.
         *  \param doulbe p_coalescence     Coalescence momentum (as defined in the upper reference).
         *  \param double (*tot_pp)(double)
         *                                  Total pp CS.
         *  \param double (*inv_pp_pbar)(double, doulbe, double)
         *                                  Invariant pbar production CS in CMF.
         *  \param double (*inv_pp_nbar)(double, doulbe, double)
         *                                  Invariant nbar production CS in CMF.
         * */
        static double inv_pp_Dbar_CM( double s, double E_dbar, double pT_dbar,
                                                              double p_coalescence,
                                                              double (*tot_pp)(double),
                                                              double (*inv_pp_pbar)(double, double, double),
                                                              double (*inv_pp_nbar)(double, double, double) ){
            E_dbar = fabs(E_dbar);
            double sq_sReduced = sqrt(s) - E_dbar;
            if (sq_sReduced<0) {
                return 0;
            }
            double sReduced     = sq_sReduced*sq_sReduced;
            
            double pL_dbar = sqrt( pow(E_dbar,2) - pow(fMass_deuteron,2) - pow(pT_dbar,2) );
            if(pL_dbar!=pL_dbar){
                return 0;
            }
            double E_pbar = sqrt( pow(fMass_proton, 2) + pow(pT_dbar/2.,2) + pow(pL_dbar/2.,2) );
            double E_nbar = sqrt( pow(fMass_neutron,2) + pow(pT_dbar/2.,2) + pow(pL_dbar/2.,2) );
            
            
            double CS           = fMass_deuteron/fMass_proton/fMass_neutron *
            (4./3. * M_PI * pow(p_coalescence,3)) *
            1./2. / tot_pp(s)*
            (
                inv_pp_pbar  (s,        E_pbar, pT_dbar/2) *
                inv_pp_nbar  (sReduced, E_nbar, pT_dbar/2) +
                inv_pp_nbar  (s,        E_nbar, pT_dbar/2) *
                inv_pp_pbar  (sReduced, E_pbar, pT_dbar/2)
            );
            
            funOut( ppCSParametrizations::inv_ppbar_Dbar_CM                     )
            varOut(E_dbar                                                       )
            varOut(pT_dbar                                                      )
            varOut(E_dbar*E_dbar-4*fMass_proton*fMass_proton-pT_dbar*pT_dbar    )
            varOut(s                                                            )
            varOut(sReduced                                                     )
            varOut(inv_pp_pbar  (s,        E_dbar/2, pT_dbar/2)                 )
            varOut(inv_pp_nbar  (sReduced, E_dbar/2, pT_dbar/2)                 )
            varOut(inv_pp_nbar  (s,        E_dbar/2, pT_dbar/2)                 )
            varOut(inv_pp_pbar  (sReduced, E_dbar/2, pT_dbar/2)                 )
            
            return CS;
            
        }
        
        
        //! Parametrization of the invariant antideuteron production (Dbar) crosssection from ppbar in CMF
        /*!
         *  Calculation accoriding to:
         *                  Chardonnet, et al.; 1997;
         *                  The production of anti-matter in our galaxy;
         *                  DOI: 10.1016/S0370-2693(97)00870-8
         *
         *
         *  \param double s                 CM energy.
         *  \param doulbe E_dbar            Energy of the produced antineutron in CMF
         *  \param doulbe pT_dbar           Transverse momentum of the produced antineutron in CMF.
         *  \param doulbe p_coalescence     Coalescence momentum (as defined in the upper reference).
         *  \param double (*tot_ppbar)(double)
         *                                  Total ppbar CS.
         *  \param double (*inv_ppbar_pbar)(double, doulbe, double)
         *                                  Invariant pbar CS in CMF.
         *  \param double (*inv_ppbar_nbar)(double, doulbe, double)
         *                                  Invariant nbar production CS in CMF.
         * */
        static double inv_ppbar_Dbar_CM( double s, double E_dbar, double pT_dbar,
                                        double p_coalescence,
                                        double (*tot_ppbar)(double),
                                        double (*inv_ppbar_pbar)(double, double, double),
                                        double (*inv_ppbar_nbar)(double, double, double) ){
            E_dbar = fabs(E_dbar);
            double sq_sReduced = sqrt(s) - E_dbar;
            if (sq_sReduced<0) {
                return 0;
            }
            double sReduced     = sq_sReduced*sq_sReduced;
            
            double pL_dbar = sqrt( pow(E_dbar,2) - pow(fMass_deuteron,2) - pow(pT_dbar,2) );
            if(pL_dbar!=pL_dbar){
                return 0;
            }
            double E_pbar = sqrt( pow(fMass_proton, 2) + pow(pT_dbar/2.,2) + pow(pL_dbar/2.,2) );
            double E_nbar = sqrt( pow(fMass_neutron,2) + pow(pT_dbar/2.,2) + pow(pL_dbar/2.,2) );
            
            
            funOut( ppCSParametrizations::inv_ppbar_Dbar_CM )
            varOut(s)
            varOut(sReduced)
            
            double CS           = fMass_deuteron/fMass_proton/fMass_neutron *
            (4./3. * M_PI * pow(p_coalescence,3)) *
            1./2. / tot_ppbar(s)*
            (
             inv_ppbar_pbar  (s,        E_pbar, pT_dbar/2) *
             inv_ppbar_nbar  (sReduced, E_nbar, pT_dbar/2) +
             inv_ppbar_nbar  (s,        E_nbar, pT_dbar/2) *
             inv_ppbar_pbar  (sReduced, E_pbar, pT_dbar/2)
            );
            return CS;
            
        }
        
        
        
        
        //! Parametrization of the invariant antihelium production (Hebar) crosssection from pp in CMF
        /*!
         *  Calculation accoriding to:
         *                  Chardonnet, et al.; 1997;
         *                  The production of anti-matter in our galaxy;
         *                  DOI: 10.1016/S0370-2693(97)00870-8
         *
         *
         *  \param double s                 CM energy.
         *  \param doulbe E_Hebar           Energy of the produced antineutron in CMF
         *  \param doulbe pT_Hebar          Transverse momentum of the produced antineutron in CMF.
         *  \param doulbe p_coalescence     Coalescence momentum (as defined in the upper reference).
         *  \param double (*tot_pp)(double)
         *                                  Total pp CS.
         *  \param double (*inv_pp_pbar)(double, doulbe, double)
         *                                  Invariant pbar production CS in CMF.
         *  \param double (*inv_pp_nbar)(double, doulbe, double)
         *                                  Invariant nbar production CS in CMF.
         * */
        static double inv_pp_Hebar_CM( double s, double E_Hebar, double pT_Hebar,
                                     double p_coalescence,
                                     double (*tot_pp)(double),
                                     double (*inv_pp_pbar)(double, double, double),
                                     double (*inv_pp_nbar)(double, double, double) ){
            E_Hebar = fabs(E_Hebar);
            double pL_Hebar = sqrt( pow(E_Hebar,2) - pow(fMass_helium3,2) - pow(pT_Hebar,2) );
            if(pL_Hebar!=pL_Hebar){
                return 0;
            }
            double E_pbar = sqrt( pow(fMass_proton, 2) + pow(pT_Hebar/3.,2) + pow(pL_Hebar/3.,2) );
            double E_nbar = sqrt( pow(fMass_neutron,2) + pow(pT_Hebar/3.,2) + pow(pL_Hebar/3.,2) );
            
            double sq_sReduced = sqrt(s) - 2*E_pbar;
            if (sq_sReduced<0) {
                return 0;
            }
            double sReduced     = sq_sReduced*sq_sReduced;
            
            sq_sReduced = sqrt(s) - 4*E_pbar;
            if (sq_sReduced<0) {
                return 0;
            }
            double sReduced_2     = sq_sReduced*sq_sReduced;
            
            
            double CS           = fMass_helium3/fMass_proton/fMass_neutron/fMass_neutron *
            pow(  (4./3. * M_PI * pow(p_coalescence,3)) / tot_pp(s), 2  )   * 1./3. *
            (
             inv_pp_pbar  (s,           E_pbar, pT_Hebar/3.) *
             inv_pp_nbar  (sReduced,    E_nbar, pT_Hebar/3.) *
             inv_pp_nbar  (sReduced_2,  E_nbar, pT_Hebar/3.) +
             
             inv_pp_pbar  (sReduced_2,  E_pbar, pT_Hebar/3.) *
             inv_pp_nbar  (sReduced,    E_nbar, pT_Hebar/3.) *
             inv_pp_nbar  (s,           E_nbar, pT_Hebar/3.) +
             
             inv_pp_pbar  (sReduced,    E_pbar, pT_Hebar/3.) *
             inv_pp_nbar  (s,           E_nbar, pT_Hebar/3.) *
             inv_pp_nbar  (sReduced_2,  E_nbar, pT_Hebar/3.)
             );
            
            return CS;
            
        }
        
        
        static double inv_ppbar_Hebar_CM( double s, double E_Hebar, double pT_Hebar,
                                         double p_coalescence,
                                         double (*tot_ppbar)(double),
                                         double (*inv_ppbar_pbar)(double, double, double),
                                         double (*inv_ppbar_nbar)(double, double, double) ){
            E_Hebar = fabs(E_Hebar);
            double pL_Hebar = sqrt( pow(E_Hebar,2) - pow(fMass_helium3,2) - pow(pT_Hebar,2) );
            if(pL_Hebar!=pL_Hebar){
                return 0;
            }
            double E_pbar = sqrt( pow(fMass_proton, 2) + pow(pT_Hebar/3.,2) + pow(pL_Hebar/3.,2) );
            double E_nbar = sqrt( pow(fMass_neutron,2) + pow(pT_Hebar/3.,2) + pow(pL_Hebar/3.,2) );
            
            double sq_sReduced = sqrt(s) - 2*E_pbar;
            if (sq_sReduced<0) {
                return 0;
            }
            double sReduced     = sq_sReduced*sq_sReduced;
            
            sq_sReduced = sqrt(s) - 4*E_pbar;
            if (sq_sReduced<0) {
                return 0;
            }
            double sReduced_2     = sq_sReduced*sq_sReduced;
            
            
            double CS           = fMass_helium3/fMass_proton/fMass_neutron/fMass_neutron *
            pow(  (4./3. * M_PI * pow(p_coalescence,3)) / tot_ppbar(s), 2  )   * 1./3. *
            (
             inv_ppbar_pbar  (s,           E_pbar, pT_Hebar/3.) *
             inv_ppbar_nbar  (sReduced,    E_nbar, pT_Hebar/3.) *
             inv_ppbar_nbar  (sReduced_2,  E_nbar, pT_Hebar/3.) +
             
             inv_ppbar_pbar  (sReduced_2,  E_pbar, pT_Hebar/3.) *
             inv_ppbar_nbar  (sReduced,    E_nbar, pT_Hebar/3.) *
             inv_ppbar_nbar  (s,           E_nbar, pT_Hebar/3.) +
             
             inv_ppbar_pbar  (sReduced,    E_pbar, pT_Hebar/3.) *
             inv_ppbar_nbar  (s,           E_nbar, pT_Hebar/3.) *
             inv_ppbar_nbar  (sReduced_2,  E_nbar, pT_Hebar/3.)
             );
            
            return CS;
            
        }
        
        
        
        
        
        
    };
    

}



//end PPCROSSSECTIONPARAMETRIZATIONS_H

#endif
