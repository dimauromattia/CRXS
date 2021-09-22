#ifndef CROSSSECTIONTRANSFORMATIONS_H
#define CROSSSECTIONTRANSFORMATIONS_H

#include <iostream>
#include <math.h>
#include <definitions.h>
#include <ppCrossSectionParametrizations.h>

namespace CRACS {
    
    /*! \brief This class contains the transformation functions from invariant CM cross section to the LAB frame.
     *
     *  Lab frame:     One proton is at rest and the other one has incident energy E_p
     *
     *
     */

    class CSTransformations{
        
    public:
        
        static bool convert_LAB_to_CM( const double T_p_LAB, const double T_pbar_LAB, const double eta_LAB, double &s, double &x_R, double &pT_pbar ){
            
            s                    = 4*fMass_proton*fMass_proton + 2 * T_p_LAB * fMass_proton;
            
            double p_pbar_LAB    = sqrt(  T_pbar_LAB*(T_pbar_LAB+2*fMass_proton)  );
            double E_p_LAB       = T_p_LAB+fMass_proton;
            double E_pbar_LAB    = T_pbar_LAB+fMass_proton;
            
            double beta          = sqrt(E_p_LAB - fMass_proton)/sqrt(E_p_LAB + fMass_proton);
            double gamma         = 1./sqrt(1 - beta*beta);
            double gammabeta     = gamma * beta;
            
            double E_pbar        = gamma * E_pbar_LAB - gammabeta * p_pbar_LAB*tanh(eta_LAB);
            double E_pbar_Max    = ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
            
            pT_pbar              = p_pbar_LAB/cosh(eta_LAB);
            x_R                  = E_pbar/E_pbar_Max;
           
            return true;
        };
        
        
        static bool convert_CM_to_LAB( const double s, const double xR, const double pT, double &T_p_LAB, double &T_pbar_LAB, double &eta_LAB, bool pL_pos=true ){
            funOut( CSTransformations::convert_CM_to_LAB )
            T_p_LAB             =   (s-4*fMass_proton*fMass_proton)/2/fMass_proton;
            
            if(T_p_LAB<0) return false;
            
            double gamma        =   sqrt(s)/2/fMass_proton;
            double gammabeta    =   sqrt(s-4*fMass_proton*fMass_proton)/2/fMass_proton;
            
            varOut(gamma)
            varOut(gammabeta)

            double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
            double E_pbar       =   xR*E_pbar_Max;
            double p_pbar       =   sqrt(E_pbar*E_pbar-fMass_proton*fMass_proton);
            
            if(E_pbar<fMass_proton) return false;
            
            
            varOut(E_pbar)
            varOut(p_pbar)
            
            if(p_pbar<pT) return false;
            double pL_pbar      =   sqrt(p_pbar*p_pbar-pT*pT);
            if (!pL_pos) {
                pL_pbar*=-1;
            }
            
            double E_pbar_LAB   = gamma * E_pbar + gammabeta * pL_pbar;
            double pL_pbar_LAB  = gammabeta * E_pbar + gamma * pL_pbar;
            
            varOut(E_pbar_LAB)
            varOut(pL_pbar_LAB)
            
            double p_pbar_LAB   = sqrt(E_pbar_LAB*E_pbar_LAB-fMass_proton*fMass_proton);
            
            double cos_theta_LAB= pL_pbar_LAB/p_pbar_LAB;
            
            
            
            T_pbar_LAB          = E_pbar_LAB-fMass_proton;
            eta_LAB             = atanh(cos_theta_LAB);
            
            varOut(cos_theta_LAB)
            varOut(eta_LAB)
            
            return true;
            
        };
        
        static bool convert_CMxF_to_LAB( const double s, const double xF, const double pT, double &T_p_LAB, double &T_pbar_LAB, double &eta_LAB ){
            funOut( CSTransformations::convert_CM_to_LAB )
            T_p_LAB             =   (s-4*fMass_proton*fMass_proton)/2/fMass_proton;
            
            if(T_p_LAB<0) return false;
            
            double gamma        =   sqrt(s)/2/fMass_proton;
            double gammabeta    =   sqrt(s-4*fMass_proton*fMass_proton)/2/fMass_proton;
            
            varOut(gamma)
            varOut(gammabeta)
            
            double pL           =   xF/2.*sqrt(s);
            double p_pbar       =   sqrt(pT*pT + pL*pL);
            double E_pbar       =   sqrt(p_pbar*p_pbar+fMass_proton*fMass_proton);
            
            if(E_pbar<fMass_proton) return false;
            
            
            varOut(E_pbar)
            varOut(p_pbar)
            
            if(p_pbar<pT) return false;
            
            double pL_pbar      = pL;
            double E_pbar_LAB   = gamma * E_pbar + gammabeta * pL_pbar;
            double pL_pbar_LAB  = gammabeta * E_pbar + gamma * pL_pbar;
            
            varOut(E_pbar_LAB)
            varOut(pL_pbar_LAB)
            
            double p_pbar_LAB   = sqrt(E_pbar_LAB*E_pbar_LAB-fMass_proton*fMass_proton);
            
            double cos_theta_LAB= pL_pbar_LAB/p_pbar_LAB;
            
            
            
            T_pbar_LAB          = E_pbar_LAB-fMass_proton;
            eta_LAB             = atanh(cos_theta_LAB);
            
            varOut(cos_theta_LAB)
            varOut(eta_LAB)
            
            return true;
            
        };
        
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
         *  \param double (*inv_pp_product_CM)(double, doulbe, double) Corresponding invariant particle production CS in CMF.
         * */
        
        static double inv_pp_product_LAB(double E_p_LAB, double E_product_LAB, double cos_theta_product_LAB, double m_product, double (*inv_pp_product_CM)(double, double, double) )
        {
            double s = 2*fMass_proton*fMass_proton + 2 * E_p_LAB * fMass_proton;
            
            //Lorentz Transformation for pbar in CM frame
            double beta          = sqrt(E_p_LAB - fMass_proton)/sqrt(E_p_LAB + fMass_proton);
            double gamma         = 1./sqrt(1 - beta*beta);
            double gammabeta     = gamma * beta;
            double p_product_LAB = sqrt(  pow( E_product_LAB, 2 ) - pow( m_product, 2 )  );
            
            double E_product   =       gamma * E_product_LAB - gammabeta * p_product_LAB*cos_theta_product_LAB;
            double pL_product  = - gammabeta * E_product_LAB + gamma     * p_product_LAB*cos_theta_product_LAB;
            if (pL_product<0) {
                E_product *= -1.;
            }
            double pT_product  = p_product_LAB*sqrt(1-cos_theta_product_LAB*cos_theta_product_LAB);
            
            funOut( CSTransformations::d3p_pp_product_LAB )
            varOut(s)
            varOut(beta)
            varOut(gamma)
            varOut(p_product_LAB)
            varOut(E_product)
            varOut(pT_product)
            
            return inv_pp_product_CM(s, E_product, pT_product);
        }
        
        
        //!  Production cross section in Lab Frame from Hep collison, transformed from pHe!
        /*!
         *  lab frame:     proton is at rest
         *  LAB frame:     helium is at rest
         *
         *  \f[  \frac{ d^3 \sigma_{Hep}^{(X), LAB} }{d^3 p_{X} }  (E_{He,lab}/n, E_{X,lab}, \cos \theta_{T,X,lab})  \f]
         *
         *  \param double E_He/n        Incident helium energy per nucleon in lab frame.
         *  \param doulbe E_product     Energy of the produced particle in lab frame.
         *  \param doulbe pT_product    Transverse momentum of the produced particle in lab frame.
         *  \param doulbe m_product     Mass of the produced particle.
         *  \param double (*inv_pHe_product_LAB)(double, doulbe, double) Corresponding invariant particle production CS in LAB frame.
         * */
        
        static double inv_Hep_product_lab(double En_He_lab, double E_product_lab, double cos_theta_product_lab, double m_product, double (*inv_pHe_product_LAB)(double, double, double) )
        {
            // We denote the Hep frame with lab while using LAB for the pHe frame!
            
            double p_product_lab    = sqrt(E_product_lab*E_product_lab-m_product*m_product);
            double pL_product_lab   = p_product_lab*cos_theta_product_lab;
            double pT_product_lab   = p_product_lab*sqrt(1-cos_theta_product_lab*cos_theta_product_lab);
            
            // Lorentz Transformation for pbar in LAB frame (He at rest!)
            double gamma         = En_He_lab*4./fMass_helium;
            double gammabeta     = sqrt(En_He_lab*En_He_lab*16./fMass_helium/fMass_helium - 1);
            
            double E_p_LAB = fMass_proton*gamma;
            
            double E_product_LAB   =   gamma        * E_product_lab     -   gammabeta   * pL_product_lab;
            double pL_product_LAB  = - gammabeta    * E_product_lab     +   gamma       * pL_product_lab;
            pL_product_LAB *= -1.; // Change of angle definition theta_lab between (pbar, He), theta_LAB is between (pbar, p)
            
            double p_product_LAB   = sqrt( pL_product_LAB*pL_product_LAB + pT_product_lab*pT_product_lab );
            double cos_theta_product_LAB = pL_product_LAB/p_product_LAB;
            
            return inv_pHe_product_LAB( E_p_LAB, E_product_LAB, cos_theta_product_LAB );
            
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
         *  \param double (*inv_pp_product_CM)(double, doulbe, double)
         *                              Corresponding invariant particle production CS in CMF.
         *  \param doulbe m_remainder   Minimum mass of the "remaining particle" (e.g. 3 m_p for pbar production, 5 m_p for Dbar production).
         *  \param int    precision     Number of steps for the integration over cos(theta).
         *  \param int    output        Output level: NO_OUT, WARN_OUT, ALL_OUT. Default: WARN_OUT
         * */
//        static double dE_pp_product_LAB(double E_p_LAB, double E_product_LAB, double m_product, double (*inv_pp_product_CM)(double, double, double), double m_remainder, int precision=100000, int output=WARN_OUT ){
//            
//            
//            double p_product_LAB = sqrt(  pow( E_product_LAB, 2 ) - pow( m_product, 2 )  );
//            
//            
//            double phi_integration              = 2*M_PI;
//            double Jacobian_and_conversion      = p_product_LAB*E_product_LAB;
//            //     Jacobian_determinant         = p_product_Lab*p_product_Lab;
//            //     dp_to_dE_conversion          = E_product_LAB/p_product_Lab;
//            
//            //Calculate maximal cos_theta
//            
//            double s = 2*fMass_proton*fMass_proton + 2 * E_p_LAB * fMass_proton;
//            
//            double beta          = sqrt(E_p_LAB - fMass_proton)/sqrt(E_p_LAB + fMass_proton);
//            double gamma         = 1./sqrt(1 - beta*beta);
//            double gammabeta     = gamma * beta;
//            
//            //double cos_theta_max = (  gamma * E_product_LAB -  (s + m_product*m_product - m_remainder*m_remainder)/2./sqrt(s) ) / (gammabeta * p_product_LAB);
//            double cos_theta_max = (  gamma * E_product_LAB -  (s + m_product*m_product - m_remainder*m_remainder)/2./sqrt(s) ) / (gammabeta * p_product_LAB);
//            
//            funOut( CSTransformations::dE_pp_product_LAB )
//            varOut(cos_theta_max)
//            
//            // if m_remainder < 0 we assume that there is no boundary on cos_theta
//            // we infer the maxmimal cos_theta numerically in a pre-loop
//            if (m_remainder<0){
//                cos_theta_max = -1;
//                double delta_cos_theta  = (1.-cos_theta_max)/precision*10;
//                double cos_theta        = cos_theta_max;
//                
//                for ( ; cos_theta<=1; cos_theta+=delta_cos_theta) {
//                    
//                    double cos_theta_product_LAB = cos_theta;
//                    double E_product   = gamma * E_product_LAB - gammabeta * p_product_LAB*cos_theta_product_LAB;
//                    double pT_product  = p_product_LAB*sqrt(1-cos_theta_product_LAB*cos_theta_product_LAB);
//                    double inv_CS_CM   = inv_pp_product_CM(s, E_product, pT_product);
//                    
//                    if(inv_CS_CM>0)
//                        break;
//                    
//                    cos_theta_max = cos_theta;
//                }
//            }
//            
//            //double E_product_CMS_max = (s + m_product*m_product - m_remainder*m_remainder)/2./sqrt(s);
//            
//            
//            varOut(beta)
//            varOut(gamma)
//            varOut(gammabeta)
//            varOut(s)
//            varOut(E_product_LAB)
//            varOut(p_product_LAB)
//            varOut(cos_theta_max)
//
//            
//            if (cos_theta_max!=cos_theta_max) { // FIXME
//                return 0;
//            }
//            cos_theta_max = std::max(cos_theta_max, -1.);
//            
//            
//            //cos_theta_max = 0.95;
//            
//            //Numerical integration for theta integration, cos(theta) in 'precision' steps
//            double integral_Trapez  = 0;
//            double integral_Upper   = 0;
//            double integral_Lower   = 0;
//            double delta_cos_theta  = (1.-cos_theta_max)/precision;
//            double cos_theta        = cos_theta_max;
//            
//            varOut(delta_cos_theta);
//            
//            if (delta_cos_theta<1e-15) { // FIXME
//                return 0;
//            }
//            
//            
//            double lastCS = 1.0/E_product_LAB *inv_pp_product_LAB(E_p_LAB, E_product_LAB, cos_theta, m_product, inv_pp_product_CM);
//            for ( ; cos_theta<=1; cos_theta+=delta_cos_theta) {
//                
//                double cos_theta_product_LAB = cos_theta;
//                double E_product   = gamma * E_product_LAB - gammabeta * p_product_LAB*cos_theta_product_LAB;
//                double pT_product  = p_product_LAB*sqrt(1-cos_theta_product_LAB*cos_theta_product_LAB);
//                double d3p_pp_product_LAB = 1.0/E_product_LAB * inv_pp_product_CM(s, E_product, pT_product);
//                
//                double CS = d3p_pp_product_LAB;
//                
//                //varOut(E_product)
//                //varOut(E_product_CMS_max)
//                //varOut(CS)
//                //varOut(" ")
//                
//                integral_Trapez  += (lastCS+CS)/2;
//                integral_Upper   += std::max(lastCS, CS);
//                integral_Lower   += std::min(lastCS, CS);
//                
//                lastCS = CS;
//            }
//            
//            integral_Trapez    *=  delta_cos_theta*phi_integration*Jacobian_and_conversion;
//            integral_Upper     *=  delta_cos_theta*phi_integration*Jacobian_and_conversion;
//            integral_Lower     *=  delta_cos_theta*phi_integration*Jacobian_and_conversion;
//            
//            double error = std::max( integral_Upper-integral_Trapez, integral_Trapez-integral_Lower );
//            
//            if(integral_Trapez>0 && error/integral_Trapez>0.01 && output>NO_OUT){
//                std::cout << warning_momentum_integration << std::endl;
//            }
//            if (output==ALL_OUT) {
//                std::cout << "Trapez Sum:      " << integral_Trapez << std::endl;
//                std::cout << "Upper Sum:       " << integral_Upper  << std::endl;
//                std::cout << "Lower Sum:       " << integral_Lower  << std::endl;
//                std::cout << "Relative Error:  " << error/integral_Trapez  << std::endl;
//            }
//            
//            
//            return integral_Trapez;
//        }

        
        //!  Energy differential production cross section in Lab Frame from pp collison, numerical integration in pseudorapidity
        /*!
         *  Lab frame:     One proton is at rest and the other one has incident energy E_p
         *
         *  \f[  \frac{ d \sigma_{pp}^{(X), LAB} }{d E_{X} }  (E_p, E_{X})  \f]
         *
         *  \param double E_p           Incident proton energy in LAB frame.
         *  \param doulbe E_product     Energy of the produced particle in LAB frame.
         *  \param doulbe m_product     Mass of the produced particle.
         *  \param double (*inv_pp_product_CM)(double, doulbe, double)
         *                              Corresponding invariant particle production CS in CMF.
         *  \param doulbe m_remainder   Minimum mass of the "remaining particle" (e.g. 3 m_p for pbar production, 5 m_p for Dbar production).
         *  \param int    precision     Number of steps for the integration over cos(theta).
         *  \param int    output        Output level: NO_OUT, WARN_OUT, ALL_OUT. Default: WARN_OUT
         * */
        static double dE_pp_product_LAB_intEta(double E_p_LAB, double E_product_LAB, double m_product, double (*inv_pp_product_CM)(double, double, double), int precision=100000, int output=WARN_OUT ){
            
            double eta_max=20.0;
            
            double p_product_LAB = sqrt(  pow( E_product_LAB, 2 ) - pow( m_product, 2 )  );
            if(p_product_LAB!=p_product_LAB)
            return 0;
            
            double phi_integration              = 2*M_PI;
            double Jacobian_and_conversion      = p_product_LAB;
            //     Jacobian_determinant         = p_product_Lab*p_product_Lab;
            //     dp_to_dE_conversion          = E_product_LAB/p_product_Lab;
            
            //Calculate maximal cos_theta
            
            double s = 2*fMass_proton*fMass_proton + 2 * E_p_LAB * fMass_proton;
            
            double beta          = sqrt(E_p_LAB - fMass_proton)/sqrt(E_p_LAB + fMass_proton);
            double gamma         = 1./sqrt(1 - beta*beta);
            double gammabeta     = gamma * beta;
            
            
            funOut( CSTransformations::dE_pp_product_LAB )
            
            
            //Numerical integration for theta integration, eta in 'precision' steps
            double integral             = 0;
            double integral_upper       = 0;
            double integral_lower       = 0;
            double contribution_first   = 0;
            double contribution_last    = 0;
            double contribution_max     = 0;
            double delta_eta            = (eta_max)/precision;
            double eta                  = 0;
            
            varOut(delta_eta)
           
            double last_CS              = 0;
            for ( ; eta<=eta_max+delta_eta/2; eta+=delta_eta) {
                
                double cos_theta_product_LAB = tanh(eta);
                
                double E_product   =       gamma * E_product_LAB - gammabeta * p_product_LAB*cos_theta_product_LAB;
                if(E_product_LAB<=m_product)
                continue;
                double pL_product  = - gammabeta * E_product_LAB + gamma     * p_product_LAB*cos_theta_product_LAB;
                
                if (pL_product<0) {
                    E_product *= -1.;
                }
                
                double pT_product  = p_product_LAB/cosh(eta);
                double d3p_pp_product_LAB = inv_pp_product_CM(s, E_product, pT_product);
                
                double CS = d3p_pp_product_LAB * pow(  cosh(eta), -2  );
                
//                if (CS!=CS){
//                    out("-------")
//                    out(cos_theta_product_LAB)
//                    out(E_product_LAB)
//                    out(p_product_LAB)
//                    out(E_p_LAB)
//                    out(s)
//                    out(E_product)
//                    out(pT_product)
//                    out(d3p_pp_product_LAB)
//                    
//                }
                
                integral        += CS;
                integral_upper  += std::max(last_CS, CS);
                integral_lower  += std::min(last_CS, CS);
                last_CS = CS;
                
                varOut("")
                varOut(CS)
                varOut(E_product)
                varOut(pT_product)
                
                if (eta<delta_eta/2) {
                    contribution_first = CS;
                }
                if (eta>eta_max-delta_eta/2) {
                    contribution_last  = CS;
                }
                if (CS>contribution_max) {
                    contribution_max=CS;
                }
                
                
            }
            
            double error     = std::max( (integral_upper-integral)/integral, (integral-integral_lower)/integral );
            double error_lim = std::max( contribution_first/contribution_max, contribution_last/contribution_max );
            
            if(integral>0 && error>0.01 && output>NO_OUT){
                std::cout << warning_eta_integration << " Adjust bin size." << std::endl;
            }
            if(integral>0 && error_lim>0.001 && output>NO_OUT){
                std::cout << warning_eta_integration << " Adjust limits."   << std::endl;
            }
            if (output==ALL_OUT) {
                out(integral)
                out(contribution_first)
                out(contribution_last)
                out(error_lim)
                out(error)
            }
            
            
            integral            *=  delta_eta*phi_integration*Jacobian_and_conversion;
            return integral;
        }
//        static double dE_pp_product_LAB_intEta(double E_p_LAB, double E_product_LAB, double m_product, double (*inv_pp_product_CM)(double, double, double, double), int precision=100000, int output=WARN_OUT ){
//            
//            double eta_max=20.0;
//            
//            double p_product_LAB = sqrt(  pow( E_product_LAB, 2 ) - pow( m_product, 2 )  );
//            if(p_product_LAB!=p_product_LAB)
//            return 0;
//            
//            
//            double phi_integration              = 2*M_PI;
//            double Jacobian_and_conversion      = p_product_LAB;
//            //     Jacobian_determinant         = p_product_Lab*p_product_Lab;
//            //     dp_to_dE_conversion          = E_product_LAB/p_product_Lab;
//            
//            //Calculate maximal cos_theta
//            
//            double s = 2*fMass_proton*fMass_proton + 2 * E_p_LAB * fMass_proton;
//            
//            double beta          = sqrt(E_p_LAB - fMass_proton)/sqrt(E_p_LAB + fMass_proton);
//            double gamma         = 1./sqrt(1 - beta*beta);
//            double gammabeta     = gamma * beta;
//            
//            
//            funOut( CSTransformations::dE_pp_product_LAB )
//            
//            
//            //Numerical integration for theta integration, eta in 'precision' steps
//            double integral             = 0;
//            double integral_upper       = 0;
//            double integral_lower       = 0;
//            double contribution_first   = 0;
//            double contribution_last    = 0;
//            double contribution_max     = 0;
//            double delta_eta            = (eta_max)/precision;
//            double eta                  = 0;
//            
//            varOut(delta_eta)
//            
//            double last_CS              = 0;
//            for ( ; eta<=eta_max+delta_eta/2; eta+=delta_eta) {
//                
//                double cos_theta_product_LAB = tanh(eta);
//                
//                double E_product   =   gamma     * E_product_LAB - gammabeta * p_product_LAB*cos_theta_product_LAB;
//                if(E_product_LAB<=m_product)
//                continue;
//                double pL_product  = - gammabeta * E_product_LAB + gamma     * p_product_LAB*cos_theta_product_LAB;
//                double pT_product  = p_product_LAB/cosh(eta);
//                double d3p_pp_product_LAB = inv_pp_product_CM(s, E_product, pT_product, pL_product);
//                
//                double CS = d3p_pp_product_LAB * pow(  cosh(eta), -2  );
//                
////                if (CS!=CS){
////                    out("-------")
////                    out(cos_theta_product_LAB)
////                    out(E_product_LAB)
////                    out(p_product_LAB)
////                    out(E_p_LAB)
////                    out(s)
////                    out(E_product)
////                    out(pT_product)
////                    out(d3p_pp_product_LAB)
////                    
////                }
//                
//                integral        += CS;
//                integral_upper  += std::max(last_CS, CS);
//                integral_lower  += std::min(last_CS, CS);
//                last_CS = CS;
//                
//                varOut("")
//                varOut(CS)
//                varOut(E_product)
//                varOut(pT_product)
//                
//                
//                if (eta<delta_eta/2) {
//                    contribution_first = CS;
//                }
//                if (eta>eta_max-delta_eta/2) {
//                    contribution_last  = CS;
//                }
//                if (CS>contribution_max) {
//                    contribution_max=CS;
//                }
//                
//                
//            }
//            
//            double error     = std::max( (integral_upper-integral)/integral, (integral-integral_lower)/integral );
//            double error_lim = std::max( contribution_first/contribution_max, contribution_last/contribution_max );
//            
//            if(integral>0 && error>0.01 && output>NO_OUT){
//                std::cout << warning_eta_integration << " Adjust bin size." << std::endl;
//            }
//            if(integral>0 && error_lim>0.001 && output>NO_OUT){
//                std::cout << warning_eta_integration << " Adjust limits."   << std::endl;
//            }
//            if (output==ALL_OUT) {
//                out(integral)
//                out(contribution_first)
//                out(contribution_last)
//                out(error_lim)
//                out(error)
//            }
//            
////            if (integral!=integral){
////                out(integral)
////            }
//            
//            integral            *=  delta_eta*phi_integration*Jacobian_and_conversion;
//            return integral;
//        }
        
        
        
        //!  Energy differential production cross section in lab Frame from Hep collison, numerical integration in pseudorapidity
        /*!
         *  lab frame:     proton is at rest
         *  LAB frame:     helium is at rest
         *
         *  \f[  \frac{ d \sigma_{Hep}^{(X), lab} }{d E_{X,lab} }   ( E_{He,lab}/n, E_{X,lab} )  \f]
         *
         *
         *  \param double E_He/n        Incident helium energy per nucleon in lab frame.
         *  \param doulbe E_product     Energy of the produced particle in lab frame.
         *  \param doulbe m_product     Mass of the produced particle.
         *  \param double (*inv_pHe_product_CM)(double, doulbe, double)
         *                              Corresponding invariant particle production CS in CM frame (nucleon, nucleon).
         *  \param doulbe m_remainder   Minimum mass of the "remaining particle" (e.g. 3 m_p for pbar production, 5 m_p for Dbar production).
         *  \param int    precision     Number of steps for the integration over cos(theta).
         *  \param int    output        Output level: NO_OUT, WARN_OUT, ALL_OUT. Default: WARN_OUT
         * */
        static double dE_Hep_product_LAB_intEta(double En_He_lab, double E_product_lab, double m_product, double (*inv_pA_product_CM)(double, double, double), int precision=100000, int output=WARN_OUT ){
            
            
            double phi_integration              = 2*M_PI;
            double eta_max=20.0;
            
            
            
            //
            //  Transformation form Hep to pHe
            //

            // We denote the Hep frame with lab while using LAB for the pHe frame!
            
            double p_product_lab = sqrt(  pow( E_product_lab, 2 ) - pow( m_product, 2 )  );
            if(p_product_lab!=p_product_lab)
            return 0;
            
            double Jacobian_and_conversion      = p_product_lab;
            //     Jacobian_determinant         = p_product_lab*p_product_lab;
            //     dp_to_dE_conversion          = E_product_lab/p_product_lab;
            
            
            // Lorentz Transformation for pbar in LAB frame (He at rest!)
            double gamma_step1         = En_He_lab*4./fMass_helium;
            double gammabeta_step1     = sqrt(En_He_lab*En_He_lab*16./fMass_helium/fMass_helium - 1);
            
            double E_p_LAB = fMass_proton*gamma_step1;
            
            
            
            //
            //  Transformation form pHe to CM
            //
            
            double s = 2*fMass_proton*fMass_proton + 2 * E_p_LAB * fMass_proton;
            
            
            double beta          = sqrt(E_p_LAB - fMass_proton)/sqrt(E_p_LAB + fMass_proton);
            double gamma         = 1./sqrt(1 - beta*beta);
            double gammabeta     = gamma * beta;
            
            
            
            //Numerical integration for theta integration, eta in 'precision' steps
            double integral             = 0;
            double integral_upper       = 0;
            double integral_lower       = 0;
            double contribution_first   = 0;
            double contribution_last    = 0;
            double contribution_max     = 0;
            double delta_eta            = (eta_max)/precision;
            double eta                  = 0;
            
            varOut(delta_eta)
            
            double last_CS              = 0;
            for ( ; eta<=eta_max+delta_eta/2; eta+=delta_eta) {
                
                double cos_theta_product_lab    = tanh(eta);
                double pT_product_lab           = p_product_lab/cosh(eta);
                double pL_product_lab           = p_product_lab*tanh(eta);
                
                
                double E_product_LAB   =   gamma_step1        * E_product_lab     -   gammabeta_step1   * pL_product_lab;
                double pL_product_LAB  = - gammabeta_step1    * E_product_lab     +   gamma_step1       * pL_product_lab;
                
                
                if(E_product_LAB<=m_product)
                continue;
                
                
                pL_product_LAB *= -1.; // Change of angle definition theta_lab between (pbar, He), theta_LAB is between (pbar, p)
                
                double p_product_LAB   = sqrt( E_product_LAB*E_product_LAB - m_product*m_product );
                double cos_theta_product_LAB = pL_product_LAB/p_product_LAB;
                if (cos_theta_product_LAB!=cos_theta_product_LAB)
                continue;
                
                
                
                
                double E_product   =       gamma * E_product_LAB - gammabeta * pL_product_LAB;
                double pT_product  = pT_product_lab;
                double pL_product  = - gammabeta * E_product_LAB + gamma     * pL_product_LAB;
                
                if (pL_product<0) {
                    E_product *= -1.;
                }
                
                double d3p_pp_product_LAB = inv_pA_product_CM(s, E_product, pT_product);
                
                double CS = d3p_pp_product_LAB * pow(  cosh(eta), -2  )*Jacobian_and_conversion;
                
//                if (CS!=CS){
//                    out("-------")
//                    out(E_product_lab)
//                    out(pT_product_lab)
//                    out(cos_theta_product_lab)
//                    out(cos_theta_product_LAB)
//                    out(E_product_LAB)
//                    out(p_product_LAB)
//                    out(En_He_lab)
//                    out(E_p_LAB)
//                    out(s)
//                    out(E_product)
//                    out(pT_product)
//                    out(d3p_pp_product_LAB)
//                    
//                }
//
//                    
//                    
//                    double E_pbar = E_product;
//                    double pT_pbar = pT_product;
//                    int A = 4;
//                    double m_A = fMass_helium;
//                    
//                    double C1  = 0.16990;
//                    double C2  = 10.28;
//                    
//                    double E_inc        =   (  s + m_A*m_A/A/A + fMass_proton*fMass_proton )/( 2.*m_A/A );
//                    double sigma_0      =   45*pow(A, 0.7)*(  1+0.016*sin( 5.3-2.63*log(A) )  );
//                    double sigma_in     =   sigma_0 * (  1. - 0.62 * exp( - E_inc*1000./200.  ) * sin( 10.9/pow( E_inc*1000, 0.28 ) )  );
//                    
//                    double E_inc_1      =   (  s  + 2*fMass_proton*fMass_proton )/( 2.*fMass_proton );
//                    double sigma_0_1    =   45*(  1+0.016*sin( 5.3 )  );
//                    double sigma_in_1   =   sigma_0_1 * (  1. - 0.62 * exp( - E_inc_1*1000./200.  ) * sin( 10.9/pow( E_inc_1*1000, 0.28 ) )  );
//                    
//                    double factor       =   pow( A, C1*log( sqrt(s)/C2)*pT_pbar  )*sigma_in/sigma_in_1;
//                    
//                    out(E_inc)
//                    out(E_inc_1)
//                    out(sigma_0)
//                    out(sigma_0_1)
//                    out(sigma_in)
//                    out(sigma_in_1)
//                    
//                    out(C1*log( sqrt(s)/C2)*pT_pbar)
//                    out(log( sqrt(s)/C2))
//                    out(pT_pbar)
//                    
//                    out(factor)
//                    CS = 1e-90;
//                    
//                }
                
                integral        += CS;
                integral_upper  += std::max(last_CS, CS);
                integral_lower  += std::min(last_CS, CS);
                last_CS = CS;
                
                varOut("")
                varOut(CS)
                varOut(E_product)
                varOut(pT_product)
                
                if (eta<delta_eta/2) {
                    contribution_first = CS;
                }
                if (eta>eta_max-delta_eta/2) {
                    contribution_last  = CS;
                }
                if (CS>contribution_max) {
                    contribution_max=CS;
                }
                
                
            }
            
            double error     = std::max( (integral_upper-integral)/integral, (integral-integral_lower)/integral );
            double error_lim = std::max( contribution_first/contribution_max, contribution_last/contribution_max );
            
            if(integral>0 && error>0.01 && output>NO_OUT){
                std::cout << warning_eta_integration << " Adjust bin size." << std::endl;
            }
            if(integral>0 && error_lim>0.001 && output>NO_OUT){
                std::cout << warning_eta_integration << " Adjust limits."   << std::endl;
            }
            if (output==ALL_OUT) {
                out(integral)
                out(contribution_first)
                out(contribution_last)
                out(error_lim)
                out(error)
            }
            
            
            integral            *=  delta_eta*phi_integration;
//            if (integral<1e-90)
//            integral = 1e-90;
            return integral;
        }
        
        
    };
    
    
}

//end CROSSSECTIONTRANSFORMATIONS_H

#endif
