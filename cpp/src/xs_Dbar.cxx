#include "math.h"
#include "iostream"

#include "gsl_integration.h"

#include "xs.h"
#include "xs_definitions.h"

namespace CRXS {
    
    double XS::p_coal__VonDoetinchen( double s ){
        double T = ( s-4*XS_definitions::fMass_proton*XS_definitions::fMass_proton) /2./XS_definitions::fMass_proton;
        
        double A = 0.153;
        double B = 4.5;
        double C = 1.47;
        
        double p = A/(1.+exp(B-log(T)/C));
        
        return p;
    }
    
    
    double XS::inv_AA_Dbar_CM( double s, double xF_dbar, double pT_dbar, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization, int coalescence){
        
        int signed_A_projectile = A_projectile;
        A_projectile = fabs(1.0001*A_projectile);
        
        int nucleons = 2;
        double p_coalescence;
        if       (coalescence==FIXED_P0) {
            p_coalescence = 0.080;
        }else if (coalescence==ENERGY_DEP__VAN_DOETINCHEM) {
            p_coalescence = p_coal__VonDoetinchen(s);
        }else{
            return -1;
        }
        
        
        double pL_dbar = xF_dbar/2.*sqrt(s);
        if(pL_dbar!=pL_dbar){
            return 0;
        }
        double E_pbar = sqrt( pow(XS_definitions::fMass_proton,  2) + pow(pT_dbar/nucleons,2) + pow(pL_dbar/nucleons,2) );
        //double E_nbar = sqrt( pow(XS_definitions::fMass_neutron, 2) + pow(pT_dbar/nucleons,2) + pow(pL_dbar/nucleons,2) );
        double E_dbar = sqrt( pow(XS_definitions::fMass_deuteron,2) + pow(pT_dbar,         2) + pow(pL_dbar,         2) );
        
        double sq__s_red = sqrt(s) - E_dbar;
        if (sq__s_red<0) {
            return 0;
        }
        double s_red = sq__s_red*sq__s_red;
        
        double inv_pp_pbar, inv_pp_nbar, inv_pp_pbar_reduced, inv_pp_nbar_reduced, AA, AA_reduced;
        
        double * C_array         = XS_definitions::Get_C_parameters        (parametrization);
        double * C_array_isospin = XS_definitions::Get_C_parameters_isospin(parametrization);
        double * D_array         = XS_definitions::Get_D_parameters        (parametrization);
        
        
        double XS;
        XS  = XS_definitions::fMass_deuteron/XS_definitions::fMass_proton/XS_definitions::fMass_neutron;
        XS *= (4./3. * 3.1415926536 * pow(p_coalescence,3)) / (pow(A_target*A_projectile, D_array[1]+D_array[2])*XS_definitions::tot_pp__diMauro(s));
        
        if      (  parametrization==KORSMEIER_II || parametrization==WINKLER  ){
            inv_pp_pbar         = XS_definitions::inv_pp_pbar_CM__Winkler(s,     E_pbar, pT_dbar/nucleons, C_array );
            inv_pp_pbar_reduced = XS_definitions::inv_pp_pbar_CM__Winkler(s_red, E_pbar, pT_dbar/nucleons, C_array );
        }else if(  parametrization==KORSMEIER_I  || parametrization==DI_MAURO_I || parametrization==DI_MAURO_II ){
            inv_pp_pbar         = XS_definitions::inv_pp_pbar_CM__diMauro(s,     E_pbar, pT_dbar/nucleons, C_array );
            inv_pp_pbar_reduced = XS_definitions::inv_pp_pbar_CM__diMauro(s_red, E_pbar, pT_dbar/nucleons, C_array );
        }else{
            printf( "Warning in CRXS::XS::inv_AA_pbar_CM. Parametrization %i is not known.", parametrization);
            return 0;
        }
        
        double deltaHyperon     = XS_definitions::deltaHyperon(s,     C_array_isospin);
        double deltaHyperon_red = XS_definitions::deltaHyperon(s_red, C_array_isospin);
        
        double deltaIsospin     = XS_definitions::deltaIsospin(s,     C_array_isospin);
        double deltaIsospin_red = XS_definitions::deltaIsospin(s_red, C_array_isospin);

        inv_pp_nbar         = inv_pp_pbar*        (1+deltaIsospin    +deltaHyperon    );
        inv_pp_nbar_reduced = inv_pp_pbar_reduced*(1+deltaIsospin_red+deltaHyperon_red);
        
        if (signed_A_projectile<0){
            inv_pp_pbar         = XS_definitions::inv_pp_p_CM__Anderson(s,     E_pbar, pT_dbar/nucleons);
            inv_pp_pbar_reduced = XS_definitions::inv_pp_p_CM__Anderson(s_red, E_pbar, pT_dbar/nucleons);
        }

        inv_pp_pbar         = inv_pp_pbar*        (1+deltaHyperon);
        inv_pp_pbar_reduced = inv_pp_pbar_reduced*(1+deltaHyperon);

        AA         = XS_definitions::factor__AA( s,     xF_dbar/nucleons, A_projectile, N_projectile, A_target, N_target, parametrization );
        AA_reduced = XS_definitions::factor__AA( s_red, xF_dbar/nucleons, A_projectile, N_projectile, A_target, N_target, parametrization );
        
        XS *= 0.5 * AA * AA_reduced * (  inv_pp_pbar * inv_pp_nbar_reduced   +   inv_pp_nbar * inv_pp_pbar_reduced  );
        
        return XS;
    }

    double XS::inv_AA_Dbar_LAB( double Tn_proj_LAB, double Tn_Dbar_LAB, double eta_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization, int coalescence ){
        int    nucleons = 2;
        double s, E_Dbar, pT_pbar, x_F;
        double T_Dbar_LAB = nucleons * Tn_Dbar_LAB;
        convert_LAB_to_CM( Tn_proj_LAB, T_Dbar_LAB, eta_LAB, s, E_Dbar, pT_pbar, x_F, D_BAR );
        return inv_AA_Dbar_CM(s, x_F, pT_pbar, A_projectile, N_projectile, A_target, N_target, parametrization, coalescence);
    }
    
    double XS::integrand__dE_AA_Dbar_LAB (double eta_LAB, void* parameters  ){
        
        double* par = (double * ) parameters;
        double Tn_proj_LAB      = par[0];
        double Tn_Dbar_LAB      = par[1];
        int    A_projectile     = par[2];
        int    N_projectile     = par[3];
        int    A_target         = par[4];
        int    N_target         = par[5];
        int    parametrization  = par[6];
        int    coalescence      = par[7];
        
//        std::cout << " Tn_proj_LAB     " <<  Tn_proj_LAB      << std::endl;
//        std::cout << " Tn_Dbar_LAB     " <<  Tn_Dbar_LAB      << std::endl;
//        std::cout << " A_projectile    " <<  A_projectile     << std::endl;
//        std::cout << " par[2]          " <<  par[2]           << std::endl;
//        std::cout << " N_projectile    " <<  N_projectile     << std::endl;
//        std::cout << " A_target        " <<  A_target         << std::endl;
//        std::cout << " N_target        " <<  N_target         << std::endl;
//        std::cout << " parametrization " <<  parametrization  << std::endl;
//        std::cout << " coalescence     " <<  coalescence      << std::endl;
        
        return  pow( cosh(eta_LAB), -2 ) * XS::inv_AA_Dbar_LAB( Tn_proj_LAB, Tn_Dbar_LAB, eta_LAB, A_projectile, N_projectile, A_target, N_target, parametrization, coalescence );
        
    }
    
    double XS::dEn_AA_Dbar_LAB( double Tn_proj_LAB, double Tn_Dbar_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization, int coalescence ){
        int nucleons = 2;
        //
        //  Integrate over all solid angle and transform to enery differential (d sigma / d E)
        //
        double E_Dbar_LAB = nucleons * Tn_Dbar_LAB + XS_definitions::fMass_deuteron;
        double p_Dbar_LAB = sqrt(  pow( E_Dbar_LAB, 2 ) - pow( XS_definitions::fMass_deuteron, 2 )  );
        if (p_Dbar_LAB!=p_Dbar_LAB){
            return 0;
        }
        double Jacobian_and_conversion      = 2*3.1415926536*p_Dbar_LAB;
        // it contains: phi_integration (2 pi), inv to d3p (1/E_pbar_LAB), Jacobian(p_pbar_Lab*p_pbar_Lab), dp to dE (E_pbar_LAB/p_pbar_Lab)
        double epsabs = 0;
        double epsrel = 1e-4;
        double res, err;
        double par[] = { Tn_proj_LAB, Tn_Dbar_LAB, 1.0001*A_projectile, 1.0001*N_projectile, 1.0001*A_target, 1.0001*N_target, 1.0001*parametrization, 1.0001*coalescence };
        
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        
        gsl_function F;
        F.function = &integrand__dE_AA_Dbar_LAB;
        F.params   = &par[0];
        
        gsl_integration_qag(&F, 0, 50, epsabs, epsrel, 1000, 2, w, &res, &err);
        
        if(err/res>epsrel){
            printf( "Warning in CRXS::XS::dE_AA_pbar_LAB. Integral accuarcy of %f is below required value of %f. \n", err/res, epsrel);
        }
        
        res *=  Jacobian_and_conversion;
        return res * nucleons;
        
    }
    
    
    
    double XS::dEn_DbarA_Dbar_LAB(  double Tn_Dbar_proj_LAB, double Tn_Dbar_prod_LAB, int A_target, int N_target, int parametrization  ){
        
        double shape      = 1.;
        double norm_shape = 0.;
        
        double log_10 = log(10);
        
        if (parametrization==APPROX_1_OVER_T) {
            norm_shape = Tn_Dbar_proj_LAB;
            if(Tn_Dbar_prod_LAB>Tn_Dbar_proj_LAB) shape=0;
        }else if (parametrization==ANDERSON) {
            shape = dE_AA_p_LAB( Tn_Dbar_proj_LAB, Tn_Dbar_prod_LAB, 1, 0, A_target, N_target, ANDERSON);
            double dlog10T    = 0.1;
            for (double log10T=-7; log10T<log10(Tn_Dbar_proj_LAB); log10T+=dlog10T) {
                double T = pow(10,log10T);
                norm_shape += T* dE_AA_p_LAB( Tn_Dbar_proj_LAB, T, 1, 0, 1, 0, ANDERSON);
            }
            norm_shape *= dlog10T * log_10;
        }else{
            printf( "Warning in CRXS::XS::dEn_DbarA_Dbar_LAB. Parametrization %i is not known.", parametrization);
            return 0;
        }
        
        double T_pbar     = Tn_Dbar_proj_LAB;
        double norm_XS    = XS_definitions::nar_pbarD(T_pbar);
        
        
//        std::cout << Tn_Dbar_proj_LAB << std::endl;
//        std::cout << Tn_Dbar_prod_LAB << std::endl;
//        std::cout << norm_XS << std::endl;
//        std::cout << shape   << std::endl;
//        std::cout << norm_XS << std::endl;
//        std::cout << " " << std::endl;
        
        return norm_XS * shape / norm_shape;
    };
    
}
