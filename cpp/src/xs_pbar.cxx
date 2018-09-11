#include "math.h"
#include "iostream"

#include "gsl_integration.h"

#include "xs.h"
#include "xs_definitions.h"

namespace CRXS {
   
    
    double XS::inv_AA_pbar_CM( double s, double xF, double pT_pbar, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization){
        
        double pL_pbar = xF*sqrt(s)/2.;
        double E_pbar  = sqrt( XS_definitions::fMass_proton*XS_definitions::fMass_proton + pL_pbar*pL_pbar + pT_pbar*pT_pbar );
        
        if      (  parametrization==KORSMEIER_II  ){
            double pp = XS_definitions::inv_pp_pbar_CM__Winkler(s, E_pbar, pT_pbar, XS_definitions::Korsmeier_II_C1_to_C16 ) ;
            double AA = XS_definitions::factor__AA( s, xF, A_projectile, N_projectile, A_target, N_target, parametrization );
            return pp * AA;
        }else if(  parametrization==KORSMEIER_I   ){
            double pp = XS_definitions::inv_pp_pbar_CM__diMauro(s, E_pbar, pT_pbar, XS_definitions::Korsmeier_I_C1_to_C11 ) ;
            double AA = XS_definitions::factor__AA( s, xF, A_projectile, N_projectile, A_target, N_target, parametrization );
            return pp * AA;
        }else if(  parametrization==WINKLER       ){
            double pp = XS_definitions::inv_pp_pbar_CM__Winkler(s, E_pbar, pT_pbar, XS_definitions::Winkler_C1_to_C16 ) ;
            double AA = XS_definitions::factor__AA( s, xF, A_projectile, N_projectile, A_target, N_target, parametrization );
            return pp * AA;
        }else if(  parametrization==DI_MAURO_I    ){
            double pp = XS_definitions::inv_pp_pbar_CM__diMauro(s, E_pbar, pT_pbar, XS_definitions::diMauro_I_C1_to_C11 ) ;
            double AA = XS_definitions::factor__AA( s, xF, A_projectile, N_projectile, A_target, N_target, parametrization );
            return pp * AA;
        }else if(  parametrization==DI_MAURO_II   ){
            double pp = XS_definitions::inv_pp_pbar_CM__diMauro(s, E_pbar, pT_pbar, XS_definitions::diMauro_II_C1_to_C11 ) ;
            double AA = XS_definitions::factor__AA( s, xF, A_projectile, N_projectile, A_target, N_target, parametrization );
            return pp * AA;
        }else{
            printf( "Warning in CRXS::XS::inv_AA_pbar_CM. Parametrization %i is not known.", parametrization);
        }
        return 0;
    }
    
    double XS::inv_AA_pbar_LAB( double Tn_proj_LAB, double T_pbar_LAB, double eta_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization){
        double s, E_pbar, pT_pbar, x_F;
        convert_LAB_to_CM( Tn_proj_LAB, T_pbar_LAB, eta_LAB, s, E_pbar, pT_pbar, x_F );
        return inv_AA_pbar_CM(s, x_F, pT_pbar, A_projectile, N_projectile, A_target, N_target, parametrization);
    }
    
    double XS::integrand__dE_AA_pbar_LAB (double eta_LAB, void* parameters  ){
        
        double* par = (double * ) parameters;
        double Tn_proj_LAB      = par[0];
        double T_pbar_LAB       = par[1];
        int    A_projectile     = par[2];
        int    N_projectile     = par[3];
        int    A_target         = par[4];
        int    N_target         = par[5];
        int    parametrization  = par[6];
        
        return  pow( cosh(eta_LAB), -2 ) * XS::inv_AA_pbar_LAB( Tn_proj_LAB, T_pbar_LAB, eta_LAB, A_projectile, N_projectile, A_target, N_target, parametrization );
        
    }
    
    double XS::dE_AA_pbar_LAB( double Tn_proj_LAB, double T_pbar_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization ){
        
        //
        //  Integrate over all solid angle and transform to enery differential (d sigma / d E)
        //
        double E_pbar_LAB = T_pbar_LAB + XS_definitions::fMass_proton;
        double p_pbar_LAB = sqrt(  pow( E_pbar_LAB, 2 ) - pow( XS_definitions::fMass_proton, 2 )  );
        if (p_pbar_LAB!=p_pbar_LAB){
            return 0;
        }
        double Jacobian_and_conversion      = 2*3.1415926536*p_pbar_LAB;
        // it contains: phi_integration (2 pi), inv to d3p (1/E_pbar_LAB), Jacobian(p_pbar_Lab*p_pbar_Lab), dp to dE (E_pbar_LAB/p_pbar_Lab)
        double epsabs = 0;
        double epsrel = 1e-4;
        double res, err;
        size_t neval;
        double par[] = { Tn_proj_LAB, T_pbar_LAB, 1.0001*A_projectile, 1.0001*N_projectile, 1.0001*A_target, 1.0001*N_target, 1.0001*parametrization };
        
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        
        gsl_function F;
        F.function = &integrand__dE_AA_pbar_LAB;
        F.params = &par[0];
        
        gsl_integration_qag(&F, 0, 50, epsabs, epsrel, 1000, 2, w, &res, &err);
        
        
        if(err/res>epsrel){
            printf( "Warning in CRXS::XS::dE_AA_pbar_LAB. Integral accuarcy of %f is below required value of %f. \n", err/res, epsrel);
        }
        
        res *=  Jacobian_and_conversion;
        return res;
        
    }
    
    
    double XS::dE_AA_pbar_LAB_incNbarAndHyperon(double Tn_proj_LAB, double T_pbar_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization){
        double s = 4*XS_definitions::fMass_proton*XS_definitions::fMass_proton + 2 * Tn_proj_LAB * XS_definitions::fMass_proton;
        double * C_array = XS_definitions::Get_C_parameters_isospin(parametrization);
        return dE_AA_pbar_LAB(Tn_proj_LAB, T_pbar_LAB, A_projectile, N_projectile, A_target, N_target, parametrization)*( 2 + 2*XS_definitions::deltaHyperon(s, C_array) + XS_definitions::deltaIsospin(s, C_array));
    }
    
    
    
    
    
    
    
    
}
