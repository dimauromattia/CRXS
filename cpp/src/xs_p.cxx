#include "math.h"
#include "iostream"

#include "gsl_integration.h"

#include "xs.h"
#include "xs_definitions.h"

namespace CRXS {
   
    
    double XS::inv_AA_p_CM( double s, double xF, double pT_p, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization){
        
        double pL_p = xF*sqrt(s)/2.;
        double E_p  = sqrt( XS_definitions::fMass_proton*XS_definitions::fMass_proton + pL_p*pL_p + pT_p*pT_p );
        
        if      (  parametrization==ANDERSON  ){
            double pp = XS_definitions::inv_pp_p_CM__Anderson(s, E_p, pT_p ) ;
            double AA = XS_definitions::factor__AA( s, xF, A_projectile, N_projectile, A_target, N_target, parametrization );
            return pp * AA;
        }else{
            printf( "Warning in CRXS::XS::inv_AA_p_CM. Parametrization %i is not known.", parametrization);
        }
        return 0;
    }
    
    double XS::inv_AA_p_LAB( double Tn_proj_LAB, double T_p_LAB, double eta_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization){
        double s, E_p, pT_p, x_F;
        convert_LAB_to_CM( Tn_proj_LAB, T_p_LAB, eta_LAB, s, E_p, pT_p, x_F );
        return inv_AA_p_CM(s, x_F, pT_p, A_projectile, N_projectile, A_target, N_target, parametrization);
    }
    
    double XS::integrand__dE_AA_p_LAB (double eta_LAB, void* parameters  ){
        
        double* par = (double * ) parameters;
        double Tn_proj_LAB      = par[0];
        double T_p_LAB          = par[1];
        int    A_projectile     = par[2];
        int    N_projectile     = par[3];
        int    A_target         = par[4];
        int    N_target         = par[5];
        int    parametrization  = par[6];
        
        return  pow( cosh(eta_LAB), -2 ) * XS::inv_AA_p_LAB( Tn_proj_LAB, T_p_LAB, eta_LAB, A_projectile, N_projectile, A_target, N_target, parametrization );
        
    }
    
    
    
    double XS::dE_AA_p_LAB( double Tn_proj_LAB, double T_p_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization ){
        
        //
        //  Integrate over all solid angle and transform to enery differential (d sigma / d E)
        //
        double E_p_LAB = T_p_LAB + XS_definitions::fMass_proton;
        double p_p_LAB = sqrt(  pow( E_p_LAB, 2 ) - pow( XS_definitions::fMass_proton, 2 )  );
        if (p_p_LAB!=p_p_LAB){
            return 0;
        }
        double Jacobian_and_conversion      = 2*3.1415926536*p_p_LAB;
        // it contains: phi_integration (2 pi), inv to d3p (1/E_p_LAB), Jacobian(p_p_Lab*p_p_Lab), dp to dE (E_p_LAB/p_p_Lab)
        double epsabs = 0;
        double epsrel = 1e-4;
        double res, err;
        size_t neval;
        double par[] = { Tn_proj_LAB, T_p_LAB, 1.0001*A_projectile, 1.0001*N_projectile, 1.0001*A_target, 1.0001*N_target, 1.0001*parametrization };
        
        gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
        
        gsl_function F;
        F.function = &integrand__dE_AA_p_LAB;
        F.params = &par[0];
        
        gsl_integration_qag(&F, 0, 50, epsabs, epsrel, 1000, 2, w, &res, &err);
        gsl_integration_workspace_free (w);
        
        if(err/res>epsrel){
            printf( "Warning in CRXS::XS::dE_AA_p_LAB. Integral accuarcy of %f is below required value of %f. \n", err/res, epsrel);
        }
        
        res *=  Jacobian_and_conversion;
        return res;
        
    }
    
    
//    double XS::dE_AA_p_LAB_incNAndHyperon(double Tn_proj_LAB, double T_p_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization){
//        double s = 4*XS_definitions::fMass_proton*XS_definitions::fMass_proton + 2 * Tn_proj_LAB * XS_definitions::fMass_proton;
//        double * C_array = XS_definitions::Get_C_parameters_isospin(parametrization);
//        return dE_AA_p_LAB(Tn_proj_LAB, T_p_LAB, A_projectile, N_projectile, A_target, N_target, parametrization)*( 2 + 2*XS_definitions::deltaHyperon(s, C_array) + XS_definitions::deltaIsospin(s, C_array));
//    }
    
    
    
    
    
    
    
    
}
