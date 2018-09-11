#include "xs_wrapper.h"
#include "xs.h"

// pbar
double inv_AA_pbar_CM( double s, double xF, double pT_pbar, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization){
    return CRXS::XS::inv_AA_pbar_CM( s, xF, pT_pbar, A_projectile, N_projectile, A_target, N_target, parametrization);
};
double inv_AA_pbar_LAB( double Tn_proj_LAB, double T_pbar_LAB, double eta_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization ){
    return CRXS::XS::inv_AA_pbar_LAB( Tn_proj_LAB, T_pbar_LAB, eta_LAB, A_projectile, N_projectile, A_target, N_target, parametrization );
};
double dE_AA_pbar_LAB( double Tn_proj_LAB, double T_pbar_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization){
    return CRXS::XS::dE_AA_pbar_LAB( Tn_proj_LAB, T_pbar_LAB, A_projectile, N_projectile, A_target, N_target, parametrization);
};
double dE_AA_pbar_LAB_incNbarAndHyperon(double Tn_proj_LAB, double T_pbar_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization){
    return CRXS::XS::dE_AA_pbar_LAB_incNbarAndHyperon( Tn_proj_LAB,  T_pbar_LAB, A_projectile, N_projectile, A_target, N_target, parametrization);
};


// Dbar
double inv_AA_Dbar_CM( double s, double xF_Dbar, double pT_Dbar, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization, int coalescence){
    return CRXS::XS::inv_AA_Dbar_CM( s, xF_Dbar, pT_Dbar, A_projectile, N_projectile, A_target, N_target, parametrization, coalescence);
};
double inv_AA_Dbar_LAB( double Tn_proj_LAB, double T_Dbar_LAB, double eta_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization, int coalescence){
    return CRXS::XS::inv_AA_Dbar_LAB( Tn_proj_LAB, T_Dbar_LAB, eta_LAB, A_projectile, N_projectile, A_target, N_target, parametrization, coalescence);
};
double dEn_AA_Dbar_LAB( double Tn_proj_LAB, double Tn_Dbar_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization, int coalescence ){
    return CRXS::XS::dEn_AA_Dbar_LAB( Tn_proj_LAB, Tn_Dbar_LAB, A_projectile, N_projectile, A_target, N_target, parametrization, coalescence );
};

