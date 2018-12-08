
// pbar
double inv_AA_pbar_CM( double s, double xF, double pT_pbar, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization);
double inv_AA_pbar_LAB( double Tn_proj_LAB, double T_pbar_LAB, double eta_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization );
double dE_AA_pbar_LAB( double Tn_proj_LAB, double T_pbar_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization);
double dE_AA_pbar_LAB_incNbarAndHyperon(double Tn_proj_LAB, double T_pbar_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization);


// Dbar
double inv_AA_Dbar_CM( double s, double xF_Dbar, double pT_Dbar, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization, int coalescence);
double inv_AA_Dbar_LAB( double Tn_proj_LAB, double T_Dbar_LAB, double eta_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization, int coalescence);
double dEn_AA_Dbar_LAB( double Tn_proj_LAB, double Tn_Dbar_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization, int coalescence );

// p
double inv_AA_p_CM( double s, double xF, double pT_p, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization);
double inv_AA_p_LAB( double Tn_proj_LAB, double T_p_LAB, double eta_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization );
double dE_AA_p_LAB( double Tn_proj_LAB, double T_p_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization);

