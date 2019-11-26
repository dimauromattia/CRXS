#include "iostream"
#include "fstream"
#include "sstream"
#include "math.h"
#include "stdexcept"

#include "xs.h"
#include "xs_definitions.h"
#include "crxs.h"



#define C_array_to_double(NAM) double C##NAM = C_array[NAM];

namespace CRXS {
    double XS_definitions::fMass_proton   = 0.9382720813;
    double XS_definitions::fMass_neutron  = 0.9395654133;
    double XS_definitions::fMass_deuteron = 1.8756;         // 0.9382720813*2;
    double XS_definitions::fMass_helium3  = 0.9382720813*3; //FIXME
    
    double XS_definitions::inv_pp_pbar_CM__Winkler(double s, double E_pbar_d, double pT_pbar, double* C_array, int len_C_array){
        
        C_array_to_double( 0);
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
            R = (1 +C9*pow(10-sqrt(s),5))  *  exp(C10*pow(10-sqrt(s),C0)*pow(x_R-fMass_proton/E_pbar_Max, 2));
        }
        double sigma_in = C11 +C12*log(sqrt(s)) + C13*pow(log(sqrt(s)), 2);
        double X = C8 * pow(log(sqrt(s)/4./fMass_proton),2);
        
        //double f0_p = R * sigma_in * C5 * pow(1-x_R, C6) * exp(-m_T/C7); // 2014 paper
        double f0_p = R * sigma_in * C5 * pow(1-x_R, C6) * pow( 1+X*(m_T-fMass_proton), -1./X/C7 );
        double invCsCM = f0_p;
        return invCsCM;
    }
    
    
    double XS_definitions::deltaHyperon(double s, double* C_array, int len_C_array){
        
        C_array_to_double( 1);
        C_array_to_double( 2);
        C_array_to_double( 3);
        C_array_to_double( 4);
        
        double factor   = 0.81;
        
        double hyperon = C1 + C2/(1+pow(C3/s,C4));
        hyperon       *= factor;  // branching ratio to Lambda and Sigma  ###### UNCERTAINTY??
        return hyperon;
    }
    
    double XS_definitions::deltaIsospin(double s,  double* C_array, int len_C_array){
        
        C_array_to_double(14);
        C_array_to_double(15);
        C_array_to_double(16);
        
        return C14/(1+pow(s/C15, C16));
        
    }
    
    
    
    
    double XS_definitions::inv_pp_pbar_CM__diMauro( double s, double E_pbar, double pT_pbar, double* C_array, int len_C_array ){
        
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
    
    
    double XS_definitions::tot_pp__diMauro(double s){
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
    
    
    
    
    double XS_definitions::el_pp__diMauro(double s){
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
    
    
    
    
    double XS_definitions::pbar_overlap_function_projectile(double x_F){
        
        int n = 25;
        
        if( x_F < -0.2499 ) return 0;
        if( x_F >  0.25   ) return 1;
        
        double xF[] = {-0.25, -0.225, -0.2,   -0.175, -0.15, -0.125, -0.1,  -0.075, -0.05, -0.025, 0.,  0.025, 0.05,  0.075, 0.1,   0.125, 0.15, 0.175,   0.2,   0.225,   0.25 };
        double F [] = {0.   , 0.0003, 0.0008, 0.0027, 0.010, 0.035,  0.110, 0.197,  0.295,  0.4,   0.5, 0.6,   0.705, 0.803, 0.890, 0.965, 0.990, 0.9973, 0.9992, 0.9997, 1.0  };
        
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
    
    double XS_definitions::pbar_overlap_function_target(double x_F){
        return 1.-pbar_overlap_function_projectile(x_F);
    }
    
    double * XS_definitions::Get_D_parameters(int parametrization){
        if(parametrization==CRXS::KORSMEIER_I){
             return Korsmeier_I_D1_to_D2;
        }else if(parametrization==CRXS::KORSMEIER_II){
            return Korsmeier_II_D1_to_D2;
        }else if(parametrization==CRXS::KORSMEIER_III){
            return Korsmeier_III_D1_to_D2;
        }else if(parametrization==CRXS::WINKLER){
            return Winkler_D1_to_D2;
        }else if(parametrization==CRXS::WINKLER_II){
            return Winkler_II_D1_to_D2;
        }else if(parametrization==CRXS::DI_MAURO_I){
            return diMauro_I_D1_to_D2;
        }else if(parametrization==CRXS::DI_MAURO_II){
            return diMauro_II_D1_to_D2;
        }else if(parametrization==CRXS::ANDERSON){
            return Korsmeier_II_D1_to_D2;
        }else if(parametrization==CRXS::WINKLER_SELF){
            return Winkler_SELF_D1_to_D2;
        }else if(parametrization==CRXS::DI_MAURO_SELF){
            return diMauro_SELF_D1_to_D2;
        }else{
            printf( "Warning in CRXS::XS_definitions::Get_D_parameters. Parametrizatino %i is not known.", parametrization);
        }
        return Dummy;
        
    }
    
    double *  XS_definitions::Get_C_parameters(int parametrization){
        if(parametrization==CRXS::KORSMEIER_I){
            return Korsmeier_I_C1_to_C11;
        }else if(parametrization==CRXS::KORSMEIER_II){
            return Korsmeier_II_C1_to_C16;
        }else if(parametrization==CRXS::KORSMEIER_III){
            return Korsmeier_III_C1_to_C16;
        }else if(parametrization==CRXS::WINKLER){
            return Winkler_C1_to_C16;
        }else if(parametrization==CRXS::WINKLER_II){
            return Winkler_II_C1_to_C16;
        }else if(parametrization==CRXS::DI_MAURO_I){
            return diMauro_I_C1_to_C11;
        }else if(parametrization==CRXS::DI_MAURO_II){
            return diMauro_II_C1_to_C11;
        }else if(parametrization==CRXS::ANDERSON){
                return Korsmeier_II_C1_to_C16;
        }else if(parametrization==CRXS::WINKLER_SELF){
            return Winkler_SELF_C1_to_C16;
        }else if(parametrization==CRXS::DI_MAURO_SELF){
            return diMauro_SELF_C1_to_C11;
        }else{
            printf( "Warning in CRXS::XS_definitions::Get_C_parameters. Parametrizatino %i is not known.", parametrization);
        }
        return Dummy;
    }
    
    double *  XS_definitions::Get_C_parameters_isospin(int parametrization){
        if(parametrization==CRXS::KORSMEIER_I){
            return Korsmeier_II_C1_to_C16; // is correct since KORSMEIER I and II use same isospin parameters
        }else if(parametrization==CRXS::KORSMEIER_II){
            return Korsmeier_II_C1_to_C16;
        }else if(parametrization==CRXS::KORSMEIER_III){
            return Korsmeier_III_C1_to_C16;
        }else if(parametrization==CRXS::WINKLER){
            return Winkler_C1_to_C16;
        }else if(parametrization==CRXS::WINKLER_II){
            return Winkler_II_C1_to_C16;
        }else if(parametrization==CRXS::DI_MAURO_I){
            return diMauro_I_C1_to_C16;
        }else if(parametrization==CRXS::DI_MAURO_II){
            return diMauro_II_C1_to_C16;
        }else if(parametrization==CRXS::ANDERSON){
            return Anderson_C1_to_C16;
        }else if(parametrization==CRXS::WINKLER_SELF){
            return Winkler_SELF_C1_to_C16;
        }else if(parametrization==CRXS::DI_MAURO_SELF){
            return diMauro_SELF_C1_to_C16;
        }else{
            printf( "Warning in CRXS::XS_definitions::Get_C_parameters_isospin. Parametrizatino %i is not known.", parametrization);
        }
        return Dummy;
    }
    
    
    double XS_definitions::factor__AA( double s, double xF, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization ){
        
        if (1000*A_projectile+100*N_projectile+10*A_target+N_target==1010) {
            return 1;
        }
        double * D_array = Get_D_parameters        (parametrization);
        double * C_array = Get_C_parameters_isospin(parametrization);

        if(parametrization==DI_MAURO_I || parametrization==DI_MAURO_II){
            return pow(A_projectile*A_target, D_array[1]);
        }

        double proj = pow(A_projectile, D_array[2])*(1+   deltaIsospin(s,&C_array[0])*N_projectile/A_projectile)*pbar_overlap_function_projectile( xF );
        double targ = pow(A_target,     D_array[2])*(1+   deltaIsospin(s,&C_array[0])*N_target    /A_target    )*pbar_overlap_function_target    ( xF );
        
        if(parametrization==WINKLER || parametrization==WINKLER_II){
            proj    = pow(A_projectile, D_array[2])*(1+0.*deltaIsospin(s,&C_array[0])*N_projectile/A_projectile)*pbar_overlap_function_projectile( xF );
            targ    = pow(A_target,     D_array[2])*(1+0.*deltaIsospin(s,&C_array[0])*N_target    /A_target    )*pbar_overlap_function_target    ( xF );
        }
        
        return pow(A_projectile*A_target, D_array[1])*( proj + targ );
    }
    
    double XS_definitions::inv_pp_p_CM__Anderson(double s, double E_p, double pT_p, double* C_array, int len_C_array){
        E_p = fabs(E_p);
        if (s<4*E_p*E_p)
            return 0;
        if ( pow(pT_p, 2.) > pow(E_p, 2.) - pow(fMass_proton, 2.) )
            return 0;
        return E_p/2./M_PI*pT_p*610*exp(-pT_p/0.166);
    }
    
    
    //  ------------------------------------------------------------- #
    //   Parameter definitions:                                       #
    //  ------------------------------------------------------------- #
    //
    //
    //  Parameters from
    //   1)  Korsmeier et al. 2018 (Param I):
    //   2)  Mauro, et al. (Eq 12); 2014 DOI: 10.1103/PhysRevD.90.085017 (recommended by KDD18):
    //   3)  Mauro, et al. (Eq 13); 2014 DOI: 10.1103/PhysRevD.90.085017:
    //
    //  The names of the parameters correspond to Mauro, et al.
    //
    //                                                {  C0,   C1,          C2,          C3,           C4,          C5,          C6,          C7,          C8,            C9,        C10,    C11   }
    //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    double XS_definitions::Korsmeier_I_C1_to_C11 [] = {  -1,   3.50193e+00, 5.58513e+00, 3.99553e-02, -2.50716e-01, 2.65053e+00, 3.78145e-02, 4.29478e-02, 2.69520e+00,   0.0,       0.0,    0.0   };
    double XS_definitions::diMauro_I_C1_to_C11   [] = {  -1,   4.499,       3.41,        0.00942,      0.445,       3.502,       0.0622,      -0.247,      2.576,         0.0,       0.0,    0.0   };
    double XS_definitions::diMauro_II_C1_to_C11  [] = {  -1,   4.448,       3.735,       0.00502,      0.708,       3.527,       0.236,       -0.729,      2.517,        -1.822e-11, 3.527,  0.384 };
    double XS_definitions::diMauro_SELF_C1_to_C11[] = {  -1,   3.50193e+00, 5.58513e+00, 3.99553e-02, -2.50716e-01, 2.65053e+00, 3.78145e-02, 4.29478e-02, 2.69520e+00,   0.0,       0.0,    0.0   };
    //
    //
    //  The names of the parameters correspond to Korsmeier et al. 2018
    //
    //                                                 {  D0,   D1,   D2    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    double XS_definitions::Korsmeier_I_D1_to_D2  []  = {  -1,  0.825, 0.167 };
    double XS_definitions::diMauro_SELF_D1_to_D2 []  = {  -1,  0.825, 0.167 };
    
    // Parameters are choosen such that the XS is simply scaled by A^0.8 for projectile and target in the case of di Mauro, et al. XSs.
    double XS_definitions::diMauro_I_D1_to_D2    []  = {  -1,  0.8,   0.    };
    double XS_definitions::diMauro_II_D1_to_D2   []  = {  -1,  0.8,   0.    };
    
    
    
    //  Parameters from
    //   1)  Winker 2017 arXiv:1701.04866:
    //   2)  Korsmeier 2018 (Param II):
    //
    //  The names of the parameters correspond to Winker 2017
    //
    //                                                [ C0,    C1,     C2,   C3,     C4,   C5,          C6,      C7,          C8,    C9,          C10,         C11,   C12,    C13,   C14,   C15,    C16  ]
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    double XS_definitions::Winkler_C1_to_C16      [] = {   1,   0.31,   0.30, 21316., 0.9,  0.047,       7.76,    0.168,       0.038, 1.0e-3,      0.7,         30.9,  -1.74,  0.71,  0.114, 20736., 0.51 };
    double XS_definitions::Winkler_II_C1_to_C16   [] = {   2,   0.31,   0.30, 21316., 0.9,  0.047,       7.76,    0.168,       0.038, 1.0e-3,      0.7,         30.9,  -1.74,  0.71,  0.114, 20736., 0.51 };
    double XS_definitions::Korsmeier_II_C1_to_C16 [] = {   1,   0.31,   0.30, 21316., 0.9,  5.01767e-02, 7.79045, 1.64809e-01, 0.038, 4.74370e-04, 3.70480e+00, 30.9,  -1.74,  0.71,  0.114, 20736., 0.51 };
    double XS_definitions::Korsmeier_III_C1_to_C16[] = {   2,   0.31,   0.30, 21316., 0.9,  5.01767e-02, 7.79045, 1.64809e-01, 0.038, 4.74370e-04, 3.70480e+00, 30.9,  -1.74,  0.71,  0.114, 20736., 0.51 };
    double XS_definitions::Winkler_SELF_C1_to_C16 [] = {   1,   0.31,   0.30, 21316., 0.9,  5.01767e-02, 7.79045, 1.64809e-01, 0.038, 4.74370e-04, 3.70480e+00, 30.9,  -1.74,  0.71,  0.114, 20736., 0.51 };
    
    // This array contains only the isospin and hyperon parameters, all others are set to 0
    double XS_definitions::diMauro_SELF_C1_to_C16[] = {  -1,   0.31,   0.30, 21316., 0.9,  0,           0,       0,           0,     0,           0,           0,      0,     0,     0.114, 20736., 0.51 };
    
    // Parameters are choosen to make deltaHyperon=0 and delta isoSpin=0.3 for diMauro XS:
    double XS_definitions::diMauro_I_C1_to_C16   [] = {  -1,   0.00,   0.00,     0., 0.0,  0.,          0.,      0.,          0.,    0.,          0,            0.,    0.,    0.,    0.6,     100., 0.0 };
    double XS_definitions::diMauro_II_C1_to_C16  [] = {  -1,   0.00,   0.00,     0., 0.0,  0.,          0.,      0.,          0.,    0.,          0,            0.,    0.,    0.,    0.6,     100., 0.0 };

    // For Anderson's pp->p use the isospin factor 0
    double XS_definitions::Anderson_C1_to_C16    [] = {  -1,   0.00,   0.00,     0., 0.0,  0.,          0.,      0.,          0.,    0.,          0,            0.,    0.,    0.,    0.0,     100., 0.0 };

    

    //
    //  The names of the parameters correspond to Korsmeier et al. 2018
    //
    //                                                [  D0,  D1,    D2     ]
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    double XS_definitions::Korsmeier_II_D1_to_D2 [] = {  -1, 0.828, 0.145  };
    double XS_definitions::Korsmeier_III_D1_to_D2[] = {  -1, 0.828, 0.145  };
    double XS_definitions::Winkler_D1_to_D2      [] = {  -1, 0.839, 0.161  };      // value of <nu_He>=1.25 is translated to D_1 and D_2 (see Korsmeier et al. 2018)
    double XS_definitions::Winkler_II_D1_to_D2   [] = {  -1, 0.839, 0.161  };      // value of <nu_He>=1.25 is translated to D_1 and D_2 (see Korsmeier et al. 2018)
    double XS_definitions::Winkler_SELF_D1_to_D2 [] = {  -1, 0.828, 0.145  };
    
    double XS_definitions::Dummy                 [] = {  -1  };

    
    
    
    //  ------------------------------------------------------------- #
    //  Interpolation of (total) cross sections from tables
    //          -> pbarp and pbarD
    //  ------------------------------------------------------------- #
  
    
    double XS_definitions::totXS_get_interpolation_loglin(double x, double array[91][4] ){
        
        double  dx = 10 * (log10(x)+2);
        int     ix = dx;
        
        if(ix<0)  ix =  0;
        if(ix>90) ix = 90;
        
        double vl = log( array[ix  ][1]  );
        double vu = log( array[ix+1][1]  );
        
        return exp( vl+(vu-vl)*(dx-ix)  );
        
    }
    
    bool XS_definitions::f_totXS_IsRead = false;
    
    void XS_definitions::totXS_TableToArray( std::string file, double array[91][4] ){
        std::ifstream   ifs;
        std::string     line;
        std::string     line2;
        ifs.open(file.c_str());
        if(!ifs){
            std::cerr << ("ERROR: File: "+file+" does not exist.").c_str() << std::endl;
            return;
        }
        int idxrow = 0;
        while(!ifs.eof()){
            std::getline(ifs, line);
            std::stringstream ss(line);
            if (line == "" ) {
                continue;
            }
            if ((line.at(0) == '#') || (line.at(0) == '*')) {
                continue;
            }
            if (idxrow>=91) {
                std::cerr << "ERROR: Table in " << file << " has too many rows." << std::endl;
                exit(0);
            }
            int idxcol = 0;
            while (std::getline(ss, line2, ' ')) {
                if (line2=="") {
                    continue;
                }
                if (idxcol>=4) {
                    std::cerr << "ERROR: Table in " << file << " has too many colums." << std::endl;
                    exit(0);
                }
                try {
                    double v = std::stod(line2);
                    if (v<1e-100) v = 1e-100;
                    array[idxrow][idxcol] = v;
                } catch(const std::invalid_argument &ia) {
                    std::cerr << "ERROR: Some value in " << file << " is not a double" << std::endl;
                    throw;
                }
                idxcol++;
            }
            idxrow++;
        }
        ifs.close();
    }
    
    
    double XS_definitions::fXS__el_pbarp [91][4];
    double XS_definitions::fXS__tot_pbarp[91][4];
    double XS_definitions::fXS__tot_pbarD[91][4];
    double XS_definitions::fXS__nar_pbarD[91][4];
    
    void XS_definitions::totXS_Read(){
        
        std::string file_ppbar_el  = CRXS_config::Get_CRXS_DataDir()+"/table_ppbar_el.txt" ;
        std::string file_ppbar_tot = CRXS_config::Get_CRXS_DataDir()+"/table_ppbar_tot.txt";
        std::string file_dpbar_tot = CRXS_config::Get_CRXS_DataDir()+"/table_dpbar_tot.txt";
        std::string file_dpbar_nar = CRXS_config::Get_CRXS_DataDir()+"/table_dpbar_nar.txt";
        
        std::cout << "Read file: " << file_ppbar_el << std::endl;
        totXS_TableToArray( file_ppbar_el, fXS__el_pbarp );
        
        std::cout << "Read file: " << file_ppbar_tot << std::endl;
        totXS_TableToArray( file_ppbar_tot, fXS__tot_pbarp );
        
        std::cout << "Read file: " << file_dpbar_tot << std::endl;
        totXS_TableToArray( file_dpbar_tot, fXS__tot_pbarD );
        
        std::cout << "Read file: " << file_dpbar_nar << std::endl;
        totXS_TableToArray( file_dpbar_nar, fXS__nar_pbarD );
        
        f_totXS_IsRead = true;
    }
    
    double XS_definitions::el_pbarp (double T_pbar){
        if(!f_totXS_IsRead) totXS_Read();
        return totXS_get_interpolation_loglin(T_pbar, fXS__el_pbarp);
    };
    double XS_definitions::tot_pbarp (double T_pbar){
        if(!f_totXS_IsRead) totXS_Read();
        return totXS_get_interpolation_loglin(T_pbar, fXS__tot_pbarp);
    };
    double XS_definitions::tot_pbarD (double T_pbar){
        if(!f_totXS_IsRead) totXS_Read();
        return totXS_get_interpolation_loglin(T_pbar, fXS__tot_pbarD);
    };
    double XS_definitions::nar_pbarD (double T_pbar){
        if(!f_totXS_IsRead) totXS_Read();
        return totXS_get_interpolation_loglin(T_pbar, fXS__nar_pbarD);
    };
    
}




