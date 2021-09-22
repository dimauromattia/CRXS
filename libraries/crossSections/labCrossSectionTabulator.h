#ifndef LABCROSSSECTIONTABULATOR_H
#define LABCROSSSECTIONTABULATOR_H

#include <iostream>
#include <math.h>
#include <string>
#include <definitions.h>
#include <crossSections.h>
#include <crossSectionTransformations.h>
#include <fileTools.h>


namespace CRACS {
    
    /*! \brief Tabulator of lab CS.
     *
     */
    class CS_lab_tab{
    public:
        
        static CS_lab_tab* GetInstance();
        double  GetCrossSection ( double Tn_proj, double Tn_prod,  int type, int parametrization);
        
        void    WriteCS(  double (*lab_CS)(double, double), std::string filename, std::string dir="." , int dLog=f_dLog );
        void    ReadCS (  std::string filename, std::string dir=".", int col=f_j_max+2, int skip_header=0  );
        double  GetInterpolation( double Tn_proj, double Tn_prod );
        
        void    WriteCS_lowProd(  double (*lab_CS)(double, double), std::string filename, std::string dir="." , int dLog=f_dLog );
        
        
        double  GetInterpolation_Transform( double Tn_proj, double Tn_prod, int Targ=0, int nProj=301, int nProd=151, double TminProj=5.629632, double TmaxProj=1e4, double TminProd=1e-1, double TmaxProd=1e2 , bool changeLoop=false);
        void    WriteCS_Transform( int col, std::string filename, std::string dir, double factor=1.0  );
        void    WriteCS_Transform( int col, std::string filename, std::string dir, int nProj, int nProd, double TminProj, double TmaxProj, double TminProd, double TmaxProd, double factor, bool changeLoop );
        void    TransformAll();
        void    TransformGalprop   (std::string filename, int parametrization, int dLog=f_dLog);
        void    TransformAllGalprop();
        void    WriteAll();
        void    ReadAll ( std::string path = "default" );
        
        void test();
        
        static double dT_pp_p_LAB_Anderson  ( double T_p, double T_pbar ){
            return CS::dE_pp_p_LAB  ( T_p+fMass_proton,   T_pbar+fMass_proton, 1000, CS::ANDERSON, NO_OUT )*1e-31;
        };
        
        static double dT_pp_pbar_LAB_diMauro12  ( double T_p, double T_pbar ){
            return CS::dE_pp_pbar_LAB  ( T_p+fMass_proton,   T_pbar+fMass_proton, 1000, CS::DI_MAURO12, NO_OUT )*1e-31;
        };
        static double dT_pHe_pbar_LAB_diMauro12 ( double T_p, double T_pbar ){
            return CS::dE_pHe_pbar_LAB ( T_p+fMass_proton,   T_pbar+fMass_proton, 1000, CS::DI_MAURO12, NO_OUT )*1e-31;
        };
        static double dT_Hep_pbar_LAB_diMauro12 ( double T_p, double T_pbar ){
            return CS::dE_Hep_pbar_LAB ( T_p+fMass_proton,   T_pbar+fMass_proton, 1000, CS::DI_MAURO12, NO_OUT )*1e-31;
        };
        static double dT_HeHe_pbar_LAB_diMauro12 ( double T_p, double T_pbar ){
            return CS::dE_HeHe_pbar_LAB ( T_p+fMass_proton,   T_pbar+fMass_proton, 1000, CS::DI_MAURO12, NO_OUT )*1e-31;
        };
        
        
        
        static double dT_pp_pbar_LAB_diMauro  ( double T_p, double T_pbar ){
            return CS::dE_pp_pbar_LAB  ( T_p+fMass_proton,   T_pbar+fMass_proton, 1000, CS::DI_MAURO,   NO_OUT )*1e-31;
        };
        static double dT_pHe_pbar_LAB_diMauro  ( double T_p, double T_pbar ){
            return CS::dE_pHe_pbar_LAB  ( T_p+fMass_proton,   T_pbar+fMass_proton,1000, CS::DI_MAURO,   NO_OUT )*1e-31;
        };
        static double dT_HeHe_pbar_LAB_diMauro  ( double T_p, double T_pbar ){
            return CS::dE_HeHe_pbar_LAB ( T_p+fMass_proton,   T_pbar+fMass_proton,1000, CS::DI_MAURO,   NO_OUT )*1e-31;
        };
        
        
        static double dT_pp_pbar_LAB_TanNg    ( double T_p, double T_pbar ){
            return CS::dE_pp_pbar_LAB  ( T_p+fMass_proton,   T_pbar+fMass_proton, 1000, CS::TAN_NG,     NO_OUT )*1e-31;
        };
        static double dT_pHe_pbar_LAB_TanNg    ( double T_p, double T_pbar ){
            return CS::dE_pHe_pbar_LAB ( T_p+fMass_proton,   T_pbar+fMass_proton, 1000, CS::TAN_NG,     NO_OUT )*1e-31;
        };
        
        
        
        static double dT_pp_pbar_LAB_Duperray ( double T_p, double T_pbar ){
            return CS::dE_pp_pbar_LAB  ( T_p+fMass_proton,   T_pbar+fMass_proton, 1000, CS::DUPERRAY,   NO_OUT )*1e-31;
        };
        static double dT_pHe_pbar_LAB_Duperray ( double T_p, double T_pbar ){
            return CS::dE_pHe_pbar_LAB ( T_p+fMass_proton,   T_pbar+fMass_proton, 1000, CS::DUPERRAY,   NO_OUT )*1e-31;
        };
        static double dT_HeHe_pbar_LAB_Duperray ( double T_p, double T_pbar ){
            return CS::dE_HeHe_pbar_LAB ( T_p+fMass_proton,   T_pbar+fMass_proton,1000, CS::DUPERRAY,   NO_OUT )*1e-31;
        };
        
        
        static double dT_pp_pbar_LAB_Korsmeier ( double T_p, double T_pbar ){
            return CS::dE_pp_pbar_LAB  ( T_p+fMass_proton,   T_pbar+fMass_proton, 1000, CS::KORSMEIER,   NO_OUT )*1e-31;
        };
        static double dT_pHe_pbar_LAB_Korsmeier ( double T_p, double T_pbar ){
            return CS::dE_pHe_pbar_LAB ( T_p+fMass_proton,   T_pbar+fMass_proton, 1000, CS::KORSMEIER,   NO_OUT )*1e-31;
        };
        
        static double dT_pp_pbar_LAB_Winkler ( double T_p, double T_pbar ){
            return CS::dE_pp_pbar_LAB  ( T_p+fMass_proton,   T_pbar+fMass_proton, 100, CS::WINKLER,     NO_OUT )*1e-31;
        };
        static double dT_pHe_pbar_LAB_Winkler ( double T_p, double T_pbar ){
            return CS::dE_pHe_pbar_LAB ( T_p+fMass_proton,   T_pbar+fMass_proton, 100, CS::WINKLER,     NO_OUT )*1e-31;
        };
        static double dT_Hep_pbar_LAB_Winkler ( double Tn_He, double T_pbar ){
            return CS::dE_Hep_pbar_LAB ( Tn_He+fMass_proton, T_pbar+fMass_proton, 100, CS::WINKLER,     NO_OUT )*1e-31;
        };
        static double dT_HeHe_pbar_LAB_Winkler ( double Tn_He, double T_pbar ){
            return CS::dE_HeHe_pbar_LAB ( Tn_He+fMass_proton, T_pbar+fMass_proton, 100, CS::WINKLER,     NO_OUT )*1e-31;
        };
        
        
        static double dT_pp_pbar_LAB_WinklerWitHypWitNbar ( double T_p, double T_pbar ){
            return CS::dE_pp_pbar_LAB  ( T_p+fMass_proton,   T_pbar+fMass_proton, 100, CS::WINKLER_withHYPwithNBAR,     NO_OUT )*1e-31;
        };
        static double dT_pHe_pbar_LAB_WinklerWitHypWitNbar ( double T_p, double T_pbar ){
            return CS::dE_pHe_pbar_LAB ( T_p+fMass_proton,   T_pbar+fMass_proton, 100, CS::WINKLER_withHYPwithNBAR,     NO_OUT )*1e-31;
        };
        static double dT_Hep_pbar_LAB_WinklerWitHypWitNbar ( double Tn_He, double T_pbar ){
            return CS::dE_Hep_pbar_LAB ( Tn_He+fMass_proton, T_pbar+fMass_proton, 100, CS::WINKLER_withHYPwithNBAR,     NO_OUT )*1e-31;
        };
        static double dT_HeHe_pbar_LAB_WinklerWitHypWitNbar ( double Tn_He, double T_pbar ){
            return CS::dE_HeHe_pbar_LAB ( Tn_He+fMass_proton, T_pbar+fMass_proton, 100, CS::WINKLER_withHYPwithNBAR,     NO_OUT )*1e-31;
        };
        
        
        static double dT_pp_Dbar_LAB  ( double T_p, double Tn_Dbar ){
            // factor 2 for dT to d(T/n)
            return CS::dE_pp_Dbar_LAB     ( T_p+fMass_proton,   Tn_Dbar*2+fMass_deuteron, 0.062, 1000, CS::DI_MAURO12, NO_OUT )*2*1e-31;
        };
        static double dT_ppbar_Dbar_LAB  ( double T_p, double Tn_Dbar ){
            // factor 2 for dT to d(T/n)
            return CS::dE_ppbar_Dbar_LAB  ( T_p+fMass_proton,   Tn_Dbar*2+fMass_deuteron, 0.062, 1000, CS::DEFAULT,    NO_OUT )*2*1e-31;
        };
        
        static double dT_pp_Hebar_LAB ( double T_p, double Tn_Hebar ){
            // factor 3 for dT to d(T/n)
            return CS::dE_pp_Hebar_LAB    ( T_p+fMass_proton,   Tn_Hebar*3+fMass_helium3, 0.078, 1000, CS::DI_MAURO12, NO_OUT )*3*1e-31;
        };
        
        static double dT_ppbar_Hebar_LAB ( double T_p, double Tn_Hebar ){
            // factor 3 for dT to d(T/n)
            return CS::dE_ppbar_Hebar_LAB ( T_p+fMass_proton,   Tn_Hebar*3+fMass_helium3, 0.078, 1000, CS::DEFAULT,    NO_OUT )*3*1e-31;
        };
        

        static void Set_dLog(int dLog){
            std::cout << warnout << "Careful! If you change dLog you are no longer able to read tables from the sofware." << warnend << std::endl;
            f_dLog  = dLog;
            f_i_max =  7*f_dLog;
            f_j_max =  8*f_dLog;
        }
        
    private:
        static CS_lab_tab* fInstance;
        CS_lab_tab();
        
        std::vector<FileTool*>   fCSfiles;
        std::vector<std::string> fCSfilesnames;
        std::vector<int>         fCStype;
        std::vector<int>         fCSparametrization;
        
        static int f_dLog;
        static int f_i_max;
        static int f_j_max;
        
        FileTool* f_file;
        
        
        
    };
    
    
}

//end LABCROSSSECTIONTABULATOR_H

#endif
