#include <iostream>

#include <configHandler.h>
#include <fileTools.h>


#include <definitions.h>
#include <crossSections.h>
#include <labCrossSectionTabulator.h>

using namespace  CRACS;

std::string fProgramDescription = "Test Antiproton production cross sections";


double dT_pp_Dupperay_lab(double T_p, double T_pbar){
    return CS::dE_pp_pbar_LAB (T_p + fMass_proton, T_pbar+fMass_proton, 10000, CS::DUPERRAY)*1e-31;
}

static double dT_pp_pbar_LAB_diMauro12  ( double T_p, double T_pbar ){
    return CS::dE_pp_pbar_LAB  ( T_p+fMass_proton,   T_pbar+fMass_proton, 10000, CS::DI_MAURO12, WARN_OUT )*1e-31;
}
static double dT_pp_pbar_LAB_diMauro  ( double T_p, double T_pbar ){
    return CS::dE_pp_pbar_LAB  ( T_p+fMass_proton,   T_pbar+fMass_proton, 10000, CS::DI_MAURO,   WARN_OUT )*1e-31;
}



int main(int argc, char *argv[])
{
    // Config handling
    ConfigHandler* config = ConfigHandler::GetInstance();
    config->SetInput(argc, argv, fProgramDescription);
    
    int step = 13;
    config->AddOptionInt     (  "step",        step,       "Step  \n\t\t\t 1: tot_pp__diMauro(s)     \n\t\t\t 1: dE_pp_p_LAB(T_p, T_p)     \n\t\t\t 3: dE_pp_pbar_LAB(T_p, T_pbar)     \n\t\t\t 4: dE_pp_Dbar_LAB(T_p, T_Dbar)     \n\t\t\t 13: Comparison table di Mauro (Eq. 13)     Default!");
    double  s;
    double  T_p;
    double  T_product;
    double  pT;
    double  cosTheta = 0;
    int     n        = 100000;
    int     o        = 1;
    config->AddOptionDouble  (  "s",            s           );
    config->AddOptionDouble  (  "T_p",          T_p         );
    config->AddOptionDouble  (  "T_product",    T_product   );
    config->AddOptionDouble  (  "pT_product",   pT          );
    config->AddOptionDouble  (  "cosTheta",     cosTheta    );
    config->AddOptionInt     (  "n",            n           );
    config->AddOptionInt     (  "o",            o           );
    config->CheckOption();
    
    varOut(step)
    
    out( CRACS::ppCSParametrizations::tot_pp__diMauro( 49e6 ) )
    
    if (step==-1) {
        CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
        //cs_lab->TransformAll();
        cs_lab->WriteAll();
    }

    if (step==1) {
        CRACS::CS_lab_tab::GetInstance()->TransformAllGalprop();
    }
   
    if (step==-100) {
        std::cout << ppCSParametrizations::inv_pp_pbar_CM__Winkler(s, T_product+fMass_proton, pT+fMass_proton) << std::endl;
    }
    
    if (step==2) {
        std::cout << CS::dE_pp_p_LAB    (T_p+fMass_proton, T_product+fMass_proton,   n, CS::ANDERSON, o) << std::endl;
    }
    if (step==3) {
        std::cout << CS::dE_pp_pbar_LAB (T_p+fMass_proton, T_product+fMass_proton,   n, CS::DI_MAURO, o) << std::endl;
    }
    if (step==4) {
        std::cout << CS::dE_pp_Dbar_LAB    (T_p+fMass_proton, T_product+fMass_deuteron, 0.078, n, CS::DI_MAURO, o) << std::endl;
    }
    if (step==5) {
        std::cout << CS::dE_ppbar_Dbar_LAB (T_p+fMass_proton, T_product+fMass_deuteron, 0.078, n, CS::DEFAULT, o) << std::endl;
    }
    if (step==6){
        //std::cout <<  CS::inv_pp_Dbar_LAB(T_p +fMass_proton, T_product+fMass_deuteron, cosTheta, 0.078, CS::DI_MAURO12) << std::endl;
    }
    
//    if (step==10) {
//        for (double dT=0; dT<=5; dT+=0.05) {
//            double T_p   = pow(10, dT);
//            for (double dx=-dT-2; dx<=0; dx+=0.05) {
//                double T_pbar    = T_p * pow(10, dx);
//                double cs_diMauro12  = CS::dE_pp_pbar_LAB(T_p+fMass_proton, T_pbar+fMass_proton, n, CS::DI_MAURO12, o);
//                double cs_diMauro13  = CS::dE_pp_pbar_LAB(T_p+fMass_proton, T_pbar+fMass_proton, n, CS::DI_MAURO,   o);
//                double cs_Tan_Ng     = CS::dE_pp_pbar_LAB(T_p+fMass_proton, T_pbar+fMass_proton, n, CS::TAN_NG,     o);
//                std:: cout << T_p << "    " << T_pbar <<  "    " << cs_diMauro12 <<  "    " << cs_diMauro13 <<  "    " << cs_Tan_Ng << std::endl;
//            }
//        }
//    }
    
    
    
//    if (step==12) {
//        CRACS::FileTool f("sigma_bestfitalldef_eq12.dat");
//        f.ExtractNumberTable(3, " ");
//        for (int i=0; i<f.GetNumberOfLines()-1; i++) {
//            double Ekin_p       = f.NumberTable(i, 0);
//            double Ekin_pbar    = f.NumberTable(i, 1);
//            double cs_serpico   = f.NumberTable(i, 2)*1e31;
//            double cs_my        = CS::dE_pp_pbar_LAB(Ekin_p+fMass_proton, Ekin_pbar+fMass_proton, n, CS::DI_MAURO12, o);
//            std:: cout << Ekin_p << "    " << Ekin_pbar << "    " << cs_serpico << "    " << cs_my << std::endl;
//        }
//    }
//    
//    
//    if (step==13) {
//        CRACS::FileTool f("sigma_bestfitalldef_eq13.dat");
//        f.ExtractNumberTable(3, " ");
//        for (int i=0; i<f.GetNumberOfLines()-1; i++) {
//            double Ekin_p       = f.NumberTable(i, 0);
//            double Ekin_pbar    = f.NumberTable(i, 1);
//            double cs_serpico   = f.NumberTable(i, 2)*1e31;
//            double cs_my        = CS::dE_pp_pbar_LAB(Ekin_p+fMass_proton, Ekin_pbar+fMass_proton, n, CS::DI_MAURO, o);
//            std:: cout << Ekin_p << "    " << Ekin_pbar << "    " << cs_serpico << "    " << cs_my << std::endl;
//        }
//    }
//    
//    if (step==20) {
//        for (double dT=1; dT<=5; dT+=0.05) {
//            double T_p   = pow(10, dT);
//            for (double dx=-dT-2; dx<=0; dx+=0.05) {
//                double T_Dbar    = T_p * pow(10, dx);
//                double cs        = CS::dE_pp_Dbar_LAB(T_p+fMass_proton, T_Dbar+fMass_deuteron, 0.078, n, CS::DI_MAURO12, o);
//                std:: cout << T_p << "    " << T_Dbar <<  "    " << cs << std::endl;
//            }
//        }
//    }
//    
//    if (step==21) {
//        
//        for (double dx=-7; dx<=0; dx+=0.01) {
//            double T_Dbar    = T_p * pow(10, dx);
//            double cs        = CS::dE_pp_Dbar_LAB(T_p+fMass_proton, T_Dbar+fMass_deuteron, 0.078, n, CS::DI_MAURO12, o);
//            std:: cout << T_p << "    " << T_Dbar << "    " << cs << std::endl;
//            
//        }
//    }
    
//    if (step==22) {
//        double dCos = 0.0001;
//        
//        std::cout << "T_D\\cosTheta" << "  ";
//        for (double cosTheta=-1; cosTheta<=1; cosTheta+=dCos) {
//            std::cout << cosTheta << "     ";
//        }
//        std::cout << std::endl;
//        for (double dT=-3; dT<=0; dT+=0.05) {
//            double T_D   = T_p*pow(10, dT);
//            std::cout << T_D << "  ";
//            for (double cosTheta=-1; cosTheta<=1; cosTheta+=dCos) {
//                double cs = CS::d3p_pp_Dbar_LAB(T_p +fMass_proton, T_D+fMass_deuteron, cosTheta, 0.078, CS::DI_MAURO12);
//                std::cout << cs << "     ";
//            }
//            std::cout << std::endl;
//        }
//    }
    
    
//    if (step==30) {
//        for (double dT=0; dT<=5; dT+=0.05) {
//            double T_pbar   = pow(10, dT);
//            for (double dx=-dT-2; dx<=0; dx+=0.05) {
//                double T_Dbar    = T_pbar * pow(10, dx);
//                //double cs        = CS::dE_ppbar_Dbar_LAB(T_pbar+fMass_proton, T_Dbar+fMass_deuteron, 0.078, n, CS::DEFAULT, o);
//                
//                double cs        = CS::dE_ppbar_Dbar_LAB(T_pbar+fMass_proton, T_Dbar+fMass_deuteron, 0.078, n, CS::DEFAULT, o);
//                std:: cout << T_pbar << "    " << T_Dbar << "    " << cs << std::endl;
//            }
//        }
//    }
//    if (step==31) {
//        double T_pbar   = T_p;
//        for (double dx=-2; dx<=log10(T_pbar); dx+=0.1) {
//            double T_Dbar    = pow(10, dx);
//            double cs        = CS::dE_ppbar_Dbar_LAB(T_p+fMass_proton, T_Dbar+fMass_deuteron, 0.078, n, CS::DEFAULT, o);
//            std:: cout << T_pbar << "    " << T_Dbar << "    " << cs << std::endl;
//            
//        }
//        
//    }
//    
//    if (step==100) {
//        CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
//        cs_lab->test();
//    }
//    
//    if (step==200) {
//        
//        CRACS::CS_lab_tab::GetInstance()->TransformAllGalprop();
//        
//        /*CRACS::FileTool f("sigma_bestfitalldef_eq12.dat");
//        f.ExtractNumberTable(3, " ");
//        for (int i=0; i<f.GetNumberOfLines()-1; i++) {
//            double Ekin_p       = f.NumberTable(i, 0);
//            double Ekin_pbar    = f.NumberTable(i, 1);
//            //double cs_serpico   = f.NumberTable(i, 2);
//            double cs_pp        = CS_lab_tab::GetInstance()->GetCrossSection(Ekin_p+fMass_proton, Ekin_pbar+fMass_proton, CS::PP_PBAR,  CS::KACHELRIESS);
//            double cs_Hep       = CS_lab_tab::GetInstance()->GetCrossSection(Ekin_p+fMass_proton, Ekin_pbar+fMass_proton, CS::HeP_PBAR, CS::KACHELRIESS);
//            double cs_pHe       = CS_lab_tab::GetInstance()->GetCrossSection(Ekin_p+fMass_proton, Ekin_pbar+fMass_proton, CS::PHe_PBAR, CS::KACHELRIESS);
//            
//            std:: cout << Ekin_p << "    " << Ekin_pbar << "    "  << cs_pp << "    " << cs_Hep << "    " << cs_pHe << std::endl; //cs_serpico*2.3 << "    "
//        }*/
//    }
//    
//    if (step==300) {
//        CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
//        cs_lab->WriteCS(dT_pp_Dupperay_lab, "dT_pp_DUPPERAY.txt", "./", 10);
//    }
//    
//    if (step==500) {
//        out(step)
//        
//        
//        double i1 = CRACS::CSTransformations::dE_pp_product_LAB       (1e6, 1e3, fMass_proton, CRACS::ppCSParametrizations::inv_pp_pbar_CM__diMauro12, 3*fMass_proton);
//        double i2 = CRACS::CSTransformations::dE_pp_product_LAB_intEta(1e6, 1e3, fMass_proton, CRACS::ppCSParametrizations::inv_pp_pbar_CM__diMauro12, 10000, ALL_OUT  );
//        
//        out(i1)
//        out(i2)
//        
//        CRACS::CS_lab_tab* lab = CRACS::CS_lab_tab::GetInstance();
//        
//        lab->WriteCS( dT_pp_pbar_LAB_diMauro12,   "dT_pp_pbar_LAB_diMauro12.txt"  );
//        lab->WriteCS( dT_pp_pbar_LAB_diMauro,     "dT_pp_pbar_LAB_diMauro.txt"    );
//        
//    }
    
    return 0;
    
}


