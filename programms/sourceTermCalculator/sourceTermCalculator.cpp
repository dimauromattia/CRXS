#include <iostream>
#include <sstream>
#include <fstream>

#include <definitions.h>
#include <configHandler.h>
#include <readFluxes.h>
#include <lisFluxes.h>
#include <crossSections.h>
#include <labCrossSectionTabulator.h>
#include <sourceTerm.h>
#include <dmEnergySpectra.h>

#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TApplication.h"

#include "QString"
#include "QStringList"

using namespace  CRACS;

std::string fProgramDescription = "Program to fit the LIS spectra of proton, helium, and antiproton fluxes";



// pbar production



// pbar production cross sections
double dT_pp_pbar_LAB_dM   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::DI_MAURO12)*2.3;
}
double dT_pHe_pbar_LAB_dM   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PHe_PBAR, CS::DI_MAURO12)*2.3;
}
double dT_Hep_pbar_LAB_dM   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::HeP_PBAR, CS::DI_MAURO12)*2.3;
}
double dT_HeHe_pbar_LAB_dM   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::DI_MAURO12)*2.3*pow(4, 0.7+0.7);
}

double dT_pp_pbar_LAB_Winkler   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::WINKLER_withHYPwithNBAR);
}
double dT_pHe_pbar_LAB_Winkler   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PHe_PBAR, CS::WINKLER_withHYPwithNBAR);
}
double dT_Hep_pbar_LAB_Winkler   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::HeP_PBAR, CS::WINKLER_withHYPwithNBAR);
}
double dT_HeHe_pbar_LAB_Winkler   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::WINKLER_withHYPwithNBAR)*pow(4, 0.7+0.7);
}


double dT_pp_pbar_LAB_TanNg   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::TAN_NG)*2.3;
}
double dT_pHe_pbar_LAB_TanNg   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PHe_PBAR, CS::TAN_NG)*2.3;
}
double dT_Hep_pbar_LAB_TanNg   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::HeP_PBAR, CS::TAN_NG)*2.3;
}
double dT_HeHe_pbar_LAB_TanNg   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::TAN_NG)*2.3*pow(4, 0.7+0.7);
}



double dT_pp_pbar_LAB_KR   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::KACHELRIESS)*2.3;
}
double dT_pHe_pbar_LAB_KR   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PHe_PBAR, CS::KACHELRIESS)*2.3;
}
double dT_Hep_pbar_LAB_KR   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::HeP_PBAR, CS::KACHELRIESS)*2.3;
}
double dT_HeHe_pbar_LAB_KR   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::HeHe_PBAR, CS::KACHELRIESS)*2.3;
}



double dT_pp_pbar_LAB_Dup   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::DUPERRAY)*2.3;
}
double dT_pHe_pbar_LAB_Dup   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PHe_PBAR, CS::DUPERRAY)*2.3;
}
double dT_Hep_pbar_LAB_Dup   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::HeP_PBAR, CS::DUPERRAY)*2.3;
}
double dT_HeHe_pbar_LAB_Dup   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::DUPERRAY)*2.3*pow(4, 0.7+0.7);
}



double dT_pp_pbar_LAB_DTUNUC   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::DTUNUC)*2.3;
}
double dT_pHe_pbar_LAB_DTUNUC   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PHe_PBAR, CS::DTUNUC)*2.3;
}
double dT_Hep_pbar_LAB_DTUNUC   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::HeP_PBAR, CS::DTUNUC)*2.3;
}
double dT_HeHe_pbar_LAB_DTUNUC   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::DTUNUC)*2.3;
}


// Dbar production (by p)

double dT_pp_Dbar_LAB_dM   ( double T_p,         double Tn_Dbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   Tn_Dbar, CS::PP_DBAR, CS::DI_MAURO12);
}
double dT_pHe_Dbar_LAB_dM   ( double T_p,         double Tn_Dbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   Tn_Dbar, CS::PP_DBAR, CS::DI_MAURO12)*pow(4, 0.7);
}
double dT_Hep_Dbar_LAB_dM   ( double T_p,         double Tn_Dbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   Tn_Dbar, CS::PP_DBAR, CS::DI_MAURO12)*pow(4, 0.7);
}
double dT_HeHe_Dbar_LAB_dM   ( double T_p,         double Tn_Dbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   Tn_Dbar, CS::PP_DBAR, CS::DI_MAURO12)*pow(4, 0.7+0.7);
}


// Dbar production (by pbar)

double dT_pbarp_Dbar_LAB   ( double T_p,         double Tn_Dbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   Tn_Dbar, CS::PPbar_DBAR, CS::DEFAULT);
}
double dT_pbarHe_Dbar_LAB   ( double T_p,         double Tn_Dbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   Tn_Dbar, CS::PPbar_DBAR, CS::DEFAULT)*pow(4, 0.7);
}


// Fluxes
double protonFlux(double T_p){
    return LISFluxes::spectrum(T_p, LISFluxes::fProtonLIS);
}
double heliumFlux(double T_He){
    return LISFluxes::spectrum(T_He, LISFluxes::fHeliumLIS);
}
double carbonFlux(double T_C){
    return LISFluxes::spectrum(T_C, LISFluxes::fCarbonLIS);
}
double oxygenFlux(double T_O){
    return LISFluxes::spectrum(T_O, LISFluxes::fOxygenLIS);
}
double antiprotonFlux(double T_pbar){
    return LISFluxes::spectrum(T_pbar, LISFluxes::fAntiprotonLIS);
}



//double proton_f(double T_p){
//    return LISFluxes::spectrum(T_p, LISFluxes::fProtonLIS);
//}

int main(int argc, char *argv[])
{
    // Config handling
    ConfigHandler* config = ConfigHandler::GetInstance();
    config->SetInput(argc, argv, fProgramDescription);
    
    int     step    = 1;
    int     o       = WARN_OUT;
    int     p       = 10000;
    config->AddOptionInt    (   "step",         step,       "Step. Default: 1"                      );
    config->AddOptionInt    (   "o",            o,          "Output. Default: 0"                    );
    config->AddOptionInt    (   "p",            p,          "Precison. Default: 10000"              );
    
    config->CheckOption();
    
    double nH   = 1e+6;
    double nHe  = 0.1  * nH;
    
    LISFluxes::fit();
    
    DM_energy_spectra* DM_spectra = DM_energy_spectra::GetInstance();
    
    //out( SourceTerms::integrationFromTo(10, dT_pp_pbar_LAB_DTUNUC,  protonFlux, nH, 1, 1e6, p, o) )
    //out( SourceTerms::integrationFromTo(10, dT_pp_pbar_LAB_TanNg,   protonFlux, nH, 1, 1e6, p, o) )
    
    //out( dT_pp_pbar_LAB_DTUNUC(100, 10) );
    //out( dT_pp_pbar_LAB_TanNg (100, 10) );
    
    //out( SourceTerms::integrationFromTo( 1, dT_pp_pbar_LAB_DTUNUC,  protonFlux, nH, 1, 1e6, p, o) )
    //out( SourceTerms::integrationFromTo( 1, dT_pp_pbar_LAB_TanNg,   protonFlux, nH, 1, 1e6, p, o) )
    
    //out( dT_pp_pbar_LAB_DTUNUC(20, 1) );
    //out( dT_pp_pbar_LAB_TanNg (20, 1) );
    
    //return 0;
    
    out("pbar   di Mauro")
    std::string diMauro = "T_pbar        total          p H           p He          He p          He He";
    std::string pbar_DM = "T_pbar        DM thermal WIMP 100 GeV";
    for (double dT=-1; dT<4; dT+=1./30) {
        double T_pbar   = pow(10, dT);
        double Tn_from  = std::min(T_pbar, 6*fMass_proton);
        
        double S_p_H    = SourceTerms::integrationFromTo(T_pbar, dT_pp_pbar_LAB_dM,   protonFlux, nH,  Tn_from, 1e6, p, o);
        double S_p_He   = SourceTerms::integrationFromTo(T_pbar, dT_pHe_pbar_LAB_dM,  protonFlux, nHe, Tn_from, 1e6, p, o);
        double S_He_H   = SourceTerms::integrationFromTo(T_pbar, dT_Hep_pbar_LAB_dM,  heliumFlux, nH,  Tn_from, 1e6, p, o);
        double S_He_He  = SourceTerms::integrationFromTo(T_pbar, dT_HeHe_pbar_LAB_dM, heliumFlux, nHe, Tn_from, 1e6, p, o);
        double S        = S_p_H + S_p_He + S_He_H + S_He_He;
        
        std::stringstream ss;
        ss  << "\n" << T_pbar << "   " << S << "   " <<  S_p_H << "   " <<  S_p_He << "   " <<  S_He_H << "   " <<  S_He_He;
        diMauro+= ss.str();
        
        double S_DM = DM_spectra->sourceTerm_pbar(T_pbar);
        std::stringstream sss;
        sss  << "\n" << T_pbar << "     " << S_DM;
        pbar_DM+= sss.str();
    }
    FileTool::WriteStringToFile(diMauro, "pbar_SourceTerm_diMauro.txt");
    FileTool::WriteStringToFile(pbar_DM, "pbar_SourceTerm_DM.txt");
    
    
    out("pbar   DTUNUC")
    std::string DTUNUC = "T_pbar        total          p H           p He          He p          He He";
    for (double dT=-1; dT<2; dT+=1./30) {
        double T_pbar   = pow(10, dT);
        double Tn_from  = std::min(T_pbar, 6*fMass_proton);
        
        double S_p_H    = SourceTerms::integrationFromTo(T_pbar, dT_pp_pbar_LAB_DTUNUC,   protonFlux, nH,  Tn_from, 1e6, p, o);
        double S_p_He   = SourceTerms::integrationFromTo(T_pbar, dT_pHe_pbar_LAB_DTUNUC,  protonFlux, nHe, Tn_from, 1e6, p, o);
        double S_He_H   = SourceTerms::integrationFromTo(T_pbar, dT_Hep_pbar_LAB_DTUNUC,  heliumFlux, nH,  Tn_from, 1e6, p, o);
        double S_He_He  = SourceTerms::integrationFromTo(T_pbar, dT_HeHe_pbar_LAB_DTUNUC, heliumFlux, nHe, Tn_from, 1e6, p, o);
        double S        = S_p_H + S_p_He + S_He_H + S_He_He;
        
        std::stringstream ss;
        ss  << "\n" << T_pbar << "   " << S << "   " <<  S_p_H << "   " <<  S_p_He << "   " <<  S_He_H << "   " <<  S_He_He;
        DTUNUC+= ss.str();
    }
    FileTool::WriteStringToFile(DTUNUC, "pbar_SourceTerm_DTUNUC.txt");
    
    out("pbar   Duperray")
    std::string Duperray = "T_pbar        total          p H           p He          He p          He He";
    for (double dT=-1; dT<4; dT+=1./30) {
        double T_pbar   = pow(10, dT);
        double Tn_from  = std::min(T_pbar, 6*fMass_proton);
        
        double S_p_H    = SourceTerms::integrationFromTo(T_pbar, dT_pp_pbar_LAB_Dup,   protonFlux, nH,  Tn_from, 1e6, p, o);
        double S_p_He   = SourceTerms::integrationFromTo(T_pbar, dT_pHe_pbar_LAB_Dup,  protonFlux, nHe, Tn_from, 1e6, p, o);
        double S_He_H   = SourceTerms::integrationFromTo(T_pbar, dT_Hep_pbar_LAB_Dup,  heliumFlux, nH,  Tn_from, 1e6, p, o);
        double S_He_He  = SourceTerms::integrationFromTo(T_pbar, dT_HeHe_pbar_LAB_Dup, heliumFlux, nHe, Tn_from, 1e6, p, o);
        double S        = S_p_H + S_p_He + S_He_H + S_He_He;
        
        std::stringstream ss;
        ss  << "\n" << T_pbar << "   " << S << "   " <<  S_p_H << "   " <<  S_p_He << "   " <<  S_He_H << "   " <<  S_He_He;
        Duperray+= ss.str();
    }
    FileTool::WriteStringToFile(Duperray, "pbar_SourceTerm_Duperray.txt");
    
    out("pbar   Tan Ng")
    std::string TanNg = "T_pbar        total          p H           p He          He p          He He";
    for (double dT=-1; dT<4; dT+=1./30) {
        double T_pbar   = pow(10, dT);
        double Tn_from  = std::min(T_pbar, 6*fMass_proton);
        
        double S_p_H    = SourceTerms::integrationFromTo(T_pbar, dT_pp_pbar_LAB_TanNg,   protonFlux, nH,  Tn_from, 1e6, p, o);
        double S_p_He   = SourceTerms::integrationFromTo(T_pbar, dT_pHe_pbar_LAB_TanNg,  protonFlux, nHe, Tn_from, 1e6, p, o);
        double S_He_H   = SourceTerms::integrationFromTo(T_pbar, dT_Hep_pbar_LAB_TanNg,  heliumFlux, nH,  Tn_from, 1e6, p, o);
        double S_He_He  = SourceTerms::integrationFromTo(T_pbar, dT_HeHe_pbar_LAB_TanNg, heliumFlux, nHe, Tn_from, 1e6, p, o);
        double S        = S_p_H + S_p_He + S_He_H + S_He_He;
        
        std::stringstream ss;
        ss  << "\n" << T_pbar << "   " << S << "   " <<  S_p_H << "   " <<  S_p_He << "   " <<  S_He_H << "   " <<  S_He_He;
        TanNg+= ss.str();
    }
    FileTool::WriteStringToFile(TanNg, "pbar_SourceTerm_TanNg.txt");
    
    
    out("pbar   Kachelriess")
    std::string KR = "T_pbar        total          p H           p He          He p          He He";
    for (double dT=-1; dT<4; dT+=1./30) {
        double T_pbar   = pow(10, dT);
        double Tn_from = std::min(T_pbar, 6*fMass_proton);
        
        double S_p_H    = SourceTerms::integrationFromTo(T_pbar, dT_pp_pbar_LAB_KR,   protonFlux, nH,  Tn_from, 1e6, p, o);
        double S_p_He   = SourceTerms::integrationFromTo(T_pbar, dT_pHe_pbar_LAB_KR,  protonFlux, nHe, Tn_from, 1e6, p, o);
        double S_He_H   = SourceTerms::integrationFromTo(T_pbar, dT_Hep_pbar_LAB_KR,  heliumFlux, nH,  Tn_from, 1e6, p, o);
        double S_He_He  = SourceTerms::integrationFromTo(T_pbar, dT_HeHe_pbar_LAB_KR, heliumFlux, nHe, Tn_from, 1e6, p, o);
        double S        = S_p_H + S_p_He + S_He_H + S_He_He;
        
        std::stringstream ss;
        ss  << "\n" << T_pbar << "   " << S << "   " <<  S_p_H << "   " <<  S_p_He << "   " <<  S_He_H << "   " <<  S_He_He;
        KR+= ss.str();
    }
    FileTool::WriteStringToFile(KR, "pbar_SourceTerm_KR.txt");
    
    
    out("pbar   Winkler")
    std::string Winkler = "T_pbar        total          p H           p He          He p          He He";
    for (double dT=-1; dT<4; dT+=1./30) {
        double T_pbar   = pow(10, dT);
        double Tn_from = std::min(T_pbar, 6*fMass_proton);
        
        double S_p_H    = SourceTerms::integrationFromTo(T_pbar, dT_pp_pbar_LAB_Winkler,   protonFlux, nH,  Tn_from, 1e6, p, o);
        double S_p_He   = SourceTerms::integrationFromTo(T_pbar, dT_pHe_pbar_LAB_Winkler,  protonFlux, nHe, Tn_from, 1e6, p, o);
        double S_He_H   = SourceTerms::integrationFromTo(T_pbar, dT_Hep_pbar_LAB_Winkler,  heliumFlux, nH,  Tn_from, 1e6, p, o);
        double S_He_He  = SourceTerms::integrationFromTo(T_pbar, dT_HeHe_pbar_LAB_Winkler, heliumFlux, nHe, Tn_from, 1e6, p, o);
        double S        = S_p_H + S_p_He + S_He_H + S_He_He;
        
        std::stringstream ss;
        ss  << "\n" << T_pbar << "   " << S << "   " <<  S_p_H << "   " <<  S_p_He << "   " <<  S_He_H << "   " <<  S_He_He;
        Winkler+= ss.str();
    }
    FileTool::WriteStringToFile(Winkler, "pbar_SourceTerm_Winkler.txt");
    
    
    
    out("Dbar   ")
    std::string dbar    = "Tn_Dbar        total          p H           p He          He p          He He";
    std::string dbar_DM = "Tn_Dbar        100 GeV, thermal WIMP";
    for (double dT=-1; dT<4; dT+=1./30) {
        double Tn_Dbar   = pow(10, dT);
        double Tn_from = std::min(Tn_Dbar, 6*fMass_proton);
        
        double S_p_H    = SourceTerms::integrationFromTo(Tn_Dbar, dT_pp_Dbar_LAB_dM,   protonFlux, nH,  Tn_from, 1e6, p, o);
        double S_p_He   = SourceTerms::integrationFromTo(Tn_Dbar, dT_pHe_Dbar_LAB_dM,  protonFlux, nHe, Tn_from, 1e6, p, o);
        double S_He_H   = SourceTerms::integrationFromTo(Tn_Dbar, dT_Hep_Dbar_LAB_dM,  heliumFlux, nH,  Tn_from, 1e6, p, o);
        double S_He_He  = SourceTerms::integrationFromTo(Tn_Dbar, dT_HeHe_Dbar_LAB_dM, heliumFlux, nHe, Tn_from, 1e6, p, o);
        
        double S_pbar_H    = SourceTerms::integrationFromTo(Tn_Dbar, dT_pbarp_Dbar_LAB,   antiprotonFlux, nH,  Tn_from, 1e6, p, o);
        double S_pbar_He   = SourceTerms::integrationFromTo(Tn_Dbar, dT_pbarHe_Dbar_LAB,  antiprotonFlux, nHe, Tn_from, 1e6, p, o);
        
        double S        = S_p_H + S_p_He + S_He_H + S_He_He + S_pbar_H + S_pbar_He;
        
        std::stringstream ss;
        ss  << "\n" << Tn_Dbar << "   " << S << "   " <<  S_p_H << "   " <<  S_p_He << "   " <<  S_He_H << "   " <<  S_He_He << "   " <<  S_pbar_H << "   " <<  S_pbar_He;
        dbar+= ss.str();
        
        double S_DM = DM_spectra->sourceTerm_Dbar(Tn_Dbar);
        std::stringstream sss;
        sss  << "\n" << Tn_Dbar << "     " << S_DM;
        dbar_DM+= sss.str();
    }
    
    FileTool::WriteStringToFile(dbar,    "Dbar_SourceTerm.txt");
    FileTool::WriteStringToFile(dbar_DM, "Dbar_SourceTerm_DM.txt");
    
    
    
    return 0;
    
    
}


