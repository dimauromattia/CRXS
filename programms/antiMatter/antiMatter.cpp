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


// Dbar production (by p)

double dT_pp_Dbar_LAB       ( double T_p,         double Tn_Dbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   Tn_Dbar, CS::PP_DBAR, CS::DI_MAURO12);
}
double dT_pHe_Dbar_LAB      ( double T_p,         double Tn_Dbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   Tn_Dbar, CS::PP_DBAR, CS::DI_MAURO12)*pow(4, 0.7);
}
double dT_Cp_Dbar_LAB      ( double T_p,         double Tn_Dbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   Tn_Dbar, CS::PP_DBAR, CS::DI_MAURO12)*pow(12,0.7);
}
double dT_Op_Dbar_LAB      ( double T_p,         double Tn_Dbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   Tn_Dbar, CS::PP_DBAR, CS::DI_MAURO12)*pow(16,0.7);
}
double dT_Hep_Dbar_LAB      ( double Tn_He,       double Tn_Dbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(Tn_He, Tn_Dbar, CS::PP_DBAR, CS::DI_MAURO12)*pow(4, 0.7);
}
double dT_HeHe_Dbar_LAB     ( double Tn_He,       double Tn_Dbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(Tn_He, Tn_Dbar, CS::PP_DBAR, CS::DI_MAURO12)*pow(4, 0.7+0.7);
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

// Hebar production (by p)

double dT_pp_Hebar_LAB       ( double T_p,         double Tn_Hebar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   Tn_Hebar, CS::PP_HEBAR, CS::DI_MAURO12);
}
double dT_pHe_Hebar_LAB      ( double T_p,         double Tn_Hebar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   Tn_Hebar, CS::PP_HEBAR, CS::DI_MAURO12)*pow(4, 0.7);
}
double dT_Cp_Hebar_LAB      ( double T_p,         double Tn_Hebar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   Tn_Hebar, CS::PP_HEBAR, CS::DI_MAURO12)*pow(12,0.7);
}
double dT_Op_Hebar_LAB      ( double T_p,         double Tn_Hebar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   Tn_Hebar, CS::PP_HEBAR, CS::DI_MAURO12)*pow(16,0.7);
}
double dT_Hep_Hebar_LAB      ( double Tn_He,       double Tn_Hebar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(Tn_He, Tn_Hebar, CS::PP_HEBAR, CS::DI_MAURO12)*pow(4, 0.7);
}
double dT_HeHe_Hebar_LAB     ( double Tn_He,       double Tn_Hebar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(Tn_He, Tn_Hebar, CS::PP_HEBAR, CS::DI_MAURO12)*pow(4, 0.7+0.7);
}


// Hebar production (by pbar)

double dT_pbarp_Hebar_LAB   ( double T_p,         double Tn_Hebar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   Tn_Hebar, CS::PPbar_HEBAR, CS::DEFAULT);
}
double dT_pbarHe_Hebar_LAB   ( double T_p,         double Tn_Hebar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   Tn_Hebar, CS::PPbar_HEBAR, CS::DEFAULT)*pow(4, 0.7);
}


// pbar production cross sections
double dT_pp_pbar_LAB_dM   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::DI_MAURO12)*2.3;
}
double dT_pHe_pbar_LAB_dM   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::DI_MAURO12)*2.3*pow(4, 0.7);;
}
double dT_Hep_pbar_LAB_dM   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::DI_MAURO12)*2.3*pow(4, 0.7);;
}
double dT_Cp_pbar_LAB_dM   ( double Tn_C,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(Tn_C,   T_pbar, CS::PP_PBAR, CS::DI_MAURO12)*2.3*pow(12, 0.7);;
}
double dT_Op_pbar_LAB_dM   ( double Tn_O,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(Tn_O,   T_pbar, CS::PP_PBAR, CS::DI_MAURO12)*2.3*pow(16, 0.7);;
}
double dT_HeHe_pbar_LAB_dM   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::DI_MAURO12)*2.3*pow(4, 0.7+0.7);
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
    
    
    if (step==-1) {
        CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
        cs_lab->WriteAll();
    }
    
    
    double nH   = 1e6;
    double nHe  = 0.1  * nH;
    
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    //cs_lab->ReadAll( "." );
    
    LISFluxes::fit();
    
    DM_energy_spectra* DM_spectra = DM_energy_spectra::GetInstance();
    

    out("Dbar   ")
    std::string dbar    = "Tn_Dbar        total          p H           p He          He p          He He         C p        O p         pbar H      pbar He";
    std::string dbar_DM = "Tn_Dbar        80 GeV, thermal WIMP";
    for (double dT=-1; dT<4; dT+=1./30) {
        double Tn_Dbar   = pow(10, dT);
        double Tn_from = std::min(Tn_Dbar, 6*fMass_proton);
        
        double S_p_H    = SourceTerms::integrationFromTo(Tn_Dbar, dT_pp_Dbar_LAB,   protonFlux, nH,  Tn_from, 1e6, p, o);
        double S_p_He   = SourceTerms::integrationFromTo(Tn_Dbar, dT_pHe_Dbar_LAB,  protonFlux, nHe, Tn_from, 1e6, p, o);
        double S_He_H   = SourceTerms::integrationFromTo(Tn_Dbar, dT_Hep_Dbar_LAB,  heliumFlux, nH,  Tn_from, 1e6, p, o);
        double S_He_He  = SourceTerms::integrationFromTo(Tn_Dbar, dT_HeHe_Dbar_LAB, heliumFlux, nHe, Tn_from, 1e6, p, o);
        double S_C_H    = 0; //SourceTerms::integrationFromTo(Tn_Dbar, dT_Cp_Dbar_LAB,   carbonFlux, nH,  Tn_from, 1e6, p, o);
        double S_O_H    = 0; //SourceTerms::integrationFromTo(Tn_Dbar, dT_Op_Dbar_LAB,   oxygenFlux, nH,  Tn_from, 1e6, p, o);
        
        
        double S_pbar_H    = SourceTerms::integrationFromTo(Tn_Dbar, dT_pbarp_Dbar_LAB,   antiprotonFlux, nH,  Tn_from, 1e6, p, o);
        double S_pbar_He   = SourceTerms::integrationFromTo(Tn_Dbar, dT_pbarHe_Dbar_LAB,  antiprotonFlux, nHe, Tn_from, 1e6, p, o);
        
        double S        = S_p_H + S_p_He + S_He_H + S_He_He + S_pbar_H + S_pbar_He + S_C_H + S_O_H;
        
        std::stringstream ss;
        ss  << "\n" << Tn_Dbar << "   " << S << "   " <<  S_p_H << "   " <<  S_p_He << "   " <<  S_He_H << "   " <<  S_He_He
        << "   " <<  S_C_H << "   " <<  S_O_H
        << "   " <<  S_pbar_H << "   " <<  S_pbar_He;
        dbar+= ss.str();

        double S_DM = DM_spectra->sourceTerm_Dbar(Tn_Dbar, 13, 80, 3e-26*1e-6, 0.43e+6, 0.062);
        std::stringstream sss;
        sss  << "\n" << Tn_Dbar << "     " << S_DM;
        dbar_DM+= sss.str();
    }
    
    
    FileTool::WriteStringToFile(dbar,    "Dbar_SourceTerm.txt");
    FileTool::WriteStringToFile(dbar_DM, "Dbar_SourceTerm_DM.txt");
    
    
    std::string toWrite = "";
    out("Dbar   energy distribution")
    toWrite    = "Tn_Dbar        dN/dE";
    for (double dT=-3; dT<4; dT+=1./30) {
        double Tn_Dbar   = pow(10, dT);
        
        
        double dTn_N = DM_spectra->dTn_N_Dbar(Tn_Dbar, 70.8);
        std::stringstream sss;
        sss  << "\n" << Tn_Dbar << "     " << dTn_N;
        toWrite+= sss.str();
    }
    
    FileTool::WriteStringToFile(toWrite,    "Dbar_dN_by_dTn.txt");
    
    
    toWrite    = "T        dN/dE";
    for (double dT=-3; dT<4; dT+=1./30) {
        double T   = pow(10, dT);
        double dT_N = DM_spectra->GetSpectrum( 70.8, T, 13 );
        std::stringstream sss;
        sss  << "\n" << T << "     " << dT_N;
        toWrite+= sss.str();
    }
    
    FileTool::WriteStringToFile(toWrite,    "pbar_dN_by_dT.txt");
    
    
    
    out("pp cross section")
    std::string XS_pp    = "T_p        pp_tot           pp_el";
    for (double dT=-3; dT<4; dT+=1./30) {
        double T_p   = pow(10, dT);
        
        double s = 4*fMass_proton*fMass_proton + 2 * fMass_proton * T_p;
        
        double tot = ppCSParametrizations::tot_pp__diMauro( s );
        double el  = ppCSParametrizations:: el_pp__diMauro( s );
        
        std::stringstream sss;
        sss  << "\n" << T_p << "     " <<  tot   << "     " <<  el  ;
        XS_pp+= sss.str();
    }
    
    FileTool::WriteStringToFile(XS_pp,    "XS_tot_pp.txt");
    
    
    out("pbar   di Mauro")
    std::string diMauro = "T_pbar        total          p H           p He          He p          He He         p C         p O";
    std::string pbar_DM = "T_pbar        DM thermal WIMP 80 GeV";
    for (double dT=-1; dT<4; dT+=1./30) {
        double T_pbar   = pow(10, dT);
        double Tn_from  = std::min(T_pbar, 6*fMass_proton);
        
        double S_p_H    = SourceTerms::integrationFromTo(T_pbar, dT_pp_pbar_LAB_dM,   protonFlux, nH,  Tn_from, 1e6, p, o);
        double S_p_He   = SourceTerms::integrationFromTo(T_pbar, dT_pHe_pbar_LAB_dM,  protonFlux, nHe, Tn_from, 1e6, p, o);
        double S_He_H   = SourceTerms::integrationFromTo(T_pbar, dT_Hep_pbar_LAB_dM,  heliumFlux, nH,  Tn_from, 1e6, p, o);
        double S_He_He  = SourceTerms::integrationFromTo(T_pbar, dT_HeHe_pbar_LAB_dM, heliumFlux, nHe, Tn_from, 1e6, p, o);
        double S_C_H    = 0; //SourceTerms::integrationFromTo(T_pbar, dT_Cp_pbar_LAB_dM,   carbonFlux, nH,  Tn_from, 1e6, p, o);
        double S_O_H    = 0; //SourceTerms::integrationFromTo(T_pbar, dT_Op_pbar_LAB_dM,   oxygenFlux, nH,  Tn_from, 1e6, p, o);
        double S        = S_p_H + S_p_He + S_He_H + S_He_He +S_C_H +S_O_H;
        
        std::stringstream ss;
        ss  << "\n" << T_pbar << "   " << S << "   " <<  S_p_H << "   " <<  S_p_He << "   " <<  S_He_H << "   " <<  S_He_He<< "   " <<  S_C_H<< "   " <<  S_O_H;
        diMauro+= ss.str();
        
        double S_DM = DM_spectra->sourceTerm_pbar(T_pbar, 13, 80, 3e-26*1e-6, 0.43e6);
        
        std::stringstream sss;
        sss  << "\n" << T_pbar << "     " << S_DM;
        pbar_DM+= sss.str();
    }
    FileTool::WriteStringToFile(diMauro, "pbar_SourceTerm_diMauro.txt");
    FileTool::WriteStringToFile(pbar_DM, "pbar_SourceTerm_DM.txt");

    
    
    
    out("Hebar   ")
    std::string Hebar    = "Tn_Hebar        total          p H           p He          He p          He He         C p        O p         pbar H      pbar He";
    std::string Hebar_DM = "Tn_Hebar        80 GeV, thermal WIMP";
    for (double dT=-1; dT<4; dT+=1./30) {
        double Tn_Hebar   = pow(10, dT);
        double Tn_from = std::min(Tn_Hebar, 6*fMass_proton);
        
        double S_p_H    = SourceTerms::integrationFromTo(Tn_Hebar, dT_pp_Hebar_LAB,   protonFlux, nH,  Tn_from, 1e6, p, o);
        double S_p_He   = SourceTerms::integrationFromTo(Tn_Hebar, dT_pHe_Hebar_LAB,  protonFlux, nHe, Tn_from, 1e6, p, o);
        double S_He_H   = SourceTerms::integrationFromTo(Tn_Hebar, dT_Hep_Hebar_LAB,  heliumFlux, nH,  Tn_from, 1e6, p, o);
        double S_He_He  = SourceTerms::integrationFromTo(Tn_Hebar, dT_HeHe_Hebar_LAB, heliumFlux, nHe, Tn_from, 1e6, p, o);
        double S_C_H    = 0; //SourceTerms::integrationFromTo(Tn_Hebar, dT_Cp_Hebar_LAB,   carbonFlux, nH,  Tn_from, 1e6, p, o);
        double S_O_H    = 0; //SourceTerms::integrationFromTo(Tn_Hebar, dT_Op_Hebar_LAB,   oxygenFlux, nH,  Tn_from, 1e6, p, o);
        
        
        double S_pbar_H    = SourceTerms::integrationFromTo(Tn_Hebar, dT_pbarp_Hebar_LAB,   antiprotonFlux, nH,  Tn_from, 1e6, p, o);
        double S_pbar_He   = SourceTerms::integrationFromTo(Tn_Hebar, dT_pbarHe_Hebar_LAB,  antiprotonFlux, nHe, Tn_from, 1e6, p, o);
        
        double S        = S_p_H + S_p_He + S_He_H + S_He_He + S_pbar_H + S_pbar_He + S_C_H + S_O_H;
        
        std::stringstream ss;
        ss  << "\n" << Tn_Hebar << "   " << S << "   " <<  S_p_H << "   " <<  S_p_He << "   " <<  S_He_H << "   " <<  S_He_He
        << "   " <<  S_C_H << "   " <<  S_O_H
        << "   " <<  S_pbar_H << "   " <<  S_pbar_He;
        Hebar+= ss.str();
        
        double S_DM = DM_spectra->sourceTerm_HeBar(Tn_Hebar, 13, 80, 3e-26*1e-6, 0.43e+6, 0.062);
        std::stringstream sss;
        sss  << "\n" << Tn_Hebar << "     " << S_DM;
        Hebar_DM+= sss.str();
    }
    
    
    FileTool::WriteStringToFile(Hebar,    "Hebar_SourceTerm.txt");
    FileTool::WriteStringToFile(Hebar_DM, "Hebar_SourceTerm_DM.txt");
    
    
    out("Hebar   energy distribution")
    toWrite    = "Tn_Hebar        dN/dE";
    for (double dT=-3; dT<4; dT+=1./30) {
        double Tn_Hebar   = pow(10, dT);
        
        
        double dTn_N = DM_spectra->dTn_N_Hebar(Tn_Hebar, 70.8);
        std::stringstream sss;
        sss  << "\n" << Tn_Hebar << "     " << dTn_N;
        toWrite+= sss.str();
    }
    
    FileTool::WriteStringToFile(toWrite,    "Hebar_dN_by_dTn.txt");
    
    out("He Spectrum")
    toWrite    = "Tn_He        Flux";
    for (double dT=-3; dT<4; dT+=1./30) {
        double Tn_He   = pow(10, dT);
        double Flux = LISFluxes::spectrum(Tn_He, LISFluxes::fHeliumLIS);
        
        std::stringstream sss;
        sss  << "\n" << Tn_He << "     " << Flux;
        toWrite+= sss.str();
    }
    
    FileTool::WriteStringToFile(toWrite,    "He_flux.txt");
    
    return 0;
    
    
}


