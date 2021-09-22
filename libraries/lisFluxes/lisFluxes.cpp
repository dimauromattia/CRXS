#include "lisFluxes.h"

#include <iostream>
#include <sstream>

#include <definitions.h>
#include <configHandler.h>
#include <fileTools.h>
#include <readFluxes.h>
#include <graph.h>

#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TMinuit.h"

double CRACS::LISFluxes::fSM      = 0.6;
//double CRACS::LISFluxes::fSM_lower     = 0.4;
//double CRACS::LISFluxes::fSM_upper     = 0.7;

CRACS::Graph CRACS::LISFluxes::fProton_LIS;
CRACS::Graph CRACS::LISFluxes::fHelium_LIS;
CRACS::Graph CRACS::LISFluxes::fCarbon_LIS;
CRACS::Graph CRACS::LISFluxes::fOxygen_LIS;
CRACS::Graph CRACS::LISFluxes::fNitrogen_LIS;

CRACS::Graph CRACS::LISFluxes::fLithium_LIS;
CRACS::Graph CRACS::LISFluxes::fBeryllium_LIS;
CRACS::Graph CRACS::LISFluxes::fBoron_LIS;

CRACS::Graph CRACS::LISFluxes::fAntiproton_LIS;

double       CRACS::LISFluxes::fParameter[100];
int          CRACS::LISFluxes::fNParam;
double       CRACS::LISFluxes::fStart[100];
double       CRACS::LISFluxes::fStep [100];

void CRACS::LISFluxes::applySolarModulation(CRACS::Graph& flux, double SM_potential, int A, int Z, double power, int type){
    
    if (type!=RIGIDITY && type!=EKINPERN)
    return;
    if (type==RIGIDITY)
    flux.ConvertSpectrumFromRigidityToEkinPerN(Z, A, power);
    
    for (int i=0; i<flux.Size(); i++) {
        
        double E_LIS;
        double phi_LIS;
        double E_LIS_err;
        double phi_LIS_err;
        
        flux.GetPoint(i, E_LIS, phi_LIS, E_LIS_err, phi_LIS_err);
        
        
        phi_LIS_err  *=  pow(E_LIS, -power);
        phi_LIS      *=  pow(E_LIS, -power);
        
        //Force field approximation
        double E        =   E_LIS - fabs(1.*Z)/A*SM_potential;
        double E_err    =   E_LIS_err;
        double factor   =   (E*E+2*E*fMass_proton)/(E_LIS*E_LIS+2*E_LIS*fMass_proton);
        double phi      =   phi_LIS*factor;
        double phi_err  =   phi_LIS_err*factor;
        if (E<=0)
        continue;
        
        phi      *=  pow(E, power);
        phi_err  *=  pow(E, power);
        
        
        flux.SetPoint(i, E, phi, E_err, phi_err);
        
    }
    
    if (type==RIGIDITY)
    flux.ConvertSpectrumFromRigidityToEkinPerN(Z, A, power);
    
    return;
}

double CRACS::LISFluxes::brokenPowerLaw(double E, double norm, double gamma1, double gamma2, double gamma3, double break1, double break2, double s1){
    
    if(E<=break2){
        
        double delta21 = gamma2-gamma1;
        return norm  *  pow(  E,                                -gamma1      )   /   pow( break1,            -gamma1  ) *
        pow(  pow(E, 1./s1)+pow(break1, 1./s1), -s1*delta21  )   /   pow( pow(2, s1)*break1, -delta21 );
        
    }
    if(E!=E)
    return NAN;
    return     brokenPowerLaw(break2, norm, gamma1, gamma2, gamma3, break1, break2, s1)  *
    pow(  E,                                -gamma3     )   /   pow( break2,            -gamma3 );
    
}

double CRACS::LISFluxes::brokenPowerLaw_delta(double E, double norm, double gamma, double delta1, double delta2, double break1, double break2, double s1){
    
    return brokenPowerLaw(E, norm, gamma-delta1, gamma, gamma+delta2, break1, break2, s1);
    
}


double CRACS::LISFluxes::brokenPowerLaw(double *x, double *par)
{
    double E        = x  [0];
    double norm     = par[0];
    double gamma1   = par[1];
    double gamma2   = par[2];
    double gamma3   = par[3];
    double break1   = par[4];
    double break2   = par[5];
    double s1       = par[6];
    
    return brokenPowerLaw(E, norm, gamma1, gamma2, gamma3, break1, break2, s1);
}

double CRACS::LISFluxes::spectrum(double Ekin_per_n, int type){
    return brokenPowerLaw(&Ekin_per_n,&fLISParameters[type][0])*pow(Ekin_per_n, -2.7);;
}

void CRACS::LISFluxes::fit(bool save){
    
    fit_simultaneously_pHeCNO       ( save );
    fit_simultaneously_pbar         ( save );
    
    fit_simultaneously_LiBeB        ( save );
    
    
    if(save){
        
        CRACS::ReadFluxes fluxes;
        fluxes.GetSpectrum("proton",   "ams", CRACS::EKINPERN).WriteToFile("proton.txt");
        fluxes.GetSpectrum("helium",   "ams", CRACS::EKINPERN).WriteToFile("helium.txt");
        fluxes.GetSpectrum("carbon",   "ams", CRACS::EKINPERN).WriteToFile("carbon.txt");
        fluxes.GetSpectrum("oxygen",   "ams", CRACS::EKINPERN).WriteToFile("oxygen.txt");
        fluxes.GetSpectrum("nitrogen", "ams", CRACS::EKINPERN).WriteToFile("nitrogen.txt");
        
        fluxes.GetSpectrum("lithium",    "ams", CRACS::EKINPERN).WriteToFile("lithium.txt");
        fluxes.GetSpectrum("beryllium",  "ams", CRACS::EKINPERN).WriteToFile("beryllium.txt");
        fluxes.GetSpectrum("boron",      "ams", CRACS::EKINPERN).WriteToFile("boron.txt");
        fluxes.GetSpectrum("pbar",       "ams", CRACS::EKINPERN).WriteToFile("antiproton.txt");
        
        std::stringstream ss;
        ss << "Parametrizations, species in row 1:proton  2:helium  3:carbon  4:oxygen  5:antiproton" << std::endl;
        for (int i = 0; i<10; i++) {
            for (int j=0; j<7; j++) {
                ss << fLISParameters[i][j] << "       ";
            }
            ss << std::endl;
        }
        FileTool::WriteStringToFile(ss.str(), "power_law_param.txt");
    }
    
}

//void CRACS::LISFluxes::fit(bool save)
//{
//    
//    CRACS::ReadFluxes fluxes;
//    
//    CRACS::Graph proton        = fluxes.GetSpectrum("proton",      "ams", CRACS::EKINPERN);
//    CRACS::Graph helium        = fluxes.GetSpectrum("helium",      "ams", CRACS::EKINPERN);
//    CRACS::Graph carbon        = fluxes.GetSpectrum("carbon",      "ams", CRACS::EKINPERN);
//    CRACS::Graph oxygen        = fluxes.GetSpectrum("oxygen",      "ams", CRACS::EKINPERN);
//    CRACS::Graph antiproton    = fluxes.GetSpectrum("pbar",        "ams", CRACS::EKINPERN);
//    
//    CRACS::Graph proton_LIS    = proton.       Copy(); proton_LIS.SetName("proton_LIS");
//    CRACS::Graph helium_LIS    = helium.       Copy(); proton_LIS.SetName("helium_LIS");
//    CRACS::Graph carbon_LIS    = carbon.       Copy(); proton_LIS.SetName("carbon_LIS");
//    CRACS::Graph oxygen_LIS    = oxygen.       Copy(); proton_LIS.SetName("oxygen_LIS");
//    CRACS::Graph antiproton_LIS= antiproton.   Copy(); proton_LIS.SetName("antiproton_LIS");
//    
//    TF1 *brokenPowerLawFun    = new TF1("brokenPowerLawFun",CRACS::LISFluxes::brokenPowerLaw, 0.1,1e4,7);
//    brokenPowerLawFun->SetParNames("N","#gamma_{1}","#gamma_{2}","#gamma_{3}","R_{br,1}","R_{br,2}","s");
//    brokenPowerLawFun->SetParameters(8.5e3, -1.5, 0.13, 0, 3.7, 270, 0.6);
//    brokenPowerLawFun->SetParLimits(5, 250, 400);
//    brokenPowerLawFun->SetLineWidth(2);
//    
//    
//    
//    TGraphErrors* pGraph    =  proton.GetAsTGraphErrors();
//    pGraph->SetMarkerStyle  (21);
//    pGraph->SetMarkerSize   (0.5);
//    pGraph->SetMarkerColor  (4);
//    pGraph->SetLineColor    (4);
//    
//    applySolarModulation(proton_LIS,   -fSM_lower,  1,  1, 2.7);
//    proton_LIS.GetAsTGraphErrors()->Fit(brokenPowerLawFun);
//    for (int i=0; i<7; i++){
//        fLISParameters_lower[fProtonLIS][i] = brokenPowerLawFun->GetParameter(i);
//    }
//    applySolarModulation(proton_LIS,    fSM_lower,  1,  1, 2.7);
//    applySolarModulation(proton_LIS,   -fSM_upper,  1,  1, 2.7);
//    proton_LIS.GetAsTGraphErrors()->Fit(brokenPowerLawFun);
//    for (int i=0; i<7; i++){
//        fLISParameters_upper[fProtonLIS][i] = brokenPowerLawFun->GetParameter(i);
//    }
//    applySolarModulation(proton_LIS,    fSM_upper,  1,  1, 2.7);
//    applySolarModulation(proton_LIS,    -fSM_mean,  1,  1, 2.7);
//    
//    TGraphErrors* pLIS      = proton_LIS.GetAsTGraphErrors();
//    pLIS->SetMarkerStyle    (20);
//    pLIS->SetMarkerSize     (0.5);
//    pLIS->SetMarkerColor    (4);
//    pLIS->SetLineColor      (4);
//    
//    brokenPowerLawFun->SetLineColor(4);
//    pLIS->Fit(brokenPowerLawFun);
//    for (int i=0; i<7; i++){
//        fLISParameters[fProtonLIS][i] = brokenPowerLawFun->GetParameter(i);
//    }
//    
//    
//    
//    TGraphErrors* heGraph   = helium.GetAsTGraphErrors();
//    heGraph->SetMarkerStyle (21);
//    heGraph->SetMarkerSize  (0.5);
//    heGraph->SetMarkerColor (2);
//    heGraph->SetLineColor   (2);
//    
//    brokenPowerLawFun->SetParameter(0, 5e2);
//    brokenPowerLawFun->SetLineColor(2);
//    
//    applySolarModulation(helium_LIS,   -fSM_lower,  4,  2, 2.7);
//    helium_LIS.GetAsTGraphErrors()->Fit(brokenPowerLawFun);
//    for (int i=0; i<7; i++){
//        fLISParameters_lower[fHeliumLIS][i] = brokenPowerLawFun->GetParameter(i);
//    }
//    applySolarModulation(helium_LIS,    fSM_lower,  4,  2, 2.7);
//    applySolarModulation(helium_LIS,   -fSM_upper,  4,  2, 2.7);
//    helium_LIS.GetAsTGraphErrors()->Fit(brokenPowerLawFun);
//    for (int i=0; i<7; i++){
//        fLISParameters_upper[fHeliumLIS][i] = brokenPowerLawFun->GetParameter(i);
//    }
//    applySolarModulation(helium_LIS,    fSM_upper,  4,  2, 2.7);
//    applySolarModulation(helium_LIS,    -fSM_mean,  4,  2, 2.7);
//    
//    
//    TGraphErrors* heLIS      = helium_LIS.GetAsTGraphErrors();
//    heLIS->SetMarkerStyle    (20);
//    heLIS->SetMarkerSize     (0.5);
//    heLIS->SetMarkerColor    (2);
//    heLIS->SetLineColor      (2);
//    
//    heLIS->Fit(brokenPowerLawFun);
//    
//    for (int i=0; i<7; i++)
//        fLISParameters[fHeliumLIS][i] = brokenPowerLawFun->GetParameter(i);
//        
//    
//    
//    TGraphErrors* CGraph    = carbon.GetAsTGraphErrors();
//    CGraph->SetMarkerStyle  (21);
//    CGraph->SetMarkerSize   (0.5);
//    CGraph->SetMarkerColor  (1);
//    CGraph->SetLineColor    (1);
//    
//    brokenPowerLawFun->SetParameter(0, 2e1);
//    brokenPowerLawFun->SetLineColor(1);
//    
//    applySolarModulation(carbon_LIS,   -fSM_lower,  12,  6, 2.7);
//    carbon_LIS.GetAsTGraphErrors()->Fit(brokenPowerLawFun);
//    for (int i=0; i<7; i++){
//        fLISParameters_lower[fCarbonLIS][i] = brokenPowerLawFun->GetParameter(i);
//    }
//    applySolarModulation(carbon_LIS,    fSM_lower,  12,  6, 2.7);
//    applySolarModulation(carbon_LIS,   -fSM_upper,  12,  6, 2.7);
//    carbon_LIS.GetAsTGraphErrors()->Fit(brokenPowerLawFun);
//    for (int i=0; i<7; i++){
//        fLISParameters_upper[fCarbonLIS][i] = brokenPowerLawFun->GetParameter(i);
//    }
//    applySolarModulation(carbon_LIS,    fSM_upper,  12,  6, 2.7);
//    applySolarModulation(carbon_LIS,    -fSM_mean,  12,  6, 2.7);
//    
//    TGraphErrors* CLIS      = carbon_LIS.GetAsTGraphErrors();
//    CLIS->SetMarkerStyle    (20);
//    CLIS->SetMarkerSize     (0.5);
//    CLIS->SetMarkerColor    (1);
//    CLIS->SetLineColor      (1);
//    
//    CLIS->Fit(brokenPowerLawFun);
//    
//    for (int i=0; i<7; i++)
//        fLISParameters[fCarbonLIS][i] = brokenPowerLawFun->GetParameter(i);
//    
//    
//    TGraphErrors* OGraph    = oxygen.GetAsTGraphErrors();
//    OGraph->SetMarkerStyle  (21);
//    OGraph->SetMarkerSize   (0.5);
//    OGraph->SetMarkerColor  (8);
//    OGraph->SetLineColor    (8);
//    
//    brokenPowerLawFun->SetParameter(0, 2e1);
//    brokenPowerLawFun->SetLineColor(8);
//    
//    applySolarModulation(oxygen_LIS,   -fSM_lower,  16,  8, 2.7);
//    oxygen_LIS.GetAsTGraphErrors()->Fit(brokenPowerLawFun);
//    for (int i=0; i<7; i++){
//        fLISParameters_lower[fOxygenLIS][i] = brokenPowerLawFun->GetParameter(i);
//    }
//    applySolarModulation(oxygen_LIS,    fSM_lower,  16,  8, 2.7);
//    applySolarModulation(oxygen_LIS,   -fSM_upper,  16,  8, 2.7);
//    oxygen_LIS.GetAsTGraphErrors()->Fit(brokenPowerLawFun);
//    for (int i=0; i<7; i++){
//        fLISParameters_upper[fOxygenLIS][i] = brokenPowerLawFun->GetParameter(i);
//    }
//    applySolarModulation(oxygen_LIS,    fSM_upper,  16,  8, 2.7);
//    applySolarModulation(oxygen_LIS,    -fSM_mean,  16,  8, 2.7);
//    
//    TGraphErrors* OLIS      = oxygen_LIS.GetAsTGraphErrors();
//    OLIS->SetMarkerStyle    (20);
//    OLIS->SetMarkerSize     (0.5);
//    OLIS->SetMarkerColor    (8);
//    OLIS->SetLineColor      (8);
//    
//    OLIS->Fit(brokenPowerLawFun);
//    
//    for (int i=0; i<7; i++)
//        fLISParameters[fOxygenLIS][i] = brokenPowerLawFun->GetParameter(i);
//    
//    
//    TGraphErrors* pbarGraph = antiproton.GetAsTGraphErrors();
//    pbarGraph->SetMarkerStyle (21);
//    pbarGraph->SetMarkerSize  (0.5);
//    pbarGraph->SetMarkerColor (6);
//    pbarGraph->SetLineColor   (6);
//    
//    brokenPowerLawFun->SetParameters(1.5e0, -2.5, 0.4, 0, 7, 270, 0.8);
//    brokenPowerLawFun->FixParameter(5, 3.52540e+02);
//    brokenPowerLawFun->FixParameter(3, 1.87906e-01);
//    
//    
//    applySolarModulation(antiproton_LIS,   -fSM_lower,  1,  -1, 2.7);
//    antiproton_LIS.GetAsTGraphErrors()->Fit(brokenPowerLawFun);
//    brokenPowerLawFun->FixParameter(3, brokenPowerLawFun->GetParameter(2));
//    for (int i=0; i<7; i++){
//        fLISParameters_lower[fAntiprotonLIS][i] = brokenPowerLawFun->GetParameter(i);
//    }
//    applySolarModulation(antiproton_LIS,    fSM_lower,  1,  -1, 2.7);
//    applySolarModulation(antiproton_LIS,   -fSM_upper,  1,  -1, 2.7);
//    antiproton_LIS.GetAsTGraphErrors()->Fit(brokenPowerLawFun);
//    brokenPowerLawFun->FixParameter(3, brokenPowerLawFun->GetParameter(2));
//    for (int i=0; i<7; i++){
//        fLISParameters_upper[fAntiprotonLIS][i] = brokenPowerLawFun->GetParameter(i);
//    }
//    applySolarModulation(antiproton_LIS,    fSM_upper,  1,  -1, 2.7);
//    applySolarModulation(antiproton_LIS,    -fSM_mean,  1,  -1, 2.7);
//    
//    
//    TGraphErrors* pbarLIS   = antiproton_LIS.GetAsTGraphErrors();
//    pbarLIS->SetMarkerStyle   (20);
//    pbarLIS->SetMarkerSize    (0.5);
//    pbarLIS->SetMarkerColor   (6);
//    pbarLIS->SetLineColor     (6);
//    
//    brokenPowerLawFun->SetLineColor(6);
//    pbarLIS->Fit(brokenPowerLawFun);
//    brokenPowerLawFun->FixParameter(3, brokenPowerLawFun->GetParameter(2));
//    
//    for (int i=0; i<7; i++)
//        fLISParameters[fAntiprotonLIS][i] = brokenPowerLawFun->GetParameter(i);
//    
//    
//    
//    if (!save) {
//        return;
//    }
//    
//    TCanvas can("can", "can", 800, 600);
//    
//    can.cd();
//    can.SetLogx();
//    can.SetLogy();
//    
//    TH1D background("background", "Fluxes AMS-02;E_{kin}/n [GeV/n]; #Phi (E_{kin}/n)^{2.7} [(GeV/n)^{2.7} m^{-2} sr^{-1} s^{-1} (GeV/n)^{-1}]", 100, 0.5, 5e3);
//    background.GetYaxis()->SetRangeUser(1e-1, 5e4);
//    background.SetStats(false);
//    background.DrawCopy();
//    
//    pLIS        ->Draw ("same p");
//    heGraph     ->Draw ("same p");
//    heLIS       ->Draw ("same p");
//    pGraph      ->Draw ("same p");
//    CGraph      ->Draw ("same p");
//    CLIS        ->Draw ("same p");
//    OGraph      ->Draw ("same p");
//    OLIS        ->Draw ("same p");
//    pbarGraph   ->Draw ("same p");
//    pbarLIS     ->Draw ("same p");
//    
//    
//    can.Print("LisFluxes.pdf", "pdf");
//    
//    std::stringstream ss;
//    ss << "Parametrizations, species in row 1:proton  2:helium  3:carbon  4:oxygen  5:antiproton" << std::endl;
//    for (int i = 0; i<5; i++) {
//        for (int j=0; j<7; j++) {
//            ss << fLISParameters[i][j] << "       ";
//        }
//        ss << std::endl;
//    }
//    FileTool::WriteStringToFile(ss.str(), "power_law_param.txt");
//    
//    std::stringstream ss_l;
//    ss_l << "Parametrizations, species in row 1:proton  2:helium  3:carbon  4:oxygen  5:antiproton" << std::endl;
//    for (int i = 0; i<5; i++) {
//        for (int j=0; j<7; j++) {
//            ss_l << fLISParameters_lower[i][j] << "       ";
//        }
//        ss_l << std::endl;
//    }
//    FileTool::WriteStringToFile(ss_l.str(), "power_law_param_lower.txt");
//    
//    std::stringstream ss_u;
//    ss_u << "Parametrizations, species in row 1:proton  2:helium  3:carbon  4:oxygen  5:antiproton" << std::endl;
//    for (int i = 0; i<5; i++) {
//        for (int j=0; j<7; j++) {
//            ss_u << fLISParameters_upper[i][j] << "       ";
//        }
//        ss_u << std::endl;
//    }
//    FileTool::WriteStringToFile(ss_u.str(), "power_law_param_upper.txt");
//    
//    
//    
//    proton.     WriteToFile( "proton.txt"       );
//    helium.     WriteToFile( "helium.txt"       );
//    carbon.     WriteToFile( "carbon.txt"       );
//    oxygen.     WriteToFile( "oxygen.txt"       );
//    antiproton. WriteToFile( "antiproton.txt"   );
//    
//    fluxes.GetSpectrum("proton", "pamela", CRACS::EKINPERN).WriteToFile("proton_pamela.txt");
//    fluxes.GetSpectrum("helium", "pamela", CRACS::EKINPERN).WriteToFile("helium_pamela.txt");
//    fluxes.GetSpectrum("pbar",   "pamela", CRACS::EKINPERN).WriteToFile("antiproton_pamela.txt");
//    fluxes.GetSpectrum("proton", "cream" , CRACS::EKINPERN).WriteToFile("proton_cream.txt");
//    fluxes.GetSpectrum("helium", "cream" , CRACS::EKINPERN).WriteToFile("helium_cream.txt");
//    
//    proton_LIS.     WriteToFile( "proton_LIS.txt"       );
//    helium_LIS.     WriteToFile( "helium_LIS.txt"       );
//    carbon_LIS.     WriteToFile( "carbon_LIS.txt"       );
//    oxygen_LIS.     WriteToFile( "oxygen_LIS.txt"       );
//    antiproton_LIS. WriteToFile( "antiproton_LIS.txt"   );
//    
//    
//}



void CRACS::LISFluxes::fcn_pHeCNO(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    
    double chisq = 0;
    
    double break2 = par[0];
    double delta2 = par[1];
    
    double norm, gamma, break1, delta1, s1;
    
    // Proton
    
    norm   = par[ 2];
    gamma  = par[ 3];
    break1 = par[ 4];
    delta1 = par[ 5];
    s1     = par[ 6];
    
    for (int i = 0; i<fProton_LIS.Size(); i++) {
        double x,y,xe,ye;
        fProton_LIS.GetPoint(i, x, y, xe, ye);
        chisq += pow( y- CRACS::LISFluxes::brokenPowerLaw_delta(x, norm, gamma, delta1, delta2, break1, break2*2, s1),2)/pow(ye, 2);
    }
    
    // Helium
    
    norm   = par[ 7];
    gamma  = par[ 8];
    break1 = par[ 9];
    delta1 = par[10];
    s1     = par[11];
    
    for (int i = 0; i<fHelium_LIS.Size(); i++) {
        double x,y,xe,ye;
        fHelium_LIS.GetPoint(i, x, y, xe, ye);
        chisq += pow( y- CRACS::LISFluxes::brokenPowerLaw_delta(x, norm, gamma, delta1, delta2, break1, break2, s1),2)/pow(ye, 2);
    }
    
    // Carbon
    
    norm   = par[12];
    gamma  = par[13];
    break1 = par[14];
    delta1 = par[15];
    s1     = par[16];
    
    for (int i = 0; i<fCarbon_LIS.Size(); i++) {
        double x,y,xe,ye;
        fCarbon_LIS.GetPoint(i, x, y, xe, ye);
        chisq += pow( y- CRACS::LISFluxes::brokenPowerLaw_delta(x, norm, gamma, delta1, delta2, break1, break2, s1),2)/pow(ye, 2);
    }
    
    // Oxygen
    
    norm   = par[17];
    gamma  = par[18];
    break1 = par[19];
    delta1 = par[20];
    s1     = par[21];
    
    for (int i = 0; i<fOxygen_LIS.Size(); i++) {
        double x,y,xe,ye;
        fOxygen_LIS.GetPoint(i, x, y, xe, ye);
        chisq += pow( y- CRACS::LISFluxes::brokenPowerLaw_delta(x, norm, gamma, delta1, delta2, break1, break2, s1),2)/pow(ye, 2);
    }
    
    // Nitrogen
    
    norm   = par[22];
    gamma  = par[23];
    break1 = par[24];
    delta1 = par[25];
    s1     = par[26];
    
    for (int i = 0; i<fNitrogen_LIS.Size(); i++) {
        double x,y,xe,ye;
        fNitrogen_LIS.GetPoint(i, x, y, xe, ye);
        chisq += pow( y- CRACS::LISFluxes::brokenPowerLaw_delta(x, norm, gamma, delta1, delta2, break1, break2, s1),2)/pow(ye, 2);
    }
    
    f = chisq;
}


void CRACS::LISFluxes::Ifit_pHeCNO(int from, int to, int printlevel)
{
    
    
    TMinuit *gMinuit = new TMinuit(CRACS::LISFluxes::fNParam);
    gMinuit->SetPrintLevel(printlevel);
    gMinuit->SetFCN( CRACS::LISFluxes::fcn_pHeCNO );
    
    Double_t arglist[100];
    Int_t ierflg = 0;
    
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
    
    // Set starting values and step sizes for parameters
    std::string prefix = "C";
    for (int i = 0; i<fNParam; i++) {
        std::stringstream ss;
        ss << i+1;
        gMinuit->mnparm(i, (prefix+ss.str()).c_str(), fStart[i], fStep[i], 0,0,ierflg);
    }
    
    for (int i = 0; i<fNParam; i++) {
        if (i<from || i>=to) {
            gMinuit->FixParameter(i);
        }
    }
    
    arglist[0] = 500;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    
    // Print results
    //Double_t amin,edm,errdef;
    //Int_t nvpar,nparx,icstat;
    //gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    
    
    for (int i=0; i<fNParam; i++) {
        double p, pErr;
        gMinuit->GetParameter(i, p, pErr);
        CRACS::LISFluxes::fParameter[i] = p;
        if (i<=from || i<to) {
            CRACS::LISFluxes::fStart    [i] = p;
        }
    }
    
    
    
    return;
    
}






void CRACS::LISFluxes::fcn_LiBeB(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    
    double chisq = 0;
    
    double break2 = par[0];
    double delta2 = par[1];
    
    double norm, gamma, break1, delta1, s1;
    
    // Lithium
    
    norm   = par[ 2];
    gamma  = par[ 3];
    break1 = par[ 4];
    delta1 = par[ 5];
    s1     = par[ 6];
    
    for (int i = 0; i<fLithium_LIS.Size(); i++) {
        double x,y,xe,ye;
        fLithium_LIS.GetPoint(i, x, y, xe, ye);
        chisq += pow( y- CRACS::LISFluxes::brokenPowerLaw_delta(x, norm, gamma, delta1, delta2, break1, break2*2, s1),2)/pow(ye, 2);
    }
    
    // Beryllium

    norm   = par[ 7];
    gamma  = par[ 8];
    break1 = par[ 9];
    delta1 = par[10];
    s1     = par[11];

    for (int i = 0; i<fBeryllium_LIS.Size(); i++) {
        double x,y,xe,ye;
        fBeryllium_LIS.GetPoint(i, x, y, xe, ye);
        chisq += pow( y- CRACS::LISFluxes::brokenPowerLaw_delta(x, norm, gamma, delta1, delta2, break1, break2, s1),2)/pow(ye, 2);
    }

    // Boron

    norm   = par[12];
    gamma  = par[13];
    break1 = par[14];
    delta1 = par[15];
    s1     = par[16];

    for (int i = 0; i<fBoron_LIS.Size(); i++) {
        double x,y,xe,ye;
        fBoron_LIS.GetPoint(i, x, y, xe, ye);
        chisq += pow( y- CRACS::LISFluxes::brokenPowerLaw_delta(x, norm, gamma, delta1, delta2, break1, break2, s1),2)/pow(ye, 2);
    }
    
    f = chisq;
}


void CRACS::LISFluxes::Ifit_LiBeB(int from, int to, int printlevel)
{
    
    
    TMinuit *gMinuit = new TMinuit(CRACS::LISFluxes::fNParam);
    gMinuit->SetPrintLevel(printlevel);
    gMinuit->SetFCN( CRACS::LISFluxes::fcn_LiBeB );
    
    Double_t arglist[100];
    Int_t ierflg = 0;
    
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
    
    // Set starting values and step sizes for parameters
    std::string prefix = "C";
    for (int i = 0; i<fNParam; i++) {
        std::stringstream ss;
        ss << i+1;
        gMinuit->mnparm(i, (prefix+ss.str()).c_str(), fStart[i], fStep[i], 0,0,ierflg);
    }
    
    for (int i = 0; i<fNParam; i++) {
        if (i<from || i>=to) {
            gMinuit->FixParameter(i);
        }
    }
    
    arglist[0] = 500;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    
    // Print results
    //Double_t amin,edm,errdef;
    //Int_t nvpar,nparx,icstat;
    //gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    
    
    for (int i=0; i<fNParam; i++) {
        double p, pErr;
        gMinuit->GetParameter(i, p, pErr);
        CRACS::LISFluxes::fParameter[i] = p;
        if (i<=from || i<to) {
            CRACS::LISFluxes::fStart    [i] = p;
        }
    }
    
    
    
    return;
    
}



void CRACS::LISFluxes::fcn_pbar(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    
    double chisq = 0;
    
    double break2 = par[0];
    double delta2 = par[1];
    
    double norm, gamma, break1, delta1, s1;
    
    // Antiproton
    
    norm   = par[ 2];
    gamma  = par[ 3];
    break1 = par[ 4];
    delta1 = par[ 5];
    s1     = par[ 6];
    
    for (int i = 0; i<fAntiproton_LIS.Size(); i++) {
        double x,y,xe,ye;
        fAntiproton_LIS.GetPoint(i, x, y, xe, ye);
        chisq += pow( y- CRACS::LISFluxes::brokenPowerLaw_delta(x, norm, gamma, delta1, delta2, break1, break2*2, s1),2)/pow(ye, 2);
    }
    
    
    f = chisq;
}


void CRACS::LISFluxes::Ifit_pbar(int from, int to, int printlevel)
{
    
    TMinuit *gMinuit = new TMinuit(CRACS::LISFluxes::fNParam);
    gMinuit->SetPrintLevel(printlevel);
    gMinuit->SetFCN( CRACS::LISFluxes::fcn_pbar );
    
    Double_t arglist[100];
    Int_t ierflg = 0;
    
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
    
    // Set starting values and step sizes for parameters
    std::string prefix = "C";
    for (int i = 0; i<fNParam; i++) {
        std::stringstream ss;
        ss << i+1;
        gMinuit->mnparm(i, (prefix+ss.str()).c_str(), fStart[i], fStep[i], 0,0,ierflg);
    }
    
    for (int i = 0; i<fNParam; i++) {
        if (i<from || i>=to) {
            gMinuit->FixParameter(i);
        }
    }
    
    arglist[0] = 500;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    
    // Print results
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    
    
    for (int i=0; i<fNParam; i++) {
        double p, pErr;
        gMinuit->GetParameter(i, p, pErr);
        CRACS::LISFluxes::fParameter[i] = p;
        if (i<=from || i<to) {
            CRACS::LISFluxes::fStart    [i] = p;
        }
    }
    
    
    return;
    
}



void CRACS::LISFluxes::fit_simultaneously_pHeCNO(bool save)
{
    
    CRACS::LISFluxes::fNParam = 2 + 5*5;//2+5*4;
    double start[100] = {150, -0.1,       1.1e4,  0.15,   3,  1.2, 0.5,       6.0e2,  0.05,   5,  1.05, 0.6,       2.1e1,  0.05,   5,  1.1, 0.6,      2.1e1, 0.05,  5,  1.1, 0.6,      2.1e0, 0.20,  5,  1.1, 0.6       };
    double step [100] = {  1,  0.01,      1e3,    0.01,   1,  0.1, 0.1,       1e2,    0.01,   1,  0.05, 0.1,       1e1,    0.01,   1,  0.1, 0.1,      1e1,   0.01,  1,  0.1, 0.1,      1e1,   0.01,  1,  0.1, 0.1      };
    
    for (int i=0; i<fNParam; i++) {
        fStart[i] = start[i];
        fStep [i] = step [i];
    }
    
    
    double norm, break1, break2, s1;
    double  gamma, delta2, delta1;
    
    CRACS::ReadFluxes fluxes;
    
    CRACS::Graph proton        = fluxes.GetSpectrum("proton",      "ams", CRACS::EKINPERN);
    CRACS::Graph helium        = fluxes.GetSpectrum("helium",      "ams", CRACS::EKINPERN);
    CRACS::Graph carbon        = fluxes.GetSpectrum("carbon",      "ams", CRACS::EKINPERN);
    CRACS::Graph oxygen        = fluxes.GetSpectrum("oxygen",      "ams", CRACS::EKINPERN);
    CRACS::Graph nitrogen      = fluxes.GetSpectrum("nitrogen",    "ams", CRACS::EKINPERN);
    
    
    fProton_LIS    = proton.       Copy(); fProton_LIS.     SetName("proton_LIS");
    fHelium_LIS    = helium.       Copy(); fHelium_LIS.     SetName("helium_LIS");
    fCarbon_LIS    = carbon.       Copy(); fCarbon_LIS.     SetName("carbon_LIS");
    fOxygen_LIS    = oxygen.       Copy(); fOxygen_LIS.     SetName("oxygen_LIS");
    fNitrogen_LIS  = nitrogen.     Copy(); fNitrogen_LIS.   SetName("nitrogen_LIS");
    
    applySolarModulation(fProton_LIS,    -fSM,  1,  1, 2.7);
    applySolarModulation(fHelium_LIS,    -fSM,  4,  2, 2.7);
    applySolarModulation(fCarbon_LIS,    -fSM, 12,  6, 2.7);
    applySolarModulation(fOxygen_LIS,    -fSM, 16,  8, 2.7);
    applySolarModulation(fNitrogen_LIS,  -fSM, 14,  7, 2.7);
    
    
    Ifit_pHeCNO( 2, 7, -1);
    Ifit_pHeCNO( 7,12, -1);
    Ifit_pHeCNO(12,17, -1);
    Ifit_pHeCNO(17,22, -1);
    Ifit_pHeCNO(22,27, -1);
    
    int print=-1;
    if(save) print = 0;
    Ifit_pHeCNO(0,100,  print);
    
    if(save){
        proton.     WriteToFile( "proton.txt"       );
        helium.     WriteToFile( "helium.txt"       );
        carbon.     WriteToFile( "carbon.txt"       );
        oxygen.     WriteToFile( "oxygen.txt"       );
        nitrogen.   WriteToFile( "nitrogen.txt"     );
        
        fluxes.GetSpectrum("proton", "pamela", CRACS::EKINPERN).WriteToFile("proton_pamela.txt");
        fluxes.GetSpectrum("helium", "pamela", CRACS::EKINPERN).WriteToFile("helium_pamela.txt");
        fluxes.GetSpectrum("proton", "cream" , CRACS::EKINPERN).WriteToFile("proton_cream.txt");
        fluxes.GetSpectrum("helium", "cream" , CRACS::EKINPERN).WriteToFile("helium_cream.txt");
        
        fProton_LIS.     WriteToFile( "proton_LIS.txt"       );
        fHelium_LIS.     WriteToFile( "helium_LIS.txt"       );
        fCarbon_LIS.     WriteToFile( "carbon_LIS.txt"       );
        fOxygen_LIS.     WriteToFile( "oxygen_LIS.txt"       );
        
    }
    
    
    break2 = fParameter[ 0]*2;
    delta2 = fParameter[ 1];
    
    norm   = fParameter[ 2];
    gamma  = fParameter[ 3];
    break1 = fParameter[ 4];
    delta1 = fParameter[ 5];
    s1     = fParameter[ 6];
    
    int i = fProtonLIS;
    fLISParameters[i][0]  = norm;
    fLISParameters[i][1]  = gamma-delta1;
    fLISParameters[i][2]  = gamma;
    fLISParameters[i][3]  = gamma+delta2;
    fLISParameters[i][4]  = break1;
    fLISParameters[i][5]  = break2;
    fLISParameters[i][6]  = s1;
    
    
    if(save){
        std::string s = "#   T/n                 flux (T/n)^2.7";
        for (double d=0; d<4; d+=0.1) {
            double Tn       = pow(10, d);
            double flux     = brokenPowerLaw_delta(Tn, norm, gamma, delta1, delta2, break1, break2, s1);
            
            s += "\n" +  QString("%1").arg(Tn  ).leftJustified(20).toStdString();
            s +=         QString("%1").arg(flux).leftJustified(20).toStdString();
            
        }
        FileTool::WriteStringToFile(s, "proton_fit.txt");
    }
    
    break2 = fParameter[0];
    delta2 = fParameter[1];
    
    norm   = fParameter[ 7];
    gamma  = fParameter[ 8];
    break1 = fParameter[ 9];
    delta1 = fParameter[10];
    s1     = fParameter[11];
    
    i = fHeliumLIS;
    fLISParameters[i][0]  = norm;
    fLISParameters[i][1]  = gamma-delta1;
    fLISParameters[i][2]  = gamma;
    fLISParameters[i][3]  = gamma+delta2;
    fLISParameters[i][4]  = break1;
    fLISParameters[i][5]  = break2;
    fLISParameters[i][6]  = s1;
    
    if(save){
        std::string s = "#   T/n                 flux (T/n)^2.7";
        for (double d=0; d<4; d+=0.1) {
            double Tn       = pow(10, d);
            double flux     = brokenPowerLaw_delta(Tn, norm, gamma, delta1, delta2, break1, break2, s1);
            
            s += "\n" +  QString("%1").arg(Tn  ).leftJustified(20).toStdString();
            s +=         QString("%1").arg(flux).leftJustified(20).toStdString();
            
        }
        FileTool::WriteStringToFile(s, "helium_fit.txt");
    }
    
    
    break2 = fParameter[0];
    delta2 = fParameter[1];
    
    norm   = fParameter[12];
    gamma  = fParameter[13];
    break1 = fParameter[14];
    delta1 = fParameter[15];
    s1     = fParameter[16];
    
    i = fCarbonLIS;
    fLISParameters[i][0]  = norm;
    fLISParameters[i][1]  = gamma-delta1;
    fLISParameters[i][2]  = gamma;
    fLISParameters[i][3]  = gamma+delta2;
    fLISParameters[i][4]  = break1;
    fLISParameters[i][5]  = break2;
    fLISParameters[i][6]  = s1;
    
    
    if(save){
        std::string s = "#   T/n                 flux (T/n)^2.7";
        for (double d=0; d<4; d+=0.1) {
            double Tn       = pow(10, d);
            double flux     = brokenPowerLaw_delta(Tn, norm, gamma, delta1, delta2, break1, break2, s1);
            
            s += "\n" +  QString("%1").arg(Tn  ).leftJustified(20).toStdString();
            s +=         QString("%1").arg(flux).leftJustified(20).toStdString();
            
        }
        FileTool::WriteStringToFile(s, "carbon_fit.txt");
    }
    
    
    break2 = fParameter[0];
    delta2 = fParameter[1];
    
    norm   = fParameter[17];
    gamma  = fParameter[18];
    break1 = fParameter[19];
    delta1 = fParameter[20];
    s1     = fParameter[21];
    
    i = fOxygenLIS;
    fLISParameters[i][0]  = norm;
    fLISParameters[i][1]  = gamma-delta1;
    fLISParameters[i][2]  = gamma;
    fLISParameters[i][3]  = gamma+delta2;
    fLISParameters[i][4]  = break1;
    fLISParameters[i][5]  = break2;
    fLISParameters[i][6]  = s1;
    
    if(save){
        std::string s = "#   T/n                 flux (T/n)^2.7";
        for (double d=0; d<4; d+=0.1) {
            double Tn       = pow(10, d);
            double flux     = brokenPowerLaw_delta(Tn, norm, gamma, delta1, delta2, break1, break2, s1);
            
            s += "\n" +  QString("%1").arg(Tn  ).leftJustified(20).toStdString();
            s +=         QString("%1").arg(flux).leftJustified(20).toStdString();
            
        }
        FileTool::WriteStringToFile(s, "oxygen_fit.txt");
    }
    
    
    
    break2 = fParameter[0];
    delta2 = fParameter[1];
    
    norm   = fParameter[22];
    gamma  = fParameter[23];
    break1 = fParameter[24];
    delta1 = fParameter[25];
    s1     = fParameter[26];
    
    i = fNitrogenLIS;
    fLISParameters[i][0]  = norm;
    fLISParameters[i][1]  = gamma-delta1;
    fLISParameters[i][2]  = gamma;
    fLISParameters[i][3]  = gamma+delta2;
    fLISParameters[i][4]  = break1;
    fLISParameters[i][5]  = break2;
    fLISParameters[i][6]  = s1;
    
    if(save){
        std::string s = "#   T/n                 flux (T/n)^2.7";
        for (double d=0; d<4; d+=0.1) {
            double Tn       = pow(10, d);
            double flux     = brokenPowerLaw_delta(Tn, norm, gamma, delta1, delta2, break1, break2, s1);
            
            s += "\n" +  QString("%1").arg(Tn  ).leftJustified(20).toStdString();
            s +=         QString("%1").arg(flux).leftJustified(20).toStdString();
            
        }
        FileTool::WriteStringToFile(s, "nitrogen_fit.txt");
    }
    
    
}



void CRACS::LISFluxes::fit_simultaneously_LiBeB(bool save)
{
    
    CRACS::LISFluxes::fNParam = 2 + 5*3;
    double start[100] = {150, -0.2,       2.1e0, 0.40,  5,  1.5, 0.6,      2.1e0, 0.40,  5,  1.5, 0.6,      2.1e0, 0.40,  5,  1.5, 0.6      };
    double step [100] = {  1,  0.01,      1e1,   0.01,  1,  0.1, 0.1,      1e1,   0.01,  1,  0.1, 0.1,      1e1,   0.01,  1,  0.1, 0.1      };
    
    for (int i=0; i<fNParam; i++) {
        fStart[i] = start[i];
        fStep [i] = step [i];
    }
    
    
    double norm, break1, break2, s1;
    double  gamma, delta2, delta1;
    
    CRACS::ReadFluxes fluxes;
    
    CRACS::Graph lithium        = fluxes.GetSpectrum("lithium",     "ams", CRACS::EKINPERN);
    CRACS::Graph beryllium      = fluxes.GetSpectrum("beryllium",   "ams", CRACS::EKINPERN);
    CRACS::Graph boron          = fluxes.GetSpectrum("boron",       "ams", CRACS::EKINPERN);
    
    
    fLithium_LIS    = lithium.     Copy(); fLithium_LIS.     SetName("lithium_LIS");
    fBeryllium_LIS  = beryllium.   Copy(); fBeryllium_LIS.   SetName("beryllium_LIS");
    fBoron_LIS      = boron.       Copy(); fBoron_LIS.       SetName("boron_LIS");
    
    applySolarModulation(fLithium_LIS,   -fSM,  7,  3, 2.7);
    applySolarModulation(fBeryllium_LIS, -fSM,  9,  4, 2.7);
    applySolarModulation(fBoron_LIS,     -fSM, 11,  5, 2.7);
    
    
    Ifit_LiBeB( 2, 7, -1);
    
    Ifit_LiBeB( 7,12, -1);
    Ifit_LiBeB(12,17, -1);

    int print=-1;
    if(save) print = 0;
    Ifit_LiBeB(0,100,  print);
    
    if(save){
        lithium.    WriteToFile( "lithium.txt"      );
        beryllium.  WriteToFile( "beryllium.txt"    );
        boron.      WriteToFile( "boron.txt"        );
        
        fLithium_LIS.     WriteToFile( "lithium_LIS.txt"    );
        fBeryllium_LIS.   WriteToFile( "beryllium_LIS.txt"  );
        fBoron_LIS.       WriteToFile( "boron_LIS.txt"      );
    
    }
    
    
    break2 = fParameter[ 0]*2;
    delta2 = fParameter[ 1];
    
    norm   = fParameter[ 2];
    gamma  = fParameter[ 3];
    break1 = fParameter[ 4];
    delta1 = fParameter[ 5];
    s1     = fParameter[ 6];
    
    int i = fLithiumLIS;
    fLISParameters[i][0]  = norm;
    fLISParameters[i][1]  = gamma-delta1;
    fLISParameters[i][2]  = gamma;
    fLISParameters[i][3]  = gamma+delta2;
    fLISParameters[i][4]  = break1;
    fLISParameters[i][5]  = break2;
    fLISParameters[i][6]  = s1;
    
    
    if(save){
        std::string s = "#   T/n                 flux (T/n)^2.7";
        for (double d=0; d<4; d+=0.1) {
            double Tn       = pow(10, d);
            double flux     = brokenPowerLaw_delta(Tn, norm, gamma, delta1, delta2, break1, break2, s1);
            
            s += "\n" +  QString("%1").arg(Tn  ).leftJustified(20).toStdString();
            s +=         QString("%1").arg(flux).leftJustified(20).toStdString();
            
        }
        FileTool::WriteStringToFile(s, "lithium_fit.txt");
    }
    
    
    break2 = fParameter[0];
    delta2 = fParameter[1];
    
    norm   = fParameter[ 7];
    gamma  = fParameter[ 8];
    break1 = fParameter[ 9];
    delta1 = fParameter[10];
    s1     = fParameter[11];
    
    i = fBerylliumLIS;
    fLISParameters[i][0]  = norm;
    fLISParameters[i][1]  = gamma-delta1;
    fLISParameters[i][2]  = gamma;
    fLISParameters[i][3]  = gamma+delta2;
    fLISParameters[i][4]  = break1;
    fLISParameters[i][5]  = break2;
    fLISParameters[i][6]  = s1;
    
    if(save){
        std::string s = "#   T/n                 flux (T/n)^2.7";
        for (double d=0; d<4; d+=0.1) {
            double Tn       = pow(10, d);
            double flux     = brokenPowerLaw_delta(Tn, norm, gamma, delta1, delta2, break1, break2, s1);
            
            s += "\n" +  QString("%1").arg(Tn  ).leftJustified(20).toStdString();
            s +=         QString("%1").arg(flux).leftJustified(20).toStdString();
            
        }
        FileTool::WriteStringToFile(s, "beryllium_fit.txt");
    }
    
    
    break2 = fParameter[0];
    delta2 = fParameter[1];
    
    norm   = fParameter[12];
    gamma  = fParameter[13];
    break1 = fParameter[14];
    delta1 = fParameter[15];
    s1     = fParameter[16];
    
    i = fBoronLIS;
    fLISParameters[i][0]  = norm;
    fLISParameters[i][1]  = gamma-delta1;
    fLISParameters[i][2]  = gamma;
    fLISParameters[i][3]  = gamma+delta2;
    fLISParameters[i][4]  = break1;
    fLISParameters[i][5]  = break2;
    fLISParameters[i][6]  = s1;
    
    
    if(save){
        std::string s = "#   T/n                 flux (T/n)^2.7";
        for (double d=0; d<4; d+=0.1) {
            double Tn       = pow(10, d);
            double flux     = brokenPowerLaw_delta(Tn, norm, gamma, delta1, delta2, break1, break2, s1);
            
            s += "\n" +  QString("%1").arg(Tn  ).leftJustified(20).toStdString();
            s +=         QString("%1").arg(flux).leftJustified(20).toStdString();
            
        }
        FileTool::WriteStringToFile(s, "boron_fit.txt");
    }
    
    
}





void CRACS::LISFluxes::fit_simultaneously_pbar(bool save)
{
    
    CRACS::LISFluxes::fNParam = 2 + 5*1;
    double start[100] = {150,  0.0,       1.5e0,  0.15,   7,  0.9, 0.8,      };
    double step [100] = {  1,  0.01,      1e3,    0.01,   1,  0.1, 0.1,      };
    
    for (int i=0; i<fNParam; i++) {
        fStart[i] = start[i];
        fStep [i] = step [i];
    }
    
    double norm, break1, break2, s1;
    double  gamma, delta2, delta1;
    
    CRACS::ReadFluxes fluxes;
    
    CRACS::Graph antiproton = fluxes.       GetSpectrum("pbar",      "ams", CRACS::EKINPERN);
    
    fAntiproton_LIS         = antiproton.   Copy(); fProton_LIS.SetName("antiproton_LIS");
    
    applySolarModulation(fAntiproton_LIS,    -fSM,  1,  -1, 2.7);
    
    
    int print=-1;
    if(save) print = 0;
    Ifit_pbar( 2, 7, print);
    //Ifit_pbar(    );
    
    if(save){
        antiproton. WriteToFile( "antiproton.txt"   );
        fAntiproton_LIS.     WriteToFile( "antiproton_LIS.txt"       );
        fluxes.GetSpectrum("pbar",   "pamela", CRACS::EKINPERN).WriteToFile("antiproton_pamela.txt");
    }
    
    break2 = fParameter[ 0];
    delta2 = fParameter[ 1];
    
    norm   = fParameter[ 2];
    gamma  = fParameter[ 3];
    break1 = fParameter[ 4];
    delta1 = fParameter[ 5];
    s1     = fParameter[ 6];
    
    int i = fAntiprotonLIS;
    fLISParameters[i][0]  = norm;
    fLISParameters[i][1]  = gamma-delta1;
    fLISParameters[i][2]  = gamma;
    fLISParameters[i][3]  = gamma+delta2;
    fLISParameters[i][4]  = break1;
    fLISParameters[i][5]  = break2;
    fLISParameters[i][6]  = s1;
    
    if(save){
        std::string s = "#   T/n                 flux (T/n)^2.7";
        for (double d=0; d<4; d+=0.1) {
            double Tn       = pow(10, d);
            double flux     = brokenPowerLaw_delta(Tn, norm, gamma, delta1, delta2, break1, break2*2, s1);
            
            s += "\n" +  QString("%1").arg(Tn  ).leftJustified(20).toStdString();
            s +=         QString("%1").arg(flux).leftJustified(20).toStdString();
            
        }
        FileTool::WriteStringToFile(s, "antiproton_fit.txt");
    }
    
}




