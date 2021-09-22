#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

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
#include "TH2D.h"
#include "TCanvas.h"
#include "TApplication.h"

#include "QString"
#include "QStringList"

using namespace  CRACS;

std::string fProgramDescription = "Program to get uncertainties.";



// pbar production cross sections
double dT_pp_pbar_LAB   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::DI_MAURO12);
}

double dT_pp_pbar_LAB_KR   ( double T_p,         double T_pbar ){
    if (T_p<6*fMass_proton) {
        return 0;
    }
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::KACHELRIESS);
}

double dT_Hep_pbar_LAB  ( double Tn_He,       double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(Tn_He, T_pbar, CS::PP_PBAR, CS::DI_MAURO12)*pow(4,0.7);
}
double dT_pHe_pbar_LAB   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::DI_MAURO12)*pow(4,0.7);
}


// Fluxes
double protonFlux(double T_p){
    return LISFluxes::spectrum(T_p, LISFluxes::fProtonLIS);
}
double heliumFlux(double T_He){
    return LISFluxes::spectrum(T_He, LISFluxes::fHeliumLIS);
}
double antiprotonFlux(double T_pbar){
    return LISFluxes::spectrum(T_pbar, LISFluxes::fAntiprotonLIS);
}

// grid
const int fN_T_pbar         =  100;
const int fN_Tn_p           =   30;
const int fN_cos            =  100;

const double fMin_T_pbar    = 1e-1;
const double fMax_T_pbar    = 1e3;
const double fMin_Tn_p      = 1e0;
const double fMax_Tn_p      = 1e5;

// Lab T
std::vector<double> fT_pbar_LAB;
std::vector<double> fTn_p_LAB;

// Source term as function of T_pbar
std::vector<double> fQ;
std::vector<double> fQ_pp;
std::vector<double> fQ_pHe;
std::vector<double> fQ_Hep;

// Source term as function of T_pbar
std::vector<double> fC_Tpbar__pp;
std::vector<double> fC_Tpbar__pHe;
std::vector<double> fC_Tpbar__Hep;

std::vector<double> fU_Tpbar__Q;

std::vector<double> fU_Tpbar__pp;
std::vector<double> fU_Tpbar__pHe;
std::vector<double> fU_Tpbar__Hep;

std::vector< std::vector<double> > fC_Tpbar_Tp__pp;
std::vector< std::vector<double> > fC_Tpbar_Tp__pHe;
std::vector< std::vector<double> > fC_Tpbar_Tp__Hep;


double f_nH        = 1e+6;
double f_nHe       = f_nH * 0.1;

int main(int argc, char *argv[])
{
    
    // Config handling
    ConfigHandler* config = ConfigHandler::GetInstance();
    config->SetInput(argc, argv, fProgramDescription);
    
    bool    showApp = false;
    int     step    = 1;
    config->AddOptionInt    (   "step",         step,       "Step. Default: 1"                    );
    config->AddOptionTrue   (   "app",          showApp                                           );
    
    config->CheckOption();
    
    LISFluxes::fit();

    //    double ff = 1.3;
//    out(SourceTerms::integrationFromTo( 0.5, dT_pp_pbar_LAB,  protonFlux, f_nH,  1./ff, 1.*ff, 350 ))
//    out(SourceTerms::integrationFromTo( 0.5, dT_pp_pbar_LAB,  protonFlux, f_nH,  2./ff, 2.*ff, 350 ))
//    
//    out(dT_pp_pbar_LAB(0.5, 0.5))
//    out(dT_pp_pbar_LAB(0.9, 0.5))
//    out(dT_pp_pbar_LAB(1.0, 0.5))
//    out(dT_pp_pbar_LAB(1.2, 0.5))
//    
//    out(SourceTerms::integrationFromTo( 5, dT_pp_pbar_LAB,  protonFlux, f_nH,  30./ff, 30.*ff, 350 ))
//    return 1;

    
    
    TApplication* app;
    
    if (showApp) {
        app = new TApplication("RootApp", &argc, argv);
    }

    {
        ReadFluxes reader;
        Graph pbar_flux = reader.GetSpectrum("pbar", "ams", EKINPERN);
        std::vector<double> Pbar_T;
        std::vector<double> Pbar_RE;
        for (int i = 0; i<pbar_flux.Size(); i++) {
            double x, y, xe, ye;
            pbar_flux.GetPoint(i, x, y, xe, ye);
            Pbar_T. push_back( x    );
            Pbar_RE.push_back( ye/y );
        }
        TGraph pbar_relative_error( Pbar_T.size(), &Pbar_T.at(0), &Pbar_RE.at(0) );
        
        TCanvas can("can", "can");
        pbar_relative_error.SetTitle(" ; T_{pbar}; relative error");
        pbar_relative_error.Draw("al");
        can.SetLogx();
        can.Print("pbar_relative_error.pdf");
    }
    
    


    
    for (int i = 0; i<fN_T_pbar+1; i++) {
        fU_Tpbar__Q.push_back(0.05);
    }
    
    double dTn_p=log10(fMax_Tn_p/fMin_Tn_p)/fN_Tn_p;
    for ( double dT=log10(fMin_Tn_p); dT<log10(fMax_Tn_p)+dTn_p/2; dT+=dTn_p ) {
        double T_p   = pow(10, dT);
        fTn_p_LAB.push_back ( T_p );
    }
    
    double dT_pbar=log10(fMax_T_pbar/fMin_T_pbar)/fN_T_pbar;
    for ( double dT=log10(fMin_T_pbar); dT<log10(fMax_T_pbar)+dT_pbar/2; dT+=dT_pbar ) {
        double T_pbar   = pow(10, dT);
        fT_pbar_LAB.push_back ( T_pbar );
        double Tn_from = std::max(T_pbar, 6*fMass_proton);
        fQ_pp.      push_back (  SourceTerms::integrationFromTo( T_pbar, dT_pp_pbar_LAB,    protonFlux, f_nH,  Tn_from, 1e6, 10000 )   );
        fQ_pHe.     push_back (  SourceTerms::integrationFromTo( T_pbar, dT_pHe_pbar_LAB,   protonFlux, f_nHe, Tn_from, 1e6, 10000 )   );
        fQ_Hep.     push_back (  SourceTerms::integrationFromTo( T_pbar, dT_Hep_pbar_LAB,   heliumFlux, f_nH,  Tn_from, 1e6, 10000 )   );
        int i = fQ_pp.size()-1;
        fQ.         push_back (  fQ_pp.at(i)+fQ_Hep.at(i)+fQ_pHe.at(i) );
    }
    
    {
        TCanvas can("can", "can");
        TGraph q;    q.SetLineWidth(3);
        TGraph pp;
        TGraph pHe;  pHe.SetLineColor(2);
        TGraph Hep;  Hep.SetLineColor(4);
        for ( int i = 0 ; i<fT_pbar_LAB.size(); i++) {
            q.  SetPoint(i, fT_pbar_LAB.at(i),(fQ_pp. at(i)+fQ_pHe. at(i)+fQ_Hep. at(i))*pow(fT_pbar_LAB.at(i), 2.7));
            pp. SetPoint(i, fT_pbar_LAB.at(i), fQ_pp. at(i)*pow(fT_pbar_LAB.at(i), 2.7));
            pHe.SetPoint(i, fT_pbar_LAB.at(i), fQ_pHe.at(i)*pow(fT_pbar_LAB.at(i), 2.7));
            Hep.SetPoint(i, fT_pbar_LAB.at(i), fQ_Hep.at(i)*pow(fT_pbar_LAB.at(i), 2.7));
        }
        q.  SetTitle(" ; T_{pbar} [GeV]; q T_{pbar}^{2.7}");
        q.  Draw("al"    );
        pp. Draw("l same");
        pHe.Draw("l same");
        Hep.Draw("l same");
        can.SetLogx();
        can.SetLogy();
        can.Print("SourceTerm.pdf");
    }
    
    {
        TCanvas can("can", "can");
        TGraph DM;
        TGraph KR;  KR.SetLineColor(2);
        for ( int i = 0 ; i<fT_pbar_LAB.size(); i++) {
            DM.SetPoint(i, fT_pbar_LAB.at(i),2.3*(SourceTerms::integrationFromTo( fT_pbar_LAB.at(i), dT_pp_pbar_LAB,    protonFlux, f_nH,  1, 1e6, 10000 ))*pow(fT_pbar_LAB.at(i), 2.7));
            KR.SetPoint(i, fT_pbar_LAB.at(i),(SourceTerms::integrationFromTo( fT_pbar_LAB.at(i), dT_pp_pbar_LAB_KR, protonFlux, f_nH,  6*fMass_proton, 1e6, 10000 ))*pow(fT_pbar_LAB.at(i), 2.7));
        }
        DM.  SetTitle(" ; T_{pbar} [GeV]; q T_{pbar}^{2.7}");
        DM.  Draw("al"    );
        KR.  Draw("l same");
        can.SetLogx();
        can.SetLogy();
        can.Print("SourceTerm_diMauro_KR.pdf");
        return 0;
    }
    
    
    for (int i = 0; i<fQ.size(); i++) {
        fC_Tpbar__pp. push_back( fQ_pp. at(i)/fQ.at(i) );
        fC_Tpbar__Hep.push_back( fQ_Hep.at(i)/fQ.at(i) );
        fC_Tpbar__pHe.push_back( fQ_pHe.at(i)/fQ.at(i) );
        
        fU_Tpbar__pp. push_back( fU_Tpbar__Q.at(i)/fC_Tpbar__pp. at(i)/3. );
        fU_Tpbar__Hep.push_back( fU_Tpbar__Q.at(i)/fC_Tpbar__Hep.at(i)/3. );
        fU_Tpbar__pHe.push_back( fU_Tpbar__Q.at(i)/fC_Tpbar__pHe.at(i)/3. );
    }
    
    {
        TCanvas can("can", "can");
        TGraph pp;
        TGraph pHe;  pHe.SetLineColor(2);
        TGraph Hep;  Hep.SetLineColor(4);
        for ( int i = 0 ; i<fT_pbar_LAB.size(); i++) {
            pp. SetPoint(i, fT_pbar_LAB.at(i), fC_Tpbar__pp. at(i)  );
            pHe.SetPoint(i, fT_pbar_LAB.at(i), fC_Tpbar__pHe.at(i)  );
            Hep.SetPoint(i, fT_pbar_LAB.at(i), fC_Tpbar__Hep.at(i)  );
        }
        TH1F back("back", " ;  T_{pbar} [GeV]; contribution", 100, 0.1, 1000);
        back.SetStats(false);
        back.GetYaxis()->SetRangeUser(0, 1.);
        back.DrawCopy();
        Hep.Draw("l same");
        pHe.Draw("l same");
        pp. Draw("l same");
        can.SetLogx();
        can.Print("contribution_channel.pdf");
        back.SetTitle(";  T_{pbar} [GeV]; Allowed uncertainty");
        TGraph pp_u (  fT_pbar_LAB.size(), &fT_pbar_LAB.at(0), &fU_Tpbar__pp. at(0));
        TGraph pHe_u(  fT_pbar_LAB.size(), &fT_pbar_LAB.at(0), &fU_Tpbar__pHe.at(0));  pHe_u.SetLineColor(2);
        TGraph Hep_u(  fT_pbar_LAB.size(), &fT_pbar_LAB.at(0), &fU_Tpbar__Hep.at(0));  Hep_u.SetLineColor(4);
        back.GetYaxis()->SetRangeUser(0, .2);
        back.DrawCopy();
        Hep_u.Draw("l same");
        pHe_u.Draw("l same");
        pp_u. Draw("l same");
        can.Print("uncertainty_channel.pdf");
    }
   
    
    double factor = sqrt(fTn_p_LAB.at(1)/fTn_p_LAB.at(0));
    for (int j=0; j<fT_pbar_LAB.size(); j++) {
        double T_pbar = fT_pbar_LAB.at(j);
        std::vector<double> vec_pp;
        std::vector<double> vec_pHe;
        std::vector<double> vec_Hep;
        for (int i = 0; i<fTn_p_LAB.size(); i++) {
            double Tn_p = fTn_p_LAB.at(i);
            vec_pp. push_back(  SourceTerms::integrationFromTo( T_pbar, dT_pp_pbar_LAB,  protonFlux, f_nH,  Tn_p/factor, Tn_p*factor, 1000 )/fQ.at(j)  );
            vec_pHe.push_back(  SourceTerms::integrationFromTo( T_pbar, dT_pHe_pbar_LAB, protonFlux, f_nHe, Tn_p/factor, Tn_p*factor, 1000 )/fQ.at(j)  );
            vec_Hep.push_back(  SourceTerms::integrationFromTo( T_pbar, dT_Hep_pbar_LAB, heliumFlux, f_nH,  Tn_p/factor, Tn_p*factor, 1000 )/fQ.at(j)  );
        }
        fC_Tpbar_Tp__pp. push_back(vec_pp);
        fC_Tpbar_Tp__pHe.push_back(vec_pHe);
        fC_Tpbar_Tp__Hep.push_back(vec_Hep);
    }
    
    {
        double dpbar =  log10(fMax_T_pbar/fMin_T_pbar)/fN_T_pbar;
        double dp    =  log10(fMax_Tn_p/fMin_Tn_p)/fN_Tn_p;
        TH2D c_Tpbar_Tp__pp ("pp", "pp; log10(T_{pbar}/GeV);log10(T_{p}/n/(GeV/n))", fN_T_pbar+1, log10(fMin_T_pbar)-dpbar/2, log10(fMax_T_pbar)+dpbar/2, fN_Tn_p+1, log10(fMin_Tn_p)-dp/2, log10(fMax_Tn_p)+dp/2);
        for (int i = 0; i<fN_T_pbar+1; i++) {
            double sum = 0;
            for (int j = 0; j<fN_Tn_p+1; j++) {
                c_Tpbar_Tp__pp.SetBinContent(i+1, j+1, fC_Tpbar_Tp__pp.at(i).at(j));
                sum += fC_Tpbar_Tp__pp.at(i).at(j);
            }
            varOut(sum/fC_Tpbar__pp.at(i))
        }
        TCanvas can("can", "can");
        c_Tpbar_Tp__pp.SetStats(false);
        c_Tpbar_Tp__pp.Draw("colz");
        can.Print("contribution_LAB.pdf");
    }
    
    int i = 49;
    double T_pbar = fT_pbar_LAB.at(i);
    
    double max_contribution = 0;
    
    for (int j=0; j<fTn_p_LAB.size(); j++) {
        double T = fTn_p_LAB.at(j);
        double s = 2*T+4*fMass_proton;
        
        double m = fMass_proton;
        double cos_t_max = (   (T_pbar+m)*(T+2*m)   - m*(T-2*m)   )/sqrt(  T*(T+2*m) * T_pbar*(T_pbar+2*m)  );
        
        std::vector<double> pT_pbar;
        std::vector<double> xR_pbar;
        std::stringstream toWrite;
        toWrite << " s     " << s << "   -";
        toWrite << "\n T_pbar     " << T_pbar << "   -";
        double sum = 0;
        for (int ic = 0; ic<=fN_cos; ic++) {
            double cos_t        = cos_t_max + (1-cos_t_max)/fN_cos*ic;
            double pT_pbar_cm   = sqrt( T_pbar*(T_pbar+2*m)*(1-cos_t*cos_t)  );
            double xR           = (  (T+2*m)*(T_pbar+m) -  sqrt(  T*(T+2*m) * T_pbar*(T_pbar+2*m)  )*cos_t  )/(  m*(T-2*m) );
            double E_pbar_cm    = 1./sqrt(2*m)*(  sqrt(T+2*m)*(T_pbar+m) - sqrt(T* T_pbar*(T_pbar+2*m))*cos_t  );
            double p_pbar       = sqrt(T_pbar*(T_pbar+2*m));
            double contribution = (1-cos_t_max)/fN_cos*2*M_PI*p_pbar*ppCSParametrizations::inv_pp_nbar_CM__diMauro12(s, E_pbar_cm, pT_pbar_cm)*1e-31/(dT_pp_pbar_LAB(T, T_pbar));
            if (contribution!=contribution) {
                continue;
            }
            sum += contribution;
            contribution*=fC_Tpbar_Tp__pp.at(i).at(j);
            max_contribution = std::max(contribution, max_contribution);
            pT_pbar.push_back(pT_pbar_cm);
            xR_pbar.push_back(xR);
            toWrite << "\n" << xR << "     " << pT_pbar_cm << "     " << contribution;
        }
        out(sum)
        std::stringstream contributions;
        contributions <<  "c__s_xR_pT__s_" << std::setw(9) << std::setfill('0')  << int(s) << "__Tpbar_" << std::setw(9) << std::setfill('0') << std::right << int(1000*T_pbar) << ".txt";
        
        FileTool::WriteStringToFile( toWrite.str().c_str(), contributions.str().c_str() );
        
        
    }
    
    std::cout << max_contribution << std::endl;
    
    return 0;
    
    
    double T = 63;
    
    double s = 2*T+4*fMass_proton;
    
    double m = fMass_proton;
    double cos_t_max = (   (T_pbar+m)*(T+2*m)   - m*(T-2*m)   )/sqrt(  T*(T+2*m) * T_pbar*(T_pbar+2*m)  );
    
    std::vector<double> pT_pbar;
    std::vector<double> xR_pbar;
    
    for (int ic = 0; ic<=100; ic++) {
        double cos_t = cos_t_max + (1-cos_t_max)/100*ic;
        double pT_cm = sqrt( T_pbar*(T_pbar+2*m)*(1-cos_t*cos_t)  );
        double xR    = (  (T+2*m)*(T_pbar+m) -  sqrt(  T*(T+2*m) * T_pbar*(T_pbar+2*m)  )*cos_t  )/(  m*(T-2*m) );
        pT_pbar.push_back(pT_cm);
        xR_pbar.push_back(xR);
    }
    
    
    {
        TGraph path_xR_pT( xR_pbar.size(), &xR_pbar.at(0), &pT_pbar.at(0) );
        
        TCanvas can("can", "can");
        path_xR_pT.SetTitle("Path; xR; p_T [GeV]");
        path_xR_pT.SetMarkerSize(5);
        path_xR_pT.Draw("ap");
        can.Print("path_xR_pT.pdf");
    }
    
    out(T)
    out(T_pbar)
    out(s)
    
    
    
    if (showApp) {
        app->Run();
    }
    return 0;
    
    
}


