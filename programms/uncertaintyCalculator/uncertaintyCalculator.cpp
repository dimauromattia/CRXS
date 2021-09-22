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

std::string fOption = "standard";


double fErr_inner=0.03;
double fErr_outer=0.30;
double fSM       =0.6;

bool   fWithPlots=true;
bool   fPPfixed  =false;

double (*f__inv_pp_pbar_LAB__parametrization) (double, double, double) = ppCSParametrizations::inv_pp_pbar_CM__diMauro12;
double Interpolate_LAB(double T_p_LAB, double T_pbar_LAB, double eta_LAB, double* array);

// inv cross sections
double inv_pp_pbar_LAB   ( double T_p,    double T_pbar,  double eta   ){
    double cosTheta = tanh(eta);
    return CSTransformations::inv_pp_product_LAB(  T_p+fMass_proton, T_pbar+fMass_proton, cosTheta, fMass_proton,  f__inv_pp_pbar_LAB__parametrization)*1e-31;
}

// Fluxes
double protonFlux(double T_p){
    return LISFluxes::spectrum(T_p, LISFluxes::fProtonLIS);
}
double heliumFlux(double T_He){
    return LISFluxes::spectrum(T_He, LISFluxes::fHeliumLIS);
}

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


// ISM densities
const double f_nH        = 1e+6;
const double f_nHe       = f_nH * 0.1;



// Kinetic variables in LAB frame, grid and arraies
const int    fN_T_pbar      =   401;
const int    fN_Tn_p        =   401;
const int    fN_eta         =   401;

const double fMin_T_pbar    = 1e-1;
const double fMax_T_pbar    = 1e3;

const double fMin_Tn_p      = 5e0;
const double fMax_Tn_p      = 5e5;

const double fMin_eta       =  0;
const double fMax_eta       = 13;

double fT_pbar_LAB  [fN_T_pbar];
double fTn_p_LAB    [fN_Tn_p  ];
double fEta_LAB     [fN_eta   ];



// Kinetic variables in CM frame, grid and arraies
const int    fN_s           = 401;
const int    fN_xR          = 401;
const int    fN_pT          = 401;

const double fMin_s         = 15;
const double fMax_s         = 2e4;

const double fMin_pT        = 1e-2;
const double fMax_pT        = 1e1;

const double fMin_xR        =  1e-2;
const double fMax_xR        =  1;

double fS     [fN_s ];
double fxR    [fN_xR];
double fxF    [fN_xR];
double fPT    [fN_pT];


// Requiered uncertainty in the source term (e.g. dictated by an experiment like AMS-02)
double fRU_Tpbar__Q  [fN_T_pbar];

// Source term as function of T_pbar
double fQ       [fN_T_pbar];


// Arrays to store distributions in various parameter spaces
// There is not enough memory too have separate arraies for ever distribution

double f2D__Tpbar_eta           [ fN_T_pbar  * fN_eta               ];
double f2D__xR_pT               [ fN_xR      * fN_pT                ];

double f2D__Tpbar_Tp__a1        [ fN_Tn_p    * fN_T_pbar            ];
double f2D__Tpbar_Tp__a2        [ fN_Tn_p    * fN_T_pbar            ];

double f3D__Tpbar_Tp_eta__a1    [ fN_T_pbar  * fN_Tn_p   * fN_eta   ];
double f3D__Tpbar_Tp_eta__a2    [ fN_T_pbar  * fN_Tn_p   * fN_eta   ];

double f3D__s_xR_pT             [ fN_s       * fN_xR     * fN_pT    ];



double fHelper                  [ fN_T_pbar  * fN_Tn_p   * fN_eta   ];




// Functions to read and write the arrays to files. The arrys are saved in the following format:
//
//      1D xarray=NULL          1D                                  2D
//
//      #DESCRIPTION            #DESCRIPTION                        #               xarray[0          ]       ...     xarray[ nx-1           ]
//      array[0  ]              xarray[0  ]   array[0]              yarray[0   ]    array [0+ 0    *nx]       ...     array [(nx-1)+ 0    *nx]
//      ...                     ...           ...                   ...             ...                               ...
//      array[n-1]              xarray[n-1]   array[n-1]            yarray[ny-1]    array [0+(ny-1)*nx]       ...     array [(nx-1)+(ny-1)*nx]
//
//
//
void read_array_from_file   (  std::string filename, int n,          double* array                                                              );
void write_array_to_file    (  std::string filename, int n,          double* array, double* xarray=NULL,            std::string description=""  );
void write_2D_array_to_file (  std::string filename, int nx, int ny, double* array, double* xarray, double* yarray                              );




// Function to determine the containence from a given contibution.
// By containence c we mean:    50% of the contribution comes from the parameter space (array entries) with c<0.5
//                              90% of the contribution comes from the parameter space (array entries) with c<0.9
void get_containence(int n, const double* contribution, double* containence);


// Function to determine a distribution in the CM from a lab frame distribution stored in array
// The function calles an interpolation function in LAB frame
// As the inversion from CM to LAB is not injective we need to provide an int for different options:
//      1: pL > 0
//      2: pL < 0
//      0: return the minmum of 1 and 2;
double GetInCM( double s, double xR, double pT, double* array, int option=0 );

double GetInCM_xF( double s, double xF, double pT, double* array );





int main(int argc, char *argv[])
{
    ConfigHandler* config = ConfigHandler::GetInstance();
    
    config->SetInput(argc, argv, fProgramDescription);
    int     step    = 0;
    int     option  = 0;
    config->AddOptionInt    (   "step",           step,       "Step.   Default: 0"                      );
    config->AddOptionInt    (   "option",         option,     "Option. Default: 0 (standard)"                      );
    config->CheckOption();
    
    
    
    // Setup LAB grid
    double dTn_p=log10(fMax_Tn_p/fMin_Tn_p)/(fN_Tn_p-1);
    int i=0;
    for ( double dT=log10(fMin_Tn_p); dT<log10(fMax_Tn_p)+dTn_p/2; dT+=dTn_p ) {
        fTn_p_LAB[i] = pow(10, dT);     i++;
    }
    double dT_pbar=log10(fMax_T_pbar/fMin_T_pbar)/(fN_T_pbar-1);
    i = 0;
    for ( double dT=log10(fMin_T_pbar); dT<log10(fMax_T_pbar)+dT_pbar/2; dT+=dT_pbar ) {
        fT_pbar_LAB[i] = pow(10, dT);   i++;
    }
    double dEta=(fMax_eta-fMin_eta)/(fN_eta-1);
    i = 0;
    for ( double eta=fMin_eta; eta<fMax_eta+dEta/2; eta+=dEta ) {
        fEta_LAB[i] = eta;        i++;
    }

    std::system("mkdir save");
    
    write_array_to_file("save/save.fT_pbar.txt", fN_T_pbar, fT_pbar_LAB );
    write_array_to_file("save/save.fTn_p.txt",   fN_Tn_p,   fTn_p_LAB   );
    write_array_to_file("save/save.fEta.txt",    fN_eta,    fEta_LAB    );
    
    // Setup CM grid
    double dS=log10(fMax_s/fMin_s)/(fN_s-1);
    i=0;
    for ( double d=log10(fMin_s); d<log10(fMax_s)+dS/2; d+=dS ) {
        fS[i] = pow(10, d);
        i++;
    }
    //    double dXR=(fMax_xR-fMin_xR)/(fN_xR-1);
    //    i=0;
    //    for ( double xR=fMin_xR; xR<fMax_xR+dXR/2; xR+=dXR ) {
    //        fxR[i] = xR;    i++;
    //    }
    double dXR=log10(fMax_xR/fMin_xR)/(fN_xR-1);
    i=0;
    for ( double d=log10(fMin_xR); d<log10(fMax_xR)+dXR/2; d+=dXR ) {
        fxR[i] = pow(10, d);    i++;
    }
    double dXF = 2./(fN_xR-1);
    i=0;
    for ( double x=-1; x<1+dXF/2; x+=dXF ) {
        fxF[i] = x;    i++;
    }
    double dpT=log10(fMax_pT/fMin_xR)/(fN_pT-1);
    i=0;
    for ( double d=log10(fMin_pT); d<log10(fMax_pT)+dpT/2; d+=dpT ) {
        fPT[i] = pow(10, d);    i++;
    }
    
    write_array_to_file("save/save.fS.txt",     fN_s,       fS      );
    write_array_to_file("save/save.fPT.txt",    fN_pT,      fPT     );
    write_array_to_file("save/save.fxR.txt",    fN_xR,      fxR     );
    
    // Set option
    if (option==0) {
        fOption = "standard";           //      0:  standard:       pp,  di Mauro        AMS-02      0.03, 0.30
    }else if (option==1){
        fOption = "pHe";                //      1:  pHe :           pHe, di Mauro        AMS-02      0.03, 0.30
        fWithPlots=false;
        f__inv_pp_pbar_LAB__parametrization = ppCSParametrizations::inv_pHe_pbar_CM__diMauro12;
    }else if (option==2){
        fOption = "Duperray";           //      2:  Duperray :      pp,  Duperray        AMS-02      0.03, 0.30
        fWithPlots=kPAReturnPixels;
        f__inv_pp_pbar_LAB__parametrization = ppCSParametrizations::inv_pp_pbar_CM__Duperray;
    }else if (option==3){
        fOption = "3_100";              //      3:  3_100 :         pp,  di Mauro        AMS-02      0.03, 1.00
        fWithPlots=false;
        fErr_inner=0.03;
        fErr_outer=1.00;
    }else if (option==4){
        fOption = "2_30";               //      4:  2_30 :          pp,  di Mauro        AMS-02      0.02, 0.30
        fWithPlots=false;
        fErr_inner=0.02;
        fErr_outer=0.30;
    }else if (option==5){
        fOption = "pp_fixed";           //      5:  pp_fixed :      pHe,  di Mauro       AMS-02      0.03, 0.30         pp fixed
        fWithPlots= false;
        fPPfixed  = true;
        f__inv_pp_pbar_LAB__parametrization = ppCSParametrizations::inv_pHe_pbar_CM__diMauro12;
    }else if (option==6){
        fOption = "pp_fixed_6_60";      //      6: pp_fixed_6_60 :  pHe,  di Mauro      AMS-02      0.06, 0.60         pp fixed
        fWithPlots= false;
        fPPfixed  = true;
        fErr_inner= 0.06;
        fErr_outer= 0.60;
        f__inv_pp_pbar_LAB__parametrization = ppCSParametrizations::inv_pHe_pbar_CM__diMauro12;
    }else if (option==7) {
        fWithPlots=false;
        fOption = "sm_700";             //      7:  sm_700:         pp,  di Mauro        AMS-02        0.03, 0.30
        fSM     = 0.7;
    }else if (option==8){
        fOption = "Winkler";            //      8:  Winkler :       pp,  Winkler         AMS-02        0.03, 0.30
        fWithPlots=true;
        f__inv_pp_pbar_LAB__parametrization = ppCSParametrizations::inv_pp_pbar_CM__WinklerWithHypWithNbar;
            
    }else if (option==100){
        fWithPlots=false;
        fOption = "GAPS";               //    100:  GAPS:           pp, di Mauro         GAPS        0.02, 1.00
    }else if (option==101) {
        fWithPlots=false;
        fOption = "2024";               //    101:  AMS_2024:       pp,  di Mauro        AMS-02('24) 0.03, 0.30
    }else if (option==102) {
        fWithPlots=false;
        fOption = "0_05";               //    102:  0_05:           pp,  di Mauro        0.05        0.03, 0.30
    }else if (option==103) {
        fWithPlots=false;
        fOption = "0_05__D";            //    102:  0_05__D:        pp,  Duperray        0.05        0.03, 0.30
        f__inv_pp_pbar_LAB__parametrization = ppCSParametrizations::inv_pp_pbar_CM__Duperray;
    }else if (option==104) {
        fWithPlots=false;
        fOption = "0_05__W";            //    102:  0_05__W:        pp,  Winkler        0.05        0.03, 0.30
        f__inv_pp_pbar_LAB__parametrization = ppCSParametrizations::inv_pp_pbar_CM__WinklerWithHypWithNbar;
    }


    LISFluxes::fit();
    
    // Parametrize pbar uncertainty
    if (step==1) {
        
        std::system( ("mkdir "+fOption                    ).c_str()         );
        std::system( ("mkdir "+fOption+"/pbar_uncertainty").c_str()         );
        
        if (option<100) {
            ReadFluxes reader;
            Graph pbar = reader.GetSpectrum("pbar", "ams", EKINPERN);
            double pbarT                  [pbar.Size()];
            double log_pbarT              [pbar.Size()];
            double pbarRelativeUncertainty[pbar.Size()];
            
            for (int i = 0; i<pbar.Size(); i++) {
                double x, xe, y, ye;
                pbar.GetPoint(i, x, y, xe, ye); x+=fSM;
                pbarT                  [i] = x;
                log_pbarT              [i] = log(x);
                pbarRelativeUncertainty[i] = ye/y;
            }
            write_array_to_file(fOption+"/pbar_uncertainty/pbar_source_uncertainty_data.txt", pbar.Size(), pbarRelativeUncertainty, pbarT,   "T_pbar      sigma_flux/flux    info: T_pbar is demodulated by 500 MeV");
            
            TGraph pbarR(pbar.Size(), &log_pbarT[0], &pbarRelativeUncertainty[0]);
            TF1 fun_pbarR("fun_pbarR", "pol8");
            pbarR.Fit(&fun_pbarR);
            
            for ( int i=0; i<fN_T_pbar; i++ ) {
                fRU_Tpbar__Q[i] = fun_pbarR.Eval( log(fT_pbar_LAB[i]) );
            }
        }
        if(option==100){
            for ( int i=0; i<fN_T_pbar; i++ ) {
                double x = log10(fT_pbar_LAB[i]);
                fRU_Tpbar__Q[i] = pow( 2*x +0.5, 8)+0.05;
            }
        }
        if(option==101){
            std::string filename = config->GetInstance()->SoftwarePath()+"/data/AMS02/pbar_ams02.txt";
            
            CRACS::Graph pbarAMS;
            CRACS::FileTool f(filename);
            f.ExtractNumberTable(14, " ", true);
            std::vector<double> rel_stat;
            std::vector<double> rel_sys;
            for (int row = 0; row<f.NumberTableGetNrows(); row++) {
                double R = CRACS::getBinExpectationValueForPowerLaw2(f.NumberTable(row, 0), f.NumberTable(row, 1), -2.7);
                double F = f.NumberTable(row, 10)*f.NumberTable(row, 13);
                double helper = 0;
                for (int col = 11; col<=12; col++) {
                    helper += pow(f.NumberTable(row, col), 2);
                }
                double sigma = sqrt(helper)*f.NumberTable(row, 13);
                
                double sigma_stat = f.NumberTable(row, 11)*f.NumberTable(row, 13);
                double sigma_sys  = f.NumberTable(row, 12)*f.NumberTable(row, 13);
                
                
                rel_sys. push_back(sigma_sys /F);
                rel_stat.push_back(sigma_stat/F);
                
                pbarAMS.AddPoint(R*1000, F*pow(R*1000, 2.7)/1000, 0, sigma*pow(R*1000, 2.7)/1000);
            }
            
            pbarAMS.ConvertSpectrumFromRigidityToEkinPerN(-1, 1, 2.7);
            
            double pbarT                  [pbarAMS.Size()];
            double log_pbarT              [pbarAMS.Size()];
            double pbarRelativeUncertainty[pbarAMS.Size()];
            for (int i = 0; i<pbarAMS.Size(); i++) {
                double x, xe, y, ye;
                pbarAMS.GetPoint(i, x, y, xe, ye); x+=600;
                pbarT                  [i] = x/1000;
                log_pbarT              [i] = log(x/1000);
                pbarRelativeUncertainty[i] = sqrt( rel_stat[i]*rel_stat[i]*4./13. + rel_sys[i]*rel_sys[i] );
            }
            
            write_array_to_file(fOption+"/pbar_uncertainty/pbar_source_uncertainty_data.txt", pbarAMS.Size(), pbarRelativeUncertainty, pbarT,   "T_pbar      sigma_flux/flux    info: T_pbar is demodulated by 500 MeV");
            
            TGraph pbarR(pbarAMS.Size(), &log_pbarT[0], &pbarRelativeUncertainty[0]);
            TF1 fun_pbarR("fun_pbarR", "pol8");
            pbarR.Fit(&fun_pbarR);
            
            for ( int i=0; i<fN_T_pbar; i++ ) {
                fRU_Tpbar__Q[i] = fun_pbarR.Eval( log(fT_pbar_LAB[i]) );
            }
        }
        if(option==102 || option==103 || option==104){
            for ( int i=0; i<fN_T_pbar; i++ ) {
                fRU_Tpbar__Q[i] = 0.05;
            }
        }
        
        if (fPPfixed) {
            for ( int i=0; i<fN_T_pbar; i++ ) {
                double T_pbar   = fT_pbar_LAB[i];
                double Tn_from  = std::min(T_pbar, 6*fMass_proton);
                
                int p= 10000;
                double S_p_H    = SourceTerms::integrationFromTo(T_pbar, dT_pp_pbar_LAB_dM,   protonFlux, f_nH,  Tn_from, 1e6, p);
                double S_p_He   = SourceTerms::integrationFromTo(T_pbar, dT_pHe_pbar_LAB_dM,  protonFlux, f_nHe, Tn_from, 1e6, p);
                double S_He_H   = SourceTerms::integrationFromTo(T_pbar, dT_Hep_pbar_LAB_dM,  heliumFlux, f_nH,  Tn_from, 1e6, p);
                double S_He_He  = SourceTerms::integrationFromTo(T_pbar, dT_HeHe_pbar_LAB_dM, heliumFlux, f_nHe, Tn_from, 1e6, p);
                double S        = S_p_H + S_p_He + S_He_H + S_He_He;
                
                fRU_Tpbar__Q[i] *= S/(S-S_p_H);
            }
        }

        
        
        write_array_to_file(  fOption+"/pbar_uncertainty/pbar_source_uncertainty.txt",    fN_T_pbar,  fRU_Tpbar__Q,   fT_pbar_LAB,    "T_pbar      rel. uncertainty required on the source term" );
        write_array_to_file(  "save/save."+fOption+".pbar_source_uncertainty.txt",        fN_T_pbar,  fRU_Tpbar__Q  );
        
        
    }
    // Unfolding in Tp and eta
    // calculate 3D contribtuion (LAB)
    if (step==2) {
        
        /////// read calculations from previous steps
        read_array_from_file("save/save."+fOption+".pbar_source_uncertainty.txt",     fN_T_pbar,      fRU_Tpbar__Q);
        
        if (fWithPlots){
            std::system(  ("mkdir "+fOption+"/contribution_LAB"             ).c_str() );
            std::system(  ("mkdir "+fOption+"/contribution_LAB/contribution").c_str() );
            std::system(  ("mkdir "+fOption+"/contribution_LAB/containence" ).c_str() );
        }
        
        for (int j=0; j<fN_T_pbar; j++) {
            double T_pbar   = fT_pbar_LAB[j];
            double sum = 0;
            for (int i=0; i<fN_Tn_p-1; i++) {
                double Tn_p = fTn_p_LAB[i];
                for (int k=0; k<fN_eta-1; k++) {
                    double eta = fEta_LAB[k];
                    double p_pbar = sqrt(T_pbar*T_pbar+2*fMass_proton*T_pbar);
                    f3D__Tpbar_Tp_eta__a1      [i + fN_Tn_p*k + fN_Tn_p*fN_eta*j] = 8*M_PI*M_PI *f_nH * p_pbar * Tn_p *protonFlux(Tn_p) / pow(cosh(eta), 2) * inv_pp_pbar_LAB(Tn_p, T_pbar, eta);
                    sum+=f3D__Tpbar_Tp_eta__a1 [i + fN_Tn_p*k + fN_Tn_p*fN_eta*j];
                }
            }
            for (int i=0; i<fN_Tn_p-1; i++) {
                for (int k=0; k<fN_eta-1; k++) {
                    f3D__Tpbar_Tp_eta__a1 [i + fN_Tn_p*k + fN_Tn_p*fN_eta*j]/=sum;
                }
            }
            out(fT_pbar_LAB[j])
            get_containence(fN_Tn_p*fN_eta, &f3D__Tpbar_Tp_eta__a1 [j*fN_eta*fN_Tn_p], &f3D__Tpbar_Tp_eta__a2 [j*fN_eta*fN_Tn_p] );
            
            if (fWithPlots){
                write_2D_array_to_file( fOption+QString("/contribution_LAB/contribution/2D_contribution_Tp_eta__%1.txt").arg(j).toStdString(),  fN_eta,  fN_Tn_p,   &f3D__Tpbar_Tp_eta__a1[j*fN_eta*fN_Tn_p],   fEta_LAB, fTn_p_LAB  );
                write_2D_array_to_file( fOption+QString("/contribution_LAB/containence//2D_containence__Tp_eta__%1.txt").arg(j).toStdString(),  fN_eta,  fN_Tn_p,   &f3D__Tpbar_Tp_eta__a2[j*fN_eta*fN_Tn_p],   fEta_LAB, fTn_p_LAB  );
            }
        }
        
        ////// save calculation to files
        write_array_to_file("save/save."+fOption+"_contribution_LAB.txt",  fN_T_pbar*fN_Tn_p*fN_eta, f3D__Tpbar_Tp_eta__a1 );
        write_array_to_file("save/save."+fOption+"_containence__LAB.txt",  fN_T_pbar*fN_Tn_p*fN_eta, f3D__Tpbar_Tp_eta__a2 );
        
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    if (step==22) {
        
        /////// read calculations from previous steps
        read_array_from_file("save/save."+fOption+"_containence__LAB.txt",          fN_T_pbar*fN_Tn_p*fN_eta,   f3D__Tpbar_Tp_eta__a1   );
        
        std::system(  ("mkdir "+fOption+"/contribution_CM/           ").c_str() );
        std::system(  ("mkdir "+fOption+"/contribution_CM/containence").c_str() );
       
        for (int j=0; j<fN_s; j++) {
            double s = fS[j];
            
            for (int i=0; i<fN_xR; i++) {
                double xF = fxF[i];
                
                for (int k=0; k<fN_pT; k++) {
                    double pT = fPT[k];
                    
                    double pL           = xF/2.*sqrt(s);
                    double E_pbar       = sqrt( pT*pT + pL*pL+ fMass_proton*fMass_proton );
                    double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
                    
                    double xR = E_pbar/E_pbar_Max;
                    
                    bool option = true;
                    if (pL<0)
                    option = false;
                    
                    double T_p_LAB, T_pbar_LAB, eta_LAB;
                    CSTransformations::convert_CM_to_LAB(s, xR, pT, T_p_LAB, T_pbar_LAB, eta_LAB, option);
                    
                    double res = 1.0;
                    //                    if (T_pbar_LAB>48 && T_pbar_LAB<52){
                    //                        out(s)
                    //                        out(T_pbar_LAB)
                    //                        out("******")
                    res = Interpolate_LAB(T_p_LAB, T_pbar_LAB, eta_LAB, f3D__Tpbar_Tp_eta__a1 );
                    //                    }
                    
                    
                    f3D__s_xR_pT [i + fN_pT*k + fN_xR*fN_pT*j] = res;
                }
            }
            out(s)
            write_2D_array_to_file( fOption+QString("/contribution_CM/containence/2D_contribution__xF_pT__%1.txt").    arg(j).toStdString(),     fN_pT,  fN_xR,   &f3D__s_xR_pT    [j*fN_pT*fN_xR],   fPT, fxF  );
            
        }
        
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    
    
    
    // calculate 3D relative uncertatinty
    if (step==3) {
        
        /////// read calculations from previous steps
        read_array_from_file("save/save."+fOption+".pbar_source_uncertainty.txt",   fN_T_pbar,                  fRU_Tpbar__Q            );
        read_array_from_file("save/save."+fOption+"_containence__LAB.txt",          fN_T_pbar*fN_Tn_p*fN_eta,   f3D__Tpbar_Tp_eta__a1   );
        
        std::system(  ("mkdir "+fOption+"/contribution_LAB/           ").c_str() );
        std::system(  ("mkdir "+fOption+"/contribution_LAB/uncertainty").c_str() );
        
        for (int j=0; j<fN_T_pbar; j++) {
            out(fT_pbar_LAB[j])
            double threshold  = (fErr_outer-fRU_Tpbar__Q[j])/(fErr_outer-fErr_inner);
            for (int i=0; i<fN_Tn_p-1; i++) {
                for (int k=0; k<fN_eta-1; k++) {
                    if (  f3D__Tpbar_Tp_eta__a1[i + fN_Tn_p*k + fN_Tn_p*fN_eta*j] < threshold  ){
                        f3D__Tpbar_Tp_eta__a2[i + fN_Tn_p*k + fN_Tn_p*fN_eta*j] = 0;
                    }else{
                        f3D__Tpbar_Tp_eta__a2[i + fN_Tn_p*k + fN_Tn_p*fN_eta*j] = 1;
                    }
                }
            }
            write_2D_array_to_file( fOption+QString("/contribution_LAB/uncertainty/2D_uncertainty__Tp_eta__%1.txt").arg(j).toStdString(),      fN_eta,  fN_Tn_p,   &f3D__Tpbar_Tp_eta__a2[j*fN_eta*fN_Tn_p],    fEta_LAB, fTn_p_LAB  );
        }
    
        write_array_to_file("save/save."+fOption+"_uncertainty__LAB.txt",  fN_T_pbar*fN_Tn_p*fN_eta, f3D__Tpbar_Tp_eta__a2 );
        
        std::system(  ("mkdir "+fOption+"/contribution_CM"            ).c_str() );
        std::system(  ("mkdir "+fOption+"/contribution_CM/uncertainty").c_str() );
        
        
        for (int j=0; j<fN_s; j++) {
            double s = fS[j];
            
            for (int i=0; i<fN_xR-1; i++) {
                double xR = fxR[i+1];
                
                for (int k=0; k<fN_pT-1; k++) {
                    double pT = fPT[k];
                    
                    f3D__s_xR_pT       [i + fN_pT*k + fN_xR*fN_pT*j] = GetInCM( s, xR, pT, f3D__Tpbar_Tp_eta__a2, 0  );
                }
            }
            out(s)
            write_2D_array_to_file( fOption+QString("/contribution_CM/uncertainty/2D_uncertainty__xR_pT__%1.txt").    arg(j).toStdString(),     fN_pT,  fN_xR,   &f3D__s_xR_pT    [j*fN_pT*fN_xR],   fPT, fxR  );
            
        }
        
    }

    //**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**
    //**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**
    //**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**
    //**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**
    
    // calculate 3D relative uncertatinty
    if (step==33) {
        
        /////// read calculations from previous steps
        read_array_from_file("save/save."+fOption+".pbar_source_uncertainty.txt",   fN_T_pbar,                  fRU_Tpbar__Q            );
        read_array_from_file("save/save."+fOption+"_containence__LAB.txt",          fN_T_pbar*fN_Tn_p*fN_eta,   f3D__Tpbar_Tp_eta__a1   );
        
        std::system(  ("mkdir "+fOption+"/contribution_LAB/           ").c_str() );
        std::system(  ("mkdir "+fOption+"/contribution_LAB/uncertainty").c_str() );
        
        for (int j=0; j<fN_T_pbar; j++) {
            out(fT_pbar_LAB[j])
            double threshold  = (fErr_outer-fRU_Tpbar__Q[j])/(fErr_outer-fErr_inner);
            for (int i=0; i<fN_Tn_p-1; i++) {
                for (int k=0; k<fN_eta-1; k++) {
                    if (  f3D__Tpbar_Tp_eta__a1[i + fN_Tn_p*k + fN_Tn_p*fN_eta*j] < threshold  ){
                        f3D__Tpbar_Tp_eta__a2[i + fN_Tn_p*k + fN_Tn_p*fN_eta*j] = 0;
                    }else{
                        f3D__Tpbar_Tp_eta__a2[i + fN_Tn_p*k + fN_Tn_p*fN_eta*j] = 1;
                    }
                }
            }
        }
        
        std::system(  ("mkdir "+fOption+"/contribution_CM"            ).c_str() );
        std::system(  ("mkdir "+fOption+"/contribution_CM/uncertainty").c_str() );
        
        
        for (int j=0; j<fN_s; j++) {
            double s = fS[j];
            
            for (int i=0; i<fN_xR-1; i++) {
                double xF = fxF[i];
                
                for (int k=0; k<fN_pT-1; k++) {
                    double pT = fPT[k];
                    f3D__s_xR_pT       [i + fN_pT*k + fN_xR*fN_pT*j] = GetInCM_xF( s, xF, pT, f3D__Tpbar_Tp_eta__a2  );
                    

//                    double pL           = xF/2.*sqrt(s);
//                    double E_pbar       = sqrt( pT*pT + pL*pL+ fMass_proton*fMass_proton );
//                    double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
//
//                    double xR = E_pbar/E_pbar_Max;
//
//                    bool option = true;
//                    if (pL<0)
//                    option = false;
//
//                    double T_p_LAB, T_pbar_LAB, eta_LAB;
//                    CSTransformations::convert_CM_to_LAB(s, xR, pT, T_p_LAB, T_pbar_LAB, eta_LAB, option);
//                    
//                    f3D__s_xR_pT       [i + fN_pT*k + fN_xR*fN_pT*j] = Interpolate_LAB(T_p_LAB, T_pbar_LAB, eta_LAB, f3D__Tpbar_Tp_eta__a2);
                    
                }
            }
            out(s)
            write_2D_array_to_file( fOption+QString("/contribution_CM/uncertainty/2D_uncertainty__xF_pT__%1.txt").    arg(j).toStdString(),     fN_pT,  fN_xR,   &f3D__s_xR_pT    [j*fN_pT*fN_xR],   fPT, fxF  );
            
        }
        
    }
    //**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**
    //**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**
    //**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**
    //**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**
    //**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**
    //**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**
    //**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**-**
    

    // calculate 3D relative uncertatinty for data comparison
    if (step==4) {
        
        /////// read calculations from previous steps
        read_array_from_file("save/save."+fOption+".pbar_source_uncertainty.txt",   fN_T_pbar,                  fRU_Tpbar__Q            );
        read_array_from_file("save/save."+fOption+"_uncertainty__LAB.txt",          fN_T_pbar*fN_Tn_p*fN_eta,   f3D__Tpbar_Tp_eta__a1   );
        
        std::system(  ("mkdir "+fOption+"/contribution_CM/           ").c_str() );
        std::system(  ("mkdir "+fOption+"/contribution_CM/pp_CS_data ").c_str() );
        
        double sqrtS_list[] = {6.2, 13.8, 17.3, 19.4, 23.5, 27.4, 30.8, 44.8, 53., 63., 200  };
        write_array_to_file(fOption+"/contribution_CM/pp_CS_data/sqrtS_values.txt", 11, sqrtS_list);
        
        for (int j=0; j<11; j++) {
            double s = sqrtS_list[j]*sqrtS_list[j];
            for (int i=0; i<fN_xR-1; i++) {
                double xR = fxR[i+1];
                for (int k=0; k<fN_pT-1; k++) {
                    double pT = fPT[k];
                    f2D__xR_pT     [i + fN_pT*k] = GetInCM( s, xR, pT, f3D__Tpbar_Tp_eta__a1, 0  );
                }
            }
            out(s)
            write_2D_array_to_file( fOption+QString("/contribution_CM/pp_CS_data/2D_uncertainty__xR_pT__%1.txt").arg(j).toStdString(),     fN_pT,  fN_xR,   f2D__xR_pT,   fPT, fxR  );
        }
    }
    
    

    // extract uncertainty profile for fixed Tp
    if (step==5) {
        
        /////// read calculations from previous steps
        read_array_from_file("save/save.pbar_relative_uncertainty_parmetrization.txt",     fN_T_pbar,      fRU_Tpbar__Q);
        read_array_from_file("save/save."+fOption+"_uncertainty__LAB.txt",          fN_T_pbar*fN_Tn_p*fN_eta,   f3D__Tpbar_Tp_eta__a1   );
        
        
        std::system(  ("mkdir "+fOption+"/contribution_LAB/uncertainty_fixedTp").c_str() );
        
        int n = 7;
        double Tp_list[] = {10, 50, 150, 450, 2e3, 4e3, 6.5e3};
        write_array_to_file(fOption+"/contribution_LAB/uncertainty_fixedTp/Tp_values.txt", n, Tp_list);
        
        
        int Tp_index[n];
        
        // find closest T_p index
        for (int i=0; i<n; i++) {
            double diff=1e90;
            for (int j=0; j<fN_Tn_p; j++) {
                double test = fabs(log(Tp_list[i]/fTn_p_LAB[j]));
                if (test<diff) {
                    diff=test;
                    Tp_index[i] = j;
                }
            }
        }
        
        
        for (int ii=0; ii<n; ii++) {
            int i = Tp_index[ii];
            for (int j=0; j<fN_T_pbar; j++) {
                for (int k=0; k<fN_eta-1; k++) {
                        f2D__Tpbar_eta[j + fN_Tn_p*k] = f3D__Tpbar_Tp_eta__a1[i + fN_Tn_p*k + fN_Tn_p*fN_eta*j];
                }
            }
            write_2D_array_to_file( fOption+QString("/contribution_LAB/uncertainty_fixedTp/2D_uncertainty_Tpbar_eta__%1.txt").arg(ii).toStdString(),      fN_eta,  fN_T_pbar,   f2D__Tpbar_eta,    fEta_LAB, fT_pbar_LAB  );
        }
        
        
        
    }
    
    
    // extract uncertainty profile for fixed Tp, in T_pbar and pT
    if (step==52) {
        
        /////// read calculations from previous steps
        read_array_from_file("save/save.pbar_relative_uncertainty_parmetrization.txt",     fN_T_pbar,      fRU_Tpbar__Q);
        read_array_from_file("save/save."+fOption+"_uncertainty__LAB.txt",          fN_T_pbar*fN_Tn_p*fN_eta,   f3D__Tpbar_Tp_eta__a1   );
        
        
        std::system(  ("mkdir "+fOption+"/contribution_LAB/uncertainty_fixedTp").c_str() );
        
        int n = 8;
        double Tp_list[] = {10,20, 50, 150, 450, 2e3, 4e3, 6.5e3};
        write_array_to_file(fOption+"/contribution_LAB/uncertainty_fixedTp/Tp_values.txt", n, Tp_list);
        
        
        int Tp_index[n];
        
        // find closest T_p index
        for (int i=0; i<n; i++) {
            double diff=1e90;
            for (int j=0; j<fN_Tn_p; j++) {
                double test = fabs(log(Tp_list[i]/fTn_p_LAB[j]));
                if (test<diff) {
                    diff=test;
                    Tp_index[i] = j;
                }
            }
        }
        
        double p_pbar_a[fN_T_pbar];
        for (int ii=0; ii<n; ii++) {
            for (int j=0; j<fN_T_pbar; j++) {
                double p_pbar = fT_pbar_LAB[j];
                p_pbar_a[j] = p_pbar;
                for (int k=0; k<fN_pT-1; k++) {
                    double pT = fPT[k];
                    
                    double T   = sqrt(p_pbar*p_pbar+fMass_proton*fMass_proton)-fMass_proton;
                    double eta = acosh(p_pbar/pT);
                    
//                    out(p_pbar)
//                    out(pT)
//                    out(eta)
//                    out(T)
                    
                    double u = 0.3;
                    if (pT<p_pbar)
                        u  = Interpolate_LAB( Tp_list[ii], T, eta, f3D__Tpbar_Tp_eta__a1);
                    f2D__Tpbar_eta[j + fN_T_pbar*k] = u;
                }
            }
            write_2D_array_to_file( fOption+QString("/contribution_LAB/uncertainty_fixedTp/2D_uncertainty_p_pbar_p_T__%1.txt").arg(ii).toStdString(),      fN_pT,  fN_T_pbar,   f2D__Tpbar_eta,    fPT, p_pbar_a  );
        }
        
        
        
    }
    
    
    
    
    
    //    if (step==3 && option==0) {
    //
    //        /////// read calculations from previous steps
    //        read_array_from_file("save/save.pbar_relative_uncertainty_parmetrization.txt",     fN_T_pbar,      fRU_Tpbar__Q);
    //        read_array_from_file("save/save.source_term_tot.txt",  fN_T_pbar,   fQ      );
    //        read_array_from_file("save/save.source_term_pp.txt",   fN_T_pbar,   fQ_pp   );
    //        read_array_from_file("save/save.source_term_pHe.txt",  fN_T_pbar,   fQ_Hep  );
    //        read_array_from_file("save/save.source_term_Hep.txt",  fN_T_pbar,   fQ_pHe  );
    //        read_array_from_file("save/save.source_term_HeHe.txt", fN_T_pbar,   fQ_HeHe );
    //        ///////
    //
    //
    //        std::system("mkdir 2D_contribution_LAB ");
    //
    //        std::system("mkdir 2D_contribution_LAB/pp  ");
    //        std::system("mkdir 2D_contribution_LAB/pHe ");
    //        std::system("mkdir 2D_contribution_LAB/Hep ");
    //
    //
    //        // contribution and commulative contribution in T_pbar T_p 2D
    //        for (int j=0; j<fN_T_pbar; j++) {
    //            double T_pbar = fT_pbar_LAB[j];
    //            //double dLog_T_p = log(fMax_Tn_p/fMin_Tn_p)/fN_Tn_p;
    //            for (int i=0; i<fN_Tn_p-1; i++) {
    //                double Tn_p = fTn_p_LAB[i+1];
    //                fC_Tpbar_Tp__pp [i + fN_Tn_p*j] = 4*M_PI *f_nH  * Tn_p * dT_pp_pbar_LAB (Tn_p, T_pbar)*protonFlux(Tn_p)/fQ_pp[j];
    //                fC_Tpbar_Tp__pHe[i + fN_Tn_p*j] = 4*M_PI *f_nHe * Tn_p * dT_pHe_pbar_LAB(Tn_p, T_pbar)*protonFlux(Tn_p)/fQ_pHe[j];
    //                fC_Tpbar_Tp__Hep[i + fN_Tn_p*j] = 4*M_PI *f_nHe * Tn_p * dT_pHe_pbar_LAB(Tn_p, T_pbar)*protonFlux(Tn_p)/fQ_Hep[j];
    //            }
    //            commulate(fN_Tn_p, &fC_Tpbar_Tp__pp [0 + fN_Tn_p*j], &fCC_Tpbar_Tp__pp [0 + fN_Tn_p*j]);
    //            commulate(fN_Tn_p, &fC_Tpbar_Tp__pHe[0 + fN_Tn_p*j], &fCC_Tpbar_Tp__pHe[0 + fN_Tn_p*j]);
    //            commulate(fN_Tn_p, &fC_Tpbar_Tp__Hep[0 + fN_Tn_p*j], &fCC_Tpbar_Tp__Hep[0 + fN_Tn_p*j]);
    //        }
    //
    //
    //        write_2D_array_to_file( "2D_contribution_LAB/pp/2D_contribution__pp.txt",       fN_T_pbar,  fN_Tn_p,   &fC_Tpbar_Tp__pp[0],    fT_pbar_LAB, fTn_p_LAB  );
    //        write_2D_array_to_file( "2D_contribution_LAB/pp/2D_ccontribution__pp.txt",      fN_T_pbar,  fN_Tn_p,   &fCC_Tpbar_Tp__pp[0],   fT_pbar_LAB, fTn_p_LAB  );
    //
    //        write_2D_array_to_file( "2D_contribution_LAB/pHe/2D_contribution__pHe.txt",     fN_T_pbar,  fN_Tn_p,   &fC_Tpbar_Tp__pHe[0],   fT_pbar_LAB, fTn_p_LAB  );
    //        write_2D_array_to_file( "2D_contribution_LAB/pHe/2D_ccontribution__pHe.txt",    fN_T_pbar,  fN_Tn_p,   &fCC_Tpbar_Tp__pHe[0],  fT_pbar_LAB, fTn_p_LAB  );
    //
    //        write_2D_array_to_file( "2D_contribution_LAB/Hep/2D_contribution__Hep.txt",     fN_T_pbar,  fN_Tn_p,   &fC_Tpbar_Tp__Hep[0],   fT_pbar_LAB, fTn_p_LAB  );
    //        write_2D_array_to_file( "2D_contribution_LAB/Hep/2D_ccontribution__Hep.txt",    fN_T_pbar,  fN_Tn_p,   &fCC_Tpbar_Tp__Hep[0],  fT_pbar_LAB, fTn_p_LAB  );
    //
    //    }
    //
    //
    //
    
    
    //    // calculate 3D contribtuion (CM)
    //    if (step==5) {
    //
    //        /////// read calculations from previous steps
    //        read_array_from_file("save/save.pbar_relative_uncertainty_parmetrization.txt",     fN_T_pbar,      fRU_Tpbar__Q);
    //        read_array_from_file("save/save.source_term_tot.txt",  fN_T_pbar,   fQ      );
    //        read_array_from_file("save/save.source_term_pp.txt",   fN_T_pbar,   fQ_pp   );
    //        read_array_from_file("save/save.source_term_pHe.txt",  fN_T_pbar,   fQ_Hep  );
    //        read_array_from_file("save/save.source_term_Hep.txt",  fN_T_pbar,   fQ_pHe  );
    //        read_array_from_file("save/save.source_term_HeHe.txt", fN_T_pbar,   fQ_HeHe );
    //        read_array_from_file("save/save.fCC_Tpbar_Tp_eta_LAB__pp.txt", fN_T_pbar*fN_Tn_p*fN_eta, fCC_Tpbar_Tp_eta_LAB__pp);
    //        ///////
    //
    //        double Tp    = 50;
    //        double Tpbar = 5;
    //        double eta   = 3;
    //
    //        out(InterpolateCC_LAB(Tp, Tpbar, eta));
    //
    //        double s;
    //        double xR;
    //        double pT;
    //
    //        out(Tp)
    //        out(Tpbar)
    //        out(eta)
    //
    //        out("to CM")
    //        CSTransformations::convert_LAB_to_CM(Tp, Tpbar, eta, s, xR, pT);
    //        out(s)
    //        out(xR)
    //        out(pT)
    //
    //        out(GetLab(s, xR, pT))
    //
    //        out("back to LAB pL>0")
    //        CSTransformations::convert_CM_to_LAB(s, xR, pT, Tp, Tpbar, eta, true);
    //        out(Tp)
    //        out(Tpbar)
    //        out(eta)
    //
    //
    //        out("back to LAB pL<0")
    //        CSTransformations::convert_CM_to_LAB(s, xR, pT, Tp, Tpbar, eta, false);
    //        out(Tp)
    //        out(Tpbar)
    //        out(eta)
    //
    //
    //
    //        std::system("mkdir 3D_s_xR_pT");
    //        std::system("mkdir 3D_s_xR_pT/pp");
    //
    //
    //        for (int j=0; j<fN_s; j++) {
    //            double s = fS[j];
    //            for (int i=0; i<fN_xR-1; i++) {
    //                double xR = fxR[i+1];
    //                for (int k=0; k<fN_pT-1; k++) {
    //                    double pT = fPT[k];
    //                    fCC_s_xR_pT__pp       [i + fN_pT*k + fN_xR*fN_pT*j] = GetLab( s, xR, pT,  1 );
    //                }
    //            }
    //            std::cout << std::endl;
    //            out(s)
    //            if (j%fMod==0){
    //                write_2D_array_to_file( QString("3D_s_xR_pT/pp/2D_ccontribution_xR_pT__%1__pp_pos.txt").arg(j).toStdString(),     fN_pT,  fN_xR,   &fCC_s_xR_pT__pp    [j*fN_pT*fN_xR],   fPT, fxR  );
    //            }
    //        }
    //        for (int j=0; j<fN_s; j++) {
    //            double s = fS[j];
    //            for (int i=0; i<fN_xR-1; i++) {
    //                double xR = fxR[i+1];
    //                for (int k=0; k<fN_pT-1; k++) {
    //                    double pT = fPT[k];
    //                    fCC_s_xR_pT__pp       [i + fN_pT*k + fN_xR*fN_pT*j] = GetLab( s, xR, pT,  2 );
    //                }
    //            }
    //            std::cout << std::endl;
    //            out(s)
    //            if (j%fMod==0){
    //                write_2D_array_to_file( QString("3D_s_xR_pT/pp/2D_ccontribution_xR_pT__%1__pp_neg.txt").arg(j).toStdString(),     fN_pT,  fN_xR,   &fCC_s_xR_pT__pp    [j*fN_pT*fN_xR],   fPT, fxR  );
    //            }
    //        }
    //        for (int j=0; j<fN_s; j++) {
    //            double s = fS[j];
    //            for (int i=0; i<fN_xR-1; i++) {
    //                double xR = fxR[i+1];
    //                for (int k=0; k<fN_pT-1; k++) {
    //                    double pT = fPT[k];
    //                    fCC_s_xR_pT__pp       [i + fN_pT*k + fN_xR*fN_pT*j] = GetLab( s, xR, pT    );
    //                }
    //            }
    //            std::cout << std::endl;
    //            out(s)
    //            if (j%fMod==0){
    //                write_2D_array_to_file( QString("3D_s_xR_pT/pp/2D_ccontribution_xR_pT__%1__pp.txt").    arg(j).toStdString(),     fN_pT,  fN_xR,   &fCC_s_xR_pT__pp    [j*fN_pT*fN_xR],   fPT, fxR  );
    //            }
    //        }
    //        
    //        // save calculation to files
    //        write_array_to_file("save/save.fCC_s_xR_pT__pp.txt", fN_s * fN_xR * fN_pT, fCC_s_xR_pT__pp);
    //        
    //    }
    

    double Tp    = 2e3;
    double Tpbar = 350;
    double eta   = 7;
    
    double s = 11.8*11.8;
    double xR=0.2;
    double pT=0.5;
    
    CSTransformations::convert_LAB_to_CM(  Tp, Tpbar, eta, s, xR, pT  );
    
    out(sqrt(s))
    out(xR)
    out(pT)
    
    
//
//    
//    double s = 11.8*11.8;
//    double xR=0.2;
//    double pT=0.5;
//    
//    
//    out("back to LAB pL>0")
//    CSTransformations::convert_CM_to_LAB(s, xR, pT, Tp, Tpbar, eta, true);
//    out(Tp)
//    out(Tpbar)
//    out(eta)
//    
//    
//    out("back to LAB pL<0")
//    CSTransformations::convert_CM_to_LAB(s, xR, pT, Tp, Tpbar, eta, false);
//    out(Tp)
//    out(Tpbar)
//    out(eta)
//    
//
//    s = 11.8*11.8;
//    xR=0.62;
//    pT=0.2;
//    
//    
//    out("back to LAB pL>0")
//    CSTransformations::convert_CM_to_LAB(s, xR, pT, Tp, Tpbar, eta, true);
//    out(Tp)
//    out(Tpbar)
//    out(eta)
//    
//    
//    out("back to LAB pL<0")
//    CSTransformations::convert_CM_to_LAB(s, xR, pT, Tp, Tpbar, eta, false);
//    out(Tp)
//    out(Tpbar)
//    out(eta)
//    
//    
//    
    
    
    return 0;
    
    
}


// Functions for data IO of arrays
void write_array_to_file(std::string filename, int n, double* array, double* xarray, std::string description){
    std::string write = "# "+description;
    for (int i=0; i<n; i++) {
        std::stringstream ss;
        ss << "\n";
        if (xarray!=NULL) {
            ss << xarray[i] << "     ";
        }
        ss << array[i];
        write += ss.str();
    }
    CRACS::FileTool::WriteStringToFile(write, filename);
};
void read_array_from_file(std::string filename, int n, double* array){
    std::ifstream infile(filename.c_str());
    std::string dummy;
    std::getline(infile, dummy);
    
    for (int i=0; i<n; i++) {
        double a;
        infile >> a;
        array[i]=a;
    }
};
void write_2D_array_to_file(std::string filename, int nx, int ny, double* array, double* xarray, double* yarray){
    std::string table = "";
    QString line("");
    line += QString("%1 ").arg( "Table" ).leftJustified(20);
    for (int j = 0; j<nx; j++) {
        line += QString("%1 ").arg( xarray[j] ).leftJustified(20);
    }
    table += line.toStdString();
    for (int i = 0; i<ny; i++) {
        table += "\n";
        QString line("");
        line += QString("%1 ").arg( yarray[i] ).leftJustified(20);
        for (int j = 0; j<nx; j++) {
            line += QString("%1 ").arg( array[i+j*ny] ).leftJustified(20);
        }
        table += line.toStdString();
    }
    FileTool::WriteStringToFile( table, filename );
};




// Funktion related to sort_quick
void change(int &a, int&b){
    int c = a;
    a = b;
    b = c;
}
void change(double &a, double&b){
    double c = a;
    a = b;
    b = c;
}
// sort_quick implements the QuickSort algorithm. It sorts the array double* array between positions n1 and n2.
// the array int* pos traces the sorting of array
// this is a recursive function
void sort_quick (double* array, int n1, int n2, int* pos, bool print=false) {
    if (print) {
        std::cout << n2-n1 << std::endl;
    }
    if (n1>=n2)
        return;
    int p=n1;
    int a=n2;
    while (abs(a-p)>0) {
        if (print) {
            std::cout <<"  " << abs(a-p) << std::endl;
        }
        if (p<a) {
            if (array[p]<array[a]) {
                a--;
            }else{
                change(a, p);
                change(array[a], array[p]);
                change(pos  [a], pos  [p]);
                a++;
            }
        }else{
            if (array[p]>array[a]) {
                a++;
            }else{
                change(a, p);
                change(array[a], array[p]);
                change(pos  [a], pos  [p]);
                a--;
            }
        }
    }
    sort_quick(array, n1, p-1, pos, print);
    sort_quick(array, p+1, n2, pos, print);
    return;
    
    
};
// this function uses the quick sort algorithm implemented above to determine the containence
void get_containence(int n, const double* contribution, double* containence){
    if (n>fN_eta*fN_Tn_p*fN_T_pbar) {
        std::cout << "Error: please use larger helper array!" << std::endl;
        return;
    }
    int position[n];
    for (int i=0; i<n; i++) {
        position[i]=i;
        fHelper [i]=contribution[i];
    }
    sort_quick(fHelper, 0, n-1, position);
    for (int i=n-2; i>=0; i--) {
        fHelper[i] = fHelper[i]+fHelper[i+1];
    }
    for (int i=0; i<n; i++) {
        containence[position[i]]=fHelper[i]/fHelper[0];
    }
}


double ttx(double x){
    x = 1-x;
    if (x<0) {
        return 0;
    }
    return x;
}
double Interpolate_LAB(double T_p_LAB, double T_pbar_LAB, double eta_LAB, double* array){
    double i_T_p_LAB        =  fN_Tn_p   * log(T_p_LAB   /fMin_Tn_p  )/log(fMax_Tn_p  /fMin_Tn_p  );
    double i_eta_LAB        =  fN_eta    * (eta_LAB-fMin_eta)/(fMax_eta-fMin_eta);
    double i_T_pbar_LAB     =  fN_T_pbar * log(T_pbar_LAB/fMin_T_pbar)/log(fMax_T_pbar/fMin_T_pbar);
    
    
    if (i_T_p_LAB   <0 || i_T_p_LAB     >=fN_Tn_p-1   ) return 5;
    if (i_T_pbar_LAB<0 || i_T_pbar_LAB  >=fN_T_pbar-1 ) return 10;
    if (i_eta_LAB   <0 || i_eta_LAB     >=fN_eta-1    ) return 20;
    
    
    int int_T_p_LAB        = i_T_p_LAB;
    int int_eta_LAB        = i_eta_LAB;
    int int_T_pbar_LAB     = i_T_pbar_LAB;
    
    varOut(T_p_LAB)
    varOut(T_pbar_LAB)
    varOut(eta_LAB)
    
    varOut(fTn_p_LAB[int_T_p_LAB])
    varOut(fT_pbar_LAB[int_T_pbar_LAB])
    varOut(fEta_LAB[int_eta_LAB])
    
    
    
    double r1  = array[ (int_T_p_LAB+0) + fN_Tn_p*(int_eta_LAB+0) + fN_Tn_p*fN_eta*(int_T_pbar_LAB+0) ];
    double r2  = array[ (int_T_p_LAB+1) + fN_Tn_p*(int_eta_LAB+0) + fN_Tn_p*fN_eta*(int_T_pbar_LAB+0) ];
    double r3  = array[ (int_T_p_LAB+0) + fN_Tn_p*(int_eta_LAB+1) + fN_Tn_p*fN_eta*(int_T_pbar_LAB+0) ];
    double r4  = array[ (int_T_p_LAB+0) + fN_Tn_p*(int_eta_LAB+0) + fN_Tn_p*fN_eta*(int_T_pbar_LAB+1) ];
    double r5  = array[ (int_T_p_LAB+1) + fN_Tn_p*(int_eta_LAB+1) + fN_Tn_p*fN_eta*(int_T_pbar_LAB+0) ];
    double r6  = array[ (int_T_p_LAB+1) + fN_Tn_p*(int_eta_LAB+0) + fN_Tn_p*fN_eta*(int_T_pbar_LAB+1) ];
    double r7  = array[ (int_T_p_LAB+0) + fN_Tn_p*(int_eta_LAB+1) + fN_Tn_p*fN_eta*(int_T_pbar_LAB+1) ];
    double r8  = array[ (int_T_p_LAB+1) + fN_Tn_p*(int_eta_LAB+1) + fN_Tn_p*fN_eta*(int_T_pbar_LAB+1) ];
    
    
    double diff_T_p     =  i_T_p_LAB    -   int_T_p_LAB;
    double diff_eta     =  i_eta_LAB    -   int_eta_LAB;
    double diff_T_pbar  =  i_T_pbar_LAB -   int_T_pbar_LAB;
    
    double d1  = sqrt(  pow(    diff_T_p  , 2  ) + pow(    diff_eta  , 2  ) + pow(    diff_T_pbar  , 2  )  )+1e-15;
    double d2  = sqrt(  pow(  1-diff_T_p  , 2  ) + pow(    diff_eta  , 2  ) + pow(    diff_T_pbar  , 2  )  )+1e-15;
    double d3  = sqrt(  pow(    diff_T_p  , 2  ) + pow(  1-diff_eta  , 2  ) + pow(    diff_T_pbar  , 2  )  )+1e-15;
    double d4  = sqrt(  pow(    diff_T_p  , 2  ) + pow(    diff_eta  , 2  ) + pow(  1-diff_T_pbar  , 2  )  )+1e-15;
    double d5  = sqrt(  pow(  1-diff_T_p  , 2  ) + pow(  1-diff_eta  , 2  ) + pow(    diff_T_pbar  , 2  )  )+1e-15;
    double d6  = sqrt(  pow(  1-diff_T_p  , 2  ) + pow(    diff_eta  , 2  ) + pow(  1-diff_T_pbar  , 2  )  )+1e-15;
    double d7  = sqrt(  pow(    diff_T_p  , 2  ) + pow(  1-diff_eta  , 2  ) + pow(  1-diff_T_pbar  , 2  )  )+1e-15;
    double d8  = sqrt(  pow(  1-diff_T_p  , 2  ) + pow(  1-diff_eta  , 2  ) + pow(  1-diff_T_pbar  , 2  )  )+1e-15;
    
    
    
    
    double sum  = r1*ttx(d1) + r2*ttx(d2) + r3*ttx(d3) + r4*ttx(d4) + r5*ttx(d5) + r6*ttx(d6) + r7*ttx(d7) + r8*ttx(d8);
    double norm = ttx(d1) + ttx(d2) + ttx(d3) + ttx(d4) + ttx(d5) + ttx(d6) + ttx(d7) + ttx(d8);
    
    
    return sum/norm;
    
}

double Interpolate_LAB_2(double T_p_LAB, double T_pbar_LAB, double eta_LAB, double* array){
    double i_T_p_LAB        =  fN_Tn_p   * log(T_p_LAB   /fMin_Tn_p  )/log(fMax_Tn_p  /fMin_Tn_p  );
    double i_eta_LAB        =  fN_eta    * (eta_LAB-fMin_eta)/(fMax_eta-fMin_eta);
    double i_T_pbar_LAB     =  fN_T_pbar * log(T_pbar_LAB/fMin_T_pbar)/log(fMax_T_pbar/fMin_T_pbar);
    
    
    if (i_T_p_LAB   <0 || i_T_p_LAB     >=fN_Tn_p-1   ) return 5;
    if (i_T_pbar_LAB<0 || i_T_pbar_LAB  >=fN_T_pbar-1 ) return 10;
    if (i_eta_LAB   <0 || i_eta_LAB     >=fN_eta-1    ) return 20;
    
    
    int int_T_p_LAB        = i_T_p_LAB;
    int int_eta_LAB        = i_eta_LAB;
    int int_T_pbar_LAB     = i_T_pbar_LAB;
    
    varOut(T_p_LAB)
    varOut(T_pbar_LAB)
    varOut(eta_LAB)
    
    varOut(fTn_p_LAB[int_T_p_LAB])
    varOut(fT_pbar_LAB[int_T_pbar_LAB])
    varOut(fEta_LAB[int_eta_LAB])
    
    
    
    double r1  = array[ (int_T_p_LAB+0) + fN_Tn_p*(int_eta_LAB+0) + fN_Tn_p*fN_eta*(int_T_pbar_LAB+0) ];
    double r2  = array[ (int_T_p_LAB+1) + fN_Tn_p*(int_eta_LAB+0) + fN_Tn_p*fN_eta*(int_T_pbar_LAB+0) ];
    double r3  = array[ (int_T_p_LAB+0) + fN_Tn_p*(int_eta_LAB+1) + fN_Tn_p*fN_eta*(int_T_pbar_LAB+0) ];
    double r4  = array[ (int_T_p_LAB+0) + fN_Tn_p*(int_eta_LAB+0) + fN_Tn_p*fN_eta*(int_T_pbar_LAB+1) ];
    double r5  = array[ (int_T_p_LAB+1) + fN_Tn_p*(int_eta_LAB+1) + fN_Tn_p*fN_eta*(int_T_pbar_LAB+0) ];
    double r6  = array[ (int_T_p_LAB+1) + fN_Tn_p*(int_eta_LAB+0) + fN_Tn_p*fN_eta*(int_T_pbar_LAB+1) ];
    double r7  = array[ (int_T_p_LAB+0) + fN_Tn_p*(int_eta_LAB+1) + fN_Tn_p*fN_eta*(int_T_pbar_LAB+1) ];
    double r8  = array[ (int_T_p_LAB+1) + fN_Tn_p*(int_eta_LAB+1) + fN_Tn_p*fN_eta*(int_T_pbar_LAB+1) ];
    
    
    double diff_T_p     =  i_T_p_LAB    -   int_T_p_LAB;
    double diff_eta     =  i_eta_LAB    -   int_eta_LAB;
    double diff_T_pbar  =  i_T_pbar_LAB -   int_T_pbar_LAB;
    
    double d1  = sqrt(  pow(    diff_T_p  , 2  ) + pow(    diff_eta  , 2  ) + pow(    diff_T_pbar  , 2  )  )+1e-15;
    double d2  = sqrt(  pow(  1-diff_T_p  , 2  ) + pow(    diff_eta  , 2  ) + pow(    diff_T_pbar  , 2  )  )+1e-15;
    double d3  = sqrt(  pow(    diff_T_p  , 2  ) + pow(  1-diff_eta  , 2  ) + pow(    diff_T_pbar  , 2  )  )+1e-15;
    double d4  = sqrt(  pow(    diff_T_p  , 2  ) + pow(    diff_eta  , 2  ) + pow(  1-diff_T_pbar  , 2  )  )+1e-15;
    double d5  = sqrt(  pow(  1-diff_T_p  , 2  ) + pow(  1-diff_eta  , 2  ) + pow(    diff_T_pbar  , 2  )  )+1e-15;
    double d6  = sqrt(  pow(  1-diff_T_p  , 2  ) + pow(    diff_eta  , 2  ) + pow(  1-diff_T_pbar  , 2  )  )+1e-15;
    double d7  = sqrt(  pow(    diff_T_p  , 2  ) + pow(  1-diff_eta  , 2  ) + pow(  1-diff_T_pbar  , 2  )  )+1e-15;
    double d8  = sqrt(  pow(  1-diff_T_p  , 2  ) + pow(  1-diff_eta  , 2  ) + pow(  1-diff_T_pbar  , 2  )  )+1e-15;
    
    double r = r1;
    double d = d1;
    if(d2<d){r=r2;d=d2;}
    if(d3<d){r=r3;d=d3;}
    if(d4<d){r=r4;d=d4;}
    if(d5<d){r=r5;d=d5;}
    if(d6<d){r=r6;d=d6;}
    if(d7<d){r=r7;d=d7;}
    if(d8<d){r=r8;d=d8;}
    
    return r;
    
}

double GetInCM( double s, double xR, double pT, double* array, int option ){
    
    double T_p_LAB;
    double T_pbar_LAB;
    double eta_LAB;
    
    double cc1 = 2;
    double cc2 = 2;
    
    if( CSTransformations::convert_CM_to_LAB(s, xR, pT, T_p_LAB, T_pbar_LAB, eta_LAB, true) )
        cc1 = Interpolate_LAB(T_p_LAB, T_pbar_LAB, eta_LAB, array);
    if( CSTransformations::convert_CM_to_LAB(s, xR, pT, T_p_LAB, T_pbar_LAB, eta_LAB, false) )
        cc2 = Interpolate_LAB(T_p_LAB, T_pbar_LAB, eta_LAB, array);
    
    if (option==1)
        return cc1;
    if (option==2)
        return cc2;
    
    return std::min(cc1,cc2);
};


double GetInCM_xF( double s, double xF, double pT, double* array ){
    
    double T_p_LAB;
    double T_pbar_LAB;
    double eta_LAB;
    
    
    if(CSTransformations::convert_CMxF_to_LAB(s, xF, pT, T_p_LAB, T_pbar_LAB, eta_LAB)){
        return Interpolate_LAB(T_p_LAB, T_pbar_LAB, eta_LAB, array);
    }
    
    return 2;
};


