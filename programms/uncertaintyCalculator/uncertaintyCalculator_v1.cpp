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



double dT_pp_pbar_LAB   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::DI_MAURO12);
}
double dT_Hep_pbar_LAB  ( double Tn_He,       double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(Tn_He, T_pbar, CS::PHe_PBAR, CS::DI_MAURO12);
}
double dT_pHe_pbar_LAB   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::HeP_PBAR, CS::DI_MAURO12);
}
double dT_HeHe_pbar_LAB   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::HeP_PBAR, CS::DI_MAURO12)*pow(4, 0.7);
}



// inv cross sections
double inv_pp_pbar_LAB   ( double T_p,    double T_pbar,  double eta   ){
    double cosTheta = tanh(eta);
    return CSTransformations::inv_pp_product_LAB(  T_p+fMass_proton, T_pbar+fMass_proton, cosTheta, fMass_proton, ppCSParametrizations::inv_pp_pbar_CM__diMauro12  )*1e-31;
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
const int    fN_T_pbar      =   401;
const int    fN_Tn_p        =   401;
const int    fN_eta         =   401;

const double fMin_T_pbar    = 1e-1;
const double fMax_T_pbar    = 1e3;

const double fMin_Tn_p      = 1e0;
const double fMax_Tn_p      = 1e5;

const double fMin_eta       =  0;
const double fMax_eta       = 13;

// grid CM
const int    fN_s           = 401;
const int    fN_xR          = 401;
const int    fN_pT          = 401;

const double fMin_s         = 15;
const double fMax_s         = 2e4;

const double fMin_pT        = 1e-2;
const double fMax_pT        = 1e2;

const double fMin_xR        =  0;
const double fMax_xR        =  1;


int fMod = 1;


// ISM density
const double f_nH        = 1e-6;
const double f_nHe       = f_nH * 0.1;


// Lab
double fT_pbar_LAB  [fN_T_pbar];
double fTn_p_LAB    [fN_Tn_p  ];
double fEta_LAB     [fN_eta   ];

// CM
double fS     [fN_s ];
double fxR    [fN_xR];
double fPT    [fN_pT];

double fRU_Tpbar__Q  [fN_T_pbar];

// Source term as function of T_pbar
double fQ       [fN_T_pbar];
double fQ_pp    [fN_T_pbar];
double fQ_pHe   [fN_T_pbar];
double fQ_Hep   [fN_T_pbar];
double fQ_HeHe  [fN_T_pbar];


double fC_Tpbar__pp  [fN_T_pbar];
double fC_Tpbar__pHe [fN_T_pbar];
double fC_Tpbar__Hep [fN_T_pbar];
double fC_Tpbar__HeHe[fN_T_pbar];


double fRU_Tpbar__pp   [fN_T_pbar];
double fRU_Tpbar__pHe  [fN_T_pbar];
double fRU_Tpbar__Hep  [fN_T_pbar];

double fC_Tpbar_Tp__pp  [fN_Tn_p * fN_T_pbar];
double fC_Tpbar_Tp__pHe [fN_Tn_p * fN_T_pbar];
double fC_Tpbar_Tp__Hep [fN_Tn_p * fN_T_pbar];

double fCC_Tpbar_Tp__pp  [fN_Tn_p * fN_T_pbar];
double fCC_Tpbar_Tp__pHe [fN_Tn_p * fN_T_pbar];
double fCC_Tpbar_Tp__Hep [fN_Tn_p * fN_T_pbar];


double fC_Tpbar_Tp_eta_LAB__pp   [fN_T_pbar * fN_Tn_p * fN_eta];
double fCC_Tpbar_Tp_eta_LAB__pp  [fN_T_pbar * fN_Tn_p * fN_eta];

double f2D__Tpbar_eta            [fN_T_pbar * fN_eta];

double fCC_s_xR_pT__pp           [fN_s * fN_xR * fN_pT ];
//double fCC_s_xR_pT__pp_pos       [fN_s * fN_xR * fN_pT ];
//double fCC_s_xR_pT__pp_neg       [fN_s * fN_xR * fN_pT ];

double f2D__xR_pT        [fN_xR * fN_pT];

double fHelper[ fN_T_pbar * fN_Tn_p * fN_eta ];


void write_array_to_file(std::string filename, int n, double* array, double* xarray=NULL, std::string description=""){
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
    //std::cout << table << std::endl;
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
// sort_quick implements the QuickSort algorithm. It sorts the array between position n1 and n2.
// the int array pos retains the original positions
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

// this function uses the quick sort algorithm implemented above to commulate an array
void commulate(int n, double* array, double* com, bool print=false){
    if (n>fN_eta*fN_Tn_p*fN_T_pbar) {
        std::cout << "Error: please use larger helper array!" << std::endl;
        return;
    }
    int position[n];
    for (int i=0; i<n; i++) {
        position[i]=i;
        fHelper [i]=array[i];
    }
    sort_quick(fHelper, 0, n-1, position, print);
    for (int i=n-2; i>=0; i--) {
        fHelper[i] = fHelper[i]+fHelper[i+1];
    }
    for (int i=0; i<n; i++) {
        com[position[i]]=fHelper[i]/fHelper[0];
    }
}



/// tranformation from LAB to CC


double ttx(double x){
    x = 1-x;
    if (x<0) {
        return 0;
    }
    return x;
}

double InterpolateCC_LAB(double T_p_LAB, double T_pbar_LAB, double eta_LAB, double* array=fCC_Tpbar_Tp_eta_LAB__pp){
    //out(eta_LAB)
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

double GetLab( double s, double xR, double pT, int option=0, double* array=fCC_Tpbar_Tp_eta_LAB__pp ){
    
    double T_p_LAB;
    double T_pbar_LAB;
    double eta_LAB;
    
    double cc1 = 2;
    double cc2 = 2;
    
    if( CSTransformations::convert_CM_to_LAB(s, xR, pT, T_p_LAB, T_pbar_LAB, eta_LAB, true) )
        cc1 = InterpolateCC_LAB(T_p_LAB, T_pbar_LAB, eta_LAB, array);
    if( CSTransformations::convert_CM_to_LAB(s, xR, pT, T_p_LAB, T_pbar_LAB, eta_LAB, false) )
        cc2 = InterpolateCC_LAB(T_p_LAB, T_pbar_LAB, eta_LAB, array);

    if (option==1)
        return cc1;
    if (option==2)
        return cc2;
    
    return std::min(cc1,cc2);
};


int main(int argc, char *argv[])
{
    ConfigHandler* config = ConfigHandler::GetInstance();
    
    
    config->SetInput(argc, argv, fProgramDescription);
    int     step    = 0;
    config->AddOptionInt    (   "step",         step,       "Step. Default: 0"                      );
    config->CheckOption();
    
    // Setup energy grid
    double dTn_p=log10(fMax_Tn_p/fMin_Tn_p)/(fN_Tn_p-1);
    int i=0;
    for ( double dT=log10(fMin_Tn_p); dT<log10(fMax_Tn_p)+dTn_p/2; dT+=dTn_p ) {
        fTn_p_LAB[i] = pow(10, dT);
        i++;
    }
    
    double dT_pbar=log10(fMax_T_pbar/fMin_T_pbar)/(fN_T_pbar-1);
    i = 0;
    for ( double dT=log10(fMin_T_pbar); dT<log10(fMax_T_pbar)+dT_pbar/2; dT+=dT_pbar ) {
        fT_pbar_LAB[i] = pow(10, dT);
        i++;
    }
    
    double dEta=(fMax_eta-fMin_eta)/(fN_eta-1);
    i = 0;
    for ( double eta=fMin_eta; eta<fMax_eta+dEta/2; eta+=dEta ) {
        fEta_LAB[i] = eta;
        i++;
    }

    
    
    
    double dS=log10(fMax_s/fMin_s)/(fN_s-1);
    i=0;
    for ( double d=log10(fMin_s); d<log10(fMax_s)+dS/2; d+=dS ) {
        fS[i] = pow(10, d);
        i++;
    }
//    double dXR=log10(fMax_xR/fMin_xR)/(fN_xR-1);
//    i=0;
//    for ( double d=log10(fMin_xR); d<log10(fMax_xR)+dXR/2; d+=dXR ) {
//        fxR[i] = pow(10, d);
//        i++;
//    }
    double dXR=(fMax_xR-fMin_xR)/(fN_xR-1);
    i=0;
    for ( double xR=fMin_xR; xR<fMax_xR+dXR/2; xR+=dXR ) {
        fxR[i] = xR;
        i++;
    }
    double dpT=log10(fMax_pT/fMin_pT)/(fN_pT-1);
    i=0;
    for ( double d=log10(fMin_pT); d<log10(fMax_pT)+dpT/2; d+=dpT ) {
        fPT[i] = pow(10, d);
        i++;
    }
    

    
    
    
    write_array_to_file("save/save.fT_pbar.txt", fN_T_pbar, fT_pbar_LAB );
    write_array_to_file("save/save.fTn_p.txt",   fN_Tn_p,   fTn_p_LAB   );
    write_array_to_file("save/save.fEta.txt",    fN_eta,    fEta_LAB    );
    
    
    write_array_to_file("save/save.fS.txt",     fN_s,       fS      );
    write_array_to_file("save/save.fPT.txt",    fN_pT,      fPT     );
    write_array_to_file("save/save.fxR.txt",    fN_xR,      fxR     );
    
    
    
    
    
    
    LISFluxes::fit();
    ReadFluxes reader;
    
    std::system("mkdir save");
    
    
    // Determine and parametrize pbar uncertainty
    if (step==1) {
        
        std::system("mkdir pbar_uncertainty");
        
        Graph pbar = reader.GetSpectrum("pbar", "ams", EKINPERN);
        double pbarT                  [pbar.Size()];
        double log_pbarT              [pbar.Size()];
        double pbarRelativeUncertainty[pbar.Size()];
        for (int i = 0; i<pbar.Size(); i++) {
            double x, xe, y, ye;
            pbar.GetPoint(i, x, y, xe, ye); x+=0.5;
            pbarT                  [i] = x;
            log_pbarT              [i] = log(x);
            pbarRelativeUncertainty[i] = ye/y;
        }
        
        write_array_to_file("pbar_uncertainty/pbar_relative_uncertainty_data.txt",               pbar.Size(),    &pbarRelativeUncertainty[0],    &pbarT[0],          "T_pbar      rel. uncertainty");
        TGraph pbarR(pbar.Size(), &log_pbarT[0], &pbarRelativeUncertainty[0]);
        TF1 fun_pbarR("fun_pbarR", "pol8");
        pbarR.Fit(&fun_pbarR);
        
        for ( int i=0; i<fN_T_pbar; i++ ) {
            fRU_Tpbar__Q[i] = fun_pbarR.Eval( log(fT_pbar_LAB[i]) );
        }
        write_array_to_file("pbar_uncertainty/pbar_relative_uncertainty_parmetrization.txt",     fN_T_pbar,      &fRU_Tpbar__Q[0],               &fT_pbar_LAB[0],    "T_pbar      rel. uncertainty" );
        
        
        // save calculation to files
        write_array_to_file("save/save.pbar_relative_uncertainty_parmetrization.txt",     fN_T_pbar,      fRU_Tpbar__Q);
        
        
    }

    // Parametrize pbar uncertainty expected for GAPS
    if (step==11) {
        
        std::system("mkdir pbar_uncertainty");
        
       
        for ( int i=0; i<fN_T_pbar; i++ ) {
            double x = log10(fT_pbar_LAB[i]);
            fRU_Tpbar__Q[i] = pow( 2*x +0.5, 8)+0.05;
        }
        write_array_to_file("pbar_uncertainty/pbar_relative_uncertainty_parmetrization.txt",     fN_T_pbar,      &fRU_Tpbar__Q[0],               &fT_pbar_LAB[0],    "T_pbar      rel. uncertainty" );
        
        
        // save calculation to files
        write_array_to_file("save/save.pbar_relative_uncertainty_parmetrization.txt",     fN_T_pbar,      fRU_Tpbar__Q);
        
        
    }

    
    // Determine source terms
    if (step==2) {
        
        /////// read calculations from previous steps
        read_array_from_file("save/save.pbar_relative_uncertainty_parmetrization.txt",     fN_T_pbar,      fRU_Tpbar__Q);
        ///////
        
        std::system("mkdir source_terms");
        for ( int i=0; i<fN_T_pbar; i++ ) {
            double T_pbar   = fT_pbar_LAB[i];
            double Tn_from  = std::min(T_pbar, 6*fMass_proton);
            
            fQ_pp   [i] =  SourceTerms::integrationFromTo( T_pbar, dT_pp_pbar_LAB,    protonFlux, f_nH,  Tn_from, fMax_Tn_p, 10000 );
            fQ_pHe  [i] =  SourceTerms::integrationFromTo( T_pbar, dT_pHe_pbar_LAB,   protonFlux, f_nHe, Tn_from, fMax_Tn_p, 10000 );
            fQ_Hep  [i] =  SourceTerms::integrationFromTo( T_pbar, dT_Hep_pbar_LAB,   heliumFlux, f_nH,  Tn_from, fMax_Tn_p, 10000 );
            fQ_HeHe [i] =  SourceTerms::integrationFromTo( T_pbar, dT_HeHe_pbar_LAB,  heliumFlux, f_nHe, Tn_from, fMax_Tn_p, 10000 );
            
            fQ      [i] =  fQ_pp[i]+fQ_Hep[i]+fQ_pHe[i]+fQ_HeHe[i];
        }
        
        write_array_to_file("source_terms/source_term_tot.txt",  fN_T_pbar,   &fQ    [0],    &fT_pbar_LAB[0],    "T_pbar      source term tot" );
        write_array_to_file("source_terms/source_term_pp.txt",   fN_T_pbar,   &fQ_pp [0],    &fT_pbar_LAB[0],    "T_pbar      source term pp " );
        write_array_to_file("source_terms/source_term_pHe.txt",  fN_T_pbar,   &fQ_Hep[0],    &fT_pbar_LAB[0],    "T_pbar      source term Hep" );
        write_array_to_file("source_terms/source_term_Hep.txt",  fN_T_pbar,   &fQ_pHe[0],    &fT_pbar_LAB[0],    "T_pbar      source term pHe" );
        write_array_to_file("source_terms/source_term_HeHe.txt", fN_T_pbar,   &fQ_HeHe[0],   &fT_pbar_LAB[0],    "T_pbar      source term HeHe");
        
        ////// save calculations for next step
        write_array_to_file("save/save.source_term_tot.txt",  fN_T_pbar,   fQ      );
        write_array_to_file("save/save.source_term_pp.txt",   fN_T_pbar,   fQ_pp   );
        write_array_to_file("save/save.source_term_pHe.txt",  fN_T_pbar,   fQ_Hep  );
        write_array_to_file("save/save.source_term_Hep.txt",  fN_T_pbar,   fQ_pHe  );
        write_array_to_file("save/save.source_term_HeHe.txt", fN_T_pbar,   fQ_HeHe );
        //////
        
        
        
        // Determine uncertainty of on pp source term
        for ( int i=0; i<fN_T_pbar; i++ ) {
            fRU_Tpbar__pp   [i] =  fQ [i]/fQ_pp [i]  *  fRU_Tpbar__Q [i];
        }
        write_array_to_file("source_terms/relative_uncertainty__T_pbar__pp.txt",  fN_T_pbar,   &fRU_Tpbar__pp[0],    &fT_pbar_LAB[0],    "T_pbar      relative uncertainty Q_pp" );
        
        
        
    }
    
    
    if (step==3) {
        
        /////// read calculations from previous steps
        read_array_from_file("save/save.pbar_relative_uncertainty_parmetrization.txt",     fN_T_pbar,      fRU_Tpbar__Q);
        read_array_from_file("save/save.source_term_tot.txt",  fN_T_pbar,   fQ      );
        read_array_from_file("save/save.source_term_pp.txt",   fN_T_pbar,   fQ_pp   );
        read_array_from_file("save/save.source_term_pHe.txt",  fN_T_pbar,   fQ_Hep  );
        read_array_from_file("save/save.source_term_Hep.txt",  fN_T_pbar,   fQ_pHe  );
        read_array_from_file("save/save.source_term_HeHe.txt", fN_T_pbar,   fQ_HeHe );
        ///////


        std::system("mkdir 2D_contribution_LAB ");
        
        std::system("mkdir 2D_contribution_LAB/pp  ");
        std::system("mkdir 2D_contribution_LAB/pHe ");
        std::system("mkdir 2D_contribution_LAB/Hep ");
        
        
        // contribution and commulative contribution in T_pbar T_p 2D
        for (int j=0; j<fN_T_pbar; j++) {
            double T_pbar = fT_pbar_LAB[j];
            //double dLog_T_p = log(fMax_Tn_p/fMin_Tn_p)/fN_Tn_p;
            for (int i=0; i<fN_Tn_p-1; i++) {
                double Tn_p = fTn_p_LAB[i+1];
                fC_Tpbar_Tp__pp [i + fN_Tn_p*j] = 4*M_PI *f_nH  * Tn_p * dT_pp_pbar_LAB (Tn_p, T_pbar)*protonFlux(Tn_p)/fQ_pp[j];
                fC_Tpbar_Tp__pHe[i + fN_Tn_p*j] = 4*M_PI *f_nHe * Tn_p * dT_pHe_pbar_LAB(Tn_p, T_pbar)*protonFlux(Tn_p)/fQ_pHe[j];
                fC_Tpbar_Tp__Hep[i + fN_Tn_p*j] = 4*M_PI *f_nHe * Tn_p * dT_pHe_pbar_LAB(Tn_p, T_pbar)*protonFlux(Tn_p)/fQ_Hep[j];
            }
            commulate(fN_Tn_p, &fC_Tpbar_Tp__pp [0 + fN_Tn_p*j], &fCC_Tpbar_Tp__pp [0 + fN_Tn_p*j]);
            commulate(fN_Tn_p, &fC_Tpbar_Tp__pHe[0 + fN_Tn_p*j], &fCC_Tpbar_Tp__pHe[0 + fN_Tn_p*j]);
            commulate(fN_Tn_p, &fC_Tpbar_Tp__Hep[0 + fN_Tn_p*j], &fCC_Tpbar_Tp__Hep[0 + fN_Tn_p*j]);
        }
        
        
        write_2D_array_to_file( "2D_contribution_LAB/pp/2D_contribution__pp.txt",       fN_T_pbar,  fN_Tn_p,   &fC_Tpbar_Tp__pp[0],    fT_pbar_LAB, fTn_p_LAB  );
        write_2D_array_to_file( "2D_contribution_LAB/pp/2D_ccontribution__pp.txt",      fN_T_pbar,  fN_Tn_p,   &fCC_Tpbar_Tp__pp[0],   fT_pbar_LAB, fTn_p_LAB  );
        
        write_2D_array_to_file( "2D_contribution_LAB/pHe/2D_contribution__pHe.txt",     fN_T_pbar,  fN_Tn_p,   &fC_Tpbar_Tp__pHe[0],   fT_pbar_LAB, fTn_p_LAB  );
        write_2D_array_to_file( "2D_contribution_LAB/pHe/2D_ccontribution__pHe.txt",    fN_T_pbar,  fN_Tn_p,   &fCC_Tpbar_Tp__pHe[0],  fT_pbar_LAB, fTn_p_LAB  );
        
        write_2D_array_to_file( "2D_contribution_LAB/Hep/2D_contribution__Hep.txt",     fN_T_pbar,  fN_Tn_p,   &fC_Tpbar_Tp__Hep[0],   fT_pbar_LAB, fTn_p_LAB  );
        write_2D_array_to_file( "2D_contribution_LAB/Hep/2D_ccontribution__Hep.txt",    fN_T_pbar,  fN_Tn_p,   &fCC_Tpbar_Tp__Hep[0],  fT_pbar_LAB, fTn_p_LAB  );
        
    }
    // Unfolding in Tp and eta
    
    
    // calculate 3D contribtuion (LAB)
    if (step==4) {
        
        /////// read calculations from previous steps
        read_array_from_file("save/save.pbar_relative_uncertainty_parmetrization.txt",     fN_T_pbar,      fRU_Tpbar__Q);
        read_array_from_file("save/save.source_term_tot.txt",  fN_T_pbar,   fQ      );
        read_array_from_file("save/save.source_term_pp.txt",   fN_T_pbar,   fQ_pp   );
        read_array_from_file("save/save.source_term_pHe.txt",  fN_T_pbar,   fQ_Hep  );
        read_array_from_file("save/save.source_term_Hep.txt",  fN_T_pbar,   fQ_pHe  );
        read_array_from_file("save/save.source_term_HeHe.txt", fN_T_pbar,   fQ_HeHe );
        ///////
        
        std::system("mkdir 3D_Tp_Tpbar_eta_LAB");
        std::system("mkdir 3D_Tp_Tpbar_eta_LAB/pp");
        
        for (int j=0; j<fN_T_pbar; j++) {
            
            double T_pbar = fT_pbar_LAB[j];
            double dLog_T_p = log(fMax_Tn_p/fMin_Tn_p)/fN_Tn_p;
            double sum = 0;
            for (int i=0; i<fN_Tn_p-1; i++) {
                double Tn_p = fTn_p_LAB[i];
                for (int k=0; k<fN_eta-1; k++) {
                    double eta = fEta_LAB[k];
                    double p_pbar = sqrt(T_pbar*T_pbar+2*fMass_proton*T_pbar);
                    fC_Tpbar_Tp_eta_LAB__pp      [i + fN_Tn_p*k + fN_Tn_p*fN_eta*j] = 8*M_PI*M_PI *f_nH * p_pbar * Tn_p *protonFlux(Tn_p) / pow(cosh(eta), 2) * inv_pp_pbar_LAB(Tn_p, T_pbar, eta);// * dLog_T_p  * dEta
                    sum+=fC_Tpbar_Tp_eta_LAB__pp [i + fN_Tn_p*k + fN_Tn_p*fN_eta*j];
                }
            }
            
            for (int i=0; i<fN_Tn_p-1; i++) {
                for (int k=0; k<fN_eta-1; k++) {
                    fC_Tpbar_Tp_eta_LAB__pp [i + fN_Tn_p*k + fN_Tn_p*fN_eta*j]/=sum;
                }
            }
            
            std::cout << std::endl;
            out(fT_pbar_LAB[j])
            out(sum*dEta*dLog_T_p)
            out(fQ_pp[j])
            
            commulate(fN_Tn_p*fN_eta, &fC_Tpbar_Tp_eta_LAB__pp [j*fN_eta*fN_Tn_p], &fCC_Tpbar_Tp_eta_LAB__pp [j*fN_eta*fN_Tn_p] );
            
            if (j%fMod==0){
                write_2D_array_to_file( QString("3D_Tp_Tpbar_eta_LAB/pp/2D_contribution_Tp_eta__%1__pp.txt").arg(j).toStdString(),      fN_eta,  fN_Tn_p,   &fC_Tpbar_Tp_eta_LAB__pp[j*fN_eta*fN_Tn_p],    fEta_LAB, fTn_p_LAB  );
                write_2D_array_to_file( QString("3D_Tp_Tpbar_eta_LAB/pp/2D_ccontribution_Tp_eta__%1__pp.txt").arg(j).toStdString(),     fN_eta,  fN_Tn_p,   &fCC_Tpbar_Tp_eta_LAB__pp[j*fN_eta*fN_Tn_p],   fEta_LAB, fTn_p_LAB  );
            }
        }
        
        // save calculation to files
        write_array_to_file("save/save.fC_Tpbar_Tp_eta_LAB__pp.txt",  fN_T_pbar*fN_Tn_p*fN_eta, fC_Tpbar_Tp_eta_LAB__pp );
        write_array_to_file("save/save.fCC_Tpbar_Tp_eta_LAB__pp.txt", fN_T_pbar*fN_Tn_p*fN_eta, fCC_Tpbar_Tp_eta_LAB__pp);
        
    }

    // calculate 3D contribtuion (CM)
    if (step==5) {
        
        /////// read calculations from previous steps
        read_array_from_file("save/save.pbar_relative_uncertainty_parmetrization.txt",     fN_T_pbar,      fRU_Tpbar__Q);
        read_array_from_file("save/save.source_term_tot.txt",  fN_T_pbar,   fQ      );
        read_array_from_file("save/save.source_term_pp.txt",   fN_T_pbar,   fQ_pp   );
        read_array_from_file("save/save.source_term_pHe.txt",  fN_T_pbar,   fQ_Hep  );
        read_array_from_file("save/save.source_term_Hep.txt",  fN_T_pbar,   fQ_pHe  );
        read_array_from_file("save/save.source_term_HeHe.txt", fN_T_pbar,   fQ_HeHe );
        read_array_from_file("save/save.fCC_Tpbar_Tp_eta_LAB__pp.txt", fN_T_pbar*fN_Tn_p*fN_eta, fCC_Tpbar_Tp_eta_LAB__pp);
        ///////
        
        double Tp    = 50;
        double Tpbar = 5;
        double eta   = 3;
        
        out(InterpolateCC_LAB(Tp, Tpbar, eta));
        
        double s;
        double xR;
        double pT;
        
        out(Tp)
        out(Tpbar)
        out(eta)
        
        out("to CM")
        CSTransformations::convert_LAB_to_CM(Tp, Tpbar, eta, s, xR, pT);
        out(s)
        out(xR)
        out(pT)
        
        out(GetLab(s, xR, pT))
        
        out("back to LAB pL>0")
        CSTransformations::convert_CM_to_LAB(s, xR, pT, Tp, Tpbar, eta, true);
        out(Tp)
        out(Tpbar)
        out(eta)
        
        
        out("back to LAB pL<0")
        CSTransformations::convert_CM_to_LAB(s, xR, pT, Tp, Tpbar, eta, false);
        out(Tp)
        out(Tpbar)
        out(eta)
        
       
       
        std::system("mkdir 3D_s_xR_pT");
        std::system("mkdir 3D_s_xR_pT/pp");
        
        
        for (int j=0; j<fN_s; j++) {
            double s = fS[j];
            for (int i=0; i<fN_xR-1; i++) {
                double xR = fxR[i+1];
                for (int k=0; k<fN_pT-1; k++) {
                    double pT = fPT[k];
                    fCC_s_xR_pT__pp       [i + fN_pT*k + fN_xR*fN_pT*j] = GetLab( s, xR, pT,  1 );
                }
            }
            std::cout << std::endl;
            out(s)
            if (j%fMod==0){
                write_2D_array_to_file( QString("3D_s_xR_pT/pp/2D_ccontribution_xR_pT__%1__pp_pos.txt").arg(j).toStdString(),     fN_pT,  fN_xR,   &fCC_s_xR_pT__pp    [j*fN_pT*fN_xR],   fPT, fxR  );
            }
        }
        for (int j=0; j<fN_s; j++) {
            double s = fS[j];
            for (int i=0; i<fN_xR-1; i++) {
                double xR = fxR[i+1];
                for (int k=0; k<fN_pT-1; k++) {
                    double pT = fPT[k];
                    fCC_s_xR_pT__pp       [i + fN_pT*k + fN_xR*fN_pT*j] = GetLab( s, xR, pT,  2 );
                }
            }
            std::cout << std::endl;
            out(s)
            if (j%fMod==0){
                write_2D_array_to_file( QString("3D_s_xR_pT/pp/2D_ccontribution_xR_pT__%1__pp_neg.txt").arg(j).toStdString(),     fN_pT,  fN_xR,   &fCC_s_xR_pT__pp    [j*fN_pT*fN_xR],   fPT, fxR  );
            }
        }
        for (int j=0; j<fN_s; j++) {
            double s = fS[j];
            for (int i=0; i<fN_xR-1; i++) {
                double xR = fxR[i+1];
                for (int k=0; k<fN_pT-1; k++) {
                    double pT = fPT[k];
                    fCC_s_xR_pT__pp       [i + fN_pT*k + fN_xR*fN_pT*j] = GetLab( s, xR, pT    );
                }
            }
            std::cout << std::endl;
            out(s)
            if (j%fMod==0){
                write_2D_array_to_file( QString("3D_s_xR_pT/pp/2D_ccontribution_xR_pT__%1__pp.txt").    arg(j).toStdString(),     fN_pT,  fN_xR,   &fCC_s_xR_pT__pp    [j*fN_pT*fN_xR],   fPT, fxR  );
            }
        }
        
        // save calculation to files
        write_array_to_file("save/save.fCC_s_xR_pT__pp.txt", fN_s * fN_xR * fN_pT, fCC_s_xR_pT__pp);
        
    }
    
    
    // calculate 3D relative uncertatinty
    if (step==6) {
        
        /////// read calculations from previous steps
        read_array_from_file("save/save.pbar_relative_uncertainty_parmetrization.txt",     fN_T_pbar,      fRU_Tpbar__Q);
        read_array_from_file("save/save.source_term_tot.txt",  fN_T_pbar,   fQ      );
        read_array_from_file("save/save.source_term_pp.txt",   fN_T_pbar,   fQ_pp   );
        read_array_from_file("save/save.source_term_pHe.txt",  fN_T_pbar,   fQ_Hep  );
        read_array_from_file("save/save.source_term_Hep.txt",  fN_T_pbar,   fQ_pHe  );
        read_array_from_file("save/save.source_term_HeHe.txt", fN_T_pbar,   fQ_HeHe );
        read_array_from_file("save/save.fCC_Tpbar_Tp_eta_LAB__pp.txt", fN_T_pbar*fN_Tn_p*fN_eta, fCC_Tpbar_Tp_eta_LAB__pp);
        ///////
    
        
        // assume uncertainty as step function: 2% and 100%
        
        for (int j=0; j<fN_T_pbar; j++) {
            
            out(fT_pbar_LAB[j])
            
            double ru_AMS     = fRU_Tpbar__Q[j]/fQ_pp[j]*fQ[j];
            double threshold  = (1-ru_AMS)/0.98;
            
            for (int i=0; i<fN_Tn_p-1; i++) {
                for (int k=0; k<fN_eta-1; k++) {
                    if (  fCC_Tpbar_Tp_eta_LAB__pp[i + fN_Tn_p*k + fN_Tn_p*fN_eta*j] < threshold  ){
                        fC_Tpbar_Tp_eta_LAB__pp[i + fN_Tn_p*k + fN_Tn_p*fN_eta*j] = 0.02;
                    }else{
                        fC_Tpbar_Tp_eta_LAB__pp[i + fN_Tn_p*k + fN_Tn_p*fN_eta*j] = 1.00;
                    }
                }
            }
            
            if (j%fMod==0)
                write_2D_array_to_file( QString("3D_Tp_Tpbar_eta_LAB/pp/2D_relativeUncertainty_Tp_eta__%1__pp.txt").arg(j).toStdString(),      fN_eta,  fN_Tn_p,   &fC_Tpbar_Tp_eta_LAB__pp[j*fN_eta*fN_Tn_p],    fEta_LAB, fTn_p_LAB  );
            
        }
        
        write_array_to_file("save/save.fC_uncertainty_Tpbar_Tp_eta_LAB__pp.txt",   fN_T_pbar*fN_Tn_p*fN_eta,  fC_Tpbar_Tp_eta_LAB__pp );
       
        
        for (int j=0; j<fN_s; j++) {
            double s = fS[j];
            
            for (int i=0; i<fN_xR-1; i++) {
                double xR = fxR[i+1];
                
                for (int k=0; k<fN_pT-1; k++) {
                    double pT = fPT[k];
                    
                    double c1 = GetLab( s, xR, pT, 1,  fC_Tpbar_Tp_eta_LAB__pp    );
                    double c2 = GetLab( s, xR, pT, 2,  fC_Tpbar_Tp_eta_LAB__pp    );
                    
                    double c  = GetLab( s, xR, pT, 0,  fC_Tpbar_Tp_eta_LAB__pp    );
                    
                    if (c1<c || c2<c) {
                        out(c1)
                        out(c2)
                        out(c)
                    }
                    
                    fCC_s_xR_pT__pp       [i + fN_pT*k + fN_xR*fN_pT*j] = GetLab( s, xR, pT, 0,  fC_Tpbar_Tp_eta_LAB__pp    );
                }
            }
            
            std::cout << std::endl;
            out(s)
            
            
            if (j%fMod==0)
                write_2D_array_to_file( QString("3D_s_xR_pT/pp/2D_relativeUncertainty_xR_pT__%1__pp.txt").    arg(j).toStdString(),     fN_pT,  fN_xR,   &fCC_s_xR_pT__pp    [j*fN_pT*fN_xR],   fPT, fxR  );
            
        }
        
    }
    
    // calculate 3D relative uncertatinty for data comparison
    if (step==7) {
        
        /////// read calculations from previous steps
        read_array_from_file("save/save.pbar_relative_uncertainty_parmetrization.txt",     fN_T_pbar,      fRU_Tpbar__Q);
        read_array_from_file("save/save.source_term_tot.txt",  fN_T_pbar,   fQ      );
        read_array_from_file("save/save.source_term_pp.txt",   fN_T_pbar,   fQ_pp   );
        read_array_from_file("save/save.source_term_pHe.txt",  fN_T_pbar,   fQ_Hep  );
        read_array_from_file("save/save.source_term_Hep.txt",  fN_T_pbar,   fQ_pHe  );
        read_array_from_file("save/save.source_term_HeHe.txt", fN_T_pbar,   fQ_HeHe );
        read_array_from_file("save/save.fC_uncertainty_Tpbar_Tp_eta_LAB__pp.txt",   fN_T_pbar*fN_Tn_p*fN_eta,  fC_Tpbar_Tp_eta_LAB__pp );
        ///////
        
        
        std::system("mkdir 3D_s_xR_pT/pp_CS_data");
        
        double sqrtS_list[] = {6.2, 13.8, 17.3, 19.4, 23.5, 27.4, 30.8, 44.8, 53., 63.  };
        
        write_array_to_file("3D_s_xR_pT/pp_CS_data/sqrtS_values.txt", 10, sqrtS_list);
        
        for (int j=0; j<10; j++) {
            double s = sqrtS_list[j]*sqrtS_list[j];
            
            for (int i=0; i<fN_xR-1; i++) {
                double xR = fxR[i+1];
                
                for (int k=0; k<fN_pT-1; k++) {
                    double pT = fPT[k];
                    f2D__xR_pT     [i + fN_pT*k] = GetLab( s, xR, pT, 0,  fC_Tpbar_Tp_eta_LAB__pp    );
                }
            }
            
            std::cout << std::endl;
            out(s)
            
            
            if (j%fMod==0)
                write_2D_array_to_file( QString("3D_s_xR_pT/pp_CS_data/2D_relativeUncertainty_xR_pT__%1__pp.txt").    arg(j).toStdString(),     fN_pT,  fN_xR,   f2D__xR_pT,   fPT, fxR  );
            
            
        }
        
        
    }

    // extract uncertainty profile for fixed Tp
    if (step==8) {
        
        /////// read calculations from previous steps
        read_array_from_file("save/save.pbar_relative_uncertainty_parmetrization.txt",     fN_T_pbar,      fRU_Tpbar__Q);
        read_array_from_file("save/save.source_term_tot.txt",  fN_T_pbar,   fQ      );
        read_array_from_file("save/save.source_term_pp.txt",   fN_T_pbar,   fQ_pp   );
        read_array_from_file("save/save.source_term_pHe.txt",  fN_T_pbar,   fQ_Hep  );
        read_array_from_file("save/save.source_term_Hep.txt",  fN_T_pbar,   fQ_pHe  );
        read_array_from_file("save/save.source_term_HeHe.txt", fN_T_pbar,   fQ_HeHe );
        read_array_from_file("save/save.fC_uncertainty_Tpbar_Tp_eta_LAB__pp.txt",   fN_T_pbar*fN_Tn_p*fN_eta,  fC_Tpbar_Tp_eta_LAB__pp );
        ///////
        
        
        std::system("mkdir 3D_Tp_Tpbar_eta_LAB/pp_fixed_Tp");
        
        double Tp_list[] = {450, 4e3, 6.5e3};
        write_array_to_file("3D_Tp_Tpbar_eta_LAB/pp_fixed_Tp/Tp_values.txt", 3, Tp_list);
        
        
        int Tp_index[3];
        
        // find closest T_p index
        for (int i=0; i<3; i++) {
            double diff=1e90;
            for (int j=0; j<fN_Tn_p; j++) {
                double test = fabs(log(Tp_list[i]/fTn_p_LAB[j]));
                if (test<diff) {
                    diff=test;
                    Tp_index[i] = j;
                }
            }
        }
        
        
        for (int ii=0; ii<3; ii++) {
            int i = Tp_index[ii];
            for (int j=0; j<fN_T_pbar; j++) {
                for (int k=0; k<fN_eta-1; k++) {
                    
                        f2D__Tpbar_eta[j + fN_Tn_p*k] = fC_Tpbar_Tp_eta_LAB__pp[i + fN_Tn_p*k + fN_Tn_p*fN_eta*j];
                    
                }
            }
            
            write_2D_array_to_file( QString("3D_Tp_Tpbar_eta_LAB/pp_fixed_Tp/2D_relativeUncertainty_Tpbar_eta__%1__pp.txt").arg(ii).toStdString(),      fN_eta,  fN_T_pbar,   f2D__Tpbar_eta,    fEta_LAB, fT_pbar_LAB  );
            
        }
        
        
    }
    
//    double Tp    = 50;
//    double Tpbar = 5;
//    double eta   = 3;
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


