#include <iostream>
#include <sstream>

#include <definitions.h>
#include <configHandler.h>
#include <readFluxes.h>
#include <lisFluxes.h>
#include <fileTools.h>
#include <crossSections.h>
#include <labCrossSectionTabulator.h>
#include <sourceTerm.h>

#include <QString>
#include <QStringList>

#include <TMinuit.h>
#include<TCanvas.h>
#include<TF1.h>
#include<TF2.h>
#include<TF3.h>
#include<TMath.h>
#include<TAxis.h>
#include <cmath>
#include <cstdlib>
#include <string>
#include <ctime>


using namespace  CRACS;

std::string fProgramDescription = "Program to fit cross section parameterizations by means of a TMinuit fit";

int fStep        = 1;
std::string fDir = ".";
int fMod         = 30;

/******************************       Struct and vectors to store experimental data       ********************************/
struct data{
    int n;
    int A;
    std::vector<double> sqrtS;
    std::vector<double> xR;
    std::vector<double> pT;
    std::vector<double> p;
    std::vector<double> var_fix;
    std::vector<double> var_plot;
    std::vector<double> CS;
    std::vector<double> CS_err;
    std::vector<double> CS_err_stat;
    std::vector<double> CS_err_sys;
    std::vector<double> err_scale;
    double          err_norm;
    std::string     experiment;
    bool            in_fit;
    double          w;
    std::vector<std::string> col_description;
};

std::vector<data>   fData_pp_pbar_all;
std::vector<data>   fData_pHe_pbar_all;
std::vector<data>   fData_pBe_pbar_all;
std::vector<data>   fData_pC_pbar_all;

std::vector<data>   fData_pA_all;

void    ReadData( std::string include_data, bool correctHyperon,
                 std::vector<data> &fData = fData_pp_pbar_all,
                 std::string type="pp",
                 QString experiments=                   ",allaby,antreasyan,brahms,capiluppi,dekkers,guettler,johnson,phenix,na49,na61,",
                 QString hyperon_inclusive_experiments= ",allaby,antreasyan,brahms,capiluppi,dekkers,guettler,johnson,phenix,",
                 int A=1  );


void    Write_results( std::vector<data> &fData );

static int          fPrint           = 0;

/******************************       Constant for the first analysis step: scan  with pp data       ********************************/
int                 fCS              = CS::DI_MAURO12;
int                 fNCparam         = 8;
int                 fNWparam         = 0;
static double       fParameters[100]; // = {4.448, 3.735, 0.00502, 0.708, 3.527, 0.236, -0.729, 2.517, -1.822e-11, 3.527, 0.384}; // Set di Mauro 12 as default, all norm. factors 1 by default
static double       fStart_pp  [100]; // = {4.448, 3.735, 0.00502, 0.708, 3.527, 0.236, -0.729, 2.517, -1.822e-11, 3.527, 0.384};
static double       fStep_pp   [100]; // = {0.01 , 0.1 , 0.001,   0.001, 0.1 ,    0.01 , 0.01,  0.01,  0.001 , 0.01,  0.01  };

double              fCovariance [40][40];
double              fChiSq_min_pp    = 0;

int                 fPrecision = 100;


// Array for 1 and 2 sigma regions. It contains the Delta chiSq in n dimensional parameter
//      e.g.  fSigma_1_region[8]   gives the 1 sigma Delta chiSq in an 8D parameter space
double              fSigma_1_region[] = {1.000043426, 2.295815161, 3.526822180, 4.719568761, 5.887700474, 7.038515492, 8.176359741, 9.304044023, 10.42350189, 11.53612748, 12.64296378, 13.74481437, 14.84231367, 15.93597259, 17.02620974, 18.11337308, 19.19775552, 20.27960639, 21.35913992, 22.43654182};
double              fSigma_2_region[] = {4.000009776, 6.180085905, 8.024894670, 9.715641142, 11.31387084, 12.84885057, 14.33712678, 15.78911024, 17.21184692, 18.61036514, 19.98840090, 21.34881933, 22.69387460, 24.02537791, 25.34481051, 26.65340209, 27.95218680, 29.24204418, 30.52372968, 31.79789790};

void CopyStart(  double* start, double* step  ){
    for (int i=0; i<fNCparam; i++) {
        fParameters[i]  = start[i];
        fStart_pp[i]    = start[i];
        fStep_pp[i]     = step[i];
    }
    
}
void SetStart(){
    if (fCS==CS::DI_MAURO12) {
        double start[]      =   {3.59482,  5.4808,  0.0330869,  -0.137454,  2.79783,  0.0477575,  -0.0390393,  2.60951}; //4.5,   3.3,  0.009,   0.45,  3.5,   0.06, -0.25, 2.6  };
        double step []      =   {0.1 ,  0.1 , 0.001,   0.01,  0.1 ,  0.01,  0.01, 0.01 };
        fNCparam = 8;
        CopyStart( start, step );
    }
    if (fCS==CS::WINKLER) {
        double start[]      =   {0.047, 7.76, 0.168, 0.038,  1.0e-3, 0.7};
        double step []      =   {0.001, 0.10, 0.010, 0.001,  1.0e-4, 0.1};
        fNCparam = 6;
        CopyStart( start, step );
    }
}

/******************************       Constant for the second analysis step: scan/fit with pA data   ********************************/
int                 fNDparam            = 4;
int                 fNWparam_pA         = 0;
static double       fParameters_pA[100] = { 8.78792e-01,   1.14248e-01, 6.36834e+00, -0.5 };

double              fCovariance_pA[40][40];

const int           CONSTANT            =  0;
const int           LINEAR              =  1;
const int           POWERLAW            =  2;
const int           PT_RESCALING        =  3;
int                 fScaling_asymm      =  POWERLAW;

double              fChiSq_min_pA    = 0;


int                 fN                  = 501; // number of iterations for the TMinuit fit!
static Double_t     fStart_pA[100]      = { 0.7,   0.00,  6.36834e+00,  2.0   };
static Double_t     fStep_pA [100]      = { 0.01,  0.01,  0.1,          0.01  };

static Double_t     fD_lower[100]       = { 0.5,   0.05,  1, -1   };
static Double_t     fD_upper[100]       = { 1.2,   0.2,  15,  0   };

void CopyStart_pA(  double* start, double* step  ){
    for (int i=0; i<fNDparam; i++) {
        fParameters_pA[i]  = start[i];
        fStart_pA[i]       = start[i];
        fStep_pA[i]        = step[i];
    }
}



static const int UPPER =  1;
static const int MEAN  =  0;
static const int LOWER = -1;

int fHyperonKind = MEAN;
int fIsospinKind = MEAN;

double deltaHyperon(double s){
    
    double C1       = 0.31;
    double C2       = 0.30;
    double C3       = 146.*146.;
    double C4       = 0.9;
    double factor   = 0.81;
    
    if( fHyperonKind == UPPER ){
        C1 = 0.31 + 0.0375;
        C2 = 0.30 + 0.0125;
        C3 = 146.*146.;
        C4 = 0.9;
        factor = 0.85;
    }
    
    if( fHyperonKind == LOWER ){
        C1 = 0.31 - 0.0375;
        C2 = 0.30 - 0.0125;
        C3 = 146.*146.;
        C4 = 0.9;
        factor = 0.77;
    }
    
    double hyperon = C1 + C2/(1+pow(C3/s,C4));
    hyperon       *= factor;  // branching ratio to Lambda and Sigma  ###### UNCERTAINTY??
    return hyperon;
    
}

double deltaIsospin(double s){
    
    double C14 = 0.114;
    double C15 = 144*144;
    double C16 = 0.51;
    
    if( fIsospinKind == UPPER ){
        C14 = 0.114  + 0.1;
        C15 = 144*144;
        C16 = 0.51;
    }
    
    if( fIsospinKind == LOWER ){
        C14 = 0.114 - 0.1;
        C15 = 144*144;
        C16 = 0.51;
    }
    
    double isospin = C14/(1+pow(s/C15, C16));
    
    return isospin;
    
}



//******************************************************************************************************************************/
//  Functions with cross section parameterization to pp, chiSq determination for pp,
//  and function to marginalize over normalization of experimental data.
//******************************************************************************************************************************/

Double_t fun_pp_pbar_CM(double sqrtS,  double pT, double xR, Double_t *par)
{
    double s            =   sqrtS*sqrtS;
    double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrtS;
    double E_pbar       =   xR*E_pbar_Max;
    double pT_pbar      =   pT;
    double ret          =   0;
    if (fCS==CS::DI_MAURO12)
        ret          =   ppCSParametrizations::inv_pp_pbar_CM__diMauro(s, E_pbar, pT_pbar, par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], 0, 0, 0);
    if (fCS==CS::WINKLER){
        
        double C5 = par[0];
        double C6 = par[1];
        double C7 = par[2];
        double C8 = par[3];
        double C9 = par[4];
        double C10= par[5];
        
        
        
        double C1 = 0.31;
        double C2 = 0.30;
        double C3 = 146.*146.;
        double C4 = 0.9;
        
        
//        double C7 = 0.168;
//        double C8 = 0.038;
        
        
        double C11 = 30.9;
        double C12 = -1.74;
        double C13 = 0.71;
        
        double C14 = 0.114;
        double C15 = 144*144;
        double C16 = 0.51;

        ret        = ppCSParametrizations::inv_pp_pbar_CM__Winkler(s, E_pbar, pT_pbar, false, false, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15, C16);
    }
    if (fCS==-1){
        ret        = ppCSParametrizations::inv_pp_pbar_CM__Winkler(s, E_pbar, pT_pbar, false, false);
    }
    return ret;
}

void fcn_pp_pbar_CM(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    double*     cpar    = par;
    double*     wpar    = &par[fNCparam];
    Double_t    chisq   = 0;
    int         ndf     = 0;
    int         ee      = 0;
    double      delta;
    for (int e=0; e<fData_pp_pbar_all.size(); e++) {
        data & d = fData_pp_pbar_all.at(e);
        if(!d.in_fit){
            continue;
        }
        double chisq_e = 0;
        for (int i=0;i<d.n; i++) {
            double sqrtS    =   d.sqrtS[i];
            double xR       =   d.xR   [i];
            double pT       =   d.pT   [i];
            delta           =   pow( (wpar[ee]*d.CS[i] - fun_pp_pbar_CM(sqrtS, pT, xR, cpar))/d.CS_err[i]/wpar[ee], 2 );
            //if (delta!=delta) continue;
            chisq          +=   delta;
            if (fPrint) {
                chisq_e    +=   delta;
                if (delta>30) {
                    std::cout << warnout << std::endl;
                    out(d.experiment)
                    int line = i;
                    out(line)
                    out(delta)
                    out(wpar[ee])
                    out(sqrtS)
                    out(xR)
                    out(pT)
                    out(d.CS[i])
                    out(d.CS_err[i])
                    out(fun_pp_pbar_CM(sqrtS, pT, xR, cpar))
                    std::cout << warnend << std::endl;
                }
            }
        }
        chisq          +=   pow( (1-wpar[ee])/d.err_norm, 2);
        if(fPrint){
            std::cout << "*************************" << std::endl;
            chisq_e     += pow( (1-wpar[ee])/d.err_norm, 2);
            out(d.experiment)
            out(chisq_e)
            out(d.n)
            out(d.w)
            double chisq_ndf = chisq_e/d.n;
            out(chisq_ndf)
            std::cout << "*************************" << std::endl;
            
        }
        ndf+= d.n;
        d.w = wpar[ee];
        ee ++;
    }
    if(fPrint){
        std::cout << "-------------------------" << std::endl;
        std::cout << "*************************" << std::endl;
        out(chisq)
        ndf = ndf-fNCparam-fNWparam;
        out(ndf)
        out(chisq/ndf)
        std::cout << "*************************" << std::endl;
        std::cout << "-------------------------" << std::endl;
        
    }
    f = chisq;
}

void Print_chiSq(data d){
    
    double*     cpar    = fParameters;
    double      chisq   = 0;
    double      delta;
    
    double chisq_e = 0;
    for (int i=0;i<d.n; i++) {
        double sqrtS    =   d.sqrtS[i];
        double xR       =   d.xR   [i];
        double pT       =   d.pT   [i];
        delta           =   pow( ( d.w*d.CS[i] - fun_pp_pbar_CM(sqrtS, pT, xR, cpar))/d.CS_err[i]/d.w, 2 );
        //if (delta!=delta) continue;
        chisq          +=   delta;
        chisq_e    +=   delta;
        if (delta>30) {
            std::cout << warnout << std::endl;
            out(d.experiment)
            int line = i;
            out(line)
            out(delta)
            out(d.w)
            out(sqrtS)
            out(xR)
            out(pT)
            out(d.CS[i])
            out(d.CS_err[i])
            out(fun_pp_pbar_CM(sqrtS, pT, xR, cpar))
            std::cout << warnend << std::endl;
        }
        
    }
    chisq          +=   pow( (1-d.w)/d.err_norm, 2);
    std::cout << "*************************" << std::endl;
    chisq_e     += pow( (1-d.w)/d.err_norm, 2);
    out(d.experiment)
    out(chisq_e)
    out(d.n)
    out(d.w)
    double chisq_ndf = chisq_e/d.n;
    out(chisq_ndf)
    std::cout << "*************************" << std::endl;
}

void Ifit_pp()
{
    
    out(fNCparam)
    out(fNWparam)
    
    fPrint = 0;
    
    TMinuit *gMinuit = new TMinuit(fNCparam+fNWparam);
    gMinuit->SetPrintLevel(-1);
    gMinuit->SetFCN( fcn_pp_pbar_CM );
    
    Double_t arglist[100];
    Int_t ierflg = 0;
    
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
    
    // Set starting values and step sizes for parameters
    std::string prefix = "C";
    for (int i = 0; i<fNCparam; i++) {
        std::stringstream ss;
        ss << i+1;
        gMinuit->mnparm(i, (prefix+ss.str()).c_str(), fStart_pp[i], fStep_pp[i], 0,0,ierflg);
    }
    prefix = "W";
    for (int i = 0; i<fNWparam; i++) {
        std::stringstream ss;
        ss << i;
        gMinuit->mnparm(i+fNCparam, (prefix+ss.str()).c_str(), 1., 0.01, 0,0,ierflg);
    }
    
    if(fCS==CS::WINKLER){
        gMinuit->FixParameter(3);
    }

    arglist[0] = 500;
    arglist[1] = 1.;
    gMinuit->SetPrintLevel(0);
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    
    // Print results
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    
    
    int n;
    double paramErrors[fNCparam];
    
    for (int i=0; i<fNCparam+fNWparam; i++) {
        double p, pErr;
        gMinuit->GetParameter(i, p, pErr);
        fParameters[i] = p;
        paramErrors[i] = pErr;
    }
    double h, chisq;
    fPrint = 1;
    fcn_pp_pbar_CM(n, &h, chisq, fParameters, 1);
    
    fPrint = 0;
    
    CopyStart(fParameters, paramErrors);
    arglist[0] = 500;
    arglist[1] = 1.;
    
    gMinuit->SetPrintLevel(0);
    gMinuit->mnexcm("HESSE", arglist ,2,ierflg);
    
    // Print results
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    
    
    for (int i=0; i<fNCparam+fNWparam; i++) {
        double p, pErr;
        gMinuit->GetParameter(i, p, pErr);
        fParameters[i] = p;
        paramErrors[i] = pErr;
    }
    fPrint = 1;
    fcn_pp_pbar_CM(n, &h, chisq, fParameters, 1);
    
    fChiSq_min_pp = chisq;
    
    double matrix [fNCparam+fNWparam][fNCparam+fNWparam];
    
    int dim =fNCparam+fNWparam;
    if (fCS==CS::WINKLER){
        dim --;
    }
    gMinuit->mnemat(&matrix[0][0],fNCparam+fNWparam);
    
    if(fCS==CS::WINKLER){
        for (int i=dim-1; i>=3; i--) {
            for (int j = 0; j<dim; j++) {
                matrix[i+1][j] = matrix[i][j];
            }
        }
        for (int i=dim-1; i>=3; i--) {
            for (int j = 0; j<dim+1; j++) {
                matrix[j][i+1] = matrix[j][i];
            }
        }
        for (int i = 0; i<fNCparam+fNWparam; i++) {
            matrix[3][i] = 0.;
            matrix[i][3] = 0.;
        }
        matrix[3][3] = pow(0.00057,2);
    }
    
    out( "Covariance Matrix" )
    std::string write = "";
    for (int i = 0; i<fNCparam+fNWparam; i++) {
        double par,err;
        gMinuit->GetParameter(i, par, err);
        write += QString("%1 ").arg(par).leftJustified(15).toStdString();
        for (int j = 0; j<fNCparam+fNWparam; j++) {
            fCovariance[i][j] = matrix[i][j];
            std::cout.width(12); std::cout << std::left <<  fCovariance[i][j] << " ";
            write += QString("%1 ").arg(fCovariance[i][j]).leftJustified(15).toStdString();
        }
        write += "\n";
        std::cout << std::endl;
    }
    std::cout <<fDir<<std::endl;
    FileTool::WriteStringToFile(write, "covariance_matrix_pp.txt", fDir+"/" );
    
    out( "Correlation Matrix" )
    for (int i = 0; i<fNCparam+fNWparam; i++) {
        for (int j = 0; j<fNCparam+fNWparam; j++) {
            std::cout.width(12); std::cout << std::left << "&" + QString("%1").arg( fCovariance[i][j]/sqrt(fCovariance[i][i]*fCovariance[j][j]) , 6, 'f', 3).toStdString();
        }
        std::cout << std::endl;
    }
    std::cout <<fDir<<std::endl;
    
    return;
    
}


//******************************************************************************************************************************/
//  END of functions needed for the pp scan
//******************************************************************************************************************************/




//******************************************************************************************************************************/
//  Functions with factor parameterization for pA
//******************************************************************************************************************************/


static double pbar_overlap_function_projectile(double x_F){
    
    int n = 25;
    
    if( x_F < -0.2499 ) return 0;
    if( x_F >  0.25   ) return 1;
    
    double xF[] = {-0.25, -0.225, -0.2, -0.175, -0.15, -0.125, -0.1, -0.075, -0.05, -0.025, 0., 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25};
    double F[]  = {0., 0.0003, 0.0008, 0.0027, 0.010, 0.035, 0.110, 0.197, 0.295, 0.4, 0.5, 0.6, 0.705, 0.803, 0.890, 0.965, 0.990, 0.9973, 0.9992, 0.9997, 1.0};
    
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

static double pbar_overlap_function_target(double x_F){
    
    return 1.-pbar_overlap_function_projectile(x_F);
    
}


double fun_pA_pbar_factor(int A, double sqrtS,  double pT, double x_R_d, Double_t *par)
{
    
    double s = sqrtS * sqrtS;
    
    double D1  = par[0];
    double D2  = par[1];
    
    double factor1       =   pow( A, D1);
    
    double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
    
    double E_pbar = fabs(x_R_d)*E_pbar_Max;
    double p  = sqrt( E_pbar*E_pbar-fMass_proton*fMass_proton );
    double pL = sqrt( p*p-pT*pT );
    
    if (pL!=pL){
        pL = 0;
    }
    
    if (x_R_d<0) pL*=-1;
    
    
    double xF = pL*2./sqrt(s);
    
    double factor_pT = 1.;
    double nu = 1;
    if (fScaling_asymm==LINEAR){
        nu = 1+(D2-1)*(A-1)/11.;
    }
    else if(fScaling_asymm==POWERLAW){
        nu = pow(A, D2);
    }
    else if(fScaling_asymm==PT_RESCALING){
        double D3  = par[2];
        nu = pow(A, D2);
        factor_pT = exp( log(A)*D3*pT );
    }
    
    double C14 = 0.114;
    double C15 = 144*144;
    double C16 = 0.51;
    double isospin = C14/(1+pow(s/C15, C16));
    
    double factor2 = factor1 * ( (2.+isospin)/2.*nu*pbar_overlap_function_target(xF) +  pbar_overlap_function_projectile(xF) )*factor_pT;

    if (fCS==-1){
        double nu=1;
        if (A== 1) isospin = 1.0;
        if (A== 4) nu = 1.25;
        if (A==12) nu = 1.60;
        factor2 = A/nu* ( (2.+isospin)/2.*nu*pbar_overlap_function_target(xF) +  pbar_overlap_function_projectile(xF) );
    }
    
    return factor2;
    
}


double fun_AA_pbar_factor(int A1, int A2, int N1, int N2, double sqrtS,  double pT, double Epbar_d, Double_t *par)
{
    
    double s = sqrtS * sqrtS;
    
    double D1  = par[0];
    double D2  = par[1];
    
    double E_pbar = fabs(Epbar_d);
    double p  = sqrt( E_pbar*E_pbar-fMass_proton*fMass_proton );
    double pL = sqrt( p*p-pT*pT );
    
    if (pL!=pL){
        pL = 0;
        return 0;
    }
    if (Epbar_d<0) pL*=-1;
    
    
    double xF = pL*2./sqrt(s);
    
    double C14 = 0.114;
    double C15 = 144*144;
    double C16 = 0.51;
    double isospin = C14/(1+pow(s/C15, C16));
    
    
    double nu  = pow(A1, D2);
    double nu2 = pow(A2, D2);
    
    
    double factor = pow( A1, D1) * pow( A2, D1) * ( (1.+N1*1./A1*isospin)*nu*pbar_overlap_function_projectile(xF) + (1.+N2*1./A2*isospin)*nu2*pbar_overlap_function_target(xF) );
    
    return factor;
    
}

double fun_AA_pbar_factor(int A1, int A2, double sqrtS,  double pT, double Epbar_d, Double_t *par)
{
    
    int N1 = 0.5*A1;
    if(A1==1){
        N1=0;
    }
    int N2 = 0.5*A2;
    if(A2==1){
        N2=0;
    }
    return fun_AA_pbar_factor( A1, A2, N1, N2, sqrtS, pT, Epbar_d, par);
    
}




// TMinuit
//***************************************************************

void fcn_pA_pbar_CM(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    double* dpar    =   par;
    double* wpar    =&  par[fNDparam];
    int ndf         =   0;
    Double_t chisq  =   0;
    double   delta;
    int     ee      = 0;
    for (int e=0; e<fData_pA_all.size(); e++) {
        data &  d       = fData_pA_all[e];
        int     A       = d.A;
        double  chisq_e = 0;
        
        if(!d.in_fit)
            continue;
        for (int i=0;i<d.n; i++) {
            double sqrtS    =   d.sqrtS[i];
            double xR       =   d.xR   [i];
            double pT       =   d.pT   [i];
            varOut(d.experiment)
            varOut(wpar[ee])
            delta           =   pow( (wpar[ee]*d.CS[i] - fun_pp_pbar_CM(sqrtS, pT, xR, fParameters)
                                      *fun_pA_pbar_factor( A, sqrtS, pT, xR, dpar))/d.CS_err[i]/wpar[ee], 2 );
            //if (delta!=delta) continue;
            chisq          +=   delta;
            chisq_e        +=   delta;
            if (fPrint) {
                if (delta>30) {
                    std::cout << warnout << std::endl;
                    out(d.experiment)
                    int line = i;
                    out(line)
                    out(delta)
                    out(wpar[ee])
                    out(sqrtS)
                    out(xR)
                    out(pT)
                    out(d.CS[i])
                    out(d.CS_err[i])
                    out(fun_pp_pbar_CM(sqrtS, pT, xR, fParameters))
                    out(fun_pA_pbar_factor( A, sqrtS, pT, xR, dpar))
                    std::cout << warnend << std::endl;
                }
            }
        }
        chisq          +=   pow( (1-wpar[ee])/d.err_norm, 2);
        if(fPrint){
            std::cout << "*************************" << std::endl;
            chisq_e     += pow( (1-wpar[ee])/d.err_norm, 2);
            out(d.experiment)
            out(chisq_e)
            out(d.n)
            out(d.w)
            double chisq_ndf = chisq_e/d.n;
            out(chisq_ndf)
            std::cout << "*************************" << std::endl;
        }
        ndf+= d.n;
        varOut(d.experiment)
        varOut(fData_pHe_pbar_all.size())
        varOut(fData_pBe_pbar_all.size())
        varOut(fData_pC_pbar_all.size() )
        varOut(e)
        varOut(A)
        varOut(chisq_e)
        ee ++;
        //        chisq          +=   pow( log(wpar[ee])/log(1+fEps[e]), 2);
    }
    
    if(fPrint){
        std::cout << "-------------------------" << std::endl;
        std::cout << "*************************" << std::endl;
        out(chisq)
        ndf = ndf-fNDparam-fNWparam_pA;
        out(ndf)
        out(chisq/ndf)
        std::cout << "*************************" << std::endl;
        std::cout << "-------------------------" << std::endl;
    }
    f = chisq;
}

void Print_chiSq_pA(data d, bool print_detail=false){
    
    double* dpar    =   fParameters_pA;
    double  chisq   =   0;
    double  delta;
    
    int     A       = d.A;
    double  chisq_e = 0;
    
    for (int i=0;i<d.n; i++) {
        double sqrtS    =   d.sqrtS[i];
        double xR       =   d.xR   [i];
        double pT       =   d.pT   [i];
        varOut(d.experiment)
        delta           =   pow( (d.w*d.CS[i] - fun_pp_pbar_CM(sqrtS, pT, xR, fParameters)
                                  *fun_pA_pbar_factor( A, sqrtS, pT, xR, dpar))/d.CS_err[i]/d.w, 2 );
        //if (delta!=delta) continue;
        chisq          +=   delta;
        chisq_e        +=   delta;
        if (delta>30 && print_detail) {
            std::cout << warnout << std::endl;
            out(d.experiment)
            int line = i;
            out(line)
            out(delta)
            out(d.w)
            out(sqrtS)
            out(xR)
            out(pT)
            out(d.CS[i])
            out(d.CS_err[i])
            out(fun_pp_pbar_CM(sqrtS, pT, xR, fParameters))
            out(fun_pA_pbar_factor( A, sqrtS, pT, xR, dpar))
            std::cout << warnend << std::endl;
        }
        
    }
    
    std::cout << "*************************" << std::endl;
    chisq_e     += pow( (1-d.w)/d.err_norm, 2);
    out(d.experiment)
    out(chisq_e)
    out(d.n)
    out(d.w)
    double chisq_ndf = chisq_e/d.n;
    out(chisq_ndf)
    std::cout << "*************************" << std::endl;
    
}




double Ifit_pA()
{
    fPrint = 0;
    
    TMinuit *gMinuit = new TMinuit(fNDparam+fNWparam_pA);
    gMinuit->SetPrintLevel(-1);
    gMinuit->SetFCN( fcn_pA_pbar_CM );
    
    Double_t arglist[10];
    Int_t ierflg = 0;
    
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
    
    // Set starting values and step sizes for parameters
    std::string prefix = "D";
    for (int i = 0; i<fNDparam; i++) {
        std::stringstream ss;
        ss << i+1;
        gMinuit->mnparm(i, (prefix+ss.str()).c_str(), fStart_pA[i], fStep_pA[i], 0,0,ierflg);
    }
    prefix = "W_pA";
    for (int i = 0; i<fNWparam_pA; i++) {
        std::stringstream ss;
        ss << i;
        gMinuit->mnparm(i+fNDparam, (prefix+ss.str()).c_str(), 1., 0.01, 0,0,ierflg);
    }
    
    
    gMinuit->SetPrintLevel(0);
    
    // Now ready for minimization step
    arglist[0] = 1000000;
    arglist[1] = 1.;
    gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
    
    // Print results
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    
    
    int n;
    double paramErrors[fNDparam];
    
    for (int i=0; i<fNDparam+fNWparam_pA; i++) {
        double p, pErr;
        gMinuit->GetParameter(i, p, pErr);
        fParameters_pA[i] = p;
        paramErrors[i]    = pErr;
    }
    double h, chisq;
    fPrint = 1;
    fcn_pA_pbar_CM(n, &h, chisq, fParameters_pA, 1);
    fPrint = 0;
    
    CopyStart_pA(fParameters_pA, paramErrors);
    arglist[0] = 500;
    arglist[1] = 1.;
    
    gMinuit->SetPrintLevel(0);
    gMinuit->mnexcm("HESSE", arglist ,2,ierflg);
    
    // Print results
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    
   
    for (int i=0; i<fNDparam+fNWparam_pA; i++) {
        double p, pErr;
        gMinuit->GetParameter(i, p, pErr);
        //out(p)
        fParameters_pA[i] = p;
    }
    fPrint = 1;
    fcn_pA_pbar_CM (n, &h, chisq, fParameters_pA, 1);
    
    int ee = 0;
    for (int e=0; e<fData_pA_all.size(); e++) {
        data &  d       = fData_pA_all[e];
        if(d.in_fit){
            d.w = fParameters_pA[fNDparam+ee];
            ee ++;
        }
    }
    
    // Print results
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
    
    fChiSq_min_pA = chisq;
    
    double matrix [fNCparam+fNWparam][fNDparam+fNWparam_pA];
    
    int dim =fNDparam+fNWparam_pA;
    
    gMinuit->mnemat(&matrix[0][0],dim);
    
    out( "Covariance Matrix" )
    std::string write = "";
    for (int i = 0; i<dim; i++) {
        double par,err;
        gMinuit->GetParameter(i, par, err);
        write += QString("%1 ").arg(par).leftJustified(15).toStdString();
        for (int j = 0; j<dim; j++) {
            fCovariance_pA[i][j] = matrix[i][j];
            std::cout.width(12); std::cout << std::left <<  fCovariance_pA[i][j] << " ";
            write += QString("%1 ").arg(fCovariance_pA[i][j]).leftJustified(15).toStdString();
        }
        write += "\n";
        std::cout << std::endl;
    }
    std::cout <<fDir<<std::endl;
    FileTool::WriteStringToFile(write, "covariance_matrix_pA.txt", fDir+"/" );
    
    out( "Correlation Matrix" )
    for (int i = 0; i<dim; i++) {
        for (int j = 0; j<dim; j++) {
            std::cout.width(12); std::cout << std::left << "&" + QString("%1").arg( fCovariance_pA[i][j]/sqrt(fCovariance_pA[i][i]*fCovariance_pA[j][j]) , 6, 'f', 3).toStdString();
        }
        std::cout << std::endl;
    }
    std::cout <<fDir<<std::endl;
    
    
    return chisq;
    
}




//******************************************************************************************************************************/
//  END Functions with factor parameterization for pA
//******************************************************************************************************************************/





//******************************************************************************************************************************/
//  Functions for the source term calculation
//******************************************************************************************************************************/

CS_lab_tab* fTab;
double dT__pbar_LAB_interpolation( double Tn_proj, double Tn_prod ){
    return fTab->GetInterpolation( Tn_proj, Tn_prod );
}

double protonFlux(double T_p){
    return LISFluxes::spectrum(T_p, LISFluxes::fProtonLIS);
}
double heliumFlux(double T_p){
    return LISFluxes::spectrum(T_p, LISFluxes::fHeliumLIS);
}
double carbonFlux(double T_p){
    return LISFluxes::spectrum(T_p, LISFluxes::fCarbonLIS);
}
double nitrogenFlux(double T_p){
    return LISFluxes::spectrum(T_p, LISFluxes::fNitrogenLIS);
}
double oxygenFlux(double T_p){
    return LISFluxes::spectrum(T_p, LISFluxes::fOxygenLIS);
}


double lithiumFlux(double T_p){
    return LISFluxes::spectrum(T_p, LISFluxes::fLithiumLIS);
}
double berylliumFlux(double T_p){
    return LISFluxes::spectrum(T_p, LISFluxes::fBerylliumLIS);
}
double boronFlux(double T_p){
    return LISFluxes::spectrum(T_p, LISFluxes::fBoronLIS);
}


double neonFlux     (double T_p){
    return LISFluxes::spectrum(T_p, LISFluxes::fOxygenLIS)*150./1000;
}
double magnesiumFlux(double T_p){
    return LISFluxes::spectrum(T_p, LISFluxes::fOxygenLIS)*200./1000;
}
double siliconFlux  (double T_p){
    return LISFluxes::spectrum(T_p, LISFluxes::fOxygenLIS)*115./1000;
}

double ironFlux(double T_p){
    return LISFluxes::spectrum(T_p, LISFluxes::fOxygenLIS)*110./1000;
}


//////////////
//   pp
//////////////

static double inv_pp_pbar_CM( double s, double E_pbar, double pT_pbar){
    double ret          =   0;
    if (fCS==CS::DI_MAURO12)
        ret          =   ppCSParametrizations::inv_pp_pbar_CM__diMauro(s, E_pbar, pT_pbar, fParameters[0], fParameters[1], fParameters[2], fParameters[3], fParameters[4], fParameters[5], fParameters[6], fParameters[7], 0, 0, 0);
    if (fCS==CS::WINKLER){
      
        double C5 = fParameters[0];
        double C6 = fParameters[1];
        double C7 = fParameters[2];
        double C8 = fParameters[3];
        double C9 = fParameters[4];
        double C10= fParameters[5];
        
        double C1 = 0.31;
        double C2 = 0.30;
        double C3 = 146.*146.;
        double C4 = 0.9;
        
        double C11 = 30.9;
        double C12 = -1.74;
        double C13 = 0.71;
        
        double C14 = 0.114;
        double C15 = 144*144;
        double C16 = 0.51;
    
        ret        = ppCSParametrizations::inv_pp_pbar_CM__Winkler(s, E_pbar, pT_pbar, false, false, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15, C16);
    }
    return ret;
}
static double dE_pp_pbar_LAB(double E_p_LAB, double E_pbar_LAB, int precision=fPrecision, int output=NO_OUT){
    return CSTransformations::dE_pp_product_LAB_intEta (E_p_LAB, E_pbar_LAB, fMass_proton,
                                                        inv_pp_pbar_CM,
                                                        precision,
                                                        output);
}
static double dT_pp_pbar_LAB  ( double T_p, double T_pbar ){
    return dE_pp_pbar_LAB( T_p+fMass_proton,   T_pbar+fMass_proton)*1e-31;
}

//  AA

// AA
static double inv_A1A2_pbar_CM(  double s, double E_pbar_d, double pT_pbar, int A1, int A2){
    
    double E_pbar   = fabs(E_pbar_d);
    double CS       = inv_pp_pbar_CM( s, E_pbar, pT_pbar );
    if (CS<1e-90)
    return 0;
    double factor = fun_AA_pbar_factor( A1, A2, sqrt(s), pT_pbar, E_pbar_d, fParameters_pA );
    return factor * CS;
}



static int  fA1= 1;
static int  fA2= 1;

static double inv_A1A2_pbar_CM( double s, double E_pbar_d, double pT_pbar){
    return inv_A1A2_pbar_CM( s, E_pbar_d, pT_pbar, fA1, fA2);
}


static double inv_A1A2_pbar_CM_inclHyp_inclNbar( double s, double E_pbar_d, double pT_pbar){
    return inv_A1A2_pbar_CM( s, E_pbar_d, pT_pbar, fA1, fA2)*(2+deltaIsospin(s)+2*deltaHyperon(s));
}


static double dE_A1A2_pbar_LAB(double E_p_LAB, double E_pbar_LAB, int precision=fPrecision, int output=NO_OUT){
    return CSTransformations::dE_pp_product_LAB_intEta (E_p_LAB, E_pbar_LAB, fMass_proton,
                                                        inv_A1A2_pbar_CM,
                                                        precision,
                                                        output);
}
static double dE_A1A2_pbar_LAB_inclHyp_inclNbar(double E_p_LAB, double E_pbar_LAB, int precision=fPrecision, int output=NO_OUT){
    return CSTransformations::dE_pp_product_LAB_intEta (E_p_LAB, E_pbar_LAB, fMass_proton,
                                                        inv_A1A2_pbar_CM_inclHyp_inclNbar,
                                                        precision,
                                                        output);
}

static double dT_A1A2_pbar_LAB  ( double T_p, double T_pbar ){
    return dE_A1A2_pbar_LAB( T_p+fMass_proton,   T_pbar+fMass_proton)*1e-31;
}
static double dT_A1A2_pbar_LAB_inclHyp_inclNbar  ( double T_p, double T_pbar ){
    return dE_A1A2_pbar_LAB_inclHyp_inclNbar( T_p+fMass_proton,   T_pbar+fMass_proton)*1e-31;
}

//*************************************

static int  fN1= 1;
static int  fN2= 1;

// AA
static double inv_A1A2_N1N2_pbar_CM(  double s, double E_pbar_d, double pT_pbar, int A1, int A2, int N1, int N2){
    
    double E_pbar   = fabs(E_pbar_d);
    double CS       = inv_pp_pbar_CM( s, E_pbar, pT_pbar );
    if (CS<1e-90)
        return 0;
    double factor = fun_AA_pbar_factor( A1, A2, N1, N2, sqrt(s), pT_pbar, E_pbar_d, fParameters_pA );
    return factor * CS;
}

static double inv_A1A2_N1N2_pbar_CM( double s, double E_pbar_d, double pT_pbar){
    return inv_A1A2_N1N2_pbar_CM( s, E_pbar_d, pT_pbar, fA1, fA2, fN1, fN2);
}


static double inv_A1A2_N1N2_pbar_CM_inclHyp_inclNbar( double s, double E_pbar_d, double pT_pbar){
    return inv_A1A2_N1N2_pbar_CM( s, E_pbar_d, pT_pbar, fA1, fA2, fN1, fN2)*(2+deltaIsospin(s)+2*deltaHyperon(s));
}


static double dE_A1A2_N1N2_pbar_LAB(double E_p_LAB, double E_pbar_LAB, int precision=fPrecision, int output=NO_OUT){
    return CSTransformations::dE_pp_product_LAB_intEta (E_p_LAB, E_pbar_LAB, fMass_proton,
                                                        inv_A1A2_N1N2_pbar_CM,
                                                        precision,
                                                        output);
}
static double dE_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar(double E_p_LAB, double E_pbar_LAB, int precision=fPrecision, int output=NO_OUT){
    return CSTransformations::dE_pp_product_LAB_intEta (E_p_LAB, E_pbar_LAB, fMass_proton,
                                                        inv_A1A2_N1N2_pbar_CM_inclHyp_inclNbar,
                                                        precision,
                                                        output);
}

static double dT_A1A2_N1N2_pbar_LAB  ( double T_p, double T_pbar ){
    return dE_A1A2_N1N2_pbar_LAB( T_p+fMass_proton,   T_pbar+fMass_proton)*1e-31;
}
static double dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar  ( double T_p, double T_pbar ){
    return dE_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar( T_p+fMass_proton,   T_pbar+fMass_proton)*1e-31;
}

//******************************************************************************************************************************/
//  END   Functions for the source term calculation
//******************************************************************************************************************************/


double dT_pp_pbar_LAB_dM   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::DI_MAURO12);
}

double dT_pp_pbar_LAB_Wi   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::WINKLER);
}

double dT_pp_pbar_LAB_KR   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PP_PBAR, CS::KACHELRIESS);
}


double dT_pHe_pbar_LAB_dM   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PHe_PBAR, CS::DI_MAURO12);
}

double dT_pHe_pbar_LAB_Wi   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PHe_PBAR, CS::WINKLER);
}

double dT_pHe_pbar_LAB_KR   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::PHe_PBAR, CS::KACHELRIESS);
}



double dT_Hep_pbar_LAB_dM   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::HeP_PBAR, CS::DI_MAURO12);
}

double dT_Hep_pbar_LAB_Wi   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::HeP_PBAR, CS::WINKLER);
}

double dT_Hep_pbar_LAB_KR   ( double T_p,         double T_pbar ){
    CS_lab_tab * cs_lab = CS_lab_tab::GetInstance();
    return cs_lab->GetCrossSection(T_p,   T_pbar, CS::HeP_PBAR, CS::KACHELRIESS);
}


void writeSourceTerm( int A1, int A2, std::string suffix, bool withComparisonOld=false );



//***********************************************************************************************************************//
//***********************************************************************************************************************//
//***********************************************************************************************************************//
//***********************************************************************************************************************//
//***********************************************************************************************************************//


int fPP                 = 1;
int fCS_frac            = 0;
std::string fExperiment = "lhcb";

double len (double* a){
    return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}
void cross  (double* v1, double* v2, double* cross_res){
    cross_res[0] = v1[1]*v2[2]-v1[2]*v2[1];
    cross_res[1] = v1[2]*v2[0]-v1[0]*v2[2];
    cross_res[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

double dot(double* v1, double* v2){
    return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}
double dist( double* a, double *b , double* c, double* p ){
    
    double v1 [3] = { a[0]-b[0], a[1]-b[1], a[2]-b[2] };
    double v2 [3] = { a[0]-c[0], a[1]-c[1], a[2]-c[2] };
    double v3 [3] = { a[0]-p[0], a[1]-p[1], a[2]-p[2] };
    
    double n  [3] = {0,0,0};
    cross(v1, v2, n);
    double dn = len(n);
    n[0]      = n[0]/dn;
    n[1]      = n[1]/dn;
    n[2]      = n[2]/dn;
    
    return dot(n, v3);
    
}
double sign(double x){
    if (x<0)
        return -1;
    return 1;
}
bool inside(double* a, double *b , double* c, double* d, double* p){
    double dd, dp;
    
    dd = dist(a, b, c, d);
    if ( fabs(dd)<1e-10 || dd!=dd){
        return false;
    }
    dp = dist(a, b, c, p);
    if ( fabs(dp)<1e-10 || dp!=dp){
        return false;
    }
    dp = dp * sign(dd);
    dd = dd * sign(dd);
    if (dp<0 || dp>dd){
        return false;
    }
    
    dd = dist(a, b, d, c);
    if ( fabs(dd)<1e-10 || dd!=dd){
        return false;
    }
    dp = dist(a, b, d, p);
    if ( fabs(dp)<1e-10 || dp!=dp){
        return false;
    }
    dp = dp * sign(dd);
    dd = dd * sign(dd);
    if (dp<0 || dp>dd){
        return false;
    }
    
    dd = dist(a, c, d, b);
    if ( fabs(dd)<1e-10 || dd!=dd){
        return false;
    }
    dp = dist(a, c, d, p);
    if ( fabs(dp)<1e-10 || dp!=dp){
        return false;
    }
    dp = dp * sign(dd);
    dd = dd * sign(dd);
    if (dp<0 || dp>dd){
        return false;
    }
    
    dd = dist(b, c, d, a);
    if ( fabs(dd)<1e-10 || dd!=dd){
        return false;
    }
    dp = dist(b, c, d, p);
    if ( fabs(dp)<1e-10 || dp!=dp){
        return false;
    }
    dp = dp * sign(dd);
    dd = dd * sign(dd);
    if (dp<0 || dp>dd){
        return false;
    }
    
    
    return true;
}

bool in_lhcb( double ss, double pT, double xR, double ss_lower = 100, double ss_upper = 120 ){
    if ( ss< ss_lower || ss>ss_upper)
        return false;
    if ( pT< 0.47 || pT>3.1)
        return false;
    if ( xR< -0.151 || xR>0.030)
        return false;
    
    double ss_a [] = {  ss_lower,       ss_lower,       ss_lower,       ss_lower,       ss_lower,       ss_lower,       ss_lower,       ss_lower,       ss_lower,       ss_lower,
        ss_lower,       ss_lower,       ss_lower,       ss_lower,       ss_lower,       ss_lower,       ss_lower,       ss_lower,       ss_lower,       ss_lower,
        ss_upper,       ss_upper,       ss_upper,       ss_upper,       ss_upper,       ss_upper,       ss_upper,       ss_upper,       ss_upper,       ss_upper,
        ss_upper,       ss_upper,       ss_upper,       ss_upper,       ss_upper,       ss_upper,       ss_upper,       ss_upper,       ss_upper,       ss_upper    };
    double pT_a [] = {  0.47,           0.47,           0.62,           0.62,           0.75,           0.75,           0.85,           0.85,           0.97,           0.97,
        1.12,           1.12,           1.32,           1.32,           1.68,           1.68,           2.23,           2.23,           3.08,           3.08,
        0.47,           0.47,           0.62,           0.62,           0.75,           0.75,           0.85,           0.85,           0.97,           0.97,
        1.12,           1.12,           1.32,           1.32,           1.68,           1.68,           2.23,           2.23,           3.08,           3.08        };
    double xR_a [] = {  -0.0265659,     -0.0213735,     -0.0538737,     -0.0214785,     -0.0611934,     -0.0226092,     -0.0677726,     -0.0230614,     -0.0767581,     -0.02449,
        -0.0896666,     -0.0264808,     -0.109784,      0.0295951,      -0.153025,      -0.0351566,     -0.24001,       -0.0464441,     -0.150847,      -0.0703485,
        -0.0265659,     -0.0213735,     -0.0538737,     -0.0214785,     -0.0611934,     -0.0226092,     -0.0677726,     -0.0230614,     -0.0767581,     -0.02449,
        -0.0896666,     -0.0264808,     -0.109784,      0.0295951,      -0.153025,      -0.0351566,     -0.24001,       -0.0464441,     -0.150847,      -0.0703485  };
    int n = 40;
    
    double p [3] = { ss, pT, xR };
    
    for (int ia = 0; ia<n; ia++) {
        for (int ib = ia+1; ib<n; ib++) {
            for (int ic = ib+1; ic<n; ic++) {
                for (int id = ic+1; id<n; id++) {
                    double a [3] =  {ss_a[ia], pT_a[ia], xR_a[ia] };
                    double b [3] =  {ss_a[ib], pT_a[ib], xR_a[ib] };
                    double c [3] =  {ss_a[ic], pT_a[ic], xR_a[ic] };
                    double d [3] =  {ss_a[id], pT_a[id], xR_a[id] };
                    
                    if (inside(a, b, c, d, p)){
                        return true;
                    };
                }
            }
        }
    }
    if (fPP==1) {
        for (int ia = 0; ia<n; ia++) {
            for (int ib = ia+1; ib<n; ib++) {
                for (int ic = ib+1; ic<n; ic++) {
                    for (int id = ic+1; id<n; id++) {
                        double a [3] =  {ss_a[ia], pT_a[ia], -xR_a[ia] };
                        double b [3] =  {ss_a[ib], pT_a[ib], -xR_a[ib] };
                        double c [3] =  {ss_a[ic], pT_a[ic], -xR_a[ic] };
                        double d [3] =  {ss_a[id], pT_a[id], -xR_a[id] };
                        
                        if (inside(a, b, c, d, p)){
                            return true;
                        };
                    }
                }
            }
        }
    }
    
    return false;
    
}


bool in_lhcb_combined( double ss, double pT, double xR, double ss_lower = 27, double ss_upper = 120 ){
    if ( ss< ss_lower || ss>ss_upper)
        return false;
    if ( pT< 0.47 || pT>3.1)
        return false;
    if ( xR< -0.151 || xR>0.030)
        return false;
    
    double E_pbar_Max___110   =   ( 110*110          -8.*fMass_proton*fMass_proton )/2./   110  ;
    
    double E_pbar_Max_upper   =   ( ss_upper*ss_upper-8.*fMass_proton*fMass_proton )/2./ss_upper;
    double E_pbar_Max_lower   =   ( ss_lower*ss_lower-8.*fMass_proton*fMass_proton )/2./ss_lower;
    
    double f = E_pbar_Max___110/E_pbar_Max_lower;
    double g = E_pbar_Max___110/E_pbar_Max_upper;
    
    double ss_a [] = {
        ss_lower,       ss_lower,       ss_lower,       ss_lower,       ss_lower,       ss_lower,       ss_lower,       ss_lower,
        ss_upper,       ss_upper,       ss_upper,       ss_upper,       ss_upper,       ss_upper,       ss_upper,       ss_upper
    };
    double pT_a [] = {
        0.47,           0.47,           0.62,           0.62,           2.23,           2.23,           3.08,           3.08,
        0.47,           0.47,           0.62,           0.62,           2.23,           2.23,           3.08,           3.08
    };
    double xR_a [] = {
        -0.0265659*f,   -0.0213735*f,   -0.0538737*f,   -0.0214785*f,   -0.24001*f,     -0.0464441*f,   -0.150847*f,   -0.0703485*f,
        -0.0265659*g,   -0.0213735*g,   -0.0538737*g,   -0.0214785*g,   -0.24001*g,     -0.0464441*g,   -0.150847*g,    -0.0703485*g
    };
    int n = 8*2;
    
    double p [3] = { ss, pT, xR };
    
    for (int ia = 0; ia<n; ia++) {
        for (int ib = ia+1; ib<n; ib++) {
            for (int ic = ib+1; ic<n; ic++) {
                for (int id = ic+1; id<n; id++) {
                    double a [3] =  {ss_a[ia], pT_a[ia], xR_a[ia] };
                    double b [3] =  {ss_a[ib], pT_a[ib], xR_a[ib] };
                    double c [3] =  {ss_a[ic], pT_a[ic], xR_a[ic] };
                    double d [3] =  {ss_a[id], pT_a[id], xR_a[id] };
                    
                    if (inside(a, b, c, d, p)){
                        return true;
                    };
                }
            }
        }
    }
    if (fPP==1) {
        for (int ia = 0; ia<n; ia++) {
            for (int ib = ia+1; ib<n; ib++) {
                for (int ic = ib+1; ic<n; ic++) {
                    for (int id = ic+1; id<n; id++) {
                        double a [3] =  {ss_a[ia], pT_a[ia], -xR_a[ia] };
                        double b [3] =  {ss_a[ib], pT_a[ib], -xR_a[ib] };
                        double c [3] =  {ss_a[ic], pT_a[ic], -xR_a[ic] };
                        double d [3] =  {ss_a[id], pT_a[id], -xR_a[id] };
                        
                        if (inside(a, b, c, d, p)){
                            return true;
                        };
                    }
                }
            }
        }
    }
    
    return false;
    
}


bool in_na61( double ss, double pT, double xR, double ss_max=-1, double ss_min=-1 ){
    
    double ss_lower = 7.74;
    double ss_upper = 17.27;
    if (ss_max>0)
        ss_upper = ss_max;
    if (ss_min>0)
        ss_lower = ss_min;
    
    
    if ( ss< ss_lower || ss>ss_upper)
        return false;
    if ( pT< 0.05 || pT>1.05)
        return false;
    if ( xR< -0.2 || xR>0.5)
        return false;
    
    //    double ss_a [] = {
    //        7.7434,         7.7434,         7.7434,         7.7434,         7.7434,         7.7434,         7.7434,         7.7434,         7.7434,
    //        7.7434,         7.7434,         7.7434,         7.7434,         7.7434,         7.7434,         7.7434,         7.7434,         7.7434,
    //        8.7660,         8.7660,         8.7660,         8.7660,         8.7660,         8.7660,         8.7660,         8.7660,         8.7660,
    //        8.7660,         8.7660,         8.7660,         8.7660,         8.7660,         8.7660,         8.7660,         8.7660,         8.7660,
    //        12.325,         12.325,         12.325,         12.325,         12.325,         12.325,         12.325,         12.325,         12.325,         12.325,         12.325,
    //        12.325,         12.325,         12.325,         12.325,         12.325,         12.325,         12.325,         12.325,         12.325,         12.325,         12.325,
    //        17.270,         17.270,         17.270,         17.270,         17.270,         17.270,         17.270,         17.270,         17.270,         17.270,
    //        17.270,         17.270,         17.270,         17.270,         17.270,         17.270,         17.270,         17.270,         17.270,         17.270
    //    };
    //    double pT_a [] = {
    //        0.05,           0.15,           0.25,           0.35,           0.45,           0.55,           0.65,           0.75,           0.85,
    //        0.05,           0.15,           0.25,           0.35,           0.45,           0.55,           0.65,           0.75,           0.85,
    //        0.05,           0.15,           0.25,           0.35,           0.45,           0.55,           0.65,           0.75,           0.9,
    //        0.05,           0.15,           0.25,           0.35,           0.45,           0.55,           0.65,           0.75,           0.9,
    //        0.05,           0.15,           0.25,           0.35,           0.45,           0.55,           0.65,           0.75,           0.85,           0.95,           1.05,
    //        0.05,           0.15,           0.25,           0.35,           0.45,           0.55,           0.65,           0.75,           0.85,           0.95,           1.05,
    //        0.05,           0.15,           0.25,           0.35,           0.45,           0.55,           0.65,           0.75,           0.85,           0.95,
    //        0.05,           0.15,           0.25,           0.35,           0.45,           0.55,           0.65,           0.75,           0.85,           0.95
    //    };
    //    double xR_a [] = {
    //        0.27636,        0.279473,       0.285597,       0.294544,       0.306067,       0.319887,       0.335721,       0.353299,       0.372373,
    //        0.394076,       0.398515,       0.407247,       0.420005,       0.436436,       0.456143,       0.478722,       0.396405,       0.387318,
    //        0.246704,       0.249483,       0.25495,        0.262937,       0.262681,       0.274542,       0.288132,       0.303218,       0.319588,
    //        0.338215,       0.342024,       0.349518,       0.360468,       0.37457,        0.391483,       0.410862,       0.432373,       0.350582,
    //        0.16069,        0.1625,         0.166061,       0.171263,       0.177963,       0.185999,       0.195206,       0.205427,       0.216517,       0.228351,       0.240819,
    //        0.376127,       0.380364,       0.388698,       0.400875,       0.416558,       0.435368,       0.456918,       0.480841,       0.424612,       0.379111,       0.300763,
    //        0.112001,       0.113263,       0.115745,       0.119371,       0.124041,       -0.129642,      -0.136059,      -0.143182,      -0.150913,      -0.159161,
    //        0.380884,       0.385174,       0.393614,       0.405945,       0.421826,       0.364841,       0.318472,       0.486922,       0.295955,       0.372547
    //    };
    //   int n = 2*(9+9+11+10);
    if(ss_max<0 && ss_min<0){
        double ss_a [] = {
            7.7434,         7.7434,         7.7434,
            7.7434,         7.7434,         7.7434,
            8.7660,         8.7660,         8.7660,
            8.7660,         8.7660,         8.7660,
            12.325,         12.325,         12.325,
            12.325,         12.325,         12.325,
            17.270,         17.270,         17.270,
            17.270,         17.270,         17.270
        };
        double pT_a [] = {
            0.05,           0.65,           0.85,
            0.05,           0.65,           0.85,
            0.05,           0.75,           0.9,
            0.05,           0.75,           0.9,
            0.05,           0.75,           1.05,
            0.05,           0.75,           1.05,
            0.05,           0.75,           0.95,
            0.05,           0.75,           0.95
        };
        double xR_a [] = {
            0.27636,        0.335721,       0.372373,
            0.394076,       0.478722,       0.387318,
            0.246704,       0.303218,       0.319588,
            0.338215,       0.432373,       0.350582,
            0.16069,        0.205427,       0.240819,
            0.376127,       0.480841,       0.300763,
            0.112001,       -0.143182,      -0.159161,
            0.380884,       0.486922,       0.372547
        };
        int n = 2*4*3;
        double p [3] = { ss, pT, xR };
        
        for (int ia = 0; ia<n; ia++) {
            for (int ib = ia+1; ib<n; ib++) {
                for (int ic = ib+1; ic<n; ic++) {
                    for (int id = ic+1; id<n; id++) {
                        double a [3] =  {ss_a[ia], pT_a[ia], xR_a[ia] };
                        double b [3] =  {ss_a[ib], pT_a[ib], xR_a[ib] };
                        double c [3] =  {ss_a[ic], pT_a[ic], xR_a[ic] };
                        double d [3] =  {ss_a[id], pT_a[id], xR_a[id] };
                        
                        if (inside(a, b, c, d, p)){
                            return true;
                        };
                    }
                }
            }
        }
        if (fPP==1) {
            for (int ia = 0; ia<n; ia++) {
                for (int ib = ia+1; ib<n; ib++) {
                    for (int ic = ib+1; ic<n; ic++) {
                        for (int id = ic+1; id<n; id++) {
                            double a [3] =  {ss_a[ia], pT_a[ia], -xR_a[ia] };
                            double b [3] =  {ss_a[ib], pT_a[ib], -xR_a[ib] };
                            double c [3] =  {ss_a[ic], pT_a[ic], -xR_a[ic] };
                            double d [3] =  {ss_a[id], pT_a[id], -xR_a[id] };
                            
                            if (inside(a, b, c, d, p)){
                                return true;
                            };
                        }
                    }
                }
            }
        }
        
    }else if( ss_max>0 && ss_min<0 ){
        double ss_a [] = {
            7.7434,         7.7434,         7.7434,
            7.7434,         7.7434,         7.7434,
            8.7660,         8.7660,         8.7660,
            8.7660,         8.7660,         8.7660,
            12.325,         12.325,         12.325,
            12.325,         12.325,         12.325,
            17.270,         17.270,         17.270,
            17.270,         17.270,         17.270,
            ss_max,         ss_max,         ss_max,
            ss_max,         ss_max,         ss_max
        };
        double pT_a [] = {
            0.05,           0.65,           0.85,
            0.05,           0.65,           0.85,
            0.05,           0.75,           0.9,
            0.05,           0.75,           0.9,
            0.05,           0.75,           1.05,
            0.05,           0.75,           1.05,
            0.05,           0.75,           0.95,
            0.05,           0.75,           0.95,
            0.05,           0.75,           0.95,
            0.05,           0.75,           0.95
        };
        double xR_a [] = {
            0.27636,        0.335721,       0.372373,
            0.394076,       0.478722,       0.387318,
            0.246704,       0.303218,       0.319588,
            0.338215,       0.432373,       0.350582,
            0.16069,        0.205427,       0.240819,
            0.376127,       0.480841,       0.300763,
            0.112001,       -0.143182,      -0.159161,
            0.380884,       0.486922,       0.372547,
            0.112001,       -0.143182,      -0.159161,
            0.380884,       0.486922,       0.372547
        };
        int n = 2*5*3;
        
        double p [3] = { ss, pT, xR };
        for (int ia = 0; ia<n; ia++) {
            for (int ib = ia+1; ib<n; ib++) {
                for (int ic = ib+1; ic<n; ic++) {
                    for (int id = ic+1; id<n; id++) {
                        double a [3] =  {ss_a[ia], pT_a[ia], xR_a[ia] };
                        double b [3] =  {ss_a[ib], pT_a[ib], xR_a[ib] };
                        double c [3] =  {ss_a[ic], pT_a[ic], xR_a[ic] };
                        double d [3] =  {ss_a[id], pT_a[id], xR_a[id] };
                        
                        if (inside(a, b, c, d, p)){
                            return true;
                        };
                    }
                }
            }
        }
        if (fPP==1) {
            for (int ia = 0; ia<n; ia++) {
                for (int ib = ia+1; ib<n; ib++) {
                    for (int ic = ib+1; ic<n; ic++) {
                        for (int id = ic+1; id<n; id++) {
                            double a [3] =  {ss_a[ia], pT_a[ia], -xR_a[ia] };
                            double b [3] =  {ss_a[ib], pT_a[ib], -xR_a[ib] };
                            double c [3] =  {ss_a[ic], pT_a[ic], -xR_a[ic] };
                            double d [3] =  {ss_a[id], pT_a[id], -xR_a[id] };
                            
                            if (inside(a, b, c, d, p)){
                                return true;
                            };
                        }
                    }
                }
            }
        }
    }else if( ss_max<0 && ss_min>0 ){
        double ss_a [] = {
            ss_min,         ss_min,         ss_min,
            ss_min,         ss_min,         ss_min,
            7.7434,         7.7434,         7.7434,
            7.7434,         7.7434,         7.7434,
            8.7660,         8.7660,         8.7660,
            8.7660,         8.7660,         8.7660,
            12.325,         12.325,         12.325,
            12.325,         12.325,         12.325,
            17.270,         17.270,         17.270,
            17.270,         17.270,         17.270
        };
        double pT_a [] = {
            0.05,           0.65,           0.85,
            0.05,           0.65,           0.85,
            0.05,           0.65,           0.85,
            0.05,           0.65,           0.85,
            0.05,           0.75,           0.9,
            0.05,           0.75,           0.9,
            0.05,           0.75,           1.05,
            0.05,           0.75,           1.05,
            0.05,           0.75,           0.95,
            0.05,           0.75,           0.95
        };
        double xR_a [] = {
            0.27636,        0.335721,       0.372373,
            0.394076,       0.478722,       0.387318,
            0.27636,        0.335721,       0.372373,
            0.394076,       0.478722,       0.387318,
            0.246704,       0.303218,       0.319588,
            0.338215,       0.432373,       0.350582,
            0.16069,        0.205427,       0.240819,
            0.376127,       0.480841,       0.300763,
            0.112001,       -0.143182,      -0.159161,
            0.380884,       0.486922,       0.372547
        };
        int n = 2*5*3;
        
        double p [3] = { ss, pT, xR };
        for (int ia = 0; ia<n; ia++) {
            for (int ib = ia+1; ib<n; ib++) {
                for (int ic = ib+1; ic<n; ic++) {
                    for (int id = ic+1; id<n; id++) {
                        double a [3] =  {ss_a[ia], pT_a[ia], xR_a[ia] };
                        double b [3] =  {ss_a[ib], pT_a[ib], xR_a[ib] };
                        double c [3] =  {ss_a[ic], pT_a[ic], xR_a[ic] };
                        double d [3] =  {ss_a[id], pT_a[id], xR_a[id] };
                        
                        if (inside(a, b, c, d, p)){
                            return true;
                        };
                    }
                }
            }
        }
        if (fPP==1) {
            for (int ia = 0; ia<n; ia++) {
                for (int ib = ia+1; ib<n; ib++) {
                    for (int ic = ib+1; ic<n; ic++) {
                        for (int id = ic+1; id<n; id++) {
                            double a [3] =  {ss_a[ia], pT_a[ia], -xR_a[ia] };
                            double b [3] =  {ss_a[ib], pT_a[ib], -xR_a[ib] };
                            double c [3] =  {ss_a[ic], pT_a[ic], -xR_a[ic] };
                            double d [3] =  {ss_a[id], pT_a[id], -xR_a[id] };
                            
                            if (inside(a, b, c, d, p)){
                                return true;
                            };
                        }
                    }
                }
            }
        }
        
    }else if( ss_max>0 && ss_min>0 ){
        double ss_a [] = {
            ss_min,         ss_min,         ss_min,
            ss_min,         ss_min,         ss_min,
            7.7434,         7.7434,         7.7434,
            7.7434,         7.7434,         7.7434,
            8.7660,         8.7660,         8.7660,
            8.7660,         8.7660,         8.7660,
            12.325,         12.325,         12.325,
            12.325,         12.325,         12.325,
            17.270,         17.270,         17.270,
            17.270,         17.270,         17.270,
            ss_max,         ss_max,         ss_max,
            ss_max,         ss_max,         ss_max
        };
        double pT_a [] = {
            0.05,           0.65,           0.85,
            0.05,           0.65,           0.85,
            0.05,           0.65,           0.85,
            0.05,           0.65,           0.85,
            0.05,           0.75,           0.9,
            0.05,           0.75,           0.9,
            0.05,           0.75,           1.05,
            0.05,           0.75,           1.05,
            0.05,           0.75,           0.95,
            0.05,           0.75,           0.95,
            0.05,           0.75,           0.95,
            0.05,           0.75,           0.95
        };
        double xR_a [] = {
            0.27636,        0.335721,       0.372373,
            0.394076,       0.478722,       0.387318,
            0.27636,        0.335721,       0.372373,
            0.394076,       0.478722,       0.387318,
            0.246704,       0.303218,       0.319588,
            0.338215,       0.432373,       0.350582,
            0.16069,        0.205427,       0.240819,
            0.376127,       0.480841,       0.300763,
            0.112001,       -0.143182,      -0.159161,
            0.380884,       0.486922,       0.372547,
            0.112001,       -0.143182,      -0.159161,
            0.380884,       0.486922,       0.372547
        };
        int n = 2*6*3;
        
        double p [3] = { ss, pT, xR };
        for (int ia = 0; ia<n; ia++) {
            for (int ib = ia+1; ib<n; ib++) {
                for (int ic = ib+1; ic<n; ic++) {
                    for (int id = ic+1; id<n; id++) {
                        double a [3] =  {ss_a[ia], pT_a[ia], xR_a[ia] };
                        double b [3] =  {ss_a[ib], pT_a[ib], xR_a[ib] };
                        double c [3] =  {ss_a[ic], pT_a[ic], xR_a[ic] };
                        double d [3] =  {ss_a[id], pT_a[id], xR_a[id] };
                        
                        if (inside(a, b, c, d, p)){
                            return true;
                        };
                    }
                }
            }
        }
        if (fPP==1) {
            for (int ia = 0; ia<n; ia++) {
                for (int ib = ia+1; ib<n; ib++) {
                    for (int ic = ib+1; ic<n; ic++) {
                        for (int id = ic+1; id<n; id++) {
                            double a [3] =  {ss_a[ia], pT_a[ia], -xR_a[ia] };
                            double b [3] =  {ss_a[ib], pT_a[ib], -xR_a[ib] };
                            double c [3] =  {ss_a[ic], pT_a[ic], -xR_a[ic] };
                            double d [3] =  {ss_a[id], pT_a[id], -xR_a[id] };
                            
                            if (inside(a, b, c, d, p)){
                                return true;
                            };
                        }
                    }
                }
            }
        }
        
    }
    
    
    
    
    
    return false;
    
}
bool in_na49( double ss, double pT, double xR, double ss_lower = 15, double ss_upper = 20 ){
    
    if ( ss< ss_lower || ss>ss_upper)
        return false;
    if ( pT< 0.05 || pT>1.05)
        return false;
    if ( xR< -0.2 || xR>0.5)
        return false;
    
    double ss_a [] = {
        ss_lower,       ss_lower,       ss_lower,
        ss_lower,       ss_lower,       ss_lower,
        ss_upper,       ss_upper,       ss_upper,
        ss_upper,       ss_upper,       ss_upper
    };
    double pT_a [] = {
        0.1,            1.1,            1.5,
        0.1,            1.1,            1.5,
        0.1,            1.1,            1.5,
        0.1,            1.1,            1.5
    };
    double xR_a [] = {
        -0.123075,      -0.178966,      -0.216007,
        0.424685,       0.444116,       0.331053,
        -0.123075,      -0.178966,      -0.216007,
        0.424685,       0.444116,       0.331053
    };
    
    int n = 2*6;
    
    
    
    double p [3] = { ss, pT, xR };
    
    for (int ia = 0; ia<n; ia++) {
        for (int ib = ia+1; ib<n; ib++) {
            for (int ic = ib+1; ic<n; ic++) {
                for (int id = ic+1; id<n; id++) {
                    double a [3] =  {ss_a[ia], pT_a[ia], xR_a[ia] };
                    double b [3] =  {ss_a[ib], pT_a[ib], xR_a[ib] };
                    double c [3] =  {ss_a[ic], pT_a[ic], xR_a[ic] };
                    double d [3] =  {ss_a[id], pT_a[id], xR_a[id] };
                    
                    if (inside(a, b, c, d, p)){
                        return true;
                    };
                }
            }
        }
    }
    if (fPP==1) {
        for (int ia = 0; ia<n; ia++) {
            for (int ib = ia+1; ib<n; ib++) {
                for (int ic = ib+1; ic<n; ic++) {
                    for (int id = ic+1; id<n; id++) {
                        double a [3] =  {ss_a[ia], pT_a[ia], -xR_a[ia] };
                        double b [3] =  {ss_a[ib], pT_a[ib], -xR_a[ib] };
                        double c [3] =  {ss_a[ic], pT_a[ic], -xR_a[ic] };
                        double d [3] =  {ss_a[id], pT_a[id], -xR_a[id] };
                        
                        if (inside(a, b, c, d, p)){
                            return true;
                        };
                    }
                }
            }
        }
    }
    
    return false;
    
}

bool in_na49_pC( double ss, double pT, double xR, double ss_lower = 15, double ss_upper = 20 ){
    
    if ( ss< ss_lower || ss>ss_upper)
        return false;
    if ( pT< 0.05 || pT>1.05)
        return false;
    if ( xR< -0.2 || xR>0.5)
        return false;
    
    double ss_a [] = {
        ss_lower,       ss_lower,       ss_lower,       ss_lower,
        ss_lower,       ss_lower,       ss_lower,       ss_lower,
        ss_upper,       ss_upper,       ss_upper,       ss_upper,
        ss_upper,       ss_upper,       ss_upper,       ss_upper
    };
    double pT_a [] = {
        0.1,            1.1,            1.3,            1.5,
        0.1,            1.1,            1.3,            1.5,
        0.1,            1.1,            1.3,            1.5,
        0.1,            1.1,            1.3,            1.5
    };
    double xR_a [] = {
        -0.233417,      -0.267141,      -0.279494,      0.209849,
        0.327003,       0.351869,       0.279494,       0.293249,
        -0.233417,      -0.267141,      -0.279494,      0.209849,
        0.327003,       0.351869,       0.279494,       0.293249
    };
    
    int n = 2*8;
    
    
    
    double p [3] = { ss, pT, xR };
    
    for (int ia = 0; ia<n; ia++) {
        for (int ib = ia+1; ib<n; ib++) {
            for (int ic = ib+1; ic<n; ic++) {
                for (int id = ic+1; id<n; id++) {
                    double a [3] =  {ss_a[ia], pT_a[ia], xR_a[ia] };
                    double b [3] =  {ss_a[ib], pT_a[ib], xR_a[ib] };
                    double c [3] =  {ss_a[ic], pT_a[ic], xR_a[ic] };
                    double d [3] =  {ss_a[id], pT_a[id], xR_a[id] };
                    
                    if (inside(a, b, c, d, p)){
                        return true;
                    };
                }
            }
        }
    }
    if (fPP==1) {
        for (int ia = 0; ia<n; ia++) {
            for (int ib = ia+1; ib<n; ib++) {
                for (int ic = ib+1; ic<n; ic++) {
                    for (int id = ic+1; id<n; id++) {
                        double a [3] =  {ss_a[ia], pT_a[ia], -xR_a[ia] };
                        double b [3] =  {ss_a[ib], pT_a[ib], -xR_a[ib] };
                        double c [3] =  {ss_a[ic], pT_a[ic], -xR_a[ic] };
                        double d [3] =  {ss_a[id], pT_a[id], -xR_a[id] };
                        
                        if (inside(a, b, c, d, p)){
                            return true;
                        };
                    }
                }
            }
        }
    }
    
    return false;
}


static double inv_pbar_CM__exp( double s, double E_pbar, double pT_pbar){
    
    double ret          =   0;
    if (fCS_frac==0)
        ret          =  inv_A1A2_N1N2_pbar_CM(s, E_pbar, pT_pbar);
    if (fCS_frac==1){
        double factor  = 1.;
        double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt(s);
        double x_R          =   E_pbar/E_pbar_Max;
        
        if (fExperiment=="lhcb_pHe") {
            if (!in_lhcb_combined(sqrt(s), pT_pbar, x_R, 100, 120)){
                factor = 0;
            }
        }else if (fExperiment=="lhcb_pHe_87") {
            if (!in_lhcb_combined(sqrt(s), pT_pbar, x_R, 80, 100)){
                factor = 0;
            }
        }else if (fExperiment=="lhcb_pHe_43") {
            if (!in_lhcb_combined(sqrt(s), pT_pbar, x_R, 40, 50)){
                factor = 0;
            }
        }else if (fExperiment=="lhcb_pHe_combined") {
            bool in = false;
            if (in_lhcb_combined(sqrt(s), pT_pbar, x_R, 40, 50)){
                in=true;
            }else if (in_lhcb_combined(sqrt(s), pT_pbar, x_R, 50, 80)){
                in=true;
            }else if (in_lhcb_combined(sqrt(s), pT_pbar, x_R, 80, 100)){
                in=true;
            }else if (in_lhcb_combined(sqrt(s), pT_pbar, x_R, 100, 120)){
                in=true;
            }
            if(!in){
                factor = 0;
            }
        }else if (fExperiment=="lhcb_Hep") {
            if (!in_lhcb_combined(sqrt(s), pT_pbar, -x_R, 100, 120)){
                factor = 0;
            }
        }else if (fExperiment=="lhcb_Hep_87") {
            if (!in_lhcb_combined(sqrt(s), pT_pbar, -x_R, 80, 100)){
                factor = 0;
            }
        }else if (fExperiment=="lhcb_Hep_43") {
            if (!in_lhcb_combined(sqrt(s), pT_pbar, -x_R, 40, 50)){
                factor = 0;
            }
        }else if (fExperiment=="lhcb_Hep_combined") {
            bool in = false;
            if (in_lhcb_combined(sqrt(s), pT_pbar, -x_R, 40, 50)){
                in=true;
            }else if (in_lhcb_combined(sqrt(s), pT_pbar, -x_R, 50, 80)){
                in=true;
            }else if (in_lhcb_combined(sqrt(s), pT_pbar, -x_R, 80, 100)){
                in=true;
            }else if (in_lhcb_combined(sqrt(s), pT_pbar, -x_R, 100, 120)){
                in=true;
            }
            if(!in){
                factor = 0;
            }
        }else if (fExperiment=="na61"){
            if (!in_na61(sqrt(s), pT_pbar, x_R)){
                factor = 0;
            }
        }else if (fExperiment=="na61_large"){
            if (!in_na61(sqrt(s), pT_pbar, x_R, 50)){
                factor = 0;
            }
        }else if (fExperiment=="na61_low"){
            if (!in_na61(sqrt(s), pT_pbar, x_R, -1, 6.3)){
                factor = 0;
            }
        }else if (fExperiment=="na61_all"){
            if (!in_na61(sqrt(s), pT_pbar, x_R, 50, 6.3)){
                factor = 0;
            }
        }else if (fExperiment=="na49"){
            if (!in_na49(sqrt(s), pT_pbar, x_R)){
                factor = 0;
            }
        }else if (fExperiment=="na49_large"){
            if (!in_na49(sqrt(s), pT_pbar, x_R, 10, 50)){
                factor = 0;
            }
        }else if (fExperiment=="na49_pC"){
            if (!in_na49_pC(sqrt(s), pT_pbar, x_R)){
                factor = 0;
            }
        }else if (fExperiment=="na49_Cp"){
            if (!in_na49_pC(sqrt(s), pT_pbar, -x_R)){
                factor = 0;
            }
        }
        
        ret          =   factor*inv_A1A2_N1N2_pbar_CM(s, E_pbar, pT_pbar);
    }
    return ret;
}
static double inv_pbar_CM_tot__exp( double s, double E_pbar_d, double pT_pbar){
    return inv_pbar_CM__exp( s, E_pbar_d, pT_pbar)*(2+ 2*deltaHyperon(s) + 2*deltaIsospin(s) );
}

static double dE_pbar_LAB__exp(double E_p_LAB, double E_pbar_LAB, int precision=fPrecision, int output=NO_OUT){
    return CSTransformations::dE_pp_product_LAB_intEta (E_p_LAB, E_pbar_LAB, fMass_proton,
                                                        inv_pbar_CM__exp,
                                                        precision,
                                                        output);
}
static double dE_pbar_LAB_tot__exp(double E_p_LAB, double E_pbar_LAB, int precision=fPrecision, int output=NO_OUT){
    return CSTransformations::dE_pp_product_LAB_intEta (E_p_LAB, E_pbar_LAB, fMass_proton,
                                                        inv_pbar_CM_tot__exp,
                                                        precision,
                                                        output);
}
static double dT_pbar_LAB__exp  ( double T_p, double T_pbar ){
    return dE_pbar_LAB__exp( T_p+fMass_proton,   T_pbar+fMass_proton)*1e-31;
}
static double dT_pbar_LAB_tot__exp  ( double T_p, double T_pbar ){
    return dE_pbar_LAB_tot__exp( T_p+fMass_proton,   T_pbar+fMass_proton)*1e-31;
}




//***********************************************************************************************************************//
//***********************************************************************************************************************//
//***********************************************************************************************************************//
//***********************************************************************************************************************//




/***********************************************************************************************************************/
/************************************           Main program        ****************************************************/
/***********************************************************************************************************************/
int main(int argc, char *argv[])
{
    
    // Config handling
    
    ConfigHandler* config = ConfigHandler::GetInstance();
    config->SetInput(argc, argv, fProgramDescription);
    
    int cs = 0;
    int n_point = 0;
    config->AddOptionInt        (  "step",          fStep,         "default: 1"                          );
    config->AddOption           (  "dir",           fDir,          "default: ."                          );
    config->AddOptionInt        (  "scaling_asymm", fScaling_asymm,"0: constant, 1: linear, 2: powerlaw. default: 0");
    config->AddOptionInt        (  "cs",            cs            ,"0: diMauro 12, 1: Winkler,           default: 0");
    config->AddOption           (  "exp",           fExperiment,   "default: lhcb"                          );
    config->AddOptionInt        (  "npoint",        n_point       ,"Point from the sigma file            default: 0");
    
    config->AddOptionInt        (  "A1",        fA1       ,"Mass number of the projectile            default: 1");
    config->AddOptionInt        (  "A2",        fA2       ,"Mass number of the trarget               default: 1");
    
    
    
    if ( cs == 0 ){
        fCS = CS::DI_MAURO12;
    } else if ( cs==1 ){
        fCS = CS::WINKLER;
    } else if ( cs==-1 ){
        fCS = -1;
    }
    
    if (fScaling_asymm==CONSTANT){
        fNDparam = 1;
    }else if (fScaling_asymm==PT_RESCALING){
        fNDparam = 3;
    }else{
        fNDparam = 2;
    }
    
    SetStart();

    
    std::string include_data = "na49,na61,brahms";
    config->AddOption("data", include_data, "Use this to specify experiments, e.g. '--data na49,phenix', possible: allaby,antreasyan,brahms,capiluppi,dekkers,guettler,johnson,na49,na61,phenix or all,  Default: na49,na61,brahms");
    bool no_hyperon;
    config->AddOptionTrue("doNotCorrectHyperon", no_hyperon);
    
    std::string include_data_pHe = "";
    config->AddOption("data_pHe", include_data_pHe, "Use this to specify experiments, e.g. '--data na49,phenix', Default: lhcb");
    bool no_hyperon_pHe;
    config->AddOptionTrue("doNotCorrectHyperon_pHe", no_hyperon_pHe);
    
    std::string include_data_pBe = "";
    config->AddOption("data_pBe", include_data_pBe, "Use this to specify experiments, e.g. '--data na49', Default: ");
    bool no_hyperon_pBe;
    config->AddOptionTrue("doNotCorrectHyperon_pBe", no_hyperon_pBe);
    
    std::string include_data_pC = "na49";
    config->AddOption("data_pC", include_data_pC, "Use this to specify experiments, e.g. '--data na49', Default: na49");
    bool no_hyperon_pC;
    config->AddOptionTrue("doNotCorrectHyperon_pC", no_hyperon_pC);
    
    config->CheckOption();
    
    out(fStep)
    out("Read")
    
    
    ReadData        (   include_data,       !no_hyperon,
                     fData_pp_pbar_all,
                     "pp",
                     ",allaby,antreasyan,brahms,capiluppi,dekkers,guettler,johnson,phenix,na49,na61,",
                     ",allaby,antreasyan,brahms,capiluppi,dekkers,guettler,johnson,phenix,",
                     1  );
    ReadData        (   include_data_pHe,   !no_hyperon_pHe,
                     fData_pHe_pbar_all,
                     "pHe",
                     ",lhcb,",
                     ",",
                     4  );
    ReadData        (   include_data_pBe,   !no_hyperon_pBe,
                     fData_pBe_pbar_all,
                     "pBe",
                     ",veronin,",
                     ",veronin,",
                     9  );
    ReadData        (   include_data_pC,    !no_hyperon_pC,
                     fData_pC_pbar_all,
                     "pC",
                     ",na49,abramov,amann,barton,sugaya,",
                     ",abramov,amann,barton,sugaya,",
                     12  );
    
    // merge all pA data
    for (int i = 0; i<fData_pHe_pbar_all.size(); i++) {
        data & d = fData_pHe_pbar_all.at(i);
        fData_pA_all.push_back(d);
    }
    for (int i = 0; i<fData_pBe_pbar_all.size(); i++) {
        data & d = fData_pBe_pbar_all.at(i);
        fData_pA_all.push_back(d);
    }
    for (int i = 0; i<fData_pC_pbar_all.size(); i++) {
        data & d = fData_pC_pbar_all.at(i);
        fData_pA_all.push_back(d);
    }
    
    std::cout << std::endl;
    std::cout << "*****************************************************************************************************" << std::endl;
    std::cout << "***                                                                                               ***" << std::endl;
    std::cout << "***                      Fit of the pp XS                                                         ***" << std::endl;
    std::cout << "***                                                                                               ***" << std::endl;
    std::cout << "*****************************************************************************************************" << std::endl;
    std::cout << std::endl;
    
    Ifit_pp();
    
    //* Make a prefit with TMinuit for pp
    //*********************************************************************************************
    if (fStep==1) {
        std::system( (config->SoftwarePath()+"/scripts/"+"CRACS_crossSectionFitter.py --step 10 --dir " + fDir).c_str() );
        
        FileTool f_parameters = FileTool(fDir+"/parameters_pp.txt");
        f_parameters.ExtractNumberTable(fNCparam+fNWparam, " ", true);
        
        std::string write_1s = "";
        std::string write_2s = "";
        for(int n_point=0; n_point<f_parameters.NumberTableGetNrows(); n_point++){
            
            std::string counter = QString("%1").arg(n_point).rightJustified(6, '0').toStdString();
            
            double par[fNCparam+fNWparam];
            for (int jj=0; jj<fNCparam+fNWparam; jj++) {
                par[jj] = f_parameters.NumberTable(n_point, jj);
            }
            int n;
            double h, chisq;
            fPrint = 0;
            fcn_pp_pbar_CM(n, &h, chisq, par, 1);
            
            int npar = fNCparam+fNWparam;
            if(fCS==1){
                npar --;
            }
            double diff = chisq-fChiSq_min_pp;
          
            if(diff<fSigma_1_region[npar-1]){
                write_1s += f_parameters.GetLine(n_point) + "\n";
            }
            else if (diff<fSigma_2_region[npar-1]){
                write_2s += f_parameters.GetLine(n_point) + "\n";
            }
            FileTool::WriteStringToFile(write_1s, "parameters_1_sigma.txt");
            FileTool::WriteStringToFile(write_2s, "parameters_2_sigma.txt");
        }
        Write_results( fData_pp_pbar_all );
        std::system( (config->SoftwarePath()+"/scripts/"+"CRACS_crossSectionFitter.py --step 1 --dir " + fDir).c_str() );
        return 1;
    }
    
    //* TMinuit to get the pA factor
    //*********************************************************************************************

    std::cout << std::endl;
    std::cout << "*****************************************************************************************************" << std::endl;
    std::cout << "***                                                                                               ***" << std::endl;
    std::cout << "***                      Fit of the pA factor                                                     ***" << std::endl;
    std::cout << "***                                                                                               ***" << std::endl;
    std::cout << "*****************************************************************************************************" << std::endl;
    std::cout << std::endl;
    Ifit_pA();

    std::cout << std::endl;
    std::cout << "*****************************************************************************************************" << std::endl;
    std::cout << "***                                                                                               ***" << std::endl;
    std::cout << "***                      Info of LHCb compatibility                                               ***" << std::endl;
    std::cout << "***                                                                                               ***" << std::endl;
    std::cout << "*****************************************************************************************************" << std::endl;
    std::cout << std::endl;
    Print_chiSq_pA( fData_pHe_pbar_all[0] );

    std::cout << std::endl;
    std::cout << "*****************************************************************************************************" << std::endl;
    std::cout << "***                                                                                               ***" << std::endl;
    std::cout << "***                      Info of NA49 compatibility                                               ***" << std::endl;
    std::cout << "***                                                                                               ***" << std::endl;
    std::cout << "*****************************************************************************************************" << std::endl;
    std::cout << std::endl;
    Print_chiSq_pA( fData_pC_pbar_all[0] );


    if (fStep==2) {

        std::system( (config->SoftwarePath()+"/scripts/"+"CRACS_crossSectionFitter.py --step 20  --dir " + fDir).c_str() );

        FileTool f_parameters = FileTool(fDir+"/parameters_pA.txt");
        f_parameters.ExtractNumberTable(fNCparam+fNWparam+fNDparam+fNWparam_pA, " ", true);

        std::string write_1s = "";
        std::string write_2s = "";
        for(int n_point=0; n_point<f_parameters.NumberTableGetNrows(); n_point++){

            std::string counter = QString("%1").arg(n_point).rightJustified(6, '0').toStdString();

            double par[fNCparam+fNWparam];
            for (int jj=0; jj<fNCparam+fNWparam; jj++) {
                par[jj] = f_parameters.NumberTable(n_point, jj);
            }
            int n;
            double h, chisq_pp, chisq_pA;
            fPrint = 0;
            fcn_pp_pbar_CM(n, &h, chisq_pp, par, 1);

            double par_pA[fNDparam+fNWparam_pA];
            for (int jj=0; jj<fNDparam+fNWparam_pA; jj++) {
                par_pA[jj] = f_parameters.NumberTable(n_point, jj+fNCparam+fNWparam);
            }
            fPrint = 0;
            fcn_pA_pbar_CM(n, &h, chisq_pA, par_pA, 1);

            int npar = fNCparam+fNWparam+fNDparam+fNWparam_pA;
            if(fCS==1){
                npar --;
            }
            double diff = chisq_pp+chisq_pA-fChiSq_min_pp -fChiSq_min_pA ;

            if(diff<fSigma_1_region[npar-1]){
                write_1s += f_parameters.GetLine(n_point) + "\n";
            }
            else if (diff<fSigma_2_region[npar-1]){
                write_2s += f_parameters.GetLine(n_point) + "\n";
            }
            FileTool::WriteStringToFile(write_1s, "parameters_pA_1_sigma.txt");
            FileTool::WriteStringToFile(write_2s, "parameters_pA_2_sigma.txt");
        }

        Write_results( fData_pA_all );
        std::system( (config->SoftwarePath()+"/scripts/"+"CRACS_crossSectionFitter.py --step 2 --dir " + fDir).c_str() );
        return 1;


        return 1;
    }


    std::cout << std::endl;
    std::cout << "*****************************************************************************************************" << std::endl;
    std::cout << "***                                                                                               ***" << std::endl;
    std::cout << "***                      XS and source term calculations                                          ***" << std::endl;
    std::cout << "***                                                                                               ***" << std::endl;
    std::cout << "*****************************************************************************************************" << std::endl;
    std::cout << std::endl;
    
    
    
    
    
    if (fStep==3) {
        
        out(fA1)
        out(fA2)
        
        std::string AA = QString("%1-%2").arg(fA1).arg(fA2).toStdString();
        std::system( ("mkdir "+fDir+"/sourceTerm_"+AA).c_str() );
        
        fTab = CS_lab_tab::GetInstance();
        fTab->WriteCS(  dT_A1A2_pbar_LAB,                     "sourceTerm_"+AA+"/dT_"+AA+"_pbar.txt",                   fDir  );
        fTab->WriteCS(  dT_A1A2_pbar_LAB_inclHyp_inclNbar,    "sourceTerm_"+AA+"/dT_"+AA+"_pbar_tot.txt",               fDir  );

        out("upper/lower isospin")
//        fIsospinKind = UPPER; fHyperonKind = MEAN;
//        fTab->WriteCS(  dT_A1A2_pbar_LAB_inclHyp_inclNbar,    "sourceTerm_"+AA+"/dT_"+AA+"_pbar_tot_upperIsoSpin.txt",  fDir  );
//        fIsospinKind = LOWER; fHyperonKind = MEAN;
//        fTab->WriteCS(  dT_A1A2_pbar_LAB_inclHyp_inclNbar,    "sourceTerm_"+AA+"/dT_"+AA+"_pbar_tot_lowerIsoSpin.txt",  fDir  );

        out("upper/lower isospin + hyperon")
        fIsospinKind = UPPER; fHyperonKind = UPPER;
        fTab->WriteCS(  dT_A1A2_pbar_LAB_inclHyp_inclNbar,    "sourceTerm_"+AA+"/dT_"+AA+"_pbar_tot_upper.txt",         fDir  );
        fIsospinKind = LOWER; fHyperonKind = LOWER;
        fTab->WriteCS(  dT_A1A2_pbar_LAB_inclHyp_inclNbar,    "sourceTerm_"+AA+"/dT_"+AA+"_pbar_tot_lower.txt",         fDir  );

        writeSourceTerm( fA1, fA2, "_pbar"                  , true );
        writeSourceTerm( fA1, fA2, "_pbar_tot"              );
//        writeSourceTerm( fA1, fA2, "_pbar_tot_upperIsoSpin" );
//        writeSourceTerm( fA1, fA2, "_pbar_tot_lowerIsoSpin" );
        writeSourceTerm( fA1, fA2, "_pbar_tot_upper"        );
        writeSourceTerm( fA1, fA2, "_pbar_tot_lower"        );
        
        
        
        std::string cmd = config->SoftwarePath()+"/scripts/"+QString("CRACS_crossSectionFitter.py --step 30 --cs %1 --data %2 --data_pHe \"%3\" --data_pC \"%4\" --A1 %5 --A2 %6 --dir %7 --sdir %8")
        .arg(cs)
        .arg(include_data.c_str())
        .arg(include_data_pHe.c_str())
        .arg(include_data_pC.c_str())
        .arg(fA1)
        .arg(fA2)
        .arg(fDir.c_str())
        .arg(config->SoftwarePath().c_str())
        .toStdString();
        
        //std::system( cmd.c_str() );
        
        return 1;
    }
    
    
    if (fStep==30) {
        std::cout << std::endl;
        std::cout << "*****************************************************************************************************" << std::endl;
        std::cout << "***                                                                                               ***" << std::endl;
        std::cout << "***                      XS and source term calculations                                          ***" << std::endl;
        std::cout << "***                                                                                               ***" << std::endl;
        std::cout << "*****************************************************************************************************" << std::endl;
        std::cout << std::endl;
        
        std::string AA = QString("%1-%2").arg(fA1).arg(fA2).toStdString();
        std::system( ("mkdir "+fDir+"/sourceTerm_"+AA).c_str() );
        std::string counter = QString("%1").arg(n_point).rightJustified(6, '0').toStdString();
        
        FileTool f_1_sigma = FileTool(fDir+"/parameters_pA_1_sigma.txt");
        f_1_sigma.ExtractNumberTable(fNCparam+fNWparam+fNDparam+fNWparam_pA, " ", true);
        for (int jj=0; jj<fNCparam+fNWparam; jj++) {
            fParameters[jj] = f_1_sigma.NumberTable(n_point, jj);
        }
        for (int jj=0; jj<fNDparam+fNWparam_pA; jj++) {
            fParameters_pA[jj] = f_1_sigma.NumberTable(n_point, jj+fNCparam+fNWparam);
        }
        fTab = CS_lab_tab::GetInstance();
        fTab->WriteCS(  dT_A1A2_pbar_LAB,                     "sourceTerm_"+AA+"/dT_"+AA+"_pbar_1s_"+counter+".txt",                   fDir  );
        
        writeSourceTerm( fA1, fA2, "_pbar_1s_"+counter, true );
        
        return 1;
    }
    
    
    if (fStep==4){
        
        fPrecision=1000;
        
        std::system( ("mkdir "+fDir+"/sourceTerms_publish").c_str() );
        
        fTab = CS_lab_tab::GetInstance();
        std::string AA="";
        
        fA1=1;
        fN1=0;
        fA2=1;
        fN2=0;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=1;
        fN1=0;
        fA2=4;
        fN2=2;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=2;
        fN1=1;
        fA2=1;
        fN2=0;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=2;
        fN1=1;
        fA2=4;
        fN2=2;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=3;
        fN1=1;
        fA2=1;
        fN2=0;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=3;
        fN1=1;
        fA2=4;
        fN2=2;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=4;
        fN1=2;
        fA2=1;
        fN2=0;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=4;
        fN1=2;
        fA2=4;
        fN2=2;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=12;
        fN1=6;
        fA2=1;
        fN2=0;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=12;
        fN1=6;
        fA2=4;
        fN2=2;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=13;
        fN1=7;
        fA2=1;
        fN2=0;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=13;
        fN1=7;
        fA2=4;
        fN2=2;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=14;
        fN1=7;
        fA2=1;
        fN2=0;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=14;
        fN1=7;
        fA2=4;
        fN2=2;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=15;
        fN1=8;
        fA2=1;
        fN2=0;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=15;
        fN1=8;
        fA2=4;
        fN2=2;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=16;
        fN1=8;
        fA2=1;
        fN2=0;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=16;
        fN1=8;
        fA2=4;
        fN2=2;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=17;
        fN1=9;
        fA2=1;
        fN2=0;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=17;
        fN1=9;
        fA2=4;
        fN2=2;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=18;
        fN1=10;
        fA2=1;
        fN2=0;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        fA1=18;
        fN1=10;
        fA2=4;
        fN2=2;
        AA = QString("%1-%2_%3-%4").arg(fA1).arg(fN1).arg(fA2).arg(fN2).toStdString();
        out(AA)
        fTab->WriteCS(  dT_A1A2_N1N2_pbar_LAB_inclHyp_inclNbar,    "sourceTerms_publish/dT_"+AA+"_pbar_tot.txt",               fDir  );
        
        
    }
    
    
    if (fStep==100){
        
        
        double nH   = 1e+6;
        double nHe  = nH*0.1;
        double nC   = nH*0.000141254;
        
        double nISM = nH;
        
        double (*flux)(double) =protonFlux;
        
        fA1=1; fA2=1; fN1=0; fN2=0;
        fPP=1;
        if( QString(fExperiment.c_str()).startsWith("lhcb_Hep") ){
            fA1=4; fA2=1; fN1=2; fN2=0;
            fPP=0;
            flux = heliumFlux;
            nISM = nH;
        }
        if( QString(fExperiment.c_str()).startsWith("na49_Cp") ){
            fA1=12; fA2=1; fN1=6; fN2=0;
            fPP=0;
            flux=carbonFlux;
            nISM = nH;
        }
        if( QString(fExperiment.c_str()).startsWith("lhcb_pHe") ){
            fA1=1; fA2=4; fN1=0; fN2=2;
            fPP=0;
            flux = protonFlux;
            nISM = nHe;
        }
        if( QString(fExperiment.c_str()).startsWith("na49_pC") ){
            fA1=1; fA2=12; fN1=0; fN2=6;
            fPP=0;
            flux = protonFlux;
            nISM = nC;
        }
        out(fA1)
        out(fA2)
        out(fN1)
        out(fN2)
        out(fPP)
        
        std::string exp = fExperiment;
        
        std::cout << std::endl;
        std::cout << "*****************************************************************************************************" << std::endl;
        std::cout << "***                                                                                               ***" << std::endl;
        std::cout << "***                      XS and source term calculations (from LCHb phase space)                  ***" << std::endl;
        std::cout << "***                                                                                               ***" << std::endl;
        std::cout << "*****************************************************************************************************" << std::endl;
        std::cout << std::endl;
        
        LISFluxes::fit();
        std::system( ("mkdir "+fDir+"/"+exp).c_str() );
        
        fTab = CS_lab_tab::GetInstance();
        
        fTab->Set_dLog(100);
        
        fCS_frac = 0;
        fTab->WriteCS(  dT_pbar_LAB_tot__exp,           exp+"/dT_pbar_LAB.txt",               fDir  );
        fCS_frac = 1;
        fTab->WriteCS(  dT_pbar_LAB_tot__exp,           exp+"/dT_pbar_LAB__exp.txt",          fDir  );
        
        
        
        std::string sourceTerm = "   ";
        fTab->ReadCS (  exp+"/dT_pbar_LAB.txt", fDir  );
        sourceTerm += QString("%1").arg("T_pbar"        ).leftJustified(20).toStdString();
        sourceTerm += QString("%1").arg("q"             ).leftJustified(20).toStdString();
        sourceTerm += "\n";
        for (double dT=-1; dT<4; dT+=1./30) {
            double T_pbar   = pow(10, dT);
            double Tn_from  = std::min(T_pbar, 6*fMass_proton);
            //out(protonFlux(T_pbar))
            double S_p_H       = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,     flux, nISM,  Tn_from, 1e6, 10000, NO_OUT);
            sourceTerm += QString("%1").arg(T_pbar      ).leftJustified(20).toStdString();
            sourceTerm += QString("%1").arg(S_p_H       ).leftJustified(20).toStdString();
            sourceTerm += "\n";
        }
        FileTool::WriteStringToFile(sourceTerm, fDir+"/"+exp+"/source.txt" );
        
        sourceTerm = "   ";
        fTab->ReadCS (  exp+"/dT_pbar_LAB__exp.txt", fDir  );
        sourceTerm += QString("%1").arg("T_pbar"        ).leftJustified(20).toStdString();
        sourceTerm += QString("%1").arg("q"             ).leftJustified(20).toStdString();
        sourceTerm += "\n";
        for (double dT=-1; dT<4; dT+=1./30) {
            double T_pbar   = pow(10, dT);
            double Tn_from  = std::min(T_pbar, 6*fMass_proton);
            //out(protonFlux(T_pbar))
            double S_p_H       = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,     flux, nISM,  Tn_from, 1e6, 10000, NO_OUT);
            sourceTerm += QString("%1").arg(T_pbar      ).leftJustified(20).toStdString();
            sourceTerm += QString("%1").arg(S_p_H       ).leftJustified(20).toStdString();
            sourceTerm += "\n";
        }
        FileTool::WriteStringToFile(sourceTerm, fDir+"/"+exp+"/source__exp.txt" );
        return 1;
        
        
    }
    
    
    fPrecision=10;
    out( fPrecision )
    out(  dT_A1A2_N1N2_pbar_LAB(20, 1)  )
    out(  dT_A1A2_N1N2_pbar_LAB(20, 5)  )
    out(  dT_A1A2_N1N2_pbar_LAB(20,10)  )
    
    out(  dT_A1A2_N1N2_pbar_LAB(100, 1)  )
    out(  dT_A1A2_N1N2_pbar_LAB(100,10)  )
    out(  dT_A1A2_N1N2_pbar_LAB(100,50)  )
    
    out(  dT_A1A2_N1N2_pbar_LAB(500, 10)  )
    out(  dT_A1A2_N1N2_pbar_LAB(500, 50)  )
    out(  dT_A1A2_N1N2_pbar_LAB(500,100)  )
    
    fPrecision=100;
    out( fPrecision )
    out(  dT_A1A2_N1N2_pbar_LAB(20, 1)  )
    out(  dT_A1A2_N1N2_pbar_LAB(20, 5)  )
    out(  dT_A1A2_N1N2_pbar_LAB(20,10)  )
    
    out(  dT_A1A2_N1N2_pbar_LAB(100, 1)  )
    out(  dT_A1A2_N1N2_pbar_LAB(100,10)  )
    out(  dT_A1A2_N1N2_pbar_LAB(100,50)  )
    
    out(  dT_A1A2_N1N2_pbar_LAB(500, 10)  )
    out(  dT_A1A2_N1N2_pbar_LAB(500, 50)  )
    out(  dT_A1A2_N1N2_pbar_LAB(500,100)  )
    
    fPrecision=1000;
    out( fPrecision )
    out(  dT_A1A2_N1N2_pbar_LAB(20, 1)  )
    out(  dT_A1A2_N1N2_pbar_LAB(20, 5)  )
    out(  dT_A1A2_N1N2_pbar_LAB(20,10)  )
    
    out(  dT_A1A2_N1N2_pbar_LAB(100, 1)  )
    out(  dT_A1A2_N1N2_pbar_LAB(100,10)  )
    out(  dT_A1A2_N1N2_pbar_LAB(100,50)  )
    
    out(  dT_A1A2_N1N2_pbar_LAB(500, 10)  )
    out(  dT_A1A2_N1N2_pbar_LAB(500, 50)  )
    out(  dT_A1A2_N1N2_pbar_LAB(500,100)  )
    
    
    return 1;
    
} // END main


void writeSourceTerm( int A1, int A2, std::string suffix, bool withComparisonOld ){
  
    LISFluxes::fit();
    std::string AA = QString("%1-%2").arg(fA1).arg(fA2).toStdString();
    double nH   = 1e+6;
    double nHe  = nH*0.1;
    double nC   = nH*0.000141254;
    double nN   = nH*6.16595E-05;
    double nO   = nH*0.000389045;
    double nNe  = nH*1.47911E-06;       // Estimated to be similar to Mg!
    double nMg  = nH*1.47911E-06;
    double nSi  = nH*2.23872E-06;
    double nFe  = nH*2.5704E-07;
    
    
    
    fTab->ReadCS (  "sourceTerm_"+AA+"/dT_"+AA+suffix + ".txt", fDir  );
    std::string sourceTerm = "   ";
    sourceTerm += QString("%1").arg("T_pbar"        ).leftJustified(20).toStdString();
    sourceTerm += QString("%1").arg(AA.c_str()      ).leftJustified(20).toStdString();
    sourceTerm += "\n";
    for (double dT=-1; dT<4; dT+=1./30) {
        double T_pbar   = pow(10, dT);
        double Tn_from  = std::min(T_pbar, 6*fMass_proton);
        
        double S_p_H       = 0;
        double n = nH;
        if(fA2==1){
            n = nH;
        }else if(fA2==4){
            n = nHe;
        }else if(fA2==12){
            n = nC;
        }else if(fA2==14){
            n = nN;
        }else if(fA2==16){
            n = nO;
        }else if(fA2==20){
            n = nNe;
        }else if(fA2==24){
            n = nMg;
        }else if(fA2==28){
            n = nSi;
        }else if(fA2==56){
            n = nFe;
        }
        
        if(fA1==1){
            S_p_H = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,     protonFlux,      n,  Tn_from, 1e6, 10000, NO_OUT);
        }else if(fA1== 4){
            S_p_H = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,     heliumFlux,      n,  Tn_from, 1e6, 10000, NO_OUT);
        }else if(fA1==12){
            S_p_H = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,     carbonFlux,      n,  Tn_from, 1e6, 10000, NO_OUT);
        }else if(fA1==14){
            S_p_H = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,     nitrogenFlux,    n,  Tn_from, 1e6, 10000, NO_OUT);
        }else if(fA1==16){
            S_p_H = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,     oxygenFlux,      n,  Tn_from, 1e6, 10000, NO_OUT);
        }else if(fA1== 7){
            S_p_H = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,     lithiumFlux,     n,  Tn_from, 1e6, 10000, NO_OUT);
        }else if(fA1== 9){
            S_p_H = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,     berylliumFlux,   n,  Tn_from, 1e6, 10000, NO_OUT);
        }else if(fA1==11){
            S_p_H = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,     boronFlux,       n,  Tn_from, 1e6, 10000, NO_OUT);
        }else if(fA1==20){
            S_p_H = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,     neonFlux,        n,  Tn_from, 1e6, 10000, NO_OUT);
        }else if(fA1==24){
            S_p_H = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,     magnesiumFlux,   n,  Tn_from, 1e6, 10000, NO_OUT);
        }else if(fA1==28){
            S_p_H = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,     siliconFlux,     n,  Tn_from, 1e6, 10000, NO_OUT);
        }else if(fA1==56){
            S_p_H = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,     ironFlux,        n,  Tn_from, 1e6, 10000, NO_OUT);
        }
        
        sourceTerm += QString("%1").arg(T_pbar      ).leftJustified(20).toStdString();
        sourceTerm += QString("%1").arg(S_p_H       ).leftJustified(20).toStdString();
        
        if (withComparisonOld){
            if(fA1==1 && fA2==1){
                double S_p_H_DM    = SourceTerms::integrationFromTo(T_pbar, dT_pp_pbar_LAB_dM,              protonFlux, nH,  Tn_from, 1e6, 10000, NO_OUT);
                double S_p_H_KR    = SourceTerms::integrationFromTo(T_pbar, dT_pp_pbar_LAB_KR,              protonFlux, nH,  Tn_from, 1e6, 10000, NO_OUT);
                double S_p_H_WI    = SourceTerms::integrationFromTo(T_pbar, dT_pp_pbar_LAB_Wi,              protonFlux, nH,  Tn_from, 1e6, 10000, NO_OUT);
                
                sourceTerm += QString("%1").arg(S_p_H_DM       ).leftJustified(20).toStdString();
                sourceTerm += QString("%1").arg(S_p_H_KR       ).leftJustified(20).toStdString();
                sourceTerm += QString("%1").arg(S_p_H_WI       ).leftJustified(20).toStdString();
            }
            if(fA1==1 && fA2==4){
                double S_p_H_DM    = SourceTerms::integrationFromTo(T_pbar, dT_pHe_pbar_LAB_dM,             protonFlux, nHe, Tn_from, 1e6, 10000, NO_OUT);
                double S_p_H_KR    = SourceTerms::integrationFromTo(T_pbar, dT_pHe_pbar_LAB_KR,             protonFlux, nHe, Tn_from, 1e6, 10000, NO_OUT);
                double S_p_H_WI    = SourceTerms::integrationFromTo(T_pbar, dT_pHe_pbar_LAB_Wi,             protonFlux, nHe, Tn_from, 1e6, 10000, NO_OUT);
                
                sourceTerm += QString("%1").arg(S_p_H_DM       ).leftJustified(20).toStdString();
                sourceTerm += QString("%1").arg(S_p_H_KR       ).leftJustified(20).toStdString();
                sourceTerm += QString("%1").arg(S_p_H_WI       ).leftJustified(20).toStdString();
            }
            if(fA1==4 && fA2==1){
                double S_p_H_DM    = SourceTerms::integrationFromTo(T_pbar, dT_Hep_pbar_LAB_dM,             heliumFlux, nH,  Tn_from, 1e6, 10000, NO_OUT);
                double S_p_H_KR    = SourceTerms::integrationFromTo(T_pbar, dT_Hep_pbar_LAB_KR,             heliumFlux, nH,  Tn_from, 1e6, 10000, NO_OUT);
                double S_p_H_WI    = SourceTerms::integrationFromTo(T_pbar, dT_Hep_pbar_LAB_Wi,             heliumFlux, nH,  Tn_from, 1e6, 10000, NO_OUT);
                
                sourceTerm += QString("%1").arg(S_p_H_DM       ).leftJustified(20).toStdString();
                sourceTerm += QString("%1").arg(S_p_H_KR       ).leftJustified(20).toStdString();
                sourceTerm += QString("%1").arg(S_p_H_WI       ).leftJustified(20).toStdString();
            }
        }
        
       sourceTerm += "\n";
    }
    FileTool::WriteStringToFile(sourceTerm, fDir+"/sourceTerm_"+AA+"/sourceTerm_"+AA+suffix+".txt" );
    
    
}


//******************************************************************************************************************************/
//  Functions definitions to read and write data
//******************************************************************************************************************************/


void ReadData( std::string include_data, bool correctHyperon,
              std::vector<data> &fData,
              std::string type,
              QString experiments,
              QString hyperon_inclusive_experiments,
              int A){
    
    include_data =  ","+include_data+",";
    
    QStringList e_list   = experiments.                     split(",", QString::SkipEmptyParts);
    QStringList eh_list  = hyperon_inclusive_experiments.   split(" ", QString::SkipEmptyParts);
    
    
    bool all = false;
    
    QString data_string = include_data.c_str();
    
    if (data_string.indexOf(",all,")+1) all = true;
    
    varOut(include_data)
    varOut(all)
    
    std::string CRACS = ConfigHandler::GetInstance()->SoftwarePath();
    for (int i=0; i<e_list.size(); i++) {
        std::cout << "****************************************************" << std::endl;
        std::string filename = e_list.at(i).toStdString()+".txt";
        bool inFit              = false;
        bool e_correctHyperon   = correctHyperon;
        if (all || data_string.                  indexOf( ","+e_list.at(i)+"," )+1 ) inFit               = true;
        if (       hyperon_inclusive_experiments.indexOf( ","+e_list.at(i)+"," )<0 ) e_correctHyperon    = false;
        if (inFit) {
            if (type=="pp")
                fNWparam++;
            else
                fNWparam_pA++;
        }
        out( e_list.at(i).toStdString() )
        out( inFit )
        out( e_correctHyperon )
        data d;
        d.n=0;
        d.A = A;
        d.w = 1;
        FileTool file( CRACS+"/data/CS_measurements/"+type+"/converted_"+filename );
        file.ExtractNumberTable( 9 , " " );
        for (int j = 0; j<file.NumberTableGetNrows(); j++) {
            d.n ++;
            double sqrtS = file.NumberTable(j, 0);
            double hyperonCorrection = 1.;
            if (e_correctHyperon)
                hyperonCorrection =  1./(1.+  0.81 * ( 0.31+0.3/(1+pow(146./sqrtS, 2*0.9)) ) );
            d.sqrtS         .push_back(sqrtS);
            d.pT            .push_back(file.NumberTable(j, 1));
            d.xR            .push_back(file.NumberTable(j, 2));
            d.var_fix       .push_back(file.NumberTable(j, 3));
            d.var_plot      .push_back(file.NumberTable(j, 4));
            d.CS            .push_back(file.NumberTable(j, 5)*hyperonCorrection);
            d.CS_err_stat   .push_back(file.NumberTable(j, 6)*hyperonCorrection);
            d.CS_err_sys    .push_back(file.NumberTable(j, 7)*hyperonCorrection);
            d.err_scale     .push_back(file.NumberTable(j, 8));
            d.CS_err        .push_back( sqrt(d.CS_err_sys.at(j)*d.CS_err_sys.at(j)+d.CS_err_stat.at(j)*d.CS_err_stat.at(j)) );
            if (d.CS_err[j]<d.CS[j]*1e-5)
                d.CS_err[j]=d.CS[j];
        }
        
        for (int j = 0; j<file.NumberTableGetNrows(); j++) {
            QString line = file.GetLine(j).c_str();
            if (line.left(2)=="#*"){
                QStringList list = line.split(" ", QString::SkipEmptyParts);
                for (int k = 0; k<list.size()-1; k++) {
                    d.col_description.push_back(list.at(k+1).toStdString());
                }
                
                break;
            }
            
        }
        
        d.in_fit        = inFit;
        d.experiment    = e_list.at(i).toStdString();
        
        d.err_norm = 0;
        for (int j = 0; j<d.err_scale.size(); j++) {
            d.err_norm += d.err_scale.at(j);
        }
        d.err_norm /= d.err_scale.size();
        
        out(d.n)
        out(d.err_norm)
        
        fData.push_back(d);
        std::cout << "****************************************************" << std::endl;
        
    }
    
}


void Write_results( std::vector<data> &fData ){
    
    for (int i = fNCparam; i<100; i++) {
        fParameters[i] = 0;
    }
    for (int i = fNDparam; i<100; i++) {
        fParameters_pA[i] = 0;
    }
    
    for (int e=0; e<fData.size(); e++) {
        std::string write = "";
        data & d= fData.at(e);
        int A = d.A;
        
        write += "*  ";
        for (int k = 0; k<d.col_description.size(); k++) {
            write += QString("%1").arg( d.col_description.at(k).c_str()  ).leftJustified(20).toStdString();
        }
        write += QString("%1").arg("CS_parametrization" ).leftJustified(20).toStdString();
        write += QString("%1").arg("CS_min_1sigma"      ).leftJustified(20).toStdString();
        write += QString("%1").arg("CS_max_1sigma"      ).leftJustified(20).toStdString();
        write += QString("%1").arg("CS_min_2sigma"      ).leftJustified(20).toStdString();
        write += QString("%1").arg("CS_max_2sigma"      ).leftJustified(20).toStdString();
        write += QString("%1").arg("data rescale"       ).leftJustified(20).toStdString();
        write += "\n";
        double w_e = d.w;
        for (int i=0;i<d.n; i++) {
            double sqrtS    =   d.sqrtS      [i];
            double xR       =   d.xR         [i];
            double pT       =   d.pT         [i];
            double v_fix    =   d.var_fix    [i];
            double v_plot   =   d.var_plot   [i];
            double CS       =   d.CS         [i];
            double CS_err   =   d.CS_err     [i];
            double CS_stat  =   d.CS_err_stat[i];
            double CS_sys   =   d.CS_err_sys [i];
            double err_scale=   d.err_scale  [i];
            double CS_p     =   fun_pp_pbar_CM(sqrtS, pT, xR, fParameters)*fun_pA_pbar_factor( A, sqrtS, pT, xR, fParameters_pA);
            
            double CS_min   =   CS_p;
            double CS_max   =   CS_p;
            double CS_min2  =   CS_p;
            double CS_max2  =   CS_p;
            
            if (A==1){
                FileTool f_1_sigma = FileTool(fDir+"/parameters_1_sigma.txt");
                f_1_sigma.ExtractNumberTable(fNCparam+fNWparam, " ", true);
                for (int j=0; j<f_1_sigma.NumberTableGetNrows(); j++) {
                    double parameters[fNCparam+fNWparam];
                    for (int jj=0; jj<fNCparam+fNWparam; jj++) {
                        parameters[jj] = f_1_sigma.NumberTable(j, jj);
                    }
                  double cs = fun_pp_pbar_CM(sqrtS, pT, xR, parameters)*fun_pA_pbar_factor( A, sqrtS, pT, xR, fParameters_pA);
                    if (cs<CS_min) {
                        CS_min = cs;
                    }
                    if (cs>CS_max) {
                        CS_max = cs;
                    }
                }
                CS_min2  =   CS_min;
                CS_max2  =   CS_max;
                FileTool f_2_sigma = FileTool(fDir+"/parameters_2_sigma.txt");
                f_2_sigma.ExtractNumberTable(fNCparam+fNWparam, " ", true);
                for (int j=0; j<f_2_sigma.NumberTableGetNrows(); j++) {
                    double parameters[fNCparam+fNWparam];
                    for (int jj=0; jj<fNCparam+fNWparam; jj++) {
                        parameters[jj] = f_2_sigma.NumberTable(j, jj);
                    }
                    double cs = fun_pp_pbar_CM(sqrtS, pT, xR, parameters)*fun_pA_pbar_factor( A, sqrtS, pT, xR, fParameters_pA);
                    if (cs<CS_min2) {
                        CS_min2 = cs;
                    }
                    if (cs>CS_max2) {
                        CS_max2 = cs;
                    }
                }
            }else{
                FileTool f_1_sigma = FileTool(fDir+"/parameters_pA_1_sigma.txt");
                f_1_sigma.ExtractNumberTable(fNCparam+fNWparam+fNDparam+fNWparam_pA, " ", true);
                for (int j=0; j<f_1_sigma.NumberTableGetNrows(); j++) {
                    double parameters[fNCparam+fNWparam];
                    for (int jj=0; jj<fNCparam+fNWparam; jj++) {
                        parameters[jj] = f_1_sigma.NumberTable(j, jj);
                    }
                    double parameters_pA[fNDparam+fNWparam_pA];
                    for (int jj=0; jj<fNDparam+fNWparam_pA; jj++) {
                        parameters_pA[jj] = f_1_sigma.NumberTable(j, jj+fNCparam+fNWparam);
                    }
                    double cs = fun_pp_pbar_CM(sqrtS, pT, xR, parameters)*fun_pA_pbar_factor( A, sqrtS, pT, xR, parameters_pA);
                    if (cs<CS_min) {
                        CS_min = cs;
                    }
                    if (cs>CS_max) {
                        CS_max = cs;
                    }
                }
                CS_min2  =   CS_min;
                CS_max2  =   CS_max;
                FileTool f_2_sigma = FileTool(fDir+"/parameters_pA_2_sigma.txt");
                f_2_sigma.ExtractNumberTable(fNCparam+fNWparam+fNDparam+fNWparam_pA, " ", true);
                for (int j=0; j<f_2_sigma.NumberTableGetNrows(); j++) {
                    double parameters[fNCparam+fNWparam];
                    for (int jj=0; jj<fNCparam+fNWparam; jj++) {
                        parameters[jj] = f_2_sigma.NumberTable(j, jj);
                    }
                    double parameters_pA[fNDparam+fNWparam_pA];
                    for (int jj=0; jj<fNDparam+fNWparam_pA; jj++) {
                        parameters_pA[jj] = f_1_sigma.NumberTable(j, jj+fNCparam+fNWparam);
                    }
                    double cs = fun_pp_pbar_CM(sqrtS, pT, xR, parameters)*fun_pA_pbar_factor( A, sqrtS, pT, xR, parameters_pA);
                    if (cs<CS_min2) {
                        CS_min2 = cs;
                    }
                    if (cs>CS_max2) {
                        CS_max2 = cs;
                    }
                }
            }
            
            
            double T_p_LAB, T_pbar_LAB, eta;
            
            CSTransformations::convert_CM_to_LAB(sqrtS*sqrtS, xR, pT, T_p_LAB, T_pbar_LAB, eta, true);
            
            write += "   ";
            write += QString("%1").arg(sqrtS  ).leftJustified(20).toStdString();
            write += QString("%1").arg(pT     ).leftJustified(20).toStdString();
            write += QString("%1").arg(xR     ).leftJustified(20).toStdString();
            
            write += QString("%1").arg(v_fix  ).leftJustified(20).toStdString();
            write += QString("%1").arg(v_plot ).leftJustified(20).toStdString();
            
            write += QString("%1").arg(CS     ).leftJustified(20).toStdString();
            write += QString("%1").arg(CS_err ).leftJustified(20).toStdString();
            write += QString("%1").arg(CS_stat).leftJustified(20).toStdString();
            write += QString("%1").arg(CS_sys ).leftJustified(20).toStdString();
            write += QString("%1").arg(err_scale).leftJustified(20).toStdString();
            
            write += QString("%1").arg(CS_p   ).leftJustified(20).toStdString();
            
            write += QString("%1").arg(CS_min ).leftJustified(20).toStdString();
            write += QString("%1").arg(CS_max ).leftJustified(20).toStdString();
            write += QString("%1").arg(CS_min2).leftJustified(20).toStdString();
            write += QString("%1").arg(CS_max2).leftJustified(20).toStdString();
            
            write += QString("%1").arg(w_e    ).leftJustified(20).toStdString();
            
            write += "\n";
        }
        std::string Astr = QString("%1").arg(A).toStdString();
        
        if( d.in_fit ){
            FileTool::WriteStringToFile(write, "bestFit_"+Astr+"_"+d.experiment+".txt");
        }else{
            FileTool::WriteStringToFile(write, "comparison_"+Astr+"_"+d.experiment+".txt");
        }
        
    }
    
    
}


