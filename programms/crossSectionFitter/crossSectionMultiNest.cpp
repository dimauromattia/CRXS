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

#include <multinest.h>

using namespace  CRACS;

std::string fProgramDescription = "Program to fit cross section parameterizations by means of a MultiNest scan";

int fStep       = 1;

int fMNPoint    = -1;
int fMod        = 30;

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
int     mode_C  = 1;
int     mode_D  = 2;
int     mode_CD = 3;
void    Write_results( std::vector<data> &fData, FileTool &multiNest, int mode );

static int          fPrint           = 0;

/******************************       Constant for the first analysis step: scan  with pp data       ********************************/
int                 fCS              = CS::DI_MAURO12;
int                 fNCparam         = 8;
int                 fNWparam         = 0;
static double       fParameters[100]; // = {4.448, 3.735, 0.00502, 0.708, 3.527, 0.236, -0.729, 2.517, -1.822e-11, 3.527, 0.384}; // Set di Mauro 12 as default, all norm. factors 1 by default
static double       fStart_pp  [100]; // = {4.448, 3.735, 0.00502, 0.708, 3.527, 0.236, -0.729, 2.517, -1.822e-11, 3.527, 0.384};
static double       fStep_pp   [100]; // = {0.01 , 0.1 , 0.001,   0.001, 0.1 ,    0.01 , 0.01,  0.01,  0.001 , 0.01,  0.01  };

static double       fC_lower[100]; //    = { 0,  0,  0, -1,  0, -1, -1,  0}; // MultiNest exploration range
static double       fC_upper[100]; //    = {10, 10,  1,  1, 10,  1,  1, 10}; // MultiNest exploration range


// Constants and default values for MultiNest (pp scan)
std::string         fMNPrefix       = "CS_fit";
std::string         fMNDir          = ".";
int                 fNlive          = 500;
double              fEfr            = 0.5;
double              fTol            = 0.1;
int                 fNdims          = -1;
int                 fNpar           = -1;
double              fChisqMin       = 1e90;
int                 fCount          = 1;

// Array for 1 and 2 sigma regions. It contains the Delta chiSq in n dimensional parameter
//      e.g.  fSigma_1_region[8]   gives the 1 sigma Delta chiSq in an 8D parameter space
double              fSigma_1_region[] = {1.000043426, 2.295815161, 3.526822180, 4.719568761, 5.887700474, 7.038515492, 8.176359741, 9.304044023, 10.42350189, 11.53612748, 12.64296378, 13.74481437, 14.84231367, 15.93597259, 17.02620974, 18.11337308, 19.19775552, 20.27960639, 21.35913992, 22.43654182};
double              fSigma_2_region[] = {4.000009776, 6.180085905, 8.024894670, 9.715641142, 11.31387084, 12.84885057, 14.33712678, 15.78911024, 17.21184692, 18.61036514, 19.98840090, 21.34881933, 22.69387460, 24.02537791, 25.34481051, 26.65340209, 27.95218680, 29.24204418, 30.52372968, 31.79789790};
void                writeParameterRange(bool pA=false);

void CopyStart(  double* start, double* step, double* lower, double* upper  ){
    for (int i=0; i<fNCparam; i++) {
        fParameters[i]  = start[i];
        fStart_pp[i]    = start[i];
        fStep_pp[i]     = step[i];
        fC_lower[i]     = lower[i];
        fC_upper[i]     = upper[i];
    }
    
}
void SetStart(){
    if (fCS==CS::DI_MAURO12){
        double start[]      =   {4.499, 3.41, 0.00942, 0.445, 3.502, 0.0622, -0.247, 2.576 };
        double step []      =   {0.01 , 0.1 , 0.001,   0.001, 0.1 ,  0.01 ,   0.01,  0.01  };
        double lower[]      =   { 0,  0,  0, -1,  0, -1, -1,  0};
        double upper[]      =   {10, 10,  1,  1, 10,  1,  1, 10};
        
        fNCparam = 8;
        CopyStart( start, step, lower, upper );
    }
    if (fCS==CS::WINKLER){
        double start[]      =   { 0.047, 7.76, 0.168, 0.038, 1.0e-3, 0.7  };
        double step []      =   { 0.01 , 0.1 , 0.001, 0.001, 0.001 , 0.01 };
        double lower[]      =   { 0,     0,   -1,    -1,    -1,      0    };
        double upper[]      =   { 1,     20,   1,     1,     1,      10   };
        
        fNCparam = 6;
        CopyStart( start, step, lower, upper );
    }
}

/******************************       Constant for the second analysis step: scan/fit with pA data   ********************************/
int                 fNDparam            = 4;
int                 fNWparam_pA         = 0;
static double       fParameters_pA[100] = { 8.78792e-01,   1.14248e-01, 6.36834e+00, -0.5 };

int                 fN                  = 501; // number of iterations for the TMinuit fit!
static Double_t     fStart_pA[100]      = { 8.78792e-01,   1.14248e-01, 6.36834e+00, -0.5   };
static Double_t     fStep_pA [100]      = { 0.01,          0.01,        0.1,          0.01  };

static Double_t     fD_lower[100]       = { 0.5,   0.05,  1, -1   };
static Double_t     fD_upper[100]       = { 1.2,   0.2,  15,  0   };

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
    if(fCS==CS::DI_MAURO12){
        ret          =   ppCSParametrizations::inv_pp_pbar_CM__diMauro(s, E_pbar, pT_pbar, par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], 0, 0, 0);
    }else if (fCS==CS::WINKLER){
        double C1 = 0.31; double C2 = 0.30; double C3 = 146.*146.; double C4 = 0.9;
        double C5 = par[0]; double C6 = par[1]; double C7 = par[2]; double C8 = par[3]; double C9 = par[4]; double C10= par[5];
        double C11 = 30.9; double C12 = -1.74; double C13 = 0.71;
        double C14 = 0.114; double C15 = 144*144; double C16 = 0.51;
        ret          =   ppCSParametrizations::inv_pp_pbar_CM__Winkler(s, E_pbar, pT_pbar, false, false, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15, C16);
//        double ret_o = ppCSParametrizations::inv_pp_pbar_CM__Winkler(s, E_pbar, pT_pbar, false, false);
//        out(ret)
//        out(ret_o)
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


void Ifit_pp()
{
    
    out(fNCparam)
    out(fNWparam)
    
    fPrint = 0;
    
    TMinuit *gMinuit = new TMinuit(fNCparam+fNWparam);
    gMinuit->SetPrintLevel(-1);
    gMinuit->SetFCN( fcn_pp_pbar_CM );
    
    Double_t arglist[10];
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
    
    for (int it=0; it<fN; it++) {
        if(it%2){
            out("***********************")
            out(it)
            for (int i=0; i<fNCparam; i++) {
                gMinuit->FixParameter(i);
            }
            for (int i=0; i<fNWparam; i++) {
                gMinuit->Release(i+fNCparam);
                //                double val, err;
                //                gMinuit->GetParameter( i+fNCparam, val, err );
                //                out(val)
                //                out(err)
                
            }
        }else{
            for (int i=0; i<fNCparam; i++) {
                gMinuit->Release(i);
            }
            for (int i=0; i<fNWparam; i++) {
                gMinuit->FixParameter(i+fNCparam);
            }
        }
        
        if (it==fN-1) {
            gMinuit->SetPrintLevel(0);
            
        }
        // Now ready for minimization step
        arglist[0] = 1000000;
        arglist[1] = 1.;
        gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
        
        int n;
        for (int i=0; i<fNCparam+fNWparam; i++) {
            double p, pErr;
            gMinuit->GetParameter(i, p, pErr);
            //out(p)
            fParameters[i] = p;
        }
        double h, chisq;
        if (it==fN-1) {
            fPrint = 1;
        }
        fcn_pp_pbar_CM(n, &h, chisq, fParameters, 1);
        
        out(chisq)
        
        
    }
    
    // Print results
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
}



void chisq_minNorm(int e, double* par, double &chisq_min, double &w_min, std::vector<data> & d_vector = fData_pp_pbar_all)
{
    double   delta, chisq;
    chisq_min  = 1e90;
    for (double w=0.5; w<2; w+=0.005) {
        chisq = 0;
        for (int i=0;i<d_vector.at(e).n; i++) {
            double sqrtS    =   d_vector.at(e).sqrtS[i];
            double xR       =   d_vector.at(e).xR   [i];
            double pT       =   d_vector.at(e).pT   [i];
            double CS       =   d_vector.at(e).CS   [i];
            double CS_err   =   d_vector.at(e).CS_err[i];
            delta           =   pow( (w*CS - fun_pp_pbar_CM(sqrtS, pT, xR, par))/CS_err/w, 2 );
            chisq          +=   delta;
        }
        double eps  = d_vector.at(e).err_norm;
        chisq      += pow( (1-w)/eps, 2);
        if (chisq<chisq_min) {
            chisq_min   = chisq;
            w_min       = w;
        }
    }
}


/************************ MultiNest loglikelihood for pp scan  ********************************/

void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
    for (int i = 0; i<ndim; i++) {
        Cube[i] = fC_lower[i] + (  fC_upper[i]-fC_lower[i] )*Cube[i];
    }
    double  chisq   = 0;
    int     ee      = 0;
    for (int e=0; e<fData_pp_pbar_all.size(); e++) {
        if (!fData_pp_pbar_all.at(e).in_fit)
            continue;
        double chisq_e;
        double w_e;
        chisq_minNorm(e, Cube, chisq_e, w_e, fData_pp_pbar_all);
        chisq+=chisq_e;
        Cube[ndim+ee] = w_e;
        ee++;
    }
    // Dump every time a new min ChiSq is found
    if (chisq<fChisqMin){
        fChisqMin = chisq;
        out(fChisqMin)
        out(fCount)
    }
    fCount++;
    lnew = -chisq/2.0;
}

// dumper (standard taken form MultiNest example)
void dumper(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context)
{
    int i, j;
    
    double postdist[nSamples][nPar + 2];
    for( i = 0; i < nPar + 2; i++ )
    for( j = 0; j < nSamples; j++ )
    postdist[j][i] = posterior[0][i * (nSamples) + j];
    
    double pLivePts[nlive][nPar + 1];
    for( i = 0; i < nPar + 1; i++ )
    for( j = 0; j < nlive; j++ )
    pLivePts[j][i] = physLive[0][i * (nlive) + j];
}

//******************************************************************************************************************************/
//  END of functions needed for the pp scan
//******************************************************************************************************************************/




//******************************************************************************************************************************/
//  Functions with factor parameterization for pA
//******************************************************************************************************************************/

double fun_pA_pbar_factor(int A, double sqrtS,  double pT, double x_R_d, Double_t *par)
{
    double s = sqrtS *sqrtS;
    
    double D1  = par[0];
    double D2  = par[1];
    double D3  = par[2];
    double D4  = par[3];
    
    double sigma_0      =   (  1+0.016*sin( 5.3-2.63*log(A) )  );
    double sigma_0_1    =   (  1+0.016*sin( 5.3-2.63*log(1) )  );
    double factor1       =   pow( A, D1 + D2*log( sqrt(s)/D3)*pT  )*sigma_0/sigma_0_1;
    
    // z-asymmm
    double x_R          =   fabs(x_R_d);
    double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
    double E_pbar       =   x_R*E_pbar_Max;
    double p            =   sqrt( E_pbar*E_pbar-fMass_proton*fMass_proton );
    
    double xF = 0;
    if( p*p-pT*pT>0. ){
        double pL = sqrt( p*p-pT*pT );
        double p_pbar_Max   =   sqrt( E_pbar_Max*E_pbar_Max-fMass_proton*fMass_proton );
        if (x_R_d<0) pL*=-1;
        xF = pL/2./sqrt(s);
        if(pL!=pL){
            out(x_R)
            out(p)
            out(pT)
            out(p*p-pT*pT)
            out(pL)
            out(p_pbar_Max)
            out(xF)
            out(D4)
        }
    }
    double factor2 =  pow( A, D4*xF );
    
    return factor1*factor2;
}

double fun_pA_pbar_factor_s(int A, double sqrtS,  double pT, double x_R_d, Double_t *par)
{
    double s = sqrtS *sqrtS;
    double D1  = par[0];
    double D2  = par[1];
    double D3  = par[2];
    
    double sigma_0      =   (  1+0.016*sin( 5.3-2.63*log(A) )  );
    double sigma_0_1    =   (  1+0.016*sin( 5.3-2.63*log(1) )  );
    double factor1       =   pow( A, D1 + D2*log( sqrt(s)/D3)*pT  )*sigma_0/sigma_0_1;

    return factor1;
}

double fun_pA_pbar_factor_a(int A, double sqrtS,  double pT, double x_R_d, Double_t *par)
{
    double s = sqrtS *sqrtS;
    double D4  = par[3];
    double x_R          =   fabs(x_R_d);
    double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrt( s );
    double E_pbar       =   x_R*E_pbar_Max;
    double p            =   sqrt( E_pbar*E_pbar-fMass_proton*fMass_proton );
    double xF           =   0;
    if( p*p-pT*pT>0. ){
        double pL = sqrt( p*p-pT*pT );
        if (x_R_d<0) pL*=-1;
        xF = pL/2./sqrt(s);
    }
    double factor2 =  pow( A, D4*xF );
    return factor2;
}



// MultiNest Scan of the best fit point from pp
//***************************************************************

void chisq_minNorm_pA(int e, double* par, double &chisq_min, double &w_min, std::vector<data> & d_vector = fData_pA_all)
{
    double   delta, chisq;
    chisq_min  = 1e90;
    for (double w=0.5; w<2; w+=0.005) {
        chisq = 0;
        for (int i=0;i<d_vector.at(e).n; i++) {
            double sqrtS    =   d_vector.at(e).sqrtS[i];
            double xR       =   d_vector.at(e).xR   [i];
            double pT       =   d_vector.at(e).pT   [i];
            double CS       =   d_vector.at(e).CS   [i];
            double CS_err   =   d_vector.at(e).CS_err[i];
            int    A        =   d_vector.at(e).A;
            delta           =   pow( (w*CS - fun_pp_pbar_CM(sqrtS, pT, xR, fParameters)*fun_pA_pbar_factor(A, sqrtS, pT, xR, par))/CS_err/w, 2 );
            chisq          +=   delta;
        }
        double eps  = d_vector.at(e).err_norm;
        chisq      += pow( (1-w)/eps, 2);
        if (chisq<chisq_min) {
            chisq_min   = chisq;
            w_min       = w;
        }
    }
}

void LogLike_pA(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
    for (int i = 0; i<ndim; i++) {
        Cube[i] = fD_lower[i] + (  fD_upper[i]-fD_lower[i] )*Cube[i];
        varOut(Cube[i])
        varOut(fStart_pA[i])
    }
    double  chisq   = 0;
    int     ee      = 0;
    for (int e=0; e<fData_pA_all.size(); e++) {
        if (!fData_pA_all.at(e).in_fit)
            continue;
        double chisq_e;
        double w_e;
        chisq_minNorm_pA(e, Cube, chisq_e, w_e, fData_pA_all);
        chisq+=chisq_e;
        varOut(e)
        varOut(chisq_e)
        varOut(w_e)
        Cube[ndim+ee] = w_e;
        ee++;
    }
    varOut(chisq)
    varOut("*************")
    // Dump every time a new min ChiSq is found
    if (chisq<fChisqMin){
        fChisqMin = chisq;
        out(fChisqMin)
        out(fCount)
    }
    fCount++;
    lnew = -chisq/2.0;
}


// dumper (standard taken form MultiNest example)
void dumper_pA(int &nSamples, int &nlive, int &nPar, double **physLive, double **posterior, double **paramConstr, double &maxLogLike, double &logZ, double &INSlogZ, double &logZerr, void *context)
{
    int i, j;
    
    double postdist[nSamples][nPar + 2];
    for( i = 0; i < nPar + 2; i++ )
        for( j = 0; j < nSamples; j++ )
            postdist[j][i] = posterior[0][i * (nSamples) + j];
    
    double pLivePts[nlive][nPar + 1];
    for( i = 0; i < nPar + 1; i++ )
        for( j = 0; j < nlive; j++ )
            pLivePts[j][i] = physLive[0][i * (nlive) + j];
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
            if (fPrint || (delta>1e10) ) {
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
        ndf = ndf-fNCparam-fNWparam;
        out(ndf)
        out(chisq/ndf)
        std::cout << "*************************" << std::endl;
        std::cout << "-------------------------" << std::endl;
    }
    f = chisq;
}

double Ifit_pA(bool quiet=false)
{
    fPrint = 0;
    
    TMinuit *gMinuit = new TMinuit(fNDparam+fNWparam);
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
    
    double chisq;
    for (int it=0; it<fN; it++) {
        if(it%2){
            for (int i=0; i<fNDparam; i++) {
                gMinuit->FixParameter(i);
            }
            for (int i=0; i<fNWparam_pA; i++) {
                gMinuit->Release(i+fNDparam);
            }
        }else{
            for (int i=0; i<fNDparam; i++) {
                gMinuit->Release(i);
            }
            for (int i=0; i<fNWparam_pA; i++) {
                gMinuit->FixParameter(i+fNDparam);
            }
        }
        
        if (it==fN-1) {
            gMinuit->SetPrintLevel(0);
            
        }
        // Now ready for minimization step
        arglist[0] = 1000000;
        arglist[1] = 1.;
        gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
        
        int n;
        for (int i=0; i<fNDparam+fNWparam_pA; i++) {
            double p, pErr;
            gMinuit->GetParameter(i, p, pErr);
            //out(p)
            fParameters_pA[i] = p;
        }
        double h;
        if (it==fN-1 && !quiet) {
            fPrint = 1;
        }
        fcn_pA_pbar_CM (n, &h, chisq, fParameters_pA, 1);
        if(!quiet){
            out(chisq)
        }
    }
    int ee = 0;
    for (int e=0; e<fData_pA_all.size(); e++) {
        data &  d       = fData_pA_all[e];
        if(d.in_fit){
            d.w = fParameters_pA[fNDparam+ee];
        }
        ee ++;
    }
    
    if(quiet){
        return chisq;
    }
    
    // Print results
    Double_t amin,edm,errdef;
    Int_t nvpar,nparx,icstat;
    gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
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
double oxygenFlux(double T_p){
    return LISFluxes::spectrum(T_p, LISFluxes::fOxygenLIS);
}


//pp
static double inv_pp_pbar_CM( double s, double E_pbar, double pT_pbar){
    double ret          =   0;
    if(fCS==CS::DI_MAURO12){
        ret          =   ppCSParametrizations::inv_pp_pbar_CM__diMauro(s, E_pbar, pT_pbar, fParameters[0], fParameters[1], fParameters[2], fParameters[3], fParameters[4], fParameters[5], fParameters[6], fParameters[7], 0, 0, 0);
    }else if (fCS==CS::WINKLER){
        double C1 = 0.31; double C2 = 0.30; double C3 = 146.*146.; double C4 = 0.9;
        double C5 = fParameters[0]; double C6 = fParameters[1]; double C7 = fParameters[2]; double C8 = fParameters[3]; double C9 = fParameters[4]; double C10= fParameters[5];
        double C11 = 30.9; double C12 = -1.74; double C13 = 0.71;
        double C14 = 0.114; double C15 = 144*144; double C16 = 0.51;
        ret          =   ppCSParametrizations::inv_pp_pbar_CM__Winkler(s, E_pbar, pT_pbar, false, false, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15, C16);
    }
    return ret;
}
static double dE_pp_pbar_LAB(double E_p_LAB, double E_pbar_LAB, int precision=10000, int output=NO_OUT){
        return CSTransformations::dE_pp_product_LAB_intEta (E_p_LAB, E_pbar_LAB, fMass_proton,
                                                            inv_pp_pbar_CM,
                                                            precision,
                                                            output);
}
static double dT_pp_pbar_LAB  ( double T_p, double T_pbar ){
    return dE_pp_pbar_LAB( T_p+fMass_proton,   T_pbar+fMass_proton)*1e-31;
}


// AA
static double inv_A1A2_pbar_CM(  double s, double E_pbar_d, double pT_pbar, int A1, int A2){
    
    double E_pbar   = fabs(E_pbar_d);
    double CS       = inv_pp_pbar_CM( s, E_pbar, pT_pbar );
    if (CS<1e-90)
        return 0;
    
    double factor1_s = fun_pA_pbar_factor_s( A1, sqrt(s), pT_pbar,  E_pbar, fParameters_pA );
    double factor1_a = fun_pA_pbar_factor_a( A1, sqrt(s), pT_pbar, -E_pbar, fParameters_pA ); // flip sign of xF for z-asymm of incient nucleus
    double factor2_s = fun_pA_pbar_factor_s( A2, sqrt(s), pT_pbar,  E_pbar, fParameters_pA );
    double factor2_a = fun_pA_pbar_factor_a( A2, sqrt(s), pT_pbar,  E_pbar, fParameters_pA );
    
    return factor1_s * factor1_a * factor2_s * factor2_a * CS;
}

//  pHe
static double inv_pHe_pbar_CM( double s, double E_pbar_d, double pT_pbar){
    return inv_A1A2_pbar_CM( s, E_pbar_d, pT_pbar, 1, 4);
}
static double dE_pHe_pbar_LAB(double E_p_LAB, double E_pbar_LAB, int precision=10000, int output=NO_OUT){
    return CSTransformations::dE_pp_product_LAB_intEta (E_p_LAB, E_pbar_LAB, fMass_proton,
                                                        inv_pHe_pbar_CM,
                                                        precision,
                                                        output);
}
static double dT_pHe_pbar_LAB  ( double T_p, double T_pbar ){
    return dE_pHe_pbar_LAB( T_p+fMass_proton,   T_pbar+fMass_proton)*1e-31;
}

//  Hep
static double dE_Hep_pbar_LAB(double En_He_LAB, double E_pbar_LAB, int precision=10000, int output=NO_OUT){
    return CSTransformations::dE_Hep_product_LAB_intEta(En_He_LAB, E_pbar_LAB, fMass_proton,
                                                    inv_pHe_pbar_CM,
                                                    precision,
                                                    output);
}
static double dT_Hep_pbar_LAB  ( double Tn_He, double T_pbar ){
    return dE_Hep_pbar_LAB( Tn_He+fMass_proton,   T_pbar+fMass_proton)*1e-31;
}

//  HeHe
static double inv_HeHe_pbar_CM( double s, double E_pbar_d, double pT_pbar){
    return inv_A1A2_pbar_CM( s, E_pbar_d, pT_pbar, 4, 4);
}
static double dE_HeHe_pbar_LAB(double En_He_LAB, double E_pbar_LAB, int precision=10000, int output=NO_OUT){
    return CSTransformations::dE_pp_product_LAB_intEta (En_He_LAB, E_pbar_LAB, fMass_proton,
                                                        inv_HeHe_pbar_CM,
                                                        precision,
                                                        output);
}
static double dT_HeHe_pbar_LAB  ( double Tn_He, double T_pbar ){
    return dE_pHe_pbar_LAB( Tn_He+fMass_proton,   T_pbar+fMass_proton)*1e-31;
}

//******************************************************************************************************************************/
//  END   Functions for the source term calculation
//******************************************************************************************************************************/




/***********************************************************************************************************************/
/************************************           Main program        ****************************************************/
/***********************************************************************************************************************/
int main(int argc, char *argv[])
{
    // Config handling
    
    ConfigHandler* config = ConfigHandler::GetInstance();
    config->SetInput(argc, argv, fProgramDescription);
    
    int cs;
    config->AddOptionInt        (  "step",          fStep,         "default: 1"                          );
    config->AddOptionInt        (  "cs",            cs,            "default: 1  (1: di Mauro 12, 2: Winkler"        );
    config->AddOptionInt        (  "point",         fMNPoint,      "default: -1  (best fit)"             );
    config->AddOptionInt        (  "mod",           fMod,          "default: 30"                         );
    config->AddOptionInt        (  "nlive",         fNlive,        "default: 500"                        );
    config->AddOptionDouble     (  "efr",           fEfr,          "default: 0.5"                        );
    config->AddOptionDouble     (  "tol",           fTol,          "default: 0.1"                        );
    config->AddOption           (  "dir",           fMNDir,        "default: ."                          );
    config->AddOption           (  "prefix",        fMNPrefix,     "default: CS_fit"                     );
    
    if       (cs==1){
        fCS = CS::DI_MAURO12;
    }else if (cs==2){
        fCS = CS::WINKLER;
    }
    SetStart();

    std::string include_data = "allaby,antreasyan,brahms,capiluppi,dekkers,guettler,johnson,na49";
    config->AddOption("data", include_data, "Use this to specify experiments, e.g. '--data na49,phenix', Default: all");
    bool no_hyperon;
    config->AddOptionTrue("doNotCorrectHyperon", no_hyperon);
    
    std::string include_data_pHe = "lhcb";
    config->AddOption("data_pHe", include_data_pHe, "Use this to specify experiments, e.g. '--data na49,phenix', Default: lhcb");
    bool no_hyperon_pHe;
    config->AddOptionTrue("doNotCorrectHyperon_pHe", no_hyperon_pHe);
    
    std::string include_data_pBe = "";
    config->AddOption("data_pBe", include_data_pBe, "Use this to specify experiments, e.g. '--data na49', Default: all");
    bool no_hyperon_pBe;
    config->AddOptionTrue("doNotCorrectHyperon_pHe", no_hyperon_pBe);
    
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
    
    
    fNdims = fNCparam;
    fNpar  = fNCparam+fNWparam;
    std::system(  ("mkdir "+fMNDir+"/MultiNest").c_str() );
    writeParameterRange();
    
    
//    for (int i = 0; i<fNCparam; i++) {
//        out(i)
//        out(fStart_pp[i])
//    }
//    return 1;
    
//* Make a prefit with TMinuit for pp
//*********************************************************************************************
    if (fStep==-1) {
        std::system( ("echo dummy >> "+fMNDir+"/MultiNest/"+fMNPrefix+".txt").c_str() );
        out("echo dummy >> "+fMNDir+"/MultiNest"+fMNPrefix+".txt")
        Ifit_pp();
        FileTool dummy;
        Write_results( fData_pp_pbar_all, dummy, -1 );
        
        LISFluxes::fit();
        std::system( ("mkdir "+fMNDir+"/sourceTerm").c_str() );
        std::string number = QString("%1").arg(fMNPoint).rightJustified(5, '0').toStdString();
        if (fMNPoint==-1)
        number="_best";
        fTab = CS_lab_tab::GetInstance();
        fTab->WriteCS(  dT_pp_pbar_LAB, "sourceTerm/dT_pp_pbar_LAB__"+number+".txt", fMNDir  );
        double nH   = 1e+6;
        fTab->ReadCS (  "sourceTerm/dT_pp_pbar_LAB__"+number+".txt", fMNDir  );
        std::string sourceTerm = "   ";
        sourceTerm += QString("%1").arg("T_pbar"  ).leftJustified(20).toStdString();
        sourceTerm += QString("%1").arg("p_H"     ).leftJustified(20).toStdString();
        sourceTerm += "\n";
        for (double dT=-1; dT<4; dT+=1./30) {
            double T_pbar   = pow(10, dT);
            double Tn_from  = std::min(T_pbar, 6*fMass_proton);
            //out(protonFlux(T_pbar))
            double S_p_H    = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,   protonFlux, nH,  Tn_from, 1e6, 10000, NO_OUT);
            sourceTerm += QString("%1").arg(T_pbar   ).leftJustified(20).toStdString();
            sourceTerm += QString("%1").arg(S_p_H    ).leftJustified(20).toStdString();
            sourceTerm += "\n";
            
        }
        FileTool::WriteStringToFile(sourceTerm, fMNDir+"/sourceTerm/sourceTerm_pp_"+number+".txt" );
        
    }
    
//* Call MultiNest for pp scan
//*********************************************************************************************
    if (fStep==1) {
        // set the MultiNest sampling parameters
        int nlive       = fNlive;       // number of live points
        double efr      = fEfr;         // set the required efficiency
        double tol      = fTol;         // tol, defines the stopping criteria
        int ndims       = fNdims;       // dimensionality (no. of free parameters)
        int nPar        = fNpar;        // total no. of parameters including free & derived parameters
        int IS          = 0;            // do Nested Importance Sampling?
        int mmodal      = 1;            // do mode separation?
        int ceff        = 0;            // run in constant efficiency mode?
        int nClsPar     = 2;            // no. of parameters to do mode separation on
        int updInt      = 10;           // after how many iterations feedback is required & the output files should be updated
                                        // note: posterior files are updated & dumper routine is called after every updInt*10 iterations
        double Ztol     = -1E90;        // all the modes with logZ < Ztol are ignored
        int maxModes    = 100;          // expected max no. of modes (used only for memory allocation)
        int pWrap[ndims];               // which parameters to have periodic boundary conditions?
        for(int i = 0; i < ndims; i++)
        pWrap[i]    = 0;
        char root[100];                 // root for output files
        strcpy(root, (fMNDir+"/MultiNest/"+fMNPrefix).c_str());
        int seed        = -1;           // random no. generator seed, if < 0 then take the seed from system clock
        int fb          = 1;            // need feedback on standard output?
        int resume      = 1;            // resume from a previous job?
        int outfile     = 1;            // write output files?
        int initMPI     = 1;            // initialize MPI routines?, relevant only if compiling with MPI
                                        // set it to F if you want your main program to handle MPI initialization
        double logZero  = -1E90;        // points with loglike < logZero will be ignored by MultiNest
        int maxiter     = 0;            // max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it
                                        // has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
        void *context   = 0;            // not required by MultiNest, any additional information user wants to pass
        
        // calling MultiNest
        out(nlive)
        out(efr)
        out(tol)
        nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LogLike, dumper, context);
    }
    
//* Write fit results for the pp scan and prepare the jobs for the pp source term calculation
//*********************************************************************************************
    if (fStep==2) {
        FileTool multiNest(fMNDir+"/MultiNest/"+fMNPrefix+".txt");
        multiNest.ExtractNumberTable(fNpar+2, " ", true);
        
        double chi;
        double chiMin   = multiNest.NumberTable(0, 1);
        int    indexBestFit = 0;
        for (int i = 0; i<multiNest.NumberTableGetNrows(); i++) {
            chi     = multiNest.NumberTable(i, 1);
            if (chi<chiMin){
                chiMin       = chi;
                indexBestFit = i;
            }
        }
        
        int ee = 0;
        for (int e=0; e<fData_pp_pbar_all.size(); e++) {
            data &  d       = fData_pp_pbar_all[e];
            if(d.in_fit){
                d.w = multiNest.NumberTable( indexBestFit, 2+fNCparam+ee );
            }
            ee ++;
        }
        
        Write_results( fData_pp_pbar_all, multiNest, mode_C );
        
        
        // prepare step 12, i.e write a file with the commands to run the source term calculation for each MN point in the 2 sigma region
        std::string runFile;
        runFile += "CRACS_crossSectionMultiNest --step 12 --point "+QString("%1").arg(-1).toStdString()+" --dir "+fMNDir+" --prefix "+fMNPrefix+" \n";
        int cout = 0;
        for (int i = multiNest.NumberTableGetNrows()-1; i>=0; i--) {
            if( i%fMod!=0 )
            continue;
            double chisq = multiNest.NumberTable(i, 1);
            if(chisq<=chiMin+fSigma_2_region[fNCparam]){
                cout++;
                runFile += "CRACS_crossSectionMultiNest --step 12 --point "+QString("%1").arg(multiNest.NumberTableGetNrows()-i).toStdString()+" --dir "+fMNDir+" --prefix "+fMNPrefix+" \n";
            }
        }
        FileTool::WriteStringToFile(runFile, "runFile_pp.txt");
    }

//* TMinuit to get the pA factor
//*********************************************************************************************
    if (fStep==3) {
        
        FileTool multiNest(fMNDir+"/MultiNest/"+fMNPrefix+".txt");
        multiNest.ExtractNumberTable(fNpar+2, " ", true);
        
        for (int i = fNCparam; i<100; i++) {
            fParameters[i] = 0;
        }
        std::string write;
        write += "#  ";
        write += QString("weight"       ).leftJustified(20).toStdString();
        write += QString("chiSq_tot"    ).leftJustified(20).toStdString();
        for (int i = 0; i<fNCparam; i++) {
            write += QString("C%1").arg(  i+1  ).leftJustified(20).toStdString();
        }
        for (int i = 0; i<fNWparam; i++) {
            write += QString("w%1").arg(  i+1  ).leftJustified(20).toStdString();
        }
        for (int i = 0; i<fNDparam; i++) {
            write += QString("D%1").arg(  i+1  ).leftJustified(20).toStdString();
        }
        for (int i = 0; i<fNWparam_pA; i++) {
            write += QString("w_pA%1").arg(  i+1  ).leftJustified(20).toStdString();
        }
        write += QString("chiSq_pp" ).leftJustified(20).toStdString();
        write += QString("chiSq_pA" ).leftJustified(20).toStdString();
        write += "\n";
        
        for (int j = multiNest.NumberTableGetNrows()-1; j>=0; j--) {
            out(1.*j/multiNest.NumberTableGetNrows())
            //if (1.*j/multiNest.NumberTableGetNrows()<0.99) break;
            for (int k = 0; k<fNCparam; k++) {
                    fParameters[k] = multiNest.NumberTable(  j, 2+k  );
            }
            double chisq = Ifit_pA(true);
            out(chisq)
            write += "   ";
            write += QString("%1").arg(  multiNest.NumberTable(  j, 0  )          ).leftJustified(20).toStdString();
            write += QString("%1").arg(  multiNest.NumberTable(  j, 1  )+chisq    ).leftJustified(20).toStdString();
            for (int k = 0; k<fNpar; k++) {
                write += QString("%1").arg(  multiNest.NumberTable(  j, 2+k  )          ).leftJustified(20).toStdString();
            }
            for (int i = 0; i<fNDparam; i++) {
                write += QString("%1").arg(  fParameters_pA[i]          ).leftJustified(20).toStdString();
            }
            for (int i = 0; i<fNWparam_pA; i++) {
                write += QString("%1").arg(  fParameters_pA[i+fNDparam] ).leftJustified(20).toStdString();
            }
            write += QString("%1").arg(  multiNest.NumberTable(  j, 1  )        ).leftJustified(20).toStdString();
            write += QString("%1").arg(  chisq                                  ).leftJustified(20).toStdString();
            write += "\n";
        }
        FileTool::WriteStringToFile(write, fMNDir+"/MultiNest/"+fMNPrefix+"_pA_minuit.txt");
        writeParameterRange(true);
        
    }
    
//* Write results for the pA scan with errors from pp
//*********************************************************************************************
    if (fStep==4) {
        FileTool multiNest_pA(fMNDir+"/MultiNest/"+fMNPrefix+"_pA_minuit.txt");
        multiNest_pA.ExtractNumberTable(2+fNCparam+fNWparam+fNDparam+fNWparam_pA+2, " ", true);
        Write_results( fData_pA_all, multiNest_pA, mode_CD );
    }

    fNdims = fNDparam;
    fNpar  = fNDparam+fNWparam_pA;
    
//* Call MultiNest for pA scan (fix best fit pp)
//*********************************************************************************************
    if (fStep==5) {
        // read best fit point from pp scan
        FileTool multiNest(fMNDir+"/MultiNest/"+fMNPrefix+".txt");
        multiNest.ExtractNumberTable(fNCparam+fNWparam+2, " ", true);
        double chi;
        double chiMin   = multiNest.NumberTable(0, 1);
        int    indexBestFit = 0;
        for (int i = 0; i<multiNest.NumberTableGetNrows(); i++) {
            chi     = multiNest.NumberTable(i, 1);
            if (chi<chiMin){
                chiMin       = chi;
                indexBestFit = i;
            }
        }
        for (int k = 0; k<fNCparam; k++) {
            fParameters[k] = multiNest.NumberTable(  indexBestFit, 2+k  );
            out(fParameters[k])
        }
        for (int i = fNCparam; i<100; i++) {
            fParameters[i] = 0;
        }
        
        //test
        
        for (int e = 0; e<fData_pA_all.size(); e++) {
            data & d = fData_pA_all.at(e);
            if (!d.in_fit)
                continue;
            out(e)
            double chisq_min, w_e;
            chisq_minNorm_pA(e, fStart_pA, chisq_min, w_e);
            out(chisq_min)
            out(w_e)
            out("***************")
        }
        
        
        
        
        // set the MultiNest sampling parameters
        int nlive       = fNlive;       // number of live points
        double efr      = fEfr;         // set the required efficiency
        double tol      = fTol;         // tol, defines the stopping criteria
        int ndims       = fNdims;       // dimensionality (no. of free parameters)
        int nPar        = fNpar;        // total no. of parameters including free & derived parameters
        int IS          = 0;            // do Nested Importance Sampling?
        int mmodal      = 1;            // do mode separation?
        int ceff        = 0;            // run in constant efficiency mode?
        int nClsPar     = 2;            // no. of parameters to do mode separation on
        int updInt      = 10;           // after how many iterations feedback is required & the output files should be updated
                                        // note: posterior files are updated & dumper routine is called after every updInt*10 iterations
        double Ztol     = -1E90;        // all the modes with logZ < Ztol are ignored
        int maxModes    = 100;          // expected max no. of modes (used only for memory allocation)
        int pWrap[ndims];               // which parameters to have periodic boundary conditions?
        for(int i = 0; i < ndims; i++)
            pWrap[i]    = 0;
        char root[100];                 // root for output files
        strcpy(root, (fMNDir+"/MultiNest/"+fMNPrefix+"_pA_scan").c_str());
        int seed        = -1;           // random no. generator seed, if < 0 then take the seed from system clock
        int fb          = 1;            // need feedback on standard output?
        int resume      = 1;            // resume from a previous job?
        int outfile     = 1;            // write output files?
        int initMPI     = 1;            // initialize MPI routines?, relevant only if compiling with MPI
                                        // set it to F if you want your main program to handle MPI initialization
        double logZero  = -1E90;        // points with loglike < logZero will be ignored by MultiNest
        int maxiter     = 0;            // max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it
                                        // has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
        void *context   = 0;            // not required by MultiNest, any additional information user wants to pass
        
        // calling MultiNest
        out(nlive)
        out(efr)
        out(tol)
        nested::run(IS, mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, logZero, maxiter, LogLike_pA, dumper_pA, context);
    }

    if (fStep==6) {
        // read best fit point from pp scan
        FileTool multiNest(fMNDir+"/MultiNest/"+fMNPrefix+".txt");
        multiNest.ExtractNumberTable(fNCparam+fNWparam+2, " ", true);
        double chi;
        double chiMin   = multiNest.NumberTable(0, 1);
        int    indexBestFit = 0;
        for (int i = 0; i<multiNest.NumberTableGetNrows(); i++) {
            chi     = multiNest.NumberTable(i, 1);
            if (chi<chiMin){
                chiMin       = chi;
                indexBestFit = i;
            }
        }
        for (int k = 0; k<fNCparam; k++) {
            fParameters[k] = multiNest.NumberTable(  indexBestFit, 2+k  );
        }
        for (int i = fNCparam; i<100; i++) {
            fParameters[i] = 0;
        }
        FileTool multiNest_pA(fMNDir+"/MultiNest/"+fMNPrefix+"_pA_scan.txt");
        multiNest_pA.ExtractNumberTable(2+fNDparam+fNWparam_pA, " ", true);
        Write_results( fData_pA_all, multiNest_pA, mode_D );
        
        // prepare step 12, i.e write a file with the commands to run the source term calculation for each MN point in the 2 sigma region
        double chi_pA;
        double chiMin_pA   = multiNest.NumberTable(0, 1);
        int    indexBestFit_pA = 0;
        for (int i = 0; i<multiNest_pA.NumberTableGetNrows(); i++) {
            chi_pA     = multiNest_pA.NumberTable(i, 1);
            if (chi_pA<chiMin_pA){
                chiMin_pA       = chi_pA;
                indexBestFit_pA = i;
            }
        }
        std::string runFile;
        runFile += "CRACS_crossSectionMultiNest --step 16 --point "+QString("%1").arg(-1).toStdString()+" --dir "+fMNDir+" --prefix "+fMNPrefix+" \n";
        int cout = 0;
        for (int i = multiNest_pA.NumberTableGetNrows()-1; i>=0; i--) {
            if( i%fMod!=0 )
                continue;
            double chisq = multiNest_pA.NumberTable(i, 1);
            if(chisq<=chiMin_pA+fSigma_2_region[fNDparam]){
                cout++;
                runFile += "CRACS_crossSectionMultiNest --step 16 --point "+QString("%1").arg(multiNest_pA.NumberTableGetNrows()-i).toStdString()+" --dir "+fMNDir+" --prefix "+fMNPrefix+" \n";
            }
        }
        FileTool::WriteStringToFile(runFile, "runFile_pA.txt");
    }
    
    
    if (fStep==12) {
       
        out(fMNPoint)
        double nH   = 1e+6;
        LISFluxes::fit();
        fTab = CS_lab_tab::GetInstance();
        std::system( ("mkdir "+fMNDir+"/sourceTerm").c_str() );
        
        //***********************************
        //    source term  pp
        //***********************************
        FileTool multiNest(fMNDir+"/MultiNest/"+fMNPrefix+".txt");
        multiNest.ExtractNumberTable(fNCparam+fNWparam+2, " ", true);
        double chi;
        double chiMin   = multiNest.NumberTable(0, 1);
        int    indexBestFit = 0;
        for (int i = 0; i<multiNest.NumberTableGetNrows(); i++) {
            chi     = multiNest.NumberTable(i, 1);
            if (chi<chiMin){
                chiMin       = chi;
                indexBestFit = i;
            }
        }
        int index = multiNest.NumberTableGetNrows()-fMNPoint;
        if (fMNPoint==-1)
            index = indexBestFit;
        for (int i = 0; i<fNCparam; i++) {
            fParameters[i] = multiNest.NumberTable( index , 2+i  );
        }
        for (int i = fNCparam; i<100; i++) {
            fParameters[i] = 0;
        }
        std::string number = QString("%1").arg(fMNPoint).rightJustified(5, '0').toStdString();
        if (fMNPoint==-1)
            number="_best";
        fTab->WriteCS(  dT_pp_pbar_LAB, "sourceTerm/dT_pp_pbar_LAB__"+number+".txt", fMNDir  );
        // Get source term
        fTab->ReadCS (  "sourceTerm/dT_pp_pbar_LAB__"+number+".txt", fMNDir  );
        std::string sourceTerm = "   ";
        sourceTerm += QString("%1").arg("T_pbar"  ).leftJustified(20).toStdString();
        sourceTerm += QString("%1").arg("p_H"     ).leftJustified(20).toStdString();
        sourceTerm += "\n";
        for (double dT=-1; dT<4; dT+=1./30) {
            double T_pbar   = pow(10, dT);
            double Tn_from  = std::min(T_pbar, 6*fMass_proton);
            //out(protonFlux(T_pbar))
            double S_p_H    = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,   protonFlux, nH,  Tn_from, 1e6, 10000, NO_OUT);
            sourceTerm += QString("%1").arg(T_pbar   ).leftJustified(20).toStdString();
            sourceTerm += QString("%1").arg(S_p_H    ).leftJustified(20).toStdString();
            sourceTerm += "\n";
            
        }
        FileTool::WriteStringToFile(sourceTerm, fMNDir+"/sourceTerm/sourceTerm_pp_"+number+".txt" );
    }
    
    if (fStep==16) {
        
        
        out(fMNPoint)
        
        double nH   = 1e+6;
        double nHe  = 1e+5;
        double nC   = 5e+2;
        
        LISFluxes::fit();
        fTab = CS_lab_tab::GetInstance();
        std::system( ("mkdir "+fMNDir+"/sourceTerm").c_str() );
        
        //***********************************
        //    set  pp parameters to best fit
        //***********************************
        // read best fit point from pp scan
        FileTool multiNest(fMNDir+"/MultiNest/"+fMNPrefix+".txt");
        multiNest.ExtractNumberTable(fNCparam+fNWparam+2, " ", true);
        double chi;
        double chiMin   = multiNest.NumberTable(0, 1);
        int    indexBestFit = 0;
        for (int i = 0; i<multiNest.NumberTableGetNrows(); i++) {
            chi     = multiNest.NumberTable(i, 1);
            if (chi<chiMin){
                chiMin       = chi;
                indexBestFit = i;
            }
        }
        for (int k = 0; k<fNCparam; k++) {
            fParameters[k] = multiNest.NumberTable(  indexBestFit, 2+k  );
        }
        for (int i = fNCparam; i<100; i++) {
            fParameters[i] = 0;
        }
        //***********************************
        //    source term  pA and AA
        //***********************************
        
        
        FileTool multiNest_pA(fMNDir+"/MultiNest/"+fMNPrefix+"_pA_scan.txt");
        multiNest_pA.ExtractNumberTable(fNDparam+fNWparam_pA+2, " ", true);
        
        double chi_pA;
        double chiMin_pA   = multiNest.NumberTable(0, 1);
        int    indexBestFit_pA = 0;
        for (int i = 0; i<multiNest_pA.NumberTableGetNrows(); i++) {
            chi_pA     = multiNest_pA.NumberTable(i, 1);
            if (chi_pA<chiMin_pA){
                chiMin_pA       = chi_pA;
                indexBestFit_pA = i;
            }
        }
        int index = multiNest_pA.NumberTableGetNrows()-fMNPoint;
        if (fMNPoint==-1)
            index = indexBestFit_pA;
        for (int i = 0; i<fNDparam; i++) {
            fParameters_pA[i] = multiNest_pA.NumberTable(  index, 2+i  );
        }
        for (int i = fNDparam; i<100; i++) {
            fParameters_pA[i] = 0;
        }
        std::string number = QString("%1").arg(fMNPoint).rightJustified(5, '0').toStdString();
        if (fMNPoint==-1)
            number="_best";
        //      pHe
        //*******************
        fTab->WriteCS(  dT_pHe_pbar_LAB, "sourceTerm/dT_pHe_pbar_LAB__"+number+".txt", fMNDir  );
        fTab->ReadCS (  "sourceTerm/dT_pHe_pbar_LAB__"+number+".txt", fMNDir  );
        std::string sourceTerm = "   ";
        sourceTerm += QString("%1").arg("T_pbar"  ).leftJustified(20).toStdString();
        sourceTerm += QString("%1").arg("p_He"    ).leftJustified(20).toStdString();
        sourceTerm += "\n";
        for (double dT=-1; dT<4; dT+=1./30) {
            double T_pbar   = pow(10, dT);
            double Tn_from  = std::min(T_pbar, 6*fMass_proton);
            //out(protonFlux(T_pbar))
            double S    = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,   protonFlux, nHe,  Tn_from, 1e6, 10000, NO_OUT);
            sourceTerm += QString("%1").arg(T_pbar   ).leftJustified(20).toStdString();
            sourceTerm += QString("%1").arg(S        ).leftJustified(20).toStdString();
            sourceTerm += "\n";
            
        }
        FileTool::WriteStringToFile(sourceTerm, fMNDir+"/sourceTerm/sourceTerm_pHe_"+number+".txt" );
        
        //      Hep
        //*******************
        fTab->WriteCS(  dT_Hep_pbar_LAB, "sourceTerm/dT_Hep_pbar_LAB__"+number+".txt", fMNDir  );
        fTab->ReadCS (  "sourceTerm/dT_Hep_pbar_LAB__"+number+".txt", fMNDir  );
        sourceTerm = "   ";
        sourceTerm += QString("%1").arg("T_pbar"  ).leftJustified(20).toStdString();
        sourceTerm += QString("%1").arg("He_H"    ).leftJustified(20).toStdString();
        sourceTerm += "\n";
        for (double dT=-1; dT<4; dT+=1./30) {
            double T_pbar   = pow(10, dT);
            double Tn_from  = std::min(T_pbar, 6*fMass_proton);
            //out(protonFlux(T_pbar))
            double S    = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,   heliumFlux, nH,  Tn_from, 1e6, 10000, NO_OUT);
            sourceTerm += QString("%1").arg(T_pbar   ).leftJustified(20).toStdString();
            sourceTerm += QString("%1").arg(S        ).leftJustified(20).toStdString();
            sourceTerm += "\n";
            
        }
        FileTool::WriteStringToFile(sourceTerm, fMNDir+"/sourceTerm/sourceTerm_Hep_"+number+".txt" );
       
        //      HeHe
        //*******************
        fTab->WriteCS(  dT_HeHe_pbar_LAB, "sourceTerm/dT_HeHe_pbar_LAB__"+number+".txt", fMNDir  );
        fTab->ReadCS (  "sourceTerm/dT_HeHe_pbar_LAB__"+number+".txt", fMNDir  );
        sourceTerm = "   ";
        sourceTerm += QString("%1").arg("T_pbar"  ).leftJustified(20).toStdString();
        sourceTerm += QString("%1").arg("He_He"   ).leftJustified(20).toStdString();
        sourceTerm += "\n";
        for (double dT=-1; dT<4; dT+=1./30) {
            double T_pbar   = pow(10, dT);
            double Tn_from  = std::min(T_pbar, 6*fMass_proton);
            //out(protonFlux(T_pbar))
            double S    = SourceTerms::integrationFromTo(T_pbar, dT__pbar_LAB_interpolation,   heliumFlux, nHe,  Tn_from, 1e6, 10000, NO_OUT);
            sourceTerm += QString("%1").arg(T_pbar   ).leftJustified(20).toStdString();
            sourceTerm += QString("%1").arg(S        ).leftJustified(20).toStdString();
            sourceTerm += "\n";
            
        }
        FileTool::WriteStringToFile(sourceTerm, fMNDir+"/sourceTerm/sourceTerm_HeHe_"+number+".txt" );
        
        
    }
    
    
    
    if (fStep==22) {
        
        FileTool multiNest(fMNDir+"/MultiNest/"+fMNPrefix+".txt");
        multiNest.ExtractNumberTable( fNCparam+fNWparam+2, " ", true );
        
        double chi;
        double chiMin   = multiNest.NumberTable(0, 1);
        int    indexBestFit = 0;
        for (int i = 0; i<multiNest.NumberTableGetNrows(); i++) {
            chi     = multiNest.NumberTable( i, 1 );
            if (chi<chiMin){
                chiMin       = chi;
                indexBestFit = i;
            }
        }
        
        std::vector<double> T_pbar;
        std::vector<double> q_best;
        std::vector<double> q_min;
        std::vector<double> q_max;
        std::vector<double> q_min2;
        std::vector<double> q_max2;
        
        std::string number = QString("%1").arg(multiNest.NumberTableGetNrows()-indexBestFit).rightJustified(5, '0').toStdString();
        FileTool f(fMNDir+"/sourceTerm/sourceTerm__"+number+".txt");
        f.ExtractNumberTable(  2, " ", true, 2 );
        
        for (int i=0; i<f.NumberTableGetNrows(); i++) {
            T_pbar.push_back(  f.NumberTable(i, 0)  );
            q_best.push_back(  f.NumberTable(i, 1)  );
            q_min .push_back(  f.NumberTable(i, 1)  );
            q_max .push_back(  f.NumberTable(i, 1)  );
            q_min2.push_back(  f.NumberTable(i, 1)  );
            q_max2.push_back(  f.NumberTable(i, 1)  );
        }
        
        std::string runFile;
        int cout = 0;
        for (int i = multiNest.NumberTableGetNrows()-1; i>=0; i--) {
            double chisq = multiNest.NumberTable(i, 1);
            if(chisq<=chiMin+fSigma_2_region[fNCparam+fNWparam]){
                cout++;
                std::string number = QString("%1").arg(multiNest.NumberTableGetNrows()-i).rightJustified(5, '0').toStdString();
                FileTool f(fMNDir+"/sourceTerm/sourceTerm__"+number+".txt");
                if (!f.IsFileOk()) {
                    continue;
                }
                f.ExtractNumberTable(  2, " ", true, 2 );
                for (int i=0; i<f.NumberTableGetNrows(); i++) {
                    
                    double q = f.NumberTable(i, 1);
                    if(i==30){
                        out(q)
                    }
                    if (q < q_min2.at(i)) q_min2.at(i) = q;
                    if (q > q_max2.at(i)) q_max2.at(i) = q;
                    if(chisq<=chiMin+fSigma_1_region[fNCparam+fNWparam]){
                        if (q < q_min.at(i)) q_min.at(i) = q;
                        if (q > q_max.at(i)) q_max.at(i) = q;
                    }
                }
            }
        }
        FileTool::WriteStringToFile(runFile, "runFile.txt");
        out(cout)
        out(chiMin)
        std::string sourceTerm = "#  T_pbar          p H           min 1s         max 1s        min 2s      max 2s";
        for (int i=0; i<f.NumberTableGetNrows(); i++) {
            std::stringstream ss;
            ss  << "\n" << T_pbar.at(i) << "      " <<  q_best.at(i)  << "      " <<  q_min.at(i)  << "      " <<  q_max.at(i)  << "      " <<  q_min2.at(i)  << "      " <<  q_max2.at(i);
            sourceTerm+= ss.str();
        }
        FileTool::WriteStringToFile(sourceTerm, fMNDir+"/sourceTerm.txt" );
    }
} // END main



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
            hyperonCorrection =  1./(1.+0.31+0.3/(1+pow(146./sqrtS, 2*0.9)) );
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


void Write_results( std::vector<data> &fData, FileTool &multiNest, int mode){
    
    
    out("WRITE")
    
    double chi;
    double chiMin;
    int    indexBestFit = 0;
    if(mode>0){
           chiMin = multiNest.NumberTable(0, 1);
        for (int i = 0; i<multiNest.NumberTableGetNrows(); i++) {
            chi     = multiNest.NumberTable(i, 1);
            if (chi<chiMin){
                chiMin       = chi;
                indexBestFit = i;
            }
        }
        out(indexBestFit)
        out(multiNest.NumberTableGetNrows())
    }
    
    
    
    
    
    for (int i = fNCparam; i<100; i++) {
        fParameters[i] = 0;
    }
    for (int i = fNDparam; i<100; i++) {
        fParameters_pA[i] = 0;
    }
    
    for (int e=0; e<fData.size(); e++) {
        std::string write;
        data & d= fData.at(e);
        int A = d.A;
        int freeP = fNCparam;
        if(mode==mode_C || mode==mode_CD){
            for (int k = 0; k<fNCparam; k++) {
                fParameters[k] = multiNest.NumberTable(  indexBestFit, 2+k  );
            }
            freeP = fNCparam;
        }
        if(mode==mode_CD){
            for (int k = 0; k<fNDparam; k++) {
                fParameters_pA[k] = multiNest.NumberTable(  indexBestFit, 2+fNCparam+fNWparam+k  );
            }
        }
        if(mode==mode_D){
            for (int k = 0; k<fNDparam; k++) {
                fParameters_pA[k] = multiNest.NumberTable(  indexBestFit, 2+k  );
            }
            freeP = fNDparam;
            

        }
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
        out(w_e)
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
            
            if (mode>0) {
            for (int j = multiNest.NumberTableGetNrows()-1; j>=0; j--) {
                double chisq = multiNest.NumberTable(j, 1);
                if(chisq<=chiMin+fSigma_2_region[freeP]){
                    if(mode==mode_C || mode==mode_CD){
                        for (int k = 0; k<fNCparam; k++) {
                            fParameters[k] = multiNest.NumberTable(  j, 2+k  );
                        }
                    }
                    if(mode==mode_CD){
                        for (int k = 0; k<fNDparam; k++) {
                            fParameters_pA[k] = multiNest.NumberTable(  j, 2+fNCparam+fNWparam+k  );
                        }
                    }
                    if(mode==mode_D){
                        for (int k = 0; k<fNDparam; k++) {
                            fParameters_pA[k] = multiNest.NumberTable(  j, 2+k  );
                        }
                    }
                    double CS_t  =   fun_pp_pbar_CM(sqrtS, pT, xR, fParameters)*fun_pA_pbar_factor( A, sqrtS, pT, xR, fParameters_pA);
                    if (CS_t < CS_min2) CS_min2 = CS_t;
                    if (CS_t > CS_max2) CS_max2 = CS_t;
                    if(chisq<=chiMin+fSigma_1_region[freeP]){
                        if (CS_t < CS_min) CS_min = CS_t;
                        if (CS_t > CS_max) CS_max = CS_t;
                    }
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
        if (mode==mode_CD) {
            Astr += "_minuit";
        }
        if (mode==mode_D) {
            Astr += "_scan";
        }
        if (mode<0) {
            Astr += "_ppMinuit";
        }
        if( d.in_fit ){
            FileTool::WriteStringToFile(write, "bestFit_"+Astr+"_"+d.experiment+".txt");
        }else{
            FileTool::WriteStringToFile(write, "comparison_"+Astr+"_"+d.experiment+".txt");
        }
        
    }
    
    
    
}


//***************************** MultiNest *********************

void writeParameterRange(bool pA){
    std::string toWrite = "";
    std::string prefix = "C";
    for (int i = 0; i<fNCparam; i++) {
        std::stringstream ss;
        ss << i+1;
        toWrite += QString("%1").arg( (prefix+ss.str()).c_str() ).leftJustified(30).toStdString();
        toWrite += QString("%1").arg( fC_lower[i]               ).leftJustified(30).toStdString();
        toWrite += QString("%1").arg( fC_upper[i]               ).leftJustified(30).toStdString();
        toWrite += "\n";
    }
    prefix = "w";
    for (int i = 0; i<fNWparam; i++) {
        std::stringstream ss;
        ss << i;
        toWrite += QString("%1").arg( (prefix+ss.str()).c_str() ).leftJustified(30).toStdString();
        toWrite += QString("%1").arg( 0.5                       ).leftJustified(30).toStdString();
        toWrite += QString("%1").arg( 2.0                       ).leftJustified(30).toStdString();
        toWrite += "\n";
    }
    if(!pA){
        FileTool::WriteStringToFile(toWrite, fMNDir+"/MultiNest/"+fMNPrefix+".ranges");
    }else{
        std::string prefix = "D";
        for (int i = 0; i<fNDparam; i++) {
            std::stringstream ss;
            ss << i+1;
            toWrite += QString("%1").arg( (prefix+ss.str()).c_str() ).leftJustified(30).toStdString();
            toWrite += QString("%1").arg( fD_lower[i]               ).leftJustified(30).toStdString();
            toWrite += QString("%1").arg( fD_upper[i]               ).leftJustified(30).toStdString();
            toWrite += "\n";
        }
        prefix = "w_pA";
        for (int i = 0; i<fNWparam_pA; i++) {
            std::stringstream ss;
            ss << i;
            toWrite += QString("%1").arg( (prefix+ss.str()).c_str() ).leftJustified(30).toStdString();
            toWrite += QString("%1").arg( 0.5                       ).leftJustified(30).toStdString();
            toWrite += QString("%1").arg( 2.0                       ).leftJustified(30).toStdString();
            toWrite += "\n";
        }
        FileTool::WriteStringToFile(toWrite, fMNDir+"/MultiNest/"+fMNPrefix+"_pA_minuit.ranges");
        
    }
};




