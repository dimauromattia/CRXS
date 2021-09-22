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

std::string fProgramDescription = "Program to fit cross section parameterizations";






struct data{
    int n;
    std::vector<double> sqrtS;
    std::vector<double> xR;
    std::vector<double> pT;
    std::vector<double> var_fix;
    std::vector<double> var_plot;
    std::vector<double> CS;
    std::vector<double> CS_err;
    std::vector<double> CS_err_stat;
    std::vector<double> CS_err_sys;
    std::vector<double> err_scale;
    std::vector<std::string> col_description;
};


const int           fNCparam         = 8;
std::vector<data>   fData_pp_pbar;

static double     fC_lower[100] = { 0,  0, -1, -1,  0, -1, -1,  0};
static double     fC_upper[100] = {10, 10,  1,  1, 10,  1,  1, 10};

static double     fCparam [100] = {4.499, 3.41, 0.00942, 0.445, 3.502, 0.0622, -0.247, 2.576};
static double     fSteps  [100] = {0.01 , 0.1 , 0.001,   0.001, 0.1 ,    0.01 , 0.01,  0.01 };


const int           fN               = 101;
static int          fPrint           = 0;
bool                fMinuit          = false;

int                 fNWparam         = 0;
std::vector<double> fEps;

std::vector<std::string> fFilenames;


std::vector<std::string> fFilenames_new;
std::vector<data>   fData_pp_pbar_new;

void setupExperiments(){
    
    
    fFilenames.push_back  ( "AnticicPbarNA492010.dat");
    fEps.      push_back  ( 0.065 );
    fNWparam++;
    
    fFilenames.push_back  ( "AllabyPbarCERN1970.dat");
    fEps.      push_back  ( 0.15 );
    fNWparam++;
    
    fFilenames.push_back  ( "AntreasyanPbarFNAL1979.dat");
    fEps.      push_back  ( 0.20 );
    fNWparam++;
    
    fFilenames.push_back  ( "GuetlerPbarCERNISR1976.dat");
    fEps.      push_back  ( 0.04 );
    fNWparam++;
    
    fFilenames.push_back  ( "CapiluppiPbarCERNISR1974.dat");
    fEps.      push_back  ( 0.10 );
    fNWparam++;
    
    fFilenames.push_back  ( "JohnsonPbarFNAL1978.dat");
    fEps.      push_back  ( 0.10 );
    fNWparam++;
    
    fFilenames.push_back  ( "DekkersCERN1965.dat");
    fEps.      push_back  ( 0.20 );
    fNWparam++;
    
    fFilenames.push_back  ( "BRAHMS2008.dat");
    fEps.      push_back  ( 0.10 );                         // FIX_ME, find this from paper
    fNWparam++;
    
    fFilenames_new.push_back("all_phenix.txt");
    fFilenames_new.push_back("all_guetler.txt");
    fFilenames_new.push_back("all_allaby.txt");
    fFilenames_new.push_back("all_na49.txt");
    fFilenames_new.push_back("all_antreasyan.txt");
    fFilenames_new.push_back("all_brahms.txt");
    fFilenames_new.push_back("all_johnson.txt");
    fFilenames_new.push_back("all_dekkers.txt");
    fFilenames_new.push_back("all_capiluppi.txt");
    
}


void ReadData(){
    
    
    std::string CRACS = ConfigHandler::GetInstance()->SoftwarePath();
    
    for (int i=0; i<fFilenames.size(); i++) {
        data d;
        d.n=0;
        FileTool file( CRACS+"/data/CS_measurements/"+fFilenames.at(i) );
        file.ExtractNumberTable(5, " ");
        for (int j = 0; j<file.NumberTableGetNrows(); j++) {
            d.n ++;
            d.sqrtS .push_back(file.NumberTable(j, 0));
            d.pT    .push_back(file.NumberTable(j, 1));
            d.xR    .push_back(file.NumberTable(j, 2));
            d.CS    .push_back(file.NumberTable(j, 3));
            d.CS_err.push_back(file.NumberTable(j, 4));
        }
        out(d.n)
        fData_pp_pbar.push_back(d);
    }
    
}

void ReadData_new(){
    
    
    std::string CRACS = ConfigHandler::GetInstance()->SoftwarePath();
    
    for (int i=0; i<fFilenames_new.size(); i++) {
        data d;
        d.n=0;
        FileTool file( fFilenames_new.at(i) );
        file.ExtractNumberTable( 9 , " " );
        for (int j = 0; j<file.NumberTableGetNrows(); j++) {
            d.n ++;
            d.sqrtS         .push_back(file.NumberTable(j, 0));
            d.pT            .push_back(file.NumberTable(j, 1));
            d.xR            .push_back(file.NumberTable(j, 2));
            d.var_fix       .push_back(file.NumberTable(j, 3));
            d.var_plot      .push_back(file.NumberTable(j, 4));
            d.CS            .push_back(file.NumberTable(j, 5));
            d.CS_err_stat   .push_back(file.NumberTable(j, 6));
            d.CS_err_sys    .push_back(file.NumberTable(j, 7));
            d.err_scale     .push_back(file.NumberTable(j, 8));
            d.CS_err        .push_back( sqrt(d.CS_err_sys.at(j)*d.CS_err_sys.at(j)+d.CS_err_stat.at(j)*d.CS_err_stat.at(j)) );
        }
        
        for (int j = 0; j<file.GetNumberOfLines(); j++) {
            QString line = file.GetLine(j).c_str();
            if (line.left(1)=="*"){
                QStringList list = line.split(" ", QString::SkipEmptyParts);
                for (int k = 0; k<list.size()-1; k++) {
                    d.col_description.push_back(list.at(k+1).toStdString());
                }

                break;
            }
            
            
            
                
        }
        
        out(d.n)
        fData_pp_pbar_new.push_back(d);
    }
    
}




#include "TMinuit.h"

TMinuit *gMinuit;
double  fArglist[10];
int     fIerflg = 0;


Double_t fun_pp_pbar_CM__diMauro12(double sqrtS,  double pT, double xR, Double_t *par)
{
    double s            =   sqrtS*sqrtS;
    double E_pbar_Max   =   ( s-8.*fMass_proton*fMass_proton )/2./sqrtS;
    double E_pbar       =   xR*E_pbar_Max;
    double pT_pbar      =   pT;
    return ppCSParametrizations::inv_pp_pbar_CM__diMauro(s, E_pbar, pT_pbar, par[0], par[1], par[2], par[3], par[4], par[5], par[6], par[7], 0, 0, 0);
}


void fcn_pp_pbar_CM__diMauro12(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    double* wpar=par;
    
    Double_t chisq  = 0;
    double   delta;
    int ndf = 0;
    
    for (int e=0; e<fData_pp_pbar.size(); e++) {
        data d= fData_pp_pbar.at(e);
        double chisq_e = 0;
        for (int i=0;i<d.n; i++) {
            double sqrtS    =   d.sqrtS[i];
            double xR       =   d.xR   [i];
            double pT       =   d.pT   [i];
            delta           =   pow( (wpar[e]*d.CS[i] - fun_pp_pbar_CM__diMauro12(sqrtS, pT, xR, fCparam))/d.CS_err[i]/wpar[e], 2 );
            chisq          +=   delta;
//            if (fPrint) {
//                chisq_e    +=   delta;
//                if (delta>30) {
//                    std::cout << warnout << std::endl;
//                    out(fFilenames.at(e))
//                    int line = i;
//                    out(line)
//                    out(delta)
//                    out(wpar[e])
//                    out(sqrtS)
//                    out(xR)
//                    out(pT)
//                    out(d.CS[i])
//                    out(d.CS_err[i])
//                    out(fun_pp_pbar_CM__diMauro12(sqrtS, pT, xR, fCparam))
//                    std::cout << warnend << std::endl;
//                }
//            }
            
        }
        chisq          +=   pow( (1-wpar[e])/fEps[0], 2);
        if(fPrint){
            for (int i=0; i<fNCparam; i++) {
                out(fCparam[i])
            }
            for (int i=0; i<fNWparam; i++) {
                out(wpar[i])
            }
            std::cout << "*************************" << std::endl;
            chisq_e     += pow( (1-wpar[e])/fEps[0], 2);
            out(fFilenames.at(e))
            out(chisq_e)
            out(d.n)
            double chisq_ndf = chisq_e/d.n;
            out(chisq_ndf)
            std::cout << "*************************" << std::endl;
        }
        
        ndf+= d.n;
//        chisq          +=   pow( log(wpar[e])/log(1+fEps[e]), 2);
        
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




void chisq_minNorm(int e, double* par, double &chisq_min, double &w_min)
{
    chisq_min  = 1e90;
    
    double   delta, chisq;
    
    for (double w=0.5; w<2; w+=0.005) {
        chisq = 0;
        
        for (int i=0;i<fData_pp_pbar.at(e).n; i++) {
            double sqrtS    =   fData_pp_pbar.at(e).sqrtS[i];
            double xR       =   fData_pp_pbar.at(e).xR   [i];
            double pT       =   fData_pp_pbar.at(e).pT   [i];
            double CS       =   fData_pp_pbar.at(e).CS   [i];
            double CS_err   =   fData_pp_pbar.at(e).CS_err[i];
            delta           =   pow( (w*CS - fun_pp_pbar_CM__diMauro12(sqrtS, pT, xR, par))/CS_err/w, 2 );
            chisq          +=   delta;
        }
        
        chisq          +=   pow( (1-w)/fEps[0], 2);
        
        if (chisq<chisq_min) {
            chisq_min   = chisq;
            w_min       = w;
        }
        
    }
    
}

void chisq_minNorm_clever(int e, double* par, double &chisq_e, double &w_e)
{
    
    double   delta;
    
    double w_min = 0.5;
    double w_max = 2.0;
    
    double chisq_min  = 1e90;
    double chisq_max  = 1e90;
    
    
    w_e = 0.5;
    chisq_e = 0;
    for (int i=0;i<fData_pp_pbar.at(e).n; i++) {
        double sqrtS    =   fData_pp_pbar.at(e).sqrtS[i];
        double xR       =   fData_pp_pbar.at(e).xR   [i];
        double pT       =   fData_pp_pbar.at(e).pT   [i];
        double CS       =   fData_pp_pbar.at(e).CS   [i];
        double CS_err   =   fData_pp_pbar.at(e).CS_err[i];
        delta           =   pow( (w_e*CS - fun_pp_pbar_CM__diMauro12(sqrtS, pT, xR, par))/CS_err/w_e, 2 );
        chisq_e          +=   delta;
    }
    chisq_e          +=   pow( (1-w_e)/fEps[0], 2);
    chisq_min = chisq_e;

    w_e = 2;
    chisq_e = 0;
    for (int i=0;i<fData_pp_pbar.at(e).n; i++) {
        double sqrtS    =   fData_pp_pbar.at(e).sqrtS[i];
        double xR       =   fData_pp_pbar.at(e).xR   [i];
        double pT       =   fData_pp_pbar.at(e).pT   [i];
        double CS       =   fData_pp_pbar.at(e).CS   [i];
        double CS_err   =   fData_pp_pbar.at(e).CS_err[i];
        delta           =   pow( (w_e*CS - fun_pp_pbar_CM__diMauro12(sqrtS, pT, xR, par))/CS_err/w_e, 2 );
        chisq_e          +=   delta;
    }
    chisq_e          +=   pow( (1-w_e)/fEps[0], 2);
    chisq_max = chisq_e;

    
    
    for (int i = 0; i<20; i++) {
        
        w_e = (w_min+w_max)/2;

        chisq_e = 0;
        
        for (int i=0;i<fData_pp_pbar.at(e).n; i++) {
            double sqrtS    =   fData_pp_pbar.at(e).sqrtS[i];
            double xR       =   fData_pp_pbar.at(e).xR   [i];
            double pT       =   fData_pp_pbar.at(e).pT   [i];
            double CS       =   fData_pp_pbar.at(e).CS   [i];
            double CS_err   =   fData_pp_pbar.at(e).CS_err[i];
            delta           =   pow( (w_e*CS - fun_pp_pbar_CM__diMauro12(sqrtS, pT, xR, par))/CS_err/w_e, 2 );
            chisq_e        +=   delta;
        }
        
        chisq_e          +=   pow( (1-w_e)/fEps[0], 2);

        out(chisq_max)
        out(chisq_min)
        out(chisq_e)
        out("")
        
        if ( chisq_max<chisq_min ) {
            chisq_min = chisq_e;
            w_min = w_e;
        }else{
            chisq_max = chisq_e;
            w_max = w_e;
        }
        
        
        
    }
    
    out(w_e)
}



//______________________________________________________________________________



#include <multinest.h>


// Multinest defaulf values
std::string fMNPrefix       = "CS_fit";
std::string fMNDir          = ".";
int         fNlive          = 500;
double      fEfr            = 0.5;
double      fTol            = 0.1;
int         fNdims          = -1;
int         fNpar           = -1;

void writeParameterRange(){
    
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
    
    FileTool::WriteStringToFile(toWrite, fMNDir+"/MultiNest/"+fMNPrefix+".ranges");

};

/***********************************************************************************************************************/

double fSigma_1_region[] = {1.000043426, 2.295815161, 3.526822180, 4.719568761, 5.887700474, 7.038515492, 8.176359741, 9.304044023, 10.42350189, 11.53612748, 12.64296378, 13.74481437, 14.84231367, 15.93597259, 17.02620974, 18.11337308, 19.19775552, 20.27960639, 21.35913992, 22.43654182};

double fSigma_2_region[] = {4.000009776, 6.180085905, 8.024894670, 9.715641142, 11.31387084, 12.84885057, 14.33712678, 15.78911024, 17.21184692, 18.61036514, 19.98840090, 21.34881933, 22.69387460, 24.02537791, 25.34481051, 26.65340209, 27.95218680, 29.24204418, 30.52372968, 31.79789790};


/******************************************** loglikelihood routine ****************************************************/

// Now an example, sample an egg box likelihood

// Input arguments
// ndim 						= dimensionality (total number of free parameters) of the problem
// npars 						= total number of free plus derived parameters
// context						void pointer, any additional information
//
// Input/Output arguments
// Cube[npars]                  = on entry has the ndim parameters in unit-hypercube on exit, the physical parameters
//                                plus copy any derived parameters you want to store with the free parameters
//
// Output arguments
// lnew 						= loglikelihood

double fChisqMin = 1e90;

int fCount = 1;
void LogLike(double *Cube, int &ndim, int &npars, double &lnew, void *context)
{
    for (int i = 0; i<ndim; i++) {
        Cube[i] = fC_lower[i] + (  fC_upper[i]-fC_lower[i] )*Cube[i];
        fCparam[i] = Cube[i];
    }

    double chisq = 0;
    int n; double h;
    if(fMinuit){
        // Fit normalization
        fArglist[0] = 1000000;
        fArglist[1] = 1.;
        gMinuit->mnexcm("MIGRAD", fArglist ,2,fIerflg);
        
        
        
        for (int i=0; i<fNWparam; i++) {
            double p, pErr;
            gMinuit->GetParameter(i, p, pErr);
            Cube[ndim+i] = p;
        }
        
        fcn_pp_pbar_CM__diMauro12(n, &h, chisq, &Cube[ndim], 1);
        
        
        if (chisq<fChisqMin){
            fChisqMin = chisq;
            out(fChisqMin)
            for(int i=0; i<npars; i++){
                out(Cube[i])
            }
            for (int i=0; i<fNWparam; i++) {
                double p, pErr;
                gMinuit->GetParameter(i, p, pErr);
                out(Cube[ndim+i])
            }
        }
    }else{
        
        
        for (int e=0; e<fNWparam; e++) {
            double chisq_e;
            double w_e;
            chisq_minNorm(e, Cube, chisq_e, w_e);
            chisq+=chisq_e;
            Cube[ndim+e] = w_e;
        }
        if (chisq<fChisqMin){
            fChisqMin = chisq;
            out(fChisqMin)
            out(fCount)
        }
    }
    
    fCount++;
    
    // Return new loglikelihood
    lnew = -chisq/2.0;

    
}

/************************************************* dumper routine ******************************************************/

// Not important here!
// Standard function taken from the multinest example

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

/***********************************************************************************************************************/



/************************************************** Main program *******************************************************/


int fStep       = 1;
int fMNPoint    = 1;
int fMod        = 30;

static double fParameters[100];


static double inv_pp_pbar_CM( double s, double E_pbar, double pT_pbar){
    return ppCSParametrizations::inv_pp_pbar_CM__diMauro(s, E_pbar, pT_pbar, fParameters[0], fParameters[1], fParameters[2], fParameters[3], fParameters[4], fParameters[5], fParameters[6], fParameters[7], fParameters[8], fParameters[9], fParameters[10]);
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

CS_lab_tab* fTab;

double dT_pp_pbar_LAB_interpolation( double Tn_proj, double Tn_prod ){
    return fTab->GetInterpolation( Tn_proj, Tn_prod );
}

double protonFlux(double T_p){
    return LISFluxes::spectrum(T_p, LISFluxes::fProtonLIS);
}

int main(int argc, char *argv[])
{
    
    
    
    // Config handling
    
    ConfigHandler* config = ConfigHandler::GetInstance();
    config->SetInput(argc, argv, fProgramDescription);
    
    config->AddOptionInt        (  "step",          fStep,         "default: 1"                          );
    config->AddOptionInt        (  "point",         fMNPoint,      "default: 1"                          );
    config->AddOptionInt        (  "mod",           fMod,          "default: 30"                         );
    
    config->AddOptionInt        (  "nlive",         fNlive,        "default: 500"                        );
    config->AddOptionDouble     (  "efr",           fEfr,          "default: 0.5"                        );
    config->AddOptionDouble     (  "tol",           fTol,          "default: 0.1"                        );
    
    config->AddOption           (  "dir",           fMNDir,        "default: ."                          );
    config->AddOption           (  "prefix",        fMNPrefix,     "default: CS_fit"                     );
    
    config->CheckOption();
    
    out(fStep)
    
    out("Read")
    setupExperiments();
    ReadData();
    ReadData_new();
    
    fNdims = fNCparam;
    fNpar  = fNCparam+fNWparam;
    
    std::system(  ("mkdir "+fMNDir+"/MultiNest").c_str() );
    writeParameterRange();
    
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
    
    
    if (fStep==2) {
        
        
        FileTool multiNest(fMNDir+"/MultiNest/"+fMNPrefix+".txt");
        multiNest.ExtractNumberTable(fNpar+2, " ", true);
        
        
        double chi;
        double chiMin   = multiNest.NumberTable(0, 1);
        int    indexBestFit = 0;
        for (int i = 0; i<multiNest.GetNumberOfLines(); i++) {
            chi     = multiNest.NumberTable(i, 1);
            if (chi<chiMin){
                chiMin       = chi;
                indexBestFit = i;
            }
        }
        
        out(indexBestFit)
        out(multiNest.GetNumberOfLines())
        
        std::string runFile;
        int cout = 0;
        for (int i = multiNest.GetNumberOfLines()-1; i>=0; i--) {
            if( i!=indexBestFit && i%fMod!=0 )
                continue;
            double chisq = multiNest.NumberTable(i, 1);
            if(chisq<=chiMin+fSigma_2_region[fNCparam+fNWparam]){
                cout++;
                runFile += "CRACS_crossSectionMultiNest --step 3 --point "+QString("%1").arg(multiNest.GetNumberOfLines()-i).toStdString()+" --dir "+fMNDir+" --prefix "+fMNPrefix+" \n";
            }
        }
        FileTool::WriteStringToFile(runFile, "runFile.txt");
        out(cout)
        
        out(chiMin)
        
        double C[fNCparam];
        double w[fNWparam];
        for (int i=0; i<fNCparam; i++) {
            C[i] = multiNest.NumberTable (indexBestFit, 2+i );
            out(C[i])
        }
        for (int i=0; i<fNWparam; i++) {
            w[i] = multiNest.NumberTable( indexBestFit, 2+fNCparam+i );
            out(w[i])
        }
        
        
        for (int i = fNdims; i<100; i++) {
            fParameters[i] = 0;
        }
            
        // write some output
        
        for (int e=0; e<fFilenames.size(); e++) {
            
            std::string write;
            out(fFilenames.at(e))
            out(fData_pp_pbar.size())
            data d= fData_pp_pbar.at(e);
            
            
            double w_e = w[e];
            
            for (int i=0;i<d.n; i++) {
                double sqrtS    =   d.sqrtS [i];
                double xR       =   d.xR    [i];
                double pT       =   d.pT    [i];
                double CS       =   d.CS    [i];
                double CS_err   =   d.CS_err[i];
                
                double CS_p     =   fun_pp_pbar_CM__diMauro12(sqrtS, pT, xR, C);
                double CS_min   =   fun_pp_pbar_CM__diMauro12(sqrtS, pT, xR, C);
                double CS_max   =   fun_pp_pbar_CM__diMauro12(sqrtS, pT, xR, C);
                double CS_min2  =   fun_pp_pbar_CM__diMauro12(sqrtS, pT, xR, C);
                double CS_max2  =   fun_pp_pbar_CM__diMauro12(sqrtS, pT, xR, C);
                
                for (int j = multiNest.GetNumberOfLines()-1; j>=0; j--) {
                    double chisq = multiNest.NumberTable(j, 1);
                    if(chisq<=chiMin+fSigma_2_region[fNCparam+fNWparam]){
                        for (int k = 0; k<fNdims; k++) {
                            fParameters[k] = multiNest.NumberTable(  j, 2+k  );
                        }
                        double CS_t  =   fun_pp_pbar_CM__diMauro12(sqrtS, pT, xR, fParameters);
                        
                        if (CS_t < CS_min2) CS_min2 = CS_t;
                        if (CS_t > CS_max2) CS_max2 = CS_t;
                        if(chisq<=chiMin+fSigma_1_region[fNCparam+fNWparam]){
                            if (CS_t < CS_min) CS_min = CS_t;
                            if (CS_t > CS_max) CS_max = CS_t;
                        }
                        
                    }
                    
                    
                }
                
                
                double T_p_LAB, T_pbar_LAB, eta;
                
                CSTransformations::convert_CM_to_LAB(sqrtS*sqrtS, xR, pT, T_p_LAB, T_pbar_LAB, eta, true);
                
                
                write += QString("%1").arg(sqrtS  ).leftJustified(20).toStdString();
                write += QString("%1").arg(pT     ).leftJustified(20).toStdString();
                write += QString("%1").arg(xR     ).leftJustified(20).toStdString();
                write += QString("%1").arg(CS     ).leftJustified(20).toStdString();
                write += QString("%1").arg(CS_err ).leftJustified(20).toStdString();
                
                CSTransformations::convert_CM_to_LAB(sqrtS*sqrtS, xR, pT, T_p_LAB, T_pbar_LAB, eta, true);
                write += QString("%1").arg(T_p_LAB      ).leftJustified(20).toStdString();
                write += QString("%1").arg(T_pbar_LAB   ).leftJustified(20).toStdString();
                write += QString("%1").arg(eta          ).leftJustified(20).toStdString();
                
                write += QString("%1").arg(CS_p   ).leftJustified(20).toStdString();
                write += QString("%1").arg(CS_min ).leftJustified(20).toStdString();
                write += QString("%1").arg(CS_max ).leftJustified(20).toStdString();
                write += QString("%1").arg(CS_min2).leftJustified(20).toStdString();
                write += QString("%1").arg(CS_max2).leftJustified(20).toStdString();
                write += QString("%1").arg(w_e    ).leftJustified(20).toStdString();
                
                write += "\n";
                
            }
            
            FileTool::WriteStringToFile(write, "bestFit_"+fFilenames.at(e));
            
        }

        // write some output
        
        for (int e=0; e<fFilenames_new.size(); e++) {
            
            std::string write;
           
            data d= fData_pp_pbar_new.at(e);
            
            write += "*  ";
            for (int k = 0; k<d.col_description.size(); k++) {
                write += QString("%1").arg( d.col_description.at(k).c_str()  ).leftJustified(20).toStdString();
            }
            
            write += QString("%1").arg("CS_parametrization").leftJustified(20).toStdString();
            write += QString("%1").arg("CS_min_1sigma" ).leftJustified(20).toStdString();
            write += QString("%1").arg("CS_max_1sigma" ).leftJustified(20).toStdString();
            write += QString("%1").arg("CS_min_2sigma").leftJustified(20).toStdString();
            write += QString("%1").arg("CS_max_2sigma").leftJustified(20).toStdString();
            write += QString("%1").arg("data rescale"    ).leftJustified(20).toStdString();
            
            write += "\n";
            
            
            double w_e = 1.0;
            
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
                
                double CS_p     =   fun_pp_pbar_CM__diMauro12(sqrtS, pT, xR, C);
                double CS_min   =   fun_pp_pbar_CM__diMauro12(sqrtS, pT, xR, C);
                double CS_max   =   fun_pp_pbar_CM__diMauro12(sqrtS, pT, xR, C);
                double CS_min2  =   fun_pp_pbar_CM__diMauro12(sqrtS, pT, xR, C);
                double CS_max2  =   fun_pp_pbar_CM__diMauro12(sqrtS, pT, xR, C);
                
                for (int j = multiNest.GetNumberOfLines()-1; j>=0; j--) {
                    double chisq = multiNest.NumberTable(j, 1);
                    if(chisq<=chiMin+fSigma_2_region[fNCparam+fNWparam]){
                        for (int k = 0; k<fNdims; k++) {
                            fParameters[k] = multiNest.NumberTable(  j, 2+k  );
                        }
                        double CS_t  =   fun_pp_pbar_CM__diMauro12(sqrtS, pT, xR, fParameters);
                        
                        if (CS_t < CS_min2) CS_min2 = CS_t;
                        if (CS_t > CS_max2) CS_max2 = CS_t;
                        if(chisq<=chiMin+fSigma_1_region[fNCparam+fNWparam]){
                            if (CS_t < CS_min) CS_min = CS_t;
                            if (CS_t > CS_max) CS_max = CS_t;
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
            
            FileTool::WriteStringToFile(write, "bestFit_"+fFilenames_new.at(e));
            
        }

        
    }
    
    
    if (fStep==3) {
        
        out(fMNPoint)
        
        FileTool multiNest(fMNDir+"/MultiNest/"+fMNPrefix+".txt");
        multiNest.ExtractNumberTable(fNpar+2, " ", true);
        
        for (int i = 0; i<fNdims; i++) {
            fParameters[i] = multiNest.NumberTable(  multiNest.GetNumberOfLines()-fMNPoint, 2+i  );
        }
        for (int i = fNdims; i<100; i++) {
            fParameters[i] = 0;
        }
        
        std::system( ("mkdir "+fMNDir+"/sourceTerm").c_str() );
     
        fTab = CS_lab_tab::GetInstance();
        std::string number = QString("%1").arg(fMNPoint).rightJustified(5, '0').toStdString();
        fTab->WriteCS(  dT_pp_pbar_LAB, "sourceTerm/dT_pp_pbar_LAB__"+number+".txt", fMNDir  );
        
        
        // Get source term
        fTab->ReadCS (  "sourceTerm/dT_pp_pbar_LAB__"+number+".txt", fMNDir  );
        
        double nH   = 1e-6;
        
        LISFluxes::fit();
        
        
        std::string sourceTerm = "T_pbar          p H ";
        
        for (double dT=-1; dT<4; dT+=1./30) {
            double T_pbar   = pow(10, dT);
            double Tn_from  = std::min(T_pbar, 6*fMass_proton);
            
            //out(protonFlux(T_pbar))
            
            double S_p_H    = SourceTerms::integrationFromTo(T_pbar, dT_pp_pbar_LAB_interpolation,   protonFlux, nH,  Tn_from, 1e6, 10000, NO_OUT);
            
            std::stringstream ss;
            ss  << "\n" << T_pbar << "      " <<  S_p_H;
            sourceTerm+= ss.str();
            
        }
        FileTool::WriteStringToFile(sourceTerm, fMNDir+"/sourceTerm/sourceTerm__"+number+".txt" );
        
        
    
    }
    
    
    
    if (fStep==4) {
        
        FileTool multiNest(fMNDir+"/MultiNest/"+fMNPrefix+".txt");
        multiNest.ExtractNumberTable(fNpar+2, " ", true);
        
        
        double chi;
        double chiMin   = multiNest.NumberTable(0, 1);
        int    indexBestFit = 0;
        for (int i = 0; i<multiNest.GetNumberOfLines(); i++) {
            chi     = multiNest.NumberTable(i, 1);
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
        
        std::string number = QString("%1").arg(multiNest.GetNumberOfLines()-indexBestFit).rightJustified(5, '0').toStdString();
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
        for (int i = multiNest.GetNumberOfLines()-1; i>=0; i--) {
            double chisq = multiNest.NumberTable(i, 1);
            if(chisq<=chiMin+fSigma_2_region[fNCparam+fNWparam]){
                cout++;
                
                std::string number = QString("%1").arg(multiNest.GetNumberOfLines()-i).rightJustified(5, '0').toStdString();
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
    
    
    
    
    
    if (fStep==100) {
        
        
        
        gMinuit = new TMinuit(fNWparam);
        gMinuit->SetPrintLevel(-1);
        gMinuit->SetFCN( fcn_pp_pbar_CM__diMauro12 );
        
        
        fArglist[0] = 1;
        gMinuit->mnexcm("SET ERR", fArglist ,1,fIerflg);
        
        // Set starting values and step sizes for parameters
        
        std::string prefix = "W";
        for (int i = 0; i<fNWparam; i++) {
            std::stringstream ss;
            ss << i;
            gMinuit->mnparm(i, (prefix+ss.str()).c_str(), 1., 0.01, 0,0,fIerflg);
        }
        //gMinuit->SetPrintLevel(0);
        
        
        fArglist[0] = 1000000;
        fArglist[1] = 1.;
        gMinuit->mnexcm("MIGRAD", fArglist ,2,fIerflg);
        
        for (int i=0; i<fNWparam; i++) {
            double p, pErr;
            gMinuit->GetParameter(i, p, pErr);
            fParameters[i] = p;
            out(p)
        }
        
        int n; double h; double chisq;
        
        fcn_pp_pbar_CM__diMauro12(n, &h, chisq, fParameters, 1);
        
        for (int e=0; e<fNWparam; e++) {
            double chisq_e;
            double w_e;
            chisq_minNorm_clever(e, fCparam, chisq_e, w_e);
            chisq+=chisq_e;
        }
        
    }
    
    
    
    
    
    
    
    
}


