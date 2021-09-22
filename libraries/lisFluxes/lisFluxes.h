#ifndef LISFLUXES_H
#define LISFLUXES_H

#include <iostream>
#include <string>
#include <vector>
#include <math.h>

#include <definitions.h>
#include "graph.h"



namespace CRACS {
    
    /*! \brief Extract LIS Fluxes for various species.
     *
     */
    static double fLISParameters[10][7];
//    static double fLISParameters_upper[5][7];
//    static double fLISParameters_lower[5][7];
    
    class LISFluxes{
        
    public:
        static const int fProtonLIS       =   0;
        static const int fHeliumLIS       =   1;
        static const int fCarbonLIS       =   2;
        static const int fNitrogenLIS     =   3;
        static const int fOxygenLIS       =   4;
        
        static const int fLithiumLIS      =   5;
        static const int fBerylliumLIS    =   6;
        static const int fBoronLIS        =   7;

        static const int fAntiprotonLIS   =   8;
        
        static double fSM;
        
        
        static double   spectrum(double Ekin_per_n, int type=fProtonLIS);
        
        static void     applySolarModulation(Graph& flux, double SM_potential, int A, int Z, double power=1, int type=EKINPERN);
        static double   brokenPowerLaw(double E, double norm, double gamma1, double gamma2, double gamma3, double break1, double break2, double s1);
        static double   brokenPowerLaw(double *x, double *par);
        
        static void     fit                 ( bool save=false );
        
        static double   brokenPowerLaw_delta(double E, double norm, double gamma, double delta1, double delta2, double break1, double break2, double s1);
        

        static CRACS::Graph fProton_LIS;
        static CRACS::Graph fHelium_LIS;
        static CRACS::Graph fCarbon_LIS;
        static CRACS::Graph fNitrogen_LIS;
        static CRACS::Graph fOxygen_LIS;
        
        static CRACS::Graph fLithium_LIS;
        static CRACS::Graph fBeryllium_LIS;
        static CRACS::Graph fBoron_LIS;
        
        static CRACS::Graph fAntiproton_LIS;
        
        static double       fParameter[100];
        static int          fNParam;
        static double       fStart[100];
        static double       fStep [100];
        
        
        static void     fcn_pHeCNO                 ( Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag );
        static void     Ifit_pHeCNO                ( int from=0, int to=100,  int printlevel=0 );
        static void     fit_simultaneously_pHeCNO  ( bool save=false );
        
        static void     fcn_LiBeB                 ( Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag );
        static void     Ifit_LiBeB                ( int from=0, int to=100,  int printlevel=0 );
        static void     fit_simultaneously_LiBeB  ( bool save=false );
        
        static void     fcn_pbar                  ( Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag );
        static void     Ifit_pbar                 ( int from=0, int to=100, int printlevel=0 );
        static void     fit_simultaneously_pbar   ( bool save=false );
        
        
    private:
        
        //static std::vector< std::vector<double> > fLISParameters;
        
    };
    

}



//end LISFLUXES_H

#endif
