#include <iostream>
#include <sstream>

#include <definitions.h>
#include <configHandler.h>
#include <readFluxes.h>
#include <lisFluxes.h>
#include <crossSections.h>
#include <fileTools.h>

#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TApplication.h"

#include "QString"
#include "QStringList"

using namespace  CRACS;



std::string fProgramDescription = "Program to fit the LIS spectra of proton, helium, and antiproton fluxes";

double oxygenFlux(double T_p){
    return LISFluxes::spectrum(T_p, LISFluxes::fOxygenLIS);
}


int main(int argc, char *argv[])
{
    // Config handling
    ConfigHandler* config = ConfigHandler::GetInstance();
    config->SetInput(argc, argv, fProgramDescription);
    
    config->CheckOption();
    
    
    LISFluxes::fit(true);
    
    
    std::string s = "#   T/n                 flux (T/n)^2.7";
    for (double d=0; d<4; d+=0.1) {
        double Tn       = pow(10, d);
        double flux     = oxygenFlux(Tn);
        
        out(flux)
        
        s += "\n" +  QString("%1").arg(Tn  ).leftJustified(20).toStdString();
        s +=         QString("%1").arg(flux).leftJustified(20).toStdString();
        
    }
    FileTool::WriteStringToFile(s, "oxygen_fit.txt");
   
    return 0;
    
    
}


