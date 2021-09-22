#include "sourceTerm.h"

#include <iostream>
#include <sstream>

#include <definitions.h>
#include <configHandler.h>
#include <readFluxes.h>
#include <graph.h>

#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"


double CRACS::SourceTerms::integration( double T_product, double(*dT_production_CS)(double, double), double (*sourceFlux_Tn)(double), double nH, double T_thresh, double d_log_T_source, double cut_order, int output ){
    
    std::stringstream ss;  std::streambuf * old_buf = std::cout.rdbuf(ss.rdbuf());
    
    double source_Trapez  = 0;
    double source_Upper   = 0;
    double source_Lower   = 0;
    
    //std::stringstream sss;
    T_thresh = std::max(T_thresh,T_product);
    
    double T_source         = T_thresh;
    double integrand_lower  = T_source*dT_production_CS(T_source, T_product)*sourceFlux_Tn(T_source);
    double integrand_upper  = 0;
    for (double dlogT=log(T_thresh)+d_log_T_source; dlogT<=log(T_thresh*pow(10, cut_order)); dlogT+=d_log_T_source) {
        T_source         = exp(dlogT);
        
        double d1 = dT_production_CS(T_source, T_product);
        double d2 = sourceFlux_Tn(T_source);
        integrand_upper  = T_source*d1*d2;
        
        //sss << "    " << T_product   << "    " << T_source   << "    " << d1   << "    " << d2  << std::endl;
        
        source_Trapez  += sqrt(integrand_lower*integrand_upper);
        source_Upper   += std::max(integrand_lower, integrand_upper);
        source_Lower   += std::min(integrand_lower, integrand_upper);
        
        integrand_lower = integrand_upper;
    }
    source_Trapez *= 4*M_PI*nH*d_log_T_source;
    source_Upper  *= 4*M_PI*nH*d_log_T_source;
    source_Lower  *= 4*M_PI*nH*d_log_T_source;
    
    std::cout.rdbuf(old_buf);
    std::string output_summary = ss.str();
    
    if ( output_summary.find(warning_momentum_integration) != std::string::npos && output>NO_OUT ) {
        std::cout << warnout << "In at least one CS evaluation:" << warnend << std::endl;
        std::cout << "   " << warning_momentum_integration << std::endl;
    }
    if ( output_summary.find(warning_parametrization_not_available) != std::string::npos && output>NO_OUT ) {
        std::cout << warnout << "In at least one CS evaluation:" << warnend << std::endl;
        std::cout << "   " << warning_parametrization_not_available << std::endl;
    }
    
    
    //std::cout << sss.str() << std::endl;
    
    double error = std::max( source_Upper-source_Trapez, source_Trapez-source_Lower );
    if(source_Trapez>0 && error/source_Trapez>0.01 && output>NO_OUT){
        std::cout << warning_sourceterm_integration << std::endl;
    }
    if (output==ALL_OUT) {
        std::cout << "Trapez Sum:      " << source_Trapez << std::endl;
        std::cout << "Upper Sum:       " << source_Upper  << std::endl;
        std::cout << "Lower Sum:       " << source_Lower  << std::endl;
        std::cout << "Relative Error:  " << error/source_Trapez  << std::endl;
    }
    
    return source_Trapez;
}

double CRACS::SourceTerms::integrationFromTo(double T_product, double (*dT_production_CS)(double, double), double (*sourceFlux_Tn)(double), double nISM, double Tn_from, double Tn_to, int n, int output ){
    
    funOut(SourceTerms::integrationFromTo)
    
    std::stringstream ss;
    //std::streambuf * old_buf = std::cout.rdbuf(ss.rdbuf());
    
    double source_Trapez  = 0;
    double source_Upper   = 0;
    double source_Lower   = 0;
    
    //std::stringstream sss;
    double Tn_thresh = std::max(Tn_from,T_product);
    double T_source         = Tn_thresh;
    double integrand_lower  = T_source*dT_production_CS(T_source, T_product)*sourceFlux_Tn(T_source);
    double integrand_upper  = 0;
    
    varOut(Tn_from)
    varOut(Tn_thresh)
    if (Tn_to<=Tn_thresh){
        //std::cout.rdbuf(old_buf);
        return 0;
    }
    
    double d_log_T_source = log(Tn_to/Tn_thresh)/n;
    varOut(d_log_T_source)
    for (double dlogT=log(Tn_thresh)+d_log_T_source; dlogT<log(Tn_to); dlogT+=d_log_T_source) {
        T_source         = exp(dlogT);
        
        double d1 = dT_production_CS(T_source, T_product);
        double d2 = sourceFlux_Tn(T_source);
        integrand_upper  = T_source*d1*d2;
        
        //sss << "    " << T_product   << "    " << T_source   << "    " << d1   << "    " << d2  << std::endl;
        
        source_Trapez  += sqrt(integrand_lower*integrand_upper);
        if(source_Trapez!=source_Trapez){
            varOut(d1)
            varOut(d2)
            varOut(T_source)
            varOut(integrand_lower)
            varOut(integrand_upper)
            return 0;
        }
        source_Upper   += std::max(integrand_lower, integrand_upper);
        source_Lower   += std::min(integrand_lower, integrand_upper);
        
        integrand_lower = integrand_upper;
    }
    source_Trapez *= 4*M_PI*nISM*d_log_T_source;
    source_Upper  *= 4*M_PI*nISM*d_log_T_source;
    source_Lower  *= 4*M_PI*nISM*d_log_T_source;
    
    //std::cout.rdbuf(old_buf);
    std::string output_summary = ss.str();
    
    if ( output_summary.find(warning_momentum_integration) != std::string::npos && output>NO_OUT ) {
        std::cout << warnout << "In at least one CS evaluation:" << warnend << std::endl;
        std::cout << "   " << warning_momentum_integration << std::endl;
    }
    if ( output_summary.find(warning_parametrization_not_available) != std::string::npos && output>NO_OUT ) {
        std::cout << warnout << "In at least one CS evaluation:" << warnend << std::endl;
        std::cout << "   " << warning_parametrization_not_available << std::endl;
    }
    
    
    //std::cout << sss.str() << std::endl;
    
    double error = std::max( source_Upper-source_Trapez, source_Trapez-source_Lower );
    if(source_Trapez>0 && error/source_Trapez>0.01 && output>NO_OUT){
        std::cout << warning_sourceterm_integration << std::endl;
    }
    if (output==ALL_OUT) {
        std::cout << "Trapez Sum:      " << source_Trapez << std::endl;
        std::cout << "Upper Sum:       " << source_Upper  << std::endl;
        std::cout << "Lower Sum:       " << source_Lower  << std::endl;
        std::cout << "Relative Error:  " << error/source_Trapez  << std::endl;
    }
    
    return source_Trapez;
    
}
