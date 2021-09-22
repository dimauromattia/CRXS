#ifndef SOURCETERMS_H
#define SOURCETERMS_H

#include <iostream>
#include <string>
#include <vector>
#include <math.h>

#include <definitions.h>



namespace CRACS {
    
    
    /*! \brief Extract LIS Fluxes for various species.
     *
     */
    
    class SourceTerms{
        
    public:
        static double integration( double T_product, double(*dT_production_CS)(double, double), double (*sourceFlux_Tn)(double), double nH,
                                  double T_thresh = 6*fMass_proton,
                                  double d_log_T_source=0.1,
                                  double cut_order=7,
                                  int output=WARN_OUT );
        static double integrationFromTo( double T_product, double(*dT_production_CS)(double, double), double (*sourceFlux_Tn)(double),
                                        double nISM,
                                        double Tn_from,
                                        double Tn_to,
                                        int precision=10000,
                                        int output=WARN_OUT );
        
        
    private:
        
        //static std::vector< std::vector<double> > fLISParameters;
        
    };
    

}



//end SOURCETERMS_H

#endif
