#ifndef READCSDATA_H
#define READCSDATA_H

#include <iostream>
#include <string>
#include <vector>
#include <functional>

#include "graph.h"

#include "definitions.h"


namespace CRACS {
    
    /*! \brief Class to read the AMS fluxes in thins.
     *
     */

    class ReadCSData{
        
    public:
        //! Constructor. Reads the fluxes.
        /*!
        */
        ReadCSData(){Read();};
        
        
        
    private:
        
        //! Read fluxes.
        /*!
         */
        void Read();
        
        
        std::vector<double>       fPP_sqrtS;
        std::vector<double>       fPP_xR;
        std::vector<double>       fPP_pT;
        
        std::vector<double>       fPP_invCS;
        std::vector<double>       fPP_invCS_err;
        
        
        
    };
    

}



//end ReadCSData_H

#endif
