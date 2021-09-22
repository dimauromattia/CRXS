#ifndef READFLUXES_H
#define READFLUXES_H

#include <iostream>
#include <string>
#include <vector>
#include <functional>

#include "graph.h"

#include "definitions.h"


namespace CRACS {
    
    double getBinExpectationValueForPowerLaw2( double binStart, double binEnd, double power  );
    
    /*! \brief Class to read the AMS fluxes in thins.
     *
     */

    class ReadFluxes{
        
    public:
        //! Constructor. Reads the fluxes.
        /*!
        */
        ReadFluxes(){Read();};
        CRACS::Graph GetSpectrum    (std::string species, std::string experiment, spectrumType type);
        
    private:
        
        //! Read fluxes.
        /*!
         */
        void Read();
        
        CRACS::Graph ReadAMS            (std::string filename);
        CRACS::Graph ReadAMS_R          (std::string filename);
        CRACS::Graph ReadAMS_Tn         (std::string filename);
        CRACS::Graph ReadAMS_pbar       (std::string filename);
        CRACS::Graph ReadAMS_pbaroverp  (std::string filename);
        CRACS::Graph ReadAMS_boverc     (std::string filename);
        CRACS::Graph ReadAMS_lelptons   (std::string filename);
        
        CRACS::Graph ReadPamela         (std::string filename);
        CRACS::Graph ReadPamela_pHe     (std::string filename);
        
        CRACS::Graph ReadCream          (std::string filename);
        CRACS::Graph ReadCream_PandHe   (std::string filename);
        
        CRACS::Graph ReadHEAO           (std::string filename);
        
        CRACS::Graph ReadVoyager        (std::string filename);
        
        CRACS::Graph ReadAtic           (std::string filename, int A=1);
        
        
        std::vector<Graph>              fFluxes;
        std::vector<int>                fZ;
        std::vector<int>                fA;
        std::vector<std::string>        fSpecies;
        std::vector<std::string>        fExperiment;
        std::vector<spectrumType>       fType;
        
        
    };
    

}



//end READFLUXES_H

#endif
