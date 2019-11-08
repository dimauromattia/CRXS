#ifndef CRXS__CRXS_H
#define CRXS__CRXS_H

#include "string"

namespace CRXS {
    
    enum IntegrationMethod{
        GSL     =  1,
        TRAPEZE =  2,
    };
    
    class CRXS_config{
    public:
        static std::string Get_CRXS_DataDir();
        static int  IntegrationMethod;
        static void SetupIntegrationMethod( int method ){
            IntegrationMethod=method;
        };
        
    };
    
    
}

#endif




