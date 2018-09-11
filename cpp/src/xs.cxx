#include "math.h"
#include "iostream"

#include "gsl_integration.h"

#include "xs.h"
#include "xs_definitions.h"

namespace CRXS {
    
    bool XS::convert_LAB_to_CM( const double T_p_LAB, const double T_product_LAB, const double eta_LAB, double &s, double &E_product, double &pT_product, double &x_F, int product ){
        
        double m_product;
        
        if (product==P_BAR) {
            m_product = XS_definitions::fMass_proton;
        }else if (product==D_BAR) {
            m_product = XS_definitions::fMass_deuteron;
        }else if (product==HE_BAR) {
            m_product = XS_definitions::fMass_helium3;
        }else{
            return false;
        }
        
        s                    = 4*XS_definitions::fMass_proton*XS_definitions::fMass_proton + 2 * T_p_LAB * XS_definitions::fMass_proton;
        
        double p_product_LAB = sqrt(  T_product_LAB*(T_product_LAB+2*m_product)  );
        double E_p_LAB       = T_p_LAB+XS_definitions::fMass_proton;
        double E_product_LAB = T_product_LAB+m_product;
        
        double beta          = sqrt(E_p_LAB - XS_definitions::fMass_proton)/sqrt(E_p_LAB + XS_definitions::fMass_proton);
        double gamma         = 1./sqrt(1 - beta*beta);
        double gammabeta     = gamma * beta;
        
        E_product            =      gamma * E_product_LAB - gammabeta * p_product_LAB*tanh(eta_LAB);
        double pL_product    = -gammabeta * E_product_LAB + gamma     * p_product_LAB*tanh(eta_LAB);
        
        pT_product              = p_product_LAB/cosh(eta_LAB);
        //x_R                  = E_product/E_product_Max;
        //double E_product_Max    = ( s-8.*XS_definitions::fMass_proton*XS_definitions::fMass_proton )/2./sqrt( s );
        
        x_F = 2. * pL_product / sqrt(s);
        
        return true;
    }
    
   
}
