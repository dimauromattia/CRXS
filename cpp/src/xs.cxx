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

        pT_product           = p_product_LAB/cosh(eta_LAB);
        //x_R                  = E_product/E_product_Max;
        //double E_product_Max    = ( s-8.*XS_definitions::fMass_proton*XS_definitions::fMass_proton )/2./sqrt( s );

        x_F = 2. * pL_product / sqrt(s);

        return true;
    }
    
    
    void CRXS::XS::Set_SELF_C_parameters_diMauro(double *C){
        CRXS::XS_definitions::diMauro_SELF_C1_to_C11[ 1] = C[ 1];
        CRXS::XS_definitions::diMauro_SELF_C1_to_C11[ 2] = C[ 2];
        CRXS::XS_definitions::diMauro_SELF_C1_to_C11[ 3] = C[ 3];
        CRXS::XS_definitions::diMauro_SELF_C1_to_C11[ 4] = C[ 4];
        CRXS::XS_definitions::diMauro_SELF_C1_to_C11[ 5] = C[ 5];
        CRXS::XS_definitions::diMauro_SELF_C1_to_C11[ 6] = C[ 6];
        CRXS::XS_definitions::diMauro_SELF_C1_to_C11[ 7] = C[ 7];
        CRXS::XS_definitions::diMauro_SELF_C1_to_C11[ 8] = C[ 8];
        CRXS::XS_definitions::diMauro_SELF_C1_to_C11[ 9] = C[ 9];
        CRXS::XS_definitions::diMauro_SELF_C1_to_C11[10] = C[10];
        CRXS::XS_definitions::diMauro_SELF_C1_to_C11[11] = C[11];
        
        CRXS::XS_definitions::diMauro_SELF_C1_to_C16[ 1] = C[12];
        CRXS::XS_definitions::diMauro_SELF_C1_to_C16[ 2] = C[13];
        CRXS::XS_definitions::diMauro_SELF_C1_to_C16[ 3] = C[14];
        CRXS::XS_definitions::diMauro_SELF_C1_to_C16[ 4] = C[15];
        CRXS::XS_definitions::diMauro_SELF_C1_to_C16[14] = C[16];
        CRXS::XS_definitions::diMauro_SELF_C1_to_C16[15] = C[17];
        CRXS::XS_definitions::diMauro_SELF_C1_to_C16[16] = C[18];
    }
    
    
    void CRXS::XS::Set_SELF_C_parameters_Winkler(double *C){
        CRXS::XS_definitions::Winkler_SELF_C1_to_C16[ 1] = C[ 1];
        CRXS::XS_definitions::Winkler_SELF_C1_to_C16[ 2] = C[ 2];
        CRXS::XS_definitions::Winkler_SELF_C1_to_C16[ 3] = C[ 3];
        CRXS::XS_definitions::Winkler_SELF_C1_to_C16[ 4] = C[ 4];
        CRXS::XS_definitions::Winkler_SELF_C1_to_C16[ 5] = C[ 5];
        CRXS::XS_definitions::Winkler_SELF_C1_to_C16[ 6] = C[ 6];
        CRXS::XS_definitions::Winkler_SELF_C1_to_C16[ 7] = C[ 7];
        CRXS::XS_definitions::Winkler_SELF_C1_to_C16[ 8] = C[ 8];
        CRXS::XS_definitions::Winkler_SELF_C1_to_C16[ 9] = C[ 9];
        CRXS::XS_definitions::Winkler_SELF_C1_to_C16[10] = C[10];
        CRXS::XS_definitions::Winkler_SELF_C1_to_C16[11] = C[11];
        CRXS::XS_definitions::Winkler_SELF_C1_to_C16[12] = C[12];
        CRXS::XS_definitions::Winkler_SELF_C1_to_C16[13] = C[13];
        CRXS::XS_definitions::Winkler_SELF_C1_to_C16[14] = C[14];
        CRXS::XS_definitions::Winkler_SELF_C1_to_C16[15] = C[15];
        CRXS::XS_definitions::Winkler_SELF_C1_to_C16[16] = C[16];
    }
    
    void CRXS::XS::Set_SELF_D_parameters_diMauro(double *D){
        CRXS::XS_definitions::diMauro_SELF_D1_to_D2[ 1] = D[ 1];
        CRXS::XS_definitions::diMauro_SELF_D1_to_D2[ 2] = D[ 2];
    }
    
    void CRXS::XS::Set_SELF_D_parameters_Winkler(double *D){
        CRXS::XS_definitions::Winkler_SELF_D1_to_D2[ 1] = D[ 1];
        CRXS::XS_definitions::Winkler_SELF_D1_to_D2[ 2] = D[ 2];
    }
    
   
}
