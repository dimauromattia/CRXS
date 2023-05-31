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
            m_product = XS_definitions::fMass_helion3;
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
    
    bool   XS::fIsRestricted_pp = false;
    
    int    XS::fRestrictedParameterSpace_LAB = 0;
    double XS::fRestrictedParameterSpace_LAB__Tp   []={0};
    double XS::fRestrictedParameterSpace_LAB__Tpbar[]={0};
    double XS::fRestrictedParameterSpace_LAB__eta   []={0};
    
    void CRXS::XS::SetRestrictedParameterSpace_LAB( double Tp, double Tpbar, double eta ){
        if(fRestrictedParameterSpace_LAB==0){
            RemoveRestrictedParameterSpace_LAB();
        }
        fRestrictedParameterSpace_LAB__Tp   [fRestrictedParameterSpace_LAB] = Tp;
        fRestrictedParameterSpace_LAB__Tpbar[fRestrictedParameterSpace_LAB] = Tpbar;
        fRestrictedParameterSpace_LAB__eta  [fRestrictedParameterSpace_LAB] = eta;
        fRestrictedParameterSpace_LAB ++;
        
        if (Tp<fRestrictedParameterSpace_LAB__Tp[101]) {
            fRestrictedParameterSpace_LAB__Tp[101] = Tp;
        }
        if (Tp>fRestrictedParameterSpace_LAB__Tp[102]){
            fRestrictedParameterSpace_LAB__Tp[102] = Tp;
        }
        if (Tpbar<fRestrictedParameterSpace_LAB__Tpbar[101]) {
            fRestrictedParameterSpace_LAB__Tpbar[101] = Tpbar;
        }
        if (Tpbar>fRestrictedParameterSpace_LAB__Tpbar[102]){
            fRestrictedParameterSpace_LAB__Tpbar[102] = Tpbar;
        }
        if (eta<fRestrictedParameterSpace_LAB__eta[101]) {
            fRestrictedParameterSpace_LAB__eta[101] = eta;
        }
        if (eta>fRestrictedParameterSpace_LAB__eta[102]){
            fRestrictedParameterSpace_LAB__eta[102] = eta;
        }
    };

    void CRXS::XS::RemoveRestrictedParameterSpace_LAB(  ){
        fRestrictedParameterSpace_LAB = 0;
        fIsRestricted_pp = false;
        
        fRestrictedParameterSpace_LAB__Tp    [101] =   1e90;
        fRestrictedParameterSpace_LAB__Tp    [102] =  -1e90;
        fRestrictedParameterSpace_LAB__Tpbar [101] =   1e90;
        fRestrictedParameterSpace_LAB__Tpbar [102] =  -1e90;
        fRestrictedParameterSpace_LAB__eta   [101] =   1e90;
        fRestrictedParameterSpace_LAB__eta   [102] =  -1e90;
        
    };

    bool CRXS::XS::isInRestricted_LAB(double Tp, double Tpbar, double eta){
        if (Tp   <fRestrictedParameterSpace_LAB__Tp    [101]) return  false;
        if (Tp   >fRestrictedParameterSpace_LAB__Tp    [102]) return  false;
        if (Tpbar<fRestrictedParameterSpace_LAB__Tpbar [101]) return  false;
        if (Tpbar>fRestrictedParameterSpace_LAB__Tpbar [102]) return  false;
        if (eta  <fRestrictedParameterSpace_LAB__eta   [101]) return  false;
        if (eta  >fRestrictedParameterSpace_LAB__eta   [102]) return  false;
        double p [3] = { Tp, Tpbar, eta };
        int n = fRestrictedParameterSpace_LAB;
        for (int ia = 0; ia<n; ia++) {
            for (int ib = ia+1; ib<n; ib++) {
                for (int ic = ib+1; ic<n; ic++) {
                    for (int id = ic+1; id<n; id++) {
                        double a [3] =  {fRestrictedParameterSpace_LAB__Tp[ia], fRestrictedParameterSpace_LAB__Tpbar[ia], fRestrictedParameterSpace_LAB__eta [ia] };
                        double b [3] =  {fRestrictedParameterSpace_LAB__Tp[ib], fRestrictedParameterSpace_LAB__Tpbar[ib], fRestrictedParameterSpace_LAB__eta [ib] };
                        double c [3] =  {fRestrictedParameterSpace_LAB__Tp[ic], fRestrictedParameterSpace_LAB__Tpbar[ic], fRestrictedParameterSpace_LAB__eta [ic] };
                        double d [3] =  {fRestrictedParameterSpace_LAB__Tp[id], fRestrictedParameterSpace_LAB__Tpbar[id], fRestrictedParameterSpace_LAB__eta [id] };
                        if ( LA::inside(a, b, c, d, p)){
                            return true;
                        };
                    }
                }
            }
        }
        return false;
    }
    
    
    
    int    XS::fRestrictedParameterSpace_CM = 0;
    double XS::fRestrictedParameterSpace_CM__s    []={0};
    double XS::fRestrictedParameterSpace_CM__xf   []={0};
    double XS::fRestrictedParameterSpace_CM__pT   []={0};
    
    
    void CRXS::XS::SetRestrictedParameterSpace_CM( double s, double xf, double pT ){
        if(fRestrictedParameterSpace_CM==0){
            RemoveRestrictedParameterSpace_CM();
        }
        fRestrictedParameterSpace_CM__s    [fRestrictedParameterSpace_CM] = s;
        fRestrictedParameterSpace_CM__xf   [fRestrictedParameterSpace_CM] = xf;
        fRestrictedParameterSpace_CM__pT   [fRestrictedParameterSpace_CM] = pT;
        fRestrictedParameterSpace_CM ++;
        if (s<fRestrictedParameterSpace_CM__s[101]) {
            fRestrictedParameterSpace_CM__s[101] = s;
        }
        if (s>fRestrictedParameterSpace_CM__s[102]){
            fRestrictedParameterSpace_CM__s[102] = s;
        }
        if (xf<fRestrictedParameterSpace_CM__xf[101]) {
            fRestrictedParameterSpace_CM__xf[101] = xf;
        }
        if (xf>fRestrictedParameterSpace_CM__xf[102]){
            fRestrictedParameterSpace_CM__xf[102] = xf;
        }
        if (pT<fRestrictedParameterSpace_CM__pT[101]) {
            fRestrictedParameterSpace_CM__pT[101] = pT;
        }
        if (pT>fRestrictedParameterSpace_CM__pT[102]){
            fRestrictedParameterSpace_CM__pT[102] = pT;
        }
        //std::cout << fRestrictedParameterSpace_CM__s[101] << std::endl;
        //std::cout << fRestrictedParameterSpace_CM__s[102] << std::endl;
        //std::cout << fRestrictedParameterSpace_CM__xf[101] << std::endl;
        //std::cout << fRestrictedParameterSpace_CM__xf[102] << std::endl;
        //std::cout << fRestrictedParameterSpace_CM__pT[101] << std::endl;
        //std::cout << fRestrictedParameterSpace_CM__pT[102] << std::endl;
        //std::cout << "----" << std::endl;
    };
    
    void CRXS::XS::RemoveRestrictedParameterSpace_CM(  ){
        fRestrictedParameterSpace_CM = 0;
        fIsRestricted_pp = false;
        
        fRestrictedParameterSpace_CM__s [101] =   1e90;
        fRestrictedParameterSpace_CM__s [102] =  -1e90;
        fRestrictedParameterSpace_CM__xf[101] =   1e90;
        fRestrictedParameterSpace_CM__xf[102] =  -1e90;
        fRestrictedParameterSpace_CM__pT[101] =   1e90;
        fRestrictedParameterSpace_CM__pT[102] =  -1e90;
        
        
    };
    
    bool CRXS::XS::isInRestricted_CM(double s, double xf, double pT){
        if ( s<fRestrictedParameterSpace_CM__s [101]) return  false;
        if ( s>fRestrictedParameterSpace_CM__s [102]) return  false;
        if (xf<fRestrictedParameterSpace_CM__xf[101]) return  false;
        if (xf>fRestrictedParameterSpace_CM__xf[102]) return  false;
        if (pT<fRestrictedParameterSpace_CM__pT[101]) return  false;
        if (pT>fRestrictedParameterSpace_CM__pT[102]) return  false;
        double p [3] = { s, xf, pT };
        int n = fRestrictedParameterSpace_CM;
        for (int ia = 0; ia<n; ia++) {
            for (int ib = ia+1; ib<n; ib++) {
                for (int ic = ib+1; ic<n; ic++) {
                    for (int id = ic+1; id<n; id++) {
                        double a [3] =  {fRestrictedParameterSpace_CM__s[ia], fRestrictedParameterSpace_CM__xf   [ia], fRestrictedParameterSpace_CM__pT[ia] };
                        double b [3] =  {fRestrictedParameterSpace_CM__s[ib], fRestrictedParameterSpace_CM__xf   [ib], fRestrictedParameterSpace_CM__pT[ib] };
                        double c [3] =  {fRestrictedParameterSpace_CM__s[ic], fRestrictedParameterSpace_CM__xf   [ic], fRestrictedParameterSpace_CM__pT[ic] };
                        double d [3] =  {fRestrictedParameterSpace_CM__s[id], fRestrictedParameterSpace_CM__xf   [id], fRestrictedParameterSpace_CM__pT[id] };
                        if ( LA::inside(a, b, c, d, p)){
                            return true;
                        };
                    }
                }
            }
        }
        return false;
    }
    
    
    
    
    
}
