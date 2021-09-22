#ifndef CROSSSECTIONS_H
#define CROSSSECTIONS_H

#include <iostream>
#include <math.h>
#include <definitions.h>
#include <ppCrossSectionParametrizations.h>
#include <crossSectionTransformations.h>


namespace CRACS {
    
    static double fP_coal = 0.062;
    
    
    /*! \brief Energy differential cross sections for collisions.
     *
     *
     */
    
    class CS{
        public:
        static const int PP_PBAR;
        static const int PHe_PBAR;
        static const int HeP_PBAR;
        static const int HeHe_PBAR;
        
        static const int PP_DBAR;
        static const int PHe_DBAR;
        static const int HeP_DBAR;
        static const int HeHe_DBAR;
        
        static const int PPbar_DBAR;
        static const int HePbar_DBAR;
        
        
        static const int PP_HEBAR;
        static const int PHe_HEBAR;
        static const int HeP_HEBAR;
        static const int HeHe_HEBAR;
        
        static const int PPbar_HEBAR;
        static const int HePbar_HEBAR;
        
        
        static const int DEFAULT;
        static const int TAN_NG;
        static const int DI_MAURO;
        static const int DI_MAURO12;
        static const int WINKLER;
        static const int ANDERSON;
        static const int DUPERRAY;
        static const int DTUNUC;
        static const int KACHELRIESS;
        static const int KORSMEIER;
        static const int WINKLER_withHYPwithNBAR;
        
        
        //!  Energy differential antiproton production cross section in Lab Frame from pp collison
        /*!
         *  Lab frame:     One proton is at rest and the other one has incident energy E_p
         *
         *  \f[  \frac{ d \sigma_{pp}^{(\bar{p}), LAB} }{d E_{\bar{p}} }  (E_p, E_{\bar{p}})  \f]
         *
         *  \param double E_p               Incident proton energy in LAB frame.
         *  \param doulbe E_pbar            Energy of the produced particle in LAB frame.
         *  \param int    precision         Number of steps for the integration over cos(theta).
         *  \param int    parametrization   Parametrization: TAN_NG, DI_MAURO, ...   Default: DI_MAURO
         *  \param int    output            Output level: NO_OUT, WARN_OUT, ALL_OUT. Default: WARN_OUT
         *
         * */
        static double dE_pp_pbar_LAB(double E_p_LAB, double E_pbar_LAB,
                                     int precision=10000,
                                     int parametrization=DI_MAURO,
                                     int output=WARN_OUT)
        {
            if(parametrization==DI_MAURO) {
                return CSTransformations::dE_pp_product_LAB_intEta (E_p_LAB, E_pbar_LAB, fMass_proton,
                                                                    ppCSParametrizations::inv_pp_pbar_CM__diMauro,
                                                                    precision,
                                                                    output);
            }else if(parametrization==DI_MAURO12) {
                return CSTransformations::dE_pp_product_LAB_intEta(E_p_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_pp_pbar_CM__diMauro12,
                                                                   precision,
                                                                   output);
            }else if(parametrization==TAN_NG) {
                return CSTransformations::dE_pp_product_LAB_intEta(E_p_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_pp_pbar_CM__Tan_Ng,
                                                                   precision,
                                                                   output);
            }else if(parametrization==DUPERRAY) {
                return CSTransformations::dE_pp_product_LAB_intEta(E_p_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_pp_pbar_CM__Duperray,
                                                                   precision,
                                                                   output);
            }else if(parametrization==KORSMEIER) {
                return CSTransformations::dE_pp_product_LAB_intEta(E_p_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_pp_pbar_CM__Korsmeier,
                                                                   precision,
                                                                   output);
            }else if(parametrization==WINKLER) {
                return CSTransformations::dE_pp_product_LAB_intEta(E_p_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_pp_pbar_CM__Winkler,
                                                                   precision,
                                                                   output);
            }else if(parametrization==WINKLER_withHYPwithNBAR) {
                return CSTransformations::dE_pp_product_LAB_intEta(E_p_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_pp_pbar_CM__WinklerWithHypWithNbar,
                                                                   precision,
                                                                   output);
            }else{
                std::cout << warning_parametrization_not_available << std::endl;
            }
            
            return -1;
        }
        
        
        //!  Energy differential antiproton production cross section in Lab Frame from pHe collison
        /*!
         *  Lab frame:     He is at rest and proton has incident energy E_p
         *
         *  \f[  \frac{ d \sigma_{pHe}^{(\bar{p}), LAB} }{d E_{\bar{p}} }  (E_p, E_{\bar{p}})  \f]
         *
         *  \param double E_p               Incident proton energy in LAB frame.
         *  \param doulbe E_pbar            Energy of the produced particle in LAB frame.
         *  \param int    precision         Number of steps for the integration over cos(theta).
         *  \param int    parametrization   Parametrization: DUPPERAY, ...   Default: DUPPERAY
         *  \param int    output            Output level: NO_OUT, WARN_OUT, ALL_OUT. Default: WARN_OUT
         *
         * */
        static double dE_pHe_pbar_LAB(double E_p_LAB, double E_pbar_LAB,
                                      int precision=10000,
                                      int parametrization=DI_MAURO,
                                      int output=WARN_OUT)
        {
            if(parametrization==DUPERRAY) {
                return CSTransformations::dE_pp_product_LAB_intEta(E_p_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_pHe_pbar_CM__Duperray,
                                                                   precision,
                                                                   output);
            }else if(parametrization==DI_MAURO12) {
                return CSTransformations::dE_pp_product_LAB_intEta(E_p_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_pHe_pbar_CM__diMauro12,
                                                                   precision,
                                                                   output);
            }else if(parametrization==DI_MAURO) {
                return CSTransformations::dE_pp_product_LAB_intEta(E_p_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_pHe_pbar_CM__diMauro,
                                                                   precision,
                                                                   output);
            }else if(parametrization==TAN_NG) {
                return CSTransformations::dE_pp_product_LAB_intEta(E_p_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_pHe_pbar_CM__Tan_Ng,
                                                                   precision,
                                                                   output);
            }else if(parametrization==KORSMEIER) {
                return CSTransformations::dE_pp_product_LAB_intEta(E_p_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_pHe_pbar_CM__Korsmeier,
                                                                   precision,
                                                                   output);
            }else if(parametrization==WINKLER) {
                return CSTransformations::dE_pp_product_LAB_intEta(E_p_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_pHe_pbar_CM__Winkler,
                                                                   precision,
                                                                   output);
            }else if(parametrization==WINKLER_withHYPwithNBAR) {
                return CSTransformations::dE_pp_product_LAB_intEta(E_p_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_pHe_pbar_CM__WinklerWithHypWithNbar,
                                                                   precision,
                                                                   output);
            }else{
                std::cout << warning_parametrization_not_available << std::endl;
            }
            
            return -1;
        }
        
        //!  Energy differential antiproton production cross section in Lab Frame from pHe collison
        /*!
         *  Lab frame:     p is at rest and helium has incident energy E_He*4
         *
         *  \f[  \frac{ d \sigma_{pHe}^{(\bar{p}), LAB} }{d E_{\bar{p}} }  (E_p, E_{\bar{p}})  \f]
         *
         *  \param double En_He              Incident helium energy PER NUCLEON in LAB frame.
         *  \param doulbe E_pbar            Energy of the produced particle in LAB frame.
         *  \param int    precision         Number of steps for the integration over cos(theta).
         *  \param int    parametrization   Parametrization: DUPPERAY, ...   Default: DUPPERAY
         *  \param int    output            Output level: NO_OUT, WARN_OUT, ALL_OUT. Default: WARN_OUT
         *
         * */
        static double dE_Hep_pbar_LAB(double En_He_LAB, double E_pbar_LAB,
                                      int precision=10000,
                                      int parametrization=DI_MAURO,
                                      int output=WARN_OUT)
        {
            if(parametrization==DUPERRAY) {
                return CSTransformations::dE_Hep_product_LAB_intEta(En_He_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_pHe_pbar_CM__Duperray,
                                                                   precision,
                                                                   output);
            }else if(parametrization==DI_MAURO12) {
                return CSTransformations::dE_Hep_product_LAB_intEta(En_He_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_pHe_pbar_CM__diMauro12,
                                                                   precision,
                                                                   output);
            }else if(parametrization==DI_MAURO) {
                return CSTransformations::dE_Hep_product_LAB_intEta(En_He_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_pHe_pbar_CM__diMauro,
                                                                   precision,
                                                                   output);
            }else if(parametrization==TAN_NG) {
                return CSTransformations::dE_Hep_product_LAB_intEta(En_He_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_pHe_pbar_CM__Tan_Ng,
                                                                   precision,
                                                                   output);
            }else if(parametrization==KORSMEIER) {
                return CSTransformations::dE_Hep_product_LAB_intEta(En_He_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_pHe_pbar_CM__Korsmeier,
                                                                   precision,
                                                                   output);
            }else if(parametrization==WINKLER) {
                return CSTransformations::dE_pp_product_LAB_intEta(En_He_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_Hep_pbar_CM__Winkler,
                                                                   precision,
                                                                   output);
            }else if(parametrization==WINKLER_withHYPwithNBAR) {
                return CSTransformations::dE_pp_product_LAB_intEta(En_He_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_Hep_pbar_CM__WinklerWithHypWithNbar,
                                                                   precision,
                                                                   output);
            }else{
                std::cout << warning_parametrization_not_available << std::endl;
            }
            
            return -1;
        }
        
        
        //!  Energy differential antiproton production cross section in Lab Frame from pHe collison
        /*!
         *  Lab frame:     p is at rest and helium has incident energy E_He*4
         *
         *  \f[  \frac{ d \sigma_{pHe}^{(\bar{p}), LAB} }{d E_{\bar{p}} }  (E_p, E_{\bar{p}})  \f]
         *
         *  \param double En_He              Incident helium energy PER NUCLEON in LAB frame.
         *  \param doulbe E_pbar            Energy of the produced particle in LAB frame.
         *  \param int    precision         Number of steps for the integration over cos(theta).
         *  \param int    parametrization   Parametrization: DUPPERAY, ...   Default: DUPPERAY
         *  \param int    output            Output level: NO_OUT, WARN_OUT, ALL_OUT. Default: WARN_OUT
         *
         * */
        static double dE_HeHe_pbar_LAB(double En_He_LAB, double E_pbar_LAB,
                                      int precision=10000,
                                      int parametrization=DI_MAURO,
                                      int output=WARN_OUT)
        {
            if(parametrization==DI_MAURO12) {
                return CSTransformations::dE_pp_product_LAB_intEta(En_He_LAB, E_pbar_LAB, fMass_proton,
                                                                    ppCSParametrizations::inv_HeHe_pbar_CM__diMauro12,
                                                                    precision,
                                                                    output);
            }else if(parametrization==DI_MAURO) {
                return CSTransformations::dE_pp_product_LAB_intEta(En_He_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_HeHe_pbar_CM__diMauro,
                                                                   precision,
                                                                   output);
            }else if(parametrization==DUPERRAY) {
                return CSTransformations::dE_pp_product_LAB_intEta(En_He_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_HeHe_pbar_CM__Duperray,
                                                                   precision,
                                                                   output);
            }else if(parametrization==WINKLER) {
                return CSTransformations::dE_pp_product_LAB_intEta(En_He_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_HeHe_pbar_CM__Winkler,
                                                                   precision,
                                                                   output);
            }else if(parametrization==WINKLER_withHYPwithNBAR) {
                return CSTransformations::dE_pp_product_LAB_intEta(En_He_LAB, E_pbar_LAB, fMass_proton,
                                                                   ppCSParametrizations::inv_HeHe_pbar_CM__WinklerWithHypWithNbar,
                                                                   precision,
                                                                   output);
            }else{
                std::cout << warning_parametrization_not_available << std::endl;
            }
            
            return -1;
        }
        
        
        
        
        
        //FIXME: Change all integrations (below) to eta!
        
        
        //!  Energy differential antideuteron production cross section in Lab Frame from pp collison
        /*!
         *  Lab frame:     One proton is at rest and the other one has incident energy E_p
         *
         *  \f[  \frac{ d \sigma_{pp}^{(\bar{D}), LAB} }{d E_{\bar{D}} }  (E_p, E_{\bar{D}})  \f]
         *
         *  \param double E_p               Incident proton energy in LAB frame.
         *  \param doulbe E_Dbar            Energy of the produced particle in LAB frame.
         *  \param int    precision         Number of steps for the integration over cos(theta).
         *  \param int    parametrization   Parametrization: TAN_NG, DI_MAURO, ...   Default: DI_MAURO
         *  \param int    output            Output level: NO_OUT, WARN_OUT, ALL_OUT. Default: WARN_OUT
         *
         * */
        static double dE_pp_Dbar_LAB(double E_p_LAB, double E_dbar_LAB,
                                     double p_coal = 0.062,
                                     int precision=10000,
                                     int parametrization=DI_MAURO12,
                                     int output=WARN_OUT)
        {
            fP_coal = p_coal;
            if(parametrization==DI_MAURO12) {
                return CSTransformations::dE_pp_product_LAB_intEta(E_p_LAB, E_dbar_LAB, fMass_deuteron,
                                                            inv_pp_Dbar_CM_diMauro,
                                                            precision,
                                                            output);
            }else{
                std::cout << warning_parametrization_not_available << std::endl;
            }
            
            return -1;
        }
        
        static double dE_pp_Hebar_LAB(double E_p_LAB, double E_Hebar_LAB,
                                     double p_coal = 0.062,
                                     int precision=10000,
                                     int parametrization=DI_MAURO12,
                                     int output=WARN_OUT)
        {
            fP_coal = p_coal;
            if(parametrization==DI_MAURO12) {
                return CSTransformations::dE_pp_product_LAB_intEta(E_p_LAB, E_Hebar_LAB, fMass_helium3,
                                                                   inv_pp_Hebar_CM_diMauro,
                                                                   precision,
                                                                   output);
            }else{
                std::cout << warning_parametrization_not_available << std::endl;
            }
            
            return -1;
        }
        
        
        //!  Energy differential antideuteron production cross section in Lab Frame from ppbar collison
        /*!
         *  Lab frame:     One proton is at rest and the other one has incident energy E_pbar
         *
         *  \f[  \frac{ d \sigma_{p\bar{p}}^{(\bar{D}), LAB} }{d E_{\bar{D}} }  (E_{\bar{p}}, E_{\bar{D}})  \f]
         *
         *  We use the following assumptions:
         *      - \f$   E d^3\sigma/dp^3 \, (\bar{p} + p \rightarrow \bar{p} + X) =
         *              E d^3\sigma/dp^3 \, (p       + p \rightarrow p       + X)       \f$
         *      - \f$   E d^3\sigma/dp^3 \, (\bar{p} + p \rightarrow \bar{n} + X) =
         *              E d^3\sigma/dp^3 \, (p       + p \rightarrow \bar{n} + X)       \f$
         *
         *  In the default scenario we use di Mauro parametrizations.
         *
         *  \param double E_pbar_LAB        Incident proton energy in LAB frame.
         *  \param doulbe E_dbar_LAB        Energy of the produced particle in LAB frame.
         *  \param int    precision         Number of steps for the integration over cos(theta).
         *  \param int    parametrization   Parametrization: TAN_NG, DI_MAURO, ...   Default: DI_MAURO
         *  \param int    output            Output level: NO_OUT, WARN_OUT, ALL_OUT. Default: WARN_OUT
         *
         * */
        static double dE_ppbar_Dbar_LAB(double E_pbar_LAB, double E_dbar_LAB,
                                        double p_coal = 0.062,
                                        int precision=10000,
                                        int parametrization=DEFAULT,
                                        int output=WARN_OUT)
        {
            fP_coal = p_coal;
            if(parametrization==DEFAULT) {
                return CSTransformations::dE_pp_product_LAB_intEta(E_pbar_LAB, E_dbar_LAB, fMass_deuteron,
                                                            inv_ppbar_Dbar_CM_default,
                                                            precision,
                                                            output);
            }else{
                std::cout << warning_parametrization_not_available << std::endl;
            }
            
            return -1;
        }
        
        
        static double dE_ppbar_Hebar_LAB(double E_pbar_LAB, double E_Hebar_LAB,
                                        double p_coal = 0.062,
                                        int precision=10000,
                                        int parametrization=DEFAULT,
                                        int output=WARN_OUT)
        {
            fP_coal = p_coal;
            if(parametrization==DEFAULT) {
                return CSTransformations::dE_pp_product_LAB_intEta(E_pbar_LAB, E_Hebar_LAB, fMass_helium3,
                                                                   inv_ppbar_Hebar_CM_default,
                                                                   precision,
                                                                   output);
            }else{
                std::cout << warning_parametrization_not_available << std::endl;
            }
            
            return -1;
        }
        
        
        
        
        //!  Energy differential p production cross section in Lab Frame from pp collison
        /*!
         *  Lab frame:     One proton is at rest and the other one has incident energy E_p
         *
         *  \f[  \frac{ d \sigma_{pp}^{(p), LAB} }{d E_{p} }  (E_p, E_{p}(product))  \f]
         *
         *  \param double E_p               Incident proton energy in LAB frame.
         *  \param doulbe E_pOut            Energy of the produced particle in LAB frame.
         *  \param int    precision         Number of steps for the integration over cos(theta).
         *  \param int    parametrization   Parametrization: TAN_NG, DI_MAURO, ...   Default: ANDERSON
         *  \param int    output            Output level: NO_OUT, WARN_OUT, ALL_OUT. Default: WARN_OUT
         *
         * */
        static double dE_pp_p_LAB(double E_p_LAB, double E_pOut_LAB,
                                  int precision=100000,
                                  int parametrization=DI_MAURO,
                                  int output=WARN_OUT)
        {
            if(parametrization==ANDERSON) {
                return CSTransformations::dE_pp_product_LAB_intEta(E_p_LAB, E_pOut_LAB, fMass_proton,
                                                            ppCSParametrizations::inv_pp_p_CM__Anderson,
                                                            precision,
                                                            output);
            }else{
                std::cout << warning_parametrization_not_available << std::endl;
            }
            
            return -1;
        }
        
        
        //        //!  Momentum differential Dbar production cross section in Lab Frame from pp collison
        //        /*!
        //         *  Lab frame:     One proton is at rest and the other one has incident energy E_p
        //         *
        //         *  \f[  \frac{ d^3 \sigma_{pp}^{(\bar{D}), LAB} }{d p_{\bar{D}}^3 }  (E_p, E_{\bar{D}}, cos(\theta) )  \f]
        //         *
        //         *  \param double E_p               Incident proton energy in LAB frame.
        //         *  \param doulbe E_pOut            Energy of the produced particle in LAB frame.
        //         *  \param doulbe cos(\theta)       cos of angle between p and Dbar.
        //         *  \param int    precision         Number of steps for the integration over cos(theta).
        //         *  \param int    parametrization   Parametrization: TAN_NG, DI_MAURO, ...   Default: ANDERSON
        //         *  \param int    output            Output level: NO_OUT, WARN_OUT, ALL_OUT. Default: WARN_OUT
        //         *
        //         * */
        //        static double d3p_pp_Dbar_LAB(double E_p_LAB, double E_Dbar_LAB, double cos_theta_product_LAB,
        //                                      double p_coal = 0.062,
        //                                      int parametrization=DI_MAURO12)
        //        {
        //            fP_coal = p_coal;
        //            if(parametrization==DI_MAURO) {
        //                return CSTransformations::d3p_pp_product_LAB(E_p_LAB, E_Dbar_LAB, cos_theta_product_LAB,
        //                                                             fMass_deuteron,
        //                                                             inv_pp_Dbar_CM_diMauro);
        //            }else{
        //                std::cout << warning_parametrization_not_available << std::endl;
        //            }
        //
        //            return -1;
        //        }
        
        
        
        private:
        // Helper functions in order to set default values
        static double inv_pp_Dbar_CM_diMauro(double s, double E_dbar, double pT_dbar){
            return ppCSParametrizations::inv_pp_Dbar_CM( s, E_dbar, pT_dbar, fP_coal, ppCSParametrizations::tot_pp__diMauro, ppCSParametrizations::inv_pp_pbar_CM__diMauro12, ppCSParametrizations::inv_pp_nbar_CM__diMauro12);
        }
        static double inv_ppbar_Dbar_CM_default(double s, double E_dbar, double pT_dbar){
            return ppCSParametrizations::inv_ppbar_Dbar_CM( s, E_dbar, pT_dbar, fP_coal, ppCSParametrizations::tot_pp__diMauro, ppCSParametrizations::inv_pp_p_CM__Anderson, ppCSParametrizations::inv_pp_nbar_CM__diMauro12);
        }
        
        
        static double inv_pp_Hebar_CM_diMauro(double s, double E_Hebar, double pT_Hebar){
            return ppCSParametrizations::inv_pp_Hebar_CM( s, E_Hebar, pT_Hebar, fP_coal, ppCSParametrizations::tot_pp__diMauro, ppCSParametrizations::inv_pp_pbar_CM__diMauro12, ppCSParametrizations::inv_pp_nbar_CM__diMauro12);
        }
        static double inv_ppbar_Hebar_CM_default(double s, double E_Hebar, double pT_Hebar){
            return ppCSParametrizations::inv_ppbar_Hebar_CM( s, E_Hebar, pT_Hebar, fP_coal, ppCSParametrizations::tot_pp__diMauro, ppCSParametrizations::inv_pp_p_CM__Anderson, ppCSParametrizations::inv_pp_nbar_CM__diMauro12);
        }
        
        
        
    };
    
    
}



//end CROSSSECTIONS_H

#endif
