#ifndef CRXS__XS_DEFINITIONS_H
#define CRXS__XS_DEFINITIONS_H

#include "xs.h"
#include "stdio.h"

#include "string"

namespace CRXS {
    
    class XS_definitions{
        
    public:
        static double fMass_proton;
        static double fMass_neutron;
        static double fMass_deuteron;
        static double fMass_helium3;
        
        //! Parametrization of the invariant antiproton production (pbar) crosssection from pp in CMF
        /*!
         *  All cross sections are given in mbarn
         *  All enegies, momenta, and masses have unit GeV.
         *
         *  Taken from:     Winkler, M. W.; 2017;
         *                  Cosmic Ray Antiprotons at High Energies;
         *                  arXiv:1701.04866
         *
         *  \f[  E_{\bar{p}} \frac{ d^3 \sigma_{pp}^{(\bar{p})} }{d^3 p_{\bar{p}} }  (s, E_{\bar{p}}, p_{T,\bar{p}})  \f]
         *
         *  \param double  s         CM energy.
         *  \param doulbe  E_pbar    Energy of the produced antiproton in CMF
         *  \param doulbe  pT_pbar   Transverse momentum of the produced antiproton in CMF.
         *  \param doulbe* C_array   Ci=(1,...,16), parameters. C1=C_array[1], C2=C_array[2], ... (C_array[0] is not used)
         *  \return double           Cross section in mbarn/GeV^2
         * */
        static double inv_pp_pbar_CM__Winkler( double s, double E_pbar_d, double pT_pbar, double* C_array, int len_C_array=-1 );
        
        //! Parametrization of the invariant antiproton production (pbar) crosssection from pp in CMF
        /*!
         *  Taken from:     di Mauro, et al.; 2014;
         *                  A new evaluation of the antiproton production cross section for cosmic ray studies;
         *                  DOI: 10.1103/PhysRevD.90.085017
         *
         *  \f[  E_{\bar{p}} \frac{ d^3 \sigma_{pp}^{(\bar{p})} }{d^3 p_{\bar{p}} }  (s, E_{\bar{p}}, p_{T,\bar{p}})  \f]
         *
         *  \param double s         CM energy.
         *  \param doulbe E_pbar    Energy of the produced antiproton in CMF
         *  \param doulbe pT_pbar   Transverse momentum of the produced antiproton in CMF.
         *  \param doulbe* C_array   Ci=(1,...,16), parameters. C1=C_array[1], C2=C_array[2], ... (C_array[0] is not used)
         *  \return double           Cross section in mbarn/GeV^2
         * */
        static double inv_pp_pbar_CM__diMauro( double s, double E_pbar, double pT_pbar, double* C_array, int len_C_array=-1 );
        
        //! Parametrization of the total pp cross section.
        /*!
         *  Taken from:     di Mauro, et al.; 2014;
         *                  A new evaluation of the antiproton production cross section for cosmic ray studies;
         *                  DOI: 10.1103/PhysRevD.90.085017
         *
         *  \param double s         CM energy.
         *  \return double          Cross section in mbarn
         *
         * */
        static double el_pp__diMauro (double s);
        
        //! Parametrization of the elastic pp cross section.
        /*!
         *  Taken from:     di Mauro, et al.; 2014;
         *                  A new evaluation of the antiproton production cross section for cosmic ray studies;
         *                  DOI: 10.1103/PhysRevD.90.085017
         *
         *  \param double s         CM energy.
         *  \return double          Cross section in mbarn
         *
         * */
        static double tot_pp__diMauro(double s);
        
        
        //! Function to read a simple text table.
        /*!
         *  \param std::string file      File name of the table.
         *  \param double array[91][4]   Array to store the table entries.
         * */
        static void   totXS_TableToArray( std::string file, double array[91][4] );
        
        static bool   f_totXS_IsRead;   /// Bool to store whether the tables are alread read.
        
        //! Function to read all the total XS tables.
        /*!
         * */
        static void   totXS_Read();
        
        
        //! Function to interpolate the the XS tables. Interpolation is linear in log-log.
        /*!
         *  To be used with the arrays fXS__*
         *  \param double x              Value in column 0 of array.
         *  \param double array[91][4]   Table which is supposed to be interpolated (use column 1).
         *  \return double               Interpolted value
         * */
        static double totXS_get_interpolation_loglin(double x, double array[91][4] );
        
        
        static double fXS__tot_pbarp[91][4];    /// Array to store the total              pbar+p XS; column 0: T_pbar in GeV, column 1: XS in mbarn
        static double fXS__el_pbarp [91][4];    /// Array to store the elastic            pbar+p XS; column 0: T_pbar in GeV, column 1: XS in mbarn
        static double fXS__tot_pbarD[91][4];    /// Array to store the total              pbar+D XS; column 0: T_pbar in GeV, column 1: XS in mbarn
        static double fXS__nar_pbarD[91][4];    /// Array to store the non-annihilation   pbar+D XS; column 0: T_pbar in GeV, column 1: XS in mbarn
        
        
        //! Interpolation of the total pbar+p cross section.
        /*!
         *  Hand-fitted curve to the data collected by the Particle Data Group (PDG).
         *  The data and tables of the hand-fitted curve are stored in <CRXS dir>/cpp/data
         *
         *  \param  double T_pbar   Kinetic energy of the antiproton (p at rest)
         *  \return double          Cross section in mbarn
         *
         * */
        static double tot_pbarp (double T_pbar);
        //! Interpolation of the elastic pbar+p cross section.
        /*!
         *  Hand-fitted curve to the data collected by the Particle Data Group (PDG).
         *  The data and tables of the hand-fitted curve are stored in <CRXS dir>/cpp/data
         *
         *  \param  double T_pbar   Kinetic energy of the antiproton (p at rest)
         *  \return double          Cross section in mbarn
         *
         * */
        static double el_pbarp  (double T_pbar);
        //! Interpolation of the total pbar+D cross section.
        /*!
         *  Hand-fitted curve to the data collected by the Particle Data Group (PDG).
         *  The data and tables of the hand-fitted curve are stored in <CRXS dir>/cpp/data
         *
         *  \param  double T_pbar   Kinetic energy of the antiproton (D at rest)
         *  \return double          Cross section in mbarn
         *
         * */
        static double tot_pbarD (double T_pbar);
        //! Interpolation of the non-annihilating pbar+D cross section.
        /*!
         *  Hand-fitted curve to the data collected by:
         *
         *      A. Baldini, V. Flaminio, W. G. Moorhead, and D. R. O. Morrison,
         *      Total Cross-Sections for Reactions of High Energy Particles,
         *      edited by S. H., Vol. 12B (SpringerMaterials, 1988).
         *
         *  The data and tables of the hand-fitted curve are stored in <CRXS dir>/cpp/data
         *
         *  \param  double T_pbar   Kinetic energy of the antiproton (D at rest)
         *  \return double          Cross section in mbarn
         *
         * */
        static double nar_pbarD (double T_pbar);
        
        
        
        
        //! Parametrization of the invariant proton scattering crosssection from pp in CMF
        /*!
         *  Taken from:     Anderson, et al.; 1967;
         *                  PROTON AND PION SPECTRA FROM PROTON-PROTON INTERACTIONS AT 10, 20, AND 30 BeV/c*;
         *                  DOI: https://doi.org/10.1103/PhysRevLett.19.198
         *
         *
         *  \param double s         CM energy.
         *  \param doulbe E_p       Energy of the scattered proton in CMF
         *  \param doulbe pT_p      Transverse momentum of the scattered proton in CMF.
         * */
        static double inv_pp_p_CM__Anderson( double s, double E_p, double pT_p, double* C_array=Dummy, int len_C_array=-1 );
    
        
        //! Parametrization of the target and projectile overlap function in pbar production
        /*!
         *  Taken from:     NA49;
         *                  Inclusive production of protons, anti-protons, neutrons, deuterons and tritons in p+C collisions at 158 GeV/c beam momentum;
         *                  arXiv:1207.6520v3
         *
         *        Fig. 69, Tab. 14
         *
         *  \param   double x_F         Feynman parameter p/p_max.
         *  \return  double             Overlap function for projectile
         **/
        static double pbar_overlap_function_projectile(double x_F);
        
        
        
        //! Parametrization of the target and target overlap function in pbar production
        /*!
         *  Taken from:     NA49;
         *                  Inclusive production of protons, anti-protons, neutrons, deuterons and tritons in p+C collisions at 158 GeV/c beam momentum;
         *                  arXiv:1207.6520v3
         *
         *        Fig. 69, Tab. 14
         *
         *  \param   double x_F         Feynman parameter p/p_max.
         *  \return  double             Overlap function for target
         **/
        static double pbar_overlap_function_target    (double x_F);
        
        
        //! Parametrization of the nuclear scaling factor
        /*!
         *    Depending on the parametrization, taken from:     Korsmeier, et al.; 2018;
         *    Production cross sections of cosmic antiprotons in the light of new data from the NA61 and LHCb experiments;
         *          DOI: 10.1103/PhysRevD.97.103019
         *
         *    Taken from:     Winkler, M. W.; 2017;
         *                  Cosmic Ray Antiprotons at High Energies;
         *                  arXiv:1701.04866
         *
         *    In the case of Di Mauro XS we use a simple scaling of A^0.8 for both, projectile and target.
         *
         *    \param double s               CM energy, squared.
         *    \param doulbe xF              Feynman scaling (2*pL/sqrt(s) in CMF)
         *    \param int    A_projectile    Mass number of the projectile
         *    \param int    N_projectile    Number of neutrons in the projectile
         *    \param int    A_target        Mass number of the target
         *    \param int    N_target        Number of neutrons in the target
         *    \param string parametrization Cross section parametrization [Korsmeier_II (default), Korsmeier_I, Winkler, diMauro_I, diMauro_II]
         *
         *    \return double factor         Scaling factor
         **/
        static double factor__AA( double s, double xF, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization );
        
        static double * Get_D_parameters        (int parametrization);
        static double * Get_C_parameters        (int parametrization);
        static double * Get_C_parameters_isospin(int parametrization);
        
        static double deltaHyperon( double s, double* C_array, int len_C_array=-1 );
        static double deltaIsospin( double s, double* C_array, int len_C_array=-1 );
        
        // Parameter definitions:
        
        static double Korsmeier_I_C1_to_C11 [12];
        static double diMauro_I_C1_to_C11   [12];
        static double diMauro_II_C1_to_C11  [12];
        
        static double Winkler_C1_to_C16     [17];
        static double Korsmeier_II_C1_to_C16[17];
        
        static double diMauro_I_C1_to_C16   [17];
        static double diMauro_II_C1_to_C16  [17];
        
        static double Korsmeier_I_D1_to_D2   [3];
        static double Korsmeier_II_D1_to_D2  [3];
        static double Winkler_D1_to_D2       [3];
        static double diMauro_I_D1_to_D2     [3];
        static double diMauro_II_D1_to_D2    [3];
        
        
        
        static double Winkler_SELF_C1_to_C16     [17];
        static double diMauro_SELF_C1_to_C11     [12];
        
        static double diMauro_SELF_C1_to_C16     [17];
        
        static double Winkler_SELF_D1_to_D2       [3];
        static double diMauro_SELF_D1_to_D2       [3];
        
        
        static double Dummy                  [1];
        
    private:

    };
}



#endif
