#ifndef CRXS__XS_H
#define CRXS__XS_H

namespace CRXS {
    
    enum parametrization{
        KORSMEIER_I    = 1,
        KORSMEIER_II   = 2,
        WINKLER        = 3,
        DI_MAURO_I     = 4,
        DI_MAURO_II    = 5,
        ANDERSON       = 6
    };
    enum product{
        P_BAR    = 1,
        D_BAR    = 2,
        HE_BAR   = 3
    };
    
    enum coalescence{
        FIXED_P0                   = 1,
        ENERGY_DEP__VAN_DOETINCHEM = 2
    };
    
    class XS{
        
    public:
        
        
        //!Convert LAB frame kinetic variable to the CM frame. (The LAB frame is the ISM rest frame.)
        /*!
         *
         *  \param double Tn_proj_LAB    Kinetic energy per nucleus of the prjectile (in the LAB frame)
         *  \param doulbe T_prod_LAB     Kinetic energy of the product (in the LAB frame)
         *  \param doulbe eta_LAB        Pseudo rapidity of the product (in the LAB frame)
         *
         *  \param doulbe &s             Returns: CM energy
         *  \param doulbe &E_prod        Returns: Antiproton (total) energy
         *  \param doulbe &pT_prod       Returns: Product transverse momentum
         *  \param doulbe &x_F           Returns: Product Feynman scaling variable
         *
         *  \param int    product        Product (from enum [P_BAR, D_BAR, HE_BAR])
         *  \return bool                 True if conversion successful
         */
        static bool convert_LAB_to_CM( const double T_p_LAB, const double T_prod_LAB, const double eta_LAB, double &s, double &E_prod, double &pT_prod, double &x_F, int product=P_BAR );
        
         //!Invariant antiproton production cross section for general projectile and target nucleus for different XS parametrization
         /*!
          *  \param double s               CM energy, squared.
          *  \param doulbe xF              Feynman scaling (2*pL_pbar/sqrt(s) in CMF)
          *  \param doulbe pT_pbar         Transverse momentum of the antiproton
          *  \param int    A_projectile    Mass number of the projectile
          *  \param int    N_projectile    Number of neutrons in the projectile
          *  \param int    A_target        Mass number of the target
          *  \param int    N_target        Number of neutrons in the target
          *  \param int    parametrization Cross section parametrization, enum from[KORSMEIER_II (default), KORSMEIER_I, WINKLER, DI_MAURO_I, DI_MAURO_II]
          *
          *  \return double XS             Cross section in mbarn/GeV^2
          */
        static double inv_AA_pbar_CM( double s, double xF, double pT_pbar, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization);

        //! Invariant antiproton production cross section for general projectile and target nucleus for different XS parametrization as function of LAB frame kinetic variables
        /*!
         *
         *  \param double Tn_proj_LAB      Kinetic energy per nucleus of the prjectile (in the LAB frame)
         *  \param doulbe T_pbar_LAB       Kinetic energy of the antiproton (in the LAB frame)
         *  \param doulbe eta_LAB          Pseudo rapidity of the antiproton (in the LAB frame)
         *  \param int    A_projectile     Mass number of the projectile
         *  \param int    N_projectile     Number of neutrons in the projectile
         *  \param int    A_target         Mass number of the target
         *  \param int    N_target         Number of neutrons in the target
         *  \param int    parametrization  Cross section parametrization [Korsmeier_II (default), Korsmeier_I, Winkler, diMauro_I, diMauro_II]
         *
         *  \return double XS              Cross section in mbarn/GeV^2
         */
        static double inv_AA_pbar_LAB( double Tn_proj_LAB, double T_pbar_LAB, double eta_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization );
        
        
        //! Energy-differential antiproton production cross section for general projectile and target nucleus for different XS parametrization as function of LAB frame kinetic variables.
        /*!
         *  This cross section is integrated over all angles.
         *
         *  \param double Tn_proj_LAB      Kinetic energy per nucleus of the prjectile (in the LAB frame)
         *  \param doulbe T_pbar_LAB       Kinetic energy of the antiproton (in the LAB frame)
         *  \param int    A_projectile     Mass number of the projectile
         *  \param int    N_projectile     Number of neutrons in the projectile
         *  \param int    A_target         Mass number of the target
         *  \param int    N_target         Number of neutrons in the target
         *  \param int    parametrization  Cross section parametrization, enum from[KORSMEIER_II (default), KORSMEIER_I, WINKLER, DI_MAURO_I, DI_MAURO_II]
         *
         *  \return double XS              Cross section in mbarn/GeV
         */
        static double dE_AA_pbar_LAB( double Tn_proj_LAB, double T_pbar_LAB, int A_projectile=1, int N_projectile=0, int A_target=1, int N_target=0, int parametrization=KORSMEIER_II);
        //! Helper function for dE_AA_pbar_LAB
        static double integrand__dE_AA_pbar_LAB (double eta_LAB, void* parameters  );
        
        
        
         //!Energy-differential antiproton production cross section including antineutrons and antihyperons for general projectile and target nucleus and for different XS parametrization as function of LAB frame kinetic variables. This cross section is integrated over all angles.
         /*! This function should only be used for the parametrizations: Winkler, Korsmeier_I and Korsmeier_II.
          *   For diMauro apply a global factor 2.3 instead.
          *
          *  \param double Tn_proj_LAB      Kinetic energy per nucleus of the prjectile (in the LAB frame)
          *  \param doulbe T_pbar_LAB       Kinetic energy of the antiproton (in the LAB frame)
          *  \param int    A_projectile     Mass number of the projectile
          *  \param int    N_projectile     Number of neutrons in the projectile
          *  \param int    A_target         Mass number of the target
          *  \param int    N_target         Number of neutrons in the target
          *  \param int    parametrization  Cross section parametrization, enum from[KORSMEIER_II (default), KORSMEIER_I, WINKLER, DI_MAURO_I, DI_MAURO_II]
          *
          *  \return double XS              Cross section in mbarn/GeV
          */
        static double dE_AA_pbar_LAB_incNbarAndHyperon(double Tn_proj_LAB, double T_pbar_LAB, int A_projectile=1, int N_projectile=0, int A_target=1, int N_target=0, int parametrization=KORSMEIER_II);
        
        
        
        //!Invariant proton production cross section for general projectile and target nucleus for different XS parametrization
        /*!
         *  \param double s               CM energy, squared.
         *  \param doulbe xF              Feynman scaling (2*pL_p/sqrt(s) in CMF)
         *  \param doulbe pT_p            Transverse momentum of the proton
         *  \param int    A_projectile    Mass number of the projectile
         *  \param int    N_projectile    Number of neutrons in the projectile
         *  \param int    A_target        Mass number of the target
         *  \param int    N_target        Number of neutrons in the target
         *  \param int    parametrization Cross section parametrization, enum from[ANDERSON]
         *
         *  \return double XS             Cross section in mbarn/GeV^2
         */
        static double inv_AA_p_CM( double s, double xF, double pT_p, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization);
        
        //! Invariant proton production cross section for general projectile and target nucleus for different XS parametrization as function of LAB frame kinetic variables
        /*!
         *
         *  \param double Tn_proj_LAB      Kinetic energy per nucleus of the prjectile (in the LAB frame)
         *  \param doulbe T_p_LAB          Kinetic energy of the proton (in the LAB frame)
         *  \param doulbe eta_LAB          Pseudo rapidity of the proton (in the LAB frame)
         *  \param int    A_projectile     Mass number of the projectile
         *  \param int    N_projectile     Number of neutrons in the projectile
         *  \param int    A_target         Mass number of the target
         *  \param int    N_target         Number of neutrons in the target
         *  \param int    parametrization  Cross section parametrization [ANDERSON]
         *
         *  \return double XS              Cross section in mbarn/GeV^2
         */
        static double inv_AA_p_LAB( double Tn_proj_LAB, double T_p_LAB, double eta_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization );
        
        
        //! Energy-differential proton production cross section for general projectile and target nucleus for different XS parametrization as function of LAB frame kinetic variables.
        /*!
         *  This cross section is integrated over all angles.
         *
         *  \param double Tn_proj_LAB      Kinetic energy per nucleus of the prjectile (in the LAB frame)
         *  \param doulbe T_p_LAB          Kinetic energy of the proton (in the LAB frame)
         *  \param int    A_projectile     Mass number of the projectile
         *  \param int    N_projectile     Number of neutrons in the projectile
         *  \param int    A_target         Mass number of the target
         *  \param int    N_target         Number of neutrons in the target
         *  \param int    parametrization  Cross section parametrization, enum from[ANDERSON]
         *
         *  \return double XS              Cross section in mbarn/GeV
         */
        static double dE_AA_p_LAB( double Tn_proj_LAB, double T_p_LAB, int A_projectile=1, int N_projectile=0, int A_target=1, int N_target=0, int parametrization=ANDERSON);
        //! Helper function for dE_AA_p_LAB
        static double integrand__dE_AA_p_LAB (double eta_LAB, void* parameters  );
        
        
        
        
        
        //!Invariant antideuteron production cross section for general projectile and target nucleus for different XS parametrization
        /*!
         *
         *      Coalesence momentum from the paper:
         *          Deuteron and Antideuteron Production Simulation in Cosmic-ray Interactions
         *          Diego-Mauricio Gomez-Coral, et al.,
         *          (DOI 10.1103/PhysRevD.98.023012)
         *
         *      Taken from Eq. (5) with parameters of Korsmeier, et al.
         *      The parameter p_coal is defined as abs(p_proton - p_neutron)/2. (Note that there is also a different notation,
         *      without the factor 2., in the literature)
         *
         *
         *  \param  double s              CM energy, squared (in the nucleon-nucleon frame).
         *  \return double                Coalescence momentum
         */
        static double p_coal__VonDoetinchen( double s );
        
        
        //!Invariant antideuteron production cross section for general projectile and target nucleus for different XS parametrization
        /*!
         *      We calculat the cross section in the analytic coalescence model. With the formula:
         *
         *
         *
         *                 3                                                   3                   3
         *                d    sigma                m               3         d  sigma            d  sigma
         *                 dbar               1      D    4pi  pcoal                  pbar                nbar
         *        E      ------------  =  -------- -----  ---  ------   E     ------------  E     ------------
         *         dbar       3           sigmaTot m  m    3      1      pbar      3         nbar      3
         *                  dk                      p  n                         dk                  dk
         *                    dbar                                                 pbar                nbar
         *
         *        with:
         *
         *         3              3                 /   3                        3                                                              \
         *        d  sigma       d  sigma           |  d  sigma                 d  sigma                                                        |
         *                pbar           nbar     1 |          pbar                     nbar                                                    |
         *        ------------   ------------  =  - |  ------------(sS, vk    ) ------------ (sS - 2E    , vk    )  +  ({pbar} < - > {nbar})    |
         *             3              3           2 |       3             pbar       3               pbar    nbar                               |
         *           dk             dk              |     dk                       dk                                                           |
         *             pbar           nbar          \       pbar                     nbar                                                       /
         *
         *
         *      The parameter p_coal is defined as abs(p_proton - p_neutron)/2. (Note that there is also a different notation,
         *      without the factor 2., in the literature)
         *
         *
         *      There are two options for the coalesence momentum (cf. to option \param int coalescence):
         *
         *          1) Fixed to 80 MeV, which is the value tuned to the aleph experiment. If you prefer a different (fixed) value
         *             you can rescale the whole XS with (p_coal/80 MeV)^3.
         *          2) The energy-dependent coalescence momentum suggested in DOI 10.1103/PhysRevD.98.023012.
         *
         *      We recommend option 2)
         *
         *      The parameter p_coal is defined as abs(p_proton - p_neutron)/2. (Note that there is also a different notation,
         *      without the factor 2., in the literature)
         *
         *  \param double s               CM energy, squared (in the nucleon-nucleon frame).
         *  \param doulbe xF              Feynman scaling (2*pL_pbar/sqrt(s) in CMF, i.e. the nucleon-nucleon frame)
         *  \param doulbe pT_Dbar         Transverse momentum of the antideuteron
         *  \param int    A_projectile    Mass number of the projectile
         *  \param int    N_projectile    Number of neutrons in the projectile
         *  \param int    A_target        Mass number of the target
         *  \param int    N_target        Number of neutrons in the target
         *  \param int    parametrization Cross section parametrization, enum from[KORSMEIER_II (default), KORSMEIER_I, WINKLER, DI_MAURO_I, DI_MAURO_II]
         *  \param int    coalescence     Coalescence model, enum from[FIXED_P0 (default), ENERGY_DEP__VAN_DOETINCHEM]
         *
         *  \return double XS             Cross section in mbarn/GeV^2
         */
        static double inv_AA_Dbar_CM( double s, double xF_Dbar, double pT_Dbar, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization, int coalescence);
        
        //! Invariant antideuteron production cross section for general projectile and target nucleus for different XS parametrization as function of LAB frame kinetic variables
        /*!
         *  \param double Tn_proj_LAB      Kinetic energy per nucleus of the prjectile (in the LAB frame)
         *  \param doulbe T_Dbar_LAB       Kinetic energy of the antideuteron (in the LAB frame)
         *  \param doulbe eta_LAB          Pseudo rapidity of the antideuteron (in the LAB frame)
         *  \param int    A_projectile     Mass number of the projectile
         *  \param int    N_projectile     Number of neutrons in the projectile
         *  \param int    A_target         Mass number of the target
         *  \param int    N_target         Number of neutrons in the target
         *  \param int    parametrization  Cross section parametrization [Korsmeier_II (default), Korsmeier_I, Winkler, diMauro_I, diMauro_II]
         *  \param int    coalescence      Coalescence model, enum from[FIXED_P0 (default), ENERGY_DEP__VAN_DOETINCHEM], cf. inv_AA_Dbar_CM
         *
         *  \return double XS              Cross section in mbarn/GeV^2
         */
        static double inv_AA_Dbar_LAB( double Tn_proj_LAB, double T_Dbar_LAB, double eta_LAB, int A_projectile, int N_projectile, int A_target, int N_target, int parametrization, int coalescence);
        //! Helper function for inv_AA_Dbar_LAB
        static double integrand__dE_AA_Dbar_LAB (double eta_LAB, void* parameters  );
        
        
        //! Energy-differential antideuteron production cross section for general projectile and target nucleus for different XS parametrization as function of LAB frame kinetic variables.
        /*!
         *  This cross section is integrated over all angles.
         *
         *  \param double Tn_proj_LAB      Kinetic energy per nucleus of the prjectile (in the LAB frame)
         *  \param doulbe Tn_Dbar_LAB      Kinetic energy per nucleus of the antideuteron (in the LAB frame)
         *  \param int    A_projectile     Mass number of the projectile
         *  \param int    N_projectile     Number of neutrons in the projectile
         *  \param int    A_target         Mass number of the target
         *  \param int    N_target         Number of neutrons in the target
         *  \param int    parametrization  Cross section parametrization, enum from[KORSMEIER_II (default), KORSMEIER_I, WINKLER, DI_MAURO_I, DI_MAURO_II]
         *  \param int    coalescence      Coalescence model, enum from[FIXED_P0 (default), ENERGY_DEP__VAN_DOETINCHEM], cf. inv_AA_Dbar_CM
         *
         *  \return double XS              Cross section in mbarn/GeV
         */
        static double dEn_AA_Dbar_LAB( double Tn_proj_LAB, double Tn_Dbar_LAB, int A_projectile=1, int N_projectile=0, int A_target=1, int N_target=0, int parametrization=KORSMEIER_II, int coalescence=ENERGY_DEP__VAN_DOETINCHEM );
        
        
        
        
        
    private:




        
    };
}

#endif




