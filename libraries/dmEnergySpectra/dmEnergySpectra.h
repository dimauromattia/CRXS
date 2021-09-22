#ifndef DM_ENERGY_SPECTRA_H
#define DM_ENERGY_SPECTRA_H


namespace CRACS {
    
    /*! \brief Get energy spectra for antiprotons as tabulated by Cirelli (PPP4DM).
     *
     */

    class DM_energy_spectra{
        
    public:
        
        static DM_energy_spectra* GetInstance();
        
        
        static const int  finalE;
        static const int  finalMU;
        static const int  finalTAU;
        static const int  finalQ;
        static const int  finalC;
        static const int  finalB;
        static const int  finalT;
        static const int  finalW;
        static const int  finalZ;
        static const int  finalG;
        static const int  finalH;
        
        
        //! Calculates dN/dE (mDM, T, SM_finalState)
        /*!
         *  \param double mDM.  Dark matter mass (in GeV).
         *  \param double T.    Kinetic Energy (in GeV).
         *  \param int t.       Final state.
         */
        double GetSpectrum(double mDM, double T, int t);
        
        //! Anti proton source term \f q = d^3N/(dV dt dT) \f
        /*!
         *  \param double T.        Kinetic Energy (in GeV).
         *  \param int type.        Final state, default: DM_energy_spectra::finalB.
         *  \param double m_DM.     Dark Matter Mass in GeV.
         *  \param double sigma_v.  Termally averaged cross section in \f m^3/s \f
         *  \param double rho.      Dark Matter Density in \f GeV/m^3 \f .
         */
        double  sourceTerm_pbar(double T_pbar, int type=13, double m_DM=100, double sigma_v=3e-26*1e6, double rho=0.3*1e-6);
        
        
        //! Anti deuteron source term \f q = d^3N/(dV dt d(T/n)) \f
        /*!
         *  \param double Tn.       Kinetic Energy per nucleon (in GeV).
         *  \param int type.        Final state, default: DM_energy_spectra::finalB.
         *  \param double m_DM.     Dark Matter Mass in GeV.
         *  \param double sigma_v.  Termally averaged cross section in \f m^3/s \f, default: \f 3\,10^{20} m^3/s \f.
         *  \param double rho.      Dark Matter Density in \f GeV/m^3 \f, default: \f 3\,10^{-6} GeV/m^3 \f .
         */
        double sourceTerm_Dbar(double Tn_Dbar, int type=13, double m_DM=100, double sigma_v=3e-26*1e6, double rho=0.3*1e-6, double p_coal=0.062);
        
        //! Anti deuteron energy distribution \f q = dN/dE \f
        /*!
         *  \param double Tn.       Kinetic Energy per nucleon (in GeV).
         *  \param int type.        Final state, default: DM_energy_spectra::finalB.
         */
        double dTn_N_Dbar(double Tn_Dbar, double m_DM=100,  int type=13, double p_coal=0.062);
        
        
        //! Anti helium (3) source term \f q = d^3N/(dV dt d(T/n))
        /*!
         *  \param double Tn.       Kinetic Energy per nucleon (in GeV).
         *  \param int type.        Final state, default: DM_energy_spectra::finalB.
         *  \param double m_DM.     Dark Matter Mass in GeV.
         *  \param double sigma_v.  Termally averaged cross section in \f m^3/s \f, default: \f 3\,10^{20} m^3/s \f.
         *  \param double rho.      Dark Matter Density in \f GeV/m^3 \f, default: \f 3\,10^{-6} GeV/m^3 \f .
         */
        double sourceTerm_HeBar(double Tn_Hebar, int type=13, double m_DM=100, double sigma_v=3e-26*1e6, double rho=0.3*1e-6, double p_coal=0.062);
        
        //! Anti helium (3) energy distribution \f q = dN/dE \f
        /*!
         *  \param double Tn.       Kinetic Energy per nucleon (in GeV).
         *  \param int type.        Final state, default: DM_energy_spectra::finalB.
         */
        double dTn_N_Hebar(double Tn_Hebar, double m_DM=100,  int type=13, double p_coal=0.062);
        
        
    private:
        static DM_energy_spectra* fInstance;
        DM_energy_spectra();
        
        static const int   fNm_cirelli     = 62;
        static const int   fDeltaM_cirelli = 179;
        float       fDM_AntiprotonData_cirelli[30][fNm_cirelli*fDeltaM_cirelli];
        
        
        //! Read the annihilation spectra
        void Read(std::string filename="");

        
        
    };
    


}



//end DM_ENERGY_SPECTRA_H

#endif
