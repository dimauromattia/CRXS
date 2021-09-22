#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

#include <definitions.h>
#include <configHandler.h>

#include <dmEnergySpectra.h>


const int  CRACS::DM_energy_spectra::finalE    =  4;
const int  CRACS::DM_energy_spectra::finalMU   =  7;
const int  CRACS::DM_energy_spectra::finalTAU  =  10;
const int  CRACS::DM_energy_spectra::finalQ    =  11;
const int  CRACS::DM_energy_spectra::finalC    =  12;
const int  CRACS::DM_energy_spectra::finalB    =  13;
const int  CRACS::DM_energy_spectra::finalT    =  14;
const int  CRACS::DM_energy_spectra::finalW    =  17;
const int  CRACS::DM_energy_spectra::finalZ    =  20;
const int  CRACS::DM_energy_spectra::finalG    =  21;
const int  CRACS::DM_energy_spectra::finalH    =  23;


CRACS::DM_energy_spectra::DM_energy_spectra(){
    Read();
}

CRACS::DM_energy_spectra* CRACS::DM_energy_spectra::fInstance = NULL;

CRACS::DM_energy_spectra* CRACS::DM_energy_spectra::GetInstance(){
    if(!fInstance)
        fInstance = new DM_energy_spectra();
    return fInstance;
}



double CRACS::DM_energy_spectra::GetSpectrum(double mDM, double T, int t){
    
    double E_gev    = T;
    double mDM_gev  = mDM;
    
    double x = log10(E_gev/mDM_gev);
    if (x>=0 || x< -8.9) return 0;
    if (mDM_gev<5 || mDM_gev>100000){
        std::cout << "Warning: dark matter mass out of range! m = " << mDM << std::endl;
    }
    
    // Deal with the fact that for some particles the spectra at mass threshold are not
    // availabel => Shift the mDM_test above first available data in the table
    double mDM_test=mDM_gev;
    if(t==finalT && mDM_gev<180.1 )
        mDM_test=180.1;
    if(t==finalT && mDM_gev<173   )
        return 0;
    if(t==finalH && mDM_gev<130.1 )
        mDM_test=130.1;
    if(t==finalH && mDM_gev<125   )
        return 0;
    if(t==finalZ && mDM_gev<100.1 )
        mDM_test=100.1;
    if(t==finalZ && mDM_gev<45    ) // Approximation for ZZ*
        return 0;
    if(t==finalW && mDM_gev<90.1  )
        mDM_test=90.1;
    if(t==finalW && mDM_gev<44    ) // Approximation for WW*
        return 0;
    
    int m = 0;
    while (fDM_AntiprotonData_cirelli[0][m*fDeltaM_cirelli]<=mDM_test) m++;
    if(m>fDeltaM_cirelli-2) m = fDeltaM_cirelli-2;
    
    int e = (x+8.9)/0.05;
    
    //Interpolation
    double v_ll =   fDM_AntiprotonData_cirelli[t][(m  )*fDeltaM_cirelli+(e  )] ;
    double v_lu =   fDM_AntiprotonData_cirelli[t][(m  )*fDeltaM_cirelli+(e+1)] ;
    double v_ul =   fDM_AntiprotonData_cirelli[t][(m+1)*fDeltaM_cirelli+(e  )] ;
    double v_uu =   fDM_AntiprotonData_cirelli[t][(m+1)*fDeltaM_cirelli+(e+1)] ;
    
    if (v_ll>0) v_ll=log(v_ll); else v_ll=-100;
    if (v_lu>0) v_lu=log(v_lu); else v_lu=-100;
    if (v_ul>0) v_ul=log(v_ul); else v_ul=-100;
    if (v_uu>0) v_uu=log(v_uu); else v_uu=-100;
    
    
    double v_l = v_ll + ((x+8.9)/0.05-e)*(v_lu-v_ll);
    double v_u = v_ul + ((x+8.9)/0.05-e)*(v_uu-v_ul);
    
    double m_l = fDM_AntiprotonData_cirelli[0][ m   *fDeltaM_cirelli];
    double m_u = fDM_AntiprotonData_cirelli[0][(m+1)*fDeltaM_cirelli];
    
    double dNdlogX =  exp(  v_l +  (log(mDM_gev)-log(m_l))/(log(m_u)-log(m_l))*(v_u-v_l)  );
    
    // return dN/dE
    double ret = dNdlogX/E_gev/log(10);
    return ret;
};

double  CRACS::DM_energy_spectra::sourceTerm_pbar(double T_pbar, int type, double m_DM, double sigma_v, double rho){
    return 0.5*pow(rho/m_DM, 2)*sigma_v * GetSpectrum( m_DM, T_pbar, type);
}

double CRACS::DM_energy_spectra::sourceTerm_Dbar(double Tn_Dbar, int type, double m_DM, double sigma_v, double rho, double p_coal){
    // ref. to          Fiorenza Donato, Nicolao Fornengo, and Pierre Salati,
    //                  Antideuterons as a signature of supersymmetric dark matter
    //                  https://doi.org/10.1103/PhysRevD.62.043003
    double T        = 2*Tn_Dbar;
    double p        = sqrt(T*T+2*CRACS::fMass_deuteron*T);
    double T_pbar   = sqrt(p*p/4.+fMass_proton*fMass_proton)-fMass_proton;
    double pref     = 4./3.*pow(p_coal, 3)/p;
    pref           *= CRACS::fMass_deuteron/CRACS::fMass_proton/CRACS::fMass_neutron;
    
    return 0.5 * 2 * pref  *  pow(rho/m_DM, 2)*sigma_v * pow(  GetSpectrum( m_DM, T_pbar, type), 2  );
}

double CRACS::DM_energy_spectra::dTn_N_Dbar(double Tn_Dbar, double m_DM, int type, double p_coal){
    // ref. to          Fiorenza Donato, Nicolao Fornengo, and Pierre Salati,
    //                  Antideuterons as a signature of supersymmetric dark matter
    //                  https://doi.org/10.1103/PhysRevD.62.043003
    double T        = 2*Tn_Dbar;
    double p        = sqrt(T*T+2*CRACS::fMass_deuteron*T);
    double T_pbar   = sqrt(p*p/4.+fMass_proton*fMass_proton)-fMass_proton;
    double pref     = 4./3.*pow(p_coal, 3)/p;
    pref           *= CRACS::fMass_deuteron/CRACS::fMass_proton/CRACS::fMass_neutron;
    
    return  2. * pref * pow(  GetSpectrum( m_DM, T_pbar, type), 2  );
}


double CRACS::DM_energy_spectra::sourceTerm_HeBar(double Tn_Hebar, int type, double m_DM, double sigma_v, double rho, double p_coal){
    double T        = 3*Tn_Hebar;
    double p        = sqrt(T*T+2*CRACS::fMass_helium3*T);
    double T_pbar   = sqrt(p*p/9.+fMass_proton*fMass_proton)-fMass_proton;
    double pref     = pow(p_coal, 3)/p;
    pref           *= CRACS::fMass_helium3/CRACS::fMass_proton/CRACS::fMass_neutron/CRACS::fMass_neutron;
    
    return 0.5 * 3 * 3* pow(pref, 2)  *  pow(rho/m_DM, 2)*sigma_v * pow(  GetSpectrum( m_DM, T_pbar, type), 3  );
}


double CRACS::DM_energy_spectra::dTn_N_Hebar(double Tn_Hebar, double m_DM, int type, double p_coal){
    // ref. to          Fiorenza Donato, Nicolao Fornengo, and Pierre Salati,
    //                  Antideuterons as a signature of supersymmetric dark matter
    //                  https://doi.org/10.1103/PhysRevD.62.043003
    double T        = 3*Tn_Hebar;
    double p        = sqrt(T*T+2*CRACS::fMass_helium3*T);
    double T_pbar   = sqrt(p*p/9.+fMass_proton*fMass_proton)-fMass_proton;
    double pref     = pow(p_coal, 3)/p;
    pref           *= CRACS::fMass_deuteron/CRACS::fMass_proton/CRACS::fMass_neutron/CRACS::fMass_neutron;
    
    return  3. * 3 * pow(pref, 2) * pow(  GetSpectrum( m_DM, T_pbar, type), 3  );
}

//! Read the annihilation spectra
void CRACS::DM_energy_spectra::Read(std::string filename){
    CRACS::ConfigHandler* config = CRACS::ConfigHandler::GetInstance();
    if (filename=="") filename = config->SoftwarePath() + "/data/DM_Spectra_Cirelli/AtProduction_antiprotons.dat";
    float a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a30;
    
    std::ifstream infile(filename.c_str());
    std::string line; std::getline(infile, line);
    
    int count = 0;
    while (infile >> a1 >> a2 >> a3 >> a4 >> a5 >> a6 >> a7 >> a8 >> a9 >> a10 >> a11 >> a12 >> a13 >> a14 >> a15 >> a16 >> a17 >> a18 >> a19 >> a20 >> a21 >> a22 >> a23 >> a24 >> a25 >> a26 >> a27 >> a28 >> a29 >> a30)
    {
        fDM_AntiprotonData_cirelli[ 0][count] = a1;
        fDM_AntiprotonData_cirelli[ 1][count] = a2;
        fDM_AntiprotonData_cirelli[ 2][count] = a3;
        fDM_AntiprotonData_cirelli[ 3][count] = a4;
        fDM_AntiprotonData_cirelli[ 4][count] = a5;
        fDM_AntiprotonData_cirelli[ 5][count] = a6;
        fDM_AntiprotonData_cirelli[ 6][count] = a7;
        fDM_AntiprotonData_cirelli[ 7][count] = a8;
        fDM_AntiprotonData_cirelli[ 8][count] = a9;
        fDM_AntiprotonData_cirelli[ 9][count] = a10;
        fDM_AntiprotonData_cirelli[10][count] = a11;
        fDM_AntiprotonData_cirelli[11][count] = a12;
        fDM_AntiprotonData_cirelli[12][count] = a13;
        fDM_AntiprotonData_cirelli[13][count] = a14;
        fDM_AntiprotonData_cirelli[14][count] = a15;
        fDM_AntiprotonData_cirelli[15][count] = a16;
        fDM_AntiprotonData_cirelli[16][count] = a17;
        fDM_AntiprotonData_cirelli[17][count] = a18;
        fDM_AntiprotonData_cirelli[18][count] = a19;
        fDM_AntiprotonData_cirelli[19][count] = a20;
        fDM_AntiprotonData_cirelli[20][count] = a21;
        fDM_AntiprotonData_cirelli[21][count] = a22;
        fDM_AntiprotonData_cirelli[22][count] = a23;
        fDM_AntiprotonData_cirelli[23][count] = a24;
        fDM_AntiprotonData_cirelli[24][count] = a25;
        fDM_AntiprotonData_cirelli[25][count] = a26;
        fDM_AntiprotonData_cirelli[26][count] = a27;
        fDM_AntiprotonData_cirelli[27][count] = a28;
        fDM_AntiprotonData_cirelli[28][count] = a29;
        fDM_AntiprotonData_cirelli[29][count] = a30;
        
        count++;
    }
    
}

