"""@package XS_tools
    Documentation of the module XS_tools.
    
    The module provides an interface to CRXS to read the Lorentz invariant XS definitions from cpp.
    Then it provieds some function to transform the XSs to the LAB frame and integrate over all angles.
    """

import numpy  as np
import scipy.integrate   as integrate
import math

import cpp.xs_tools      as cpp_xs


fMass_proton         = 0.9382720813;     # PDG Review 2016

# ------------------------------------------------------------- #
#  Cross secion definitions:                                    #
#                                                               #
#  Name convention: T_XY_Z_F_P                                  #
#                                                               #
#     T     type                                                #
#           inv :   Lorentz invariant XS, E d^3 sigma/dp^3      #
#           dE  :   Energy-differential XS, d sigma/dE          #
#     XY    projectile X and target Y                           #
#     Z     product                                             #
#     F     frame (CM or LAB=ISM)                               #
#     P     parametrizationa and description                    #
#                                                               #
# ------------------------------------------------------------- #


def inv_pp_pbar_CM__Winkler_p(s, E_pbar, pT_pbar, C_array):
    """
        Parametrization of the invariant antiproton production (pbar) crosssection from pp in CMF
        
        Taken from:     Winkler, M. W.; 2017;
        Cosmic Ray Antiprotons at High Energies;
        arXiv:1701.04866
        
        All cross sections are given in mbarn
        All enegies, momenta, and masses have unit GeV.
        
        \f[  E_{\bar{p}} \frac{ d^3 \sigma_{pp}^{(\bar{p})} }{d^3 p_{\bar{p}} }  (s, E_{\bar{p}}, p_{T,\bar{p}})  \f]
        
        \param double s         CM energy.
        \param doulbe E_pbar    Energy of the produced antiproton in CMF
        \param doulbe pT_pbar   Transverse momentum of the produced antiproton in CMF.
        \param doulbe C_array   Ci=(1,...,16), parameters. C1=C_array[1], C2=C_array[2], ... (C_array[0] is not used)
        
        \return double XS       Cross section in mbarn/GeV^2
        """
    return cpp_xs.inv_pp_pbar_CM__Winkler(s, E_pbar, pT_pbar, C_array) # read function form cpp/cpp_xs.cpp



def inv_pp_pbar_CM__diMauro_p(s, E_pbar, pT_pbar, C_array):
    """
        Parametrization of the invariant antiproton production (pbar) crosssection from pp in CMF
        
        Taken from:     di Mauro, et al.; 2014;
        A new evaluation of the antiproton production cross section for cosmic ray studies;
        DOI: 10.1103/PhysRevD.90.085017
        
        All cross sections are given in mbarn
        All enegies, momenta, and masses have unit GeV.
        
        \f[  E_{\bar{p}} \frac{ d^3 \sigma_{pp}^{(\bar{p})} }{d^3 p_{\bar{p}} }  (s, E_{\bar{p}}, p_{T,\bar{p}})  \f]
        
        \param double s         CM energy.
        \param doulbe E_pbar    Energy of the produced antiproton in CMF
        \param doulbe pT_pbar   Transverse momentum of the produced antiproton in CMF.
        \param doulbe* C_array  Ci=(1,...,11), parameters. C1=C_array[1], C2=C_array[2], ... (C_array[0] is not used)
        
        \return double XS       Cross section in mbarn/GeV^2
        """
    return cpp_xs.inv_pp_pbar_CM__diMauro(s, E_pbar, pT_pbar, C_array) # read function form cpp/cpp_xs.cpp


# ------------------------------------------------------------- #
#  Parameter definitions:                                       #
# ------------------------------------------------------------- #
#
#
# Parameters from
# 1)  Korsmeier et al. 2018 (Param I):
# 2)  Mauro, et al. (Eq 12); 2014 DOI: 10.1103/PhysRevD.90.085017 (recommended by KDD18):
# 3)  Mauro, et al. (Eq 13); 2014 DOI: 10.1103/PhysRevD.90.085017:
#
#       The names of the parameters correspond to Mauro, et al.
#
#                                    [  C0,   C1,          C2,          C3,           C4,          C5,          C6,          C7,          C8,            C9,        C10,    C11   ]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Korsmeier_I_C1_to_C11    =  np.array([  -1,   3.50193e+00, 5.58513e+00, 3.99553e-02, -2.50716e-01, 2.65053e+00, 3.78145e-02, 4.29478e-02, 2.69520e+00,   0.0,       0.0,    0.0   ])
diMauro_I_C1_to_C11      =  np.array([  -1,   4.499,       3.41,        0.00942,      0.445,       3.502,       0.0622,      -0.247,      2.576,         0.0,       0.0,    0.0   ])
diMauro_II_C1_to_C11     =  np.array([  -1,   4.448,       3.735,       0.00502,      0.708,       3.527,       0.236,       -0.729,      2.517,        -1.822e-11, 3.527,  0.384 ])
#
#
#       The names of the parameters correspond to Korsmeier et al. 2018
#
#                                    [  D0,   D1,          D2          ]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Korsmeier_I_D1_to_D2    =  np.array([  -1,  0.825, 0.167])
#
#
# Parameters from
# 1)  Winker 2017 arXiv:1701.04866:
# 2)  Korsmeier 2018 (Param II):
#
#       The names of the parameters correspond to Winker 2017
#
#                                    [ C0,    C1,     C2,   C3,     C4,   C5,          C6,      C7,          C8,    C9,          C10,         C11,   C12,    C13,   C14,   C15,    C16  ]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Winkler_C1_to_C16        =  np.array([  -1,   0.31,   0.30, 21316., 0.9,  0.047,       7.76,    0.168,       0.038, 1.0e-3,      0.7,         30.9,  -1.74,  0.71,  0.114, 20736., 0.51 ])
Korsmeier_II_C1_to_C16   =  np.array([  -1,   0.31,   0.30, 21316., 0.9,  5.01767e-02, 7.79045, 1.64809e-01, 0.038, 4.74370e-04, 3.70480e+00, 30.9,  -1.74,  0.71,  0.114, 20736., 0.51 ])
#
#
#       The names of the parameters correspond to Korsmeier et al. 2018
#
#                                    [  D0,  D1,    D2     ]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
Korsmeier_II_D1_to_D2    =  np.array([  -1, 0.828, 0.145  ])
Winkler_D1_to_D2         =  np.array([  -1, 0.839, 0.161  ])      # value of <nu_He>=1.25 is translated to D_1 and D_2 (see Korsmeier et al. 2018)
#
#
#
parameters_C     = {'WINKLER': Winkler_C1_to_C16, 'KORSMEIER_II' : Korsmeier_II_C1_to_C16, 'KORSMEIER_I':Korsmeier_I_C1_to_C11, 'DI_MAURO_I': diMauro_I_C1_to_C11, 'DI_MAURO_II': diMauro_II_C1_to_C11 }
parameters_D     = {'WINKLER': Winkler_D1_to_D2,  'KORSMEIER_II' : Korsmeier_II_D1_to_D2,  'KORSMEIER_I':Korsmeier_I_D1_to_D2    }
#
NbarAndHyperon_C = {'WINKLER': Winkler_C1_to_C16, 'KORSMEIER_II' : Korsmeier_II_C1_to_C16, 'KORSMEIER_I':Korsmeier_II_C1_to_C16  }
# 'Korsmeier_I':Korsmeier_II_C1_to_C16 is correct, we use the same nbar and hyperon parameters for both Korsmeier et al parametrizations




def factor__AA( s, xF, A_projectile, N_projectile, A_target, N_target, parametrization ):
    """
        Parametrization of the nuclear scaling factor
        
        Taken from:     Korsmeier, et al.; 2018;
        Production cross sections of cosmic antiprotons in the light of new data from the NA61 and LHCb experiments;
        DOI: 10.1103/PhysRevD.97.103019
        
        \param double s               CM energy, squared.
        \param doulbe xF              Feynman scaling (2*pL/sqrt(s) in CMF)
        \param int    A_projectile    Mass number of the projectile
        \param int    N_projectile    Number of neutrons in the projectile
        \param int    A_target        Mass number of the target
        \param int    N_target        Number of neutrons in the target
        \param string parametrization Cross section parametrization [Korsmeier_II (default), Korsmeier_I, Winkler, diMauro_I, diMauro_II]
        
        \return double factor         Scaling factor
        """
    D_array = parameters_D[parametrization]
    proj = np.power(A_projectile, D_array[2])*(1+cpp_xs.deltaIsospin(s,NbarAndHyperon_C[parametrization])*N_projectile/A_projectile)*cpp_xs.pbar_overlap_function_projectile( xF )
    targ = np.power(A_target,     D_array[2])*(1+cpp_xs.deltaIsospin(s,NbarAndHyperon_C[parametrization])*N_target    /A_target    )*cpp_xs.pbar_overlap_function_target    ( xF )
    return np.power(A_projectile*A_target, D_array[1])*( proj + targ )


# ------------------------------------------------------------- #
#  pp to AA scaling:                                            #
# ------------------------------------------------------------- #

def inv_AA_pbar_CM(s, xF, pT_pbar, A_projectile, N_projectile, A_target, N_target, parametrization='KORSMEIER_II'):
    """
        Invariant cross section for general projectile and target nucleus for different XS parametrization
        
        \param double s               CM energy, squared.
        \param doulbe xF              Feynman scaling (2*pL_pbar/sqrt(s) in CMF)
        \param doulbe pT_pbar         Transverse momentum of the antiproton
        \param int    A_projectile    Mass number of the projectile
        \param int    N_projectile    Number of neutrons in the projectile
        \param int    A_target        Mass number of the target
        \param int    N_target        Number of neutrons in the target
        \param string parametrization Cross section parametrization [KORSMEIER_II (default), KORSMEIER_I, WINKLER, DI_MAURO_I, DI_MAURO_II]
        
        \return double XS             Cross section in mbarn/GeV^2
        """
    pL_pbar = xF*np.sqrt(s)/2.
    E_pbar  = np.sqrt( fMass_proton*fMass_proton + pL_pbar*pL_pbar + pT_pbar*pT_pbar )
    
    if   parametrization == 'WINKLER':
        global Winkler_C1_to_C16
        return inv_pp_pbar_CM__Winkler_p(s, E_pbar, pT_pbar, Winkler_C1_to_C16      ) * factor__AA( s, xF, A_projectile, N_projectile, A_target, N_target, parametrization )
    elif parametrization == 'KORSMEIER_II':
        global Korsmeier_II_C1_to_C16
        return inv_pp_pbar_CM__Winkler_p(s, E_pbar, pT_pbar, Korsmeier_II_C1_to_C16 ) * factor__AA( s, xF, A_projectile, N_projectile, A_target, N_target, parametrization )
    elif parametrization == 'KORSMEIER_I':
        global Korsmeier_I_C1_to_C11
        return inv_pp_pbar_CM__diMauro_p(s, E_pbar, pT_pbar, Korsmeier_I_C1_to_C11  ) * factor__AA( s, xF, A_projectile, N_projectile, A_target, N_target, parametrization )
    elif parametrization == 'DI_MAURO_I'  and 1000*A_projectile+100*N_projectile+10*A_target+1*N_target==1010:
        global diMauro_I_C1_to_C11
        return inv_pp_pbar_CM__diMauro_p(s, E_pbar, pT_pbar, diMauro_I_C1_to_C11    )
    elif parametrization == 'DI_MAURO_II' and 1000*A_projectile+100*N_projectile+10*A_target+1*N_target==1010:
        global diMauro_II_C1_to_C11
        return inv_pp_pbar_CM__diMauro_p(s, E_pbar, pT_pbar, diMauro_II_C1_to_C11   )
    else:
        print ('Error in inv_AA_pbar_CM. parametrization "%s" unknown for (A_proj, N_proj, A_targ, N_targ) = (%i,%i,%i,%i).' % (parametrization, A_projectile, N_projectile, A_target, N_target))
        return 0



# ------------------------------- #
#  Variable transformation:       #
#   -> LAB frame to CM frame      #
# ------------------------------- #

def convert_LAB_to_CM( Tn_proj_LAB, T_pbar_LAB, eta_LAB ):
    """
        Convert LAB frame kinetic variable to the CM frame. (The LAB frame is the ISM rest frame.)
        
        \param double Tn_proj_LAB    Kinetic energy per nucleus of the prjectile (in the LAB frame)
        \param doulbe T_pbar_LAB     Kinetic energy of the antiproton (in the LAB frame)
        \param doulbe eta_LAB        Pseudo rapidity of the antiproton (in the LAB frame)
       
       \return (double s, double E_pbar, double pT_pbar, double x_F)   List with CM fram variables: CM energy, Antiproton (total) energy, Antiproton transverse momentum, Feynman scaling variable)
        """
    global fMass_proton
    
    s = 4*fMass_proton*fMass_proton + 2 * Tn_proj_LAB * fMass_proton
    
    p_pbar_LAB    = np.sqrt(  T_pbar_LAB*(T_pbar_LAB+2*fMass_proton)  )
    En_proj_LAB       = Tn_proj_LAB+fMass_proton;
    E_pbar_LAB    = T_pbar_LAB+fMass_proton;
    
    beta          = np.sqrt(En_proj_LAB - fMass_proton)/np.sqrt(En_proj_LAB + fMass_proton);
    gamma         = 1./np.sqrt(1 - beta*beta);
    gammabeta     = gamma * beta;
    
    E_pbar        =  gamma     * E_pbar_LAB - gammabeta * p_pbar_LAB*np.tanh(eta_LAB);
    pL_pbar       = -gammabeta * E_pbar_LAB + gamma     * p_pbar_LAB*np.tanh(eta_LAB);
    pT_pbar       = p_pbar_LAB/np.cosh(eta_LAB);
    
    #E_pbar_Max    = ( s-8.*fMass_proton*fMass_proton )/2./np.sqrt( s );
    #x_R           = E_pbar/E_pbar_Max;
    
    xF            = 2. * pL_pbar / np.sqrt(s)
    
    return s, E_pbar, pT_pbar, xF


# ------------------------------------------------------------- #
#  XS Conversions:                                              #
# ------------------------------------------------------------- #



def inv_AA_pbar_LAB(Tn_proj_LAB, T_pbar_LAB, eta_LAB, A_projectile, N_projectile, A_target, N_target, parametrization='KORSMEIER_II'):
    """
        Invariant cross section for general projectile and target nucleus for different XS parametrization
        as function of LAB frame kinetic variables
        
        \param double Tn_proj_LAB      Kinetic energy per nucleus of the prjectile (in the LAB frame)
        \param doulbe T_pbar_LAB       Kinetic energy of the antiproton (in the LAB frame)
        \param doulbe eta_LAB          Pseudo rapidity of the antiproton (in the LAB frame)
        \param int    A_projectile     Mass number of the projectile
        \param int    N_projectile     Number of neutrons in the projectile
        \param int    A_target         Mass number of the target
        \param int    N_target         Number of neutrons in the target
        \param string parametrization  Cross section parametrization [KORSMEIER_II (default), KORSMEIER_I, WINKLER, DI_MAURO_I, DI_MAURO_II]
        
        \return double XS              Cross section in mbarn/GeV^2
        """
    s, _, pT_pbar, xF = convert_LAB_to_CM( Tn_proj_LAB, T_pbar_LAB, eta_LAB )
    return inv_AA_pbar_CM(s, xF, pT_pbar, A_projectile, N_projectile, A_target, N_target, parametrization)


#  Help functions of the vectorized version below
def _dE_AA_pbar_LAB(Tn_proj_LAB, T_pbar_LAB, A_projectile=1, N_projectile=0, A_target=1, N_target=0, parametrization='Korsmeier_II'):
    #
    #  Integrate over all solid angle and transform to enery differential (d sigma / d E)
    #
    global fMass_proton
    E_pbar_LAB = T_pbar_LAB + fMass_proton
    p_pbar_LAB = np.sqrt(  np.power( E_pbar_LAB, 2 ) - np.power( fMass_proton, 2 )  );
    if p_pbar_LAB!=p_pbar_LAB:
        return 0
    Jacobian_and_conversion      = 2*math.pi*p_pbar_LAB
    # it contains: phi_integration (2 pi), inv to d3p (1/E_pbar_LAB), Jacobian(p_pbar_Lab*p_pbar_Lab), dp to dE (E_pbar_LAB/p_pbar_Lab)
    res = integrate.quad( lambda eta: np.power(np.cosh(eta), -2) * inv_AA_pbar_LAB(Tn_proj_LAB, T_pbar_LAB, eta, A_projectile, N_projectile, A_target, N_target, parametrization), 0, 50  )[0]
    res *=  Jacobian_and_conversion
    return res
def _dE_AA_pbar_LAB_incNbarAndHyperon(Tn_proj_LAB, T_pbar_LAB, A_projectile=1, N_projectile=0, A_target=1, N_target=0, parametrization='Korsmeier_II'):
    s = 4*fMass_proton*fMass_proton + 2 * Tn_proj_LAB * fMass_proton;
    return _dE_AA_pbar_LAB(Tn_proj_LAB, T_pbar_LAB, A_projectile, N_projectile, A_target, N_target, parametrization)*( 2 + 2*cpp_xs.deltaHyperon(s, NbarAndHyperon_C[parametrization]) + cpp_xs.deltaIsospin(s,NbarAndHyperon_C[parametrization]) )

_dE_AA_pbar_LAB_incNbarAndHyperon_v = np.vectorize(_dE_AA_pbar_LAB_incNbarAndHyperon)
_dE_AA_pbar_LAB_v = np.vectorize(_dE_AA_pbar_LAB)


def dE_AA_pbar_LAB(Tn_proj_LAB, T_pbar_LAB, A_projectile=1, N_projectile=0, A_target=1, N_target=0, parametrization='KORSMEIER_II'):
    """
        Energy-differential cross section for general projectile and target nucleus for different XS parametrization
        as function of LAB frame kinetic variables.
        This cross section is integrated over all angles.
        
        Usage: The function is vectorized in both energies. dE_AA_pbar_LAB( [array], [arry] ) returns the correct corrisponding matrix.
        
        \param double Tn_proj_LAB      Kinetic energy per nucleus of the prjectile (in the LAB frame)
        \param doulbe T_pbar_LAB       Kinetic energy of the antiproton (in the LAB frame)
        \param int    A_projectile     Mass number of the projectile
        \param int    N_projectile     Number of neutrons in the projectile
        \param int    A_target         Mass number of the target
        \param int    N_target         Number of neutrons in the target
        \param string parametrization  Cross section parametrization [KORSMEIER_II (default), KORSMEIER_I, WINKLER, DI_MAURO_I, DI_MAURO_II]
        
        \return double XS              Cross section in mbarn/GeV
        """
    Tn_proj_LAB_v, T_pbar_LAB_v = np.meshgrid(Tn_proj_LAB, T_pbar_LAB, indexing='ij')
    ret = _dE_AA_pbar_LAB_v(Tn_proj_LAB_v, T_pbar_LAB_v, A_projectile, N_projectile, A_target, N_target, parametrization)
    if len(ret[:,0])==1:
        ret=ret[0,:]
        if len(ret)==1:
            return ret[0]
        return ret
    if len(ret[0,:])==1:
        return ret[:,0]
    return ret


def dE_AA_pbar_LAB_incNbarAndHyperon(Tn_proj_LAB, T_pbar_LAB, A_projectile=1, N_projectile=0, A_target=1, N_target=0, parametrization='KORSMEIER_II'):
    """
        Energy-differential cross section including antineutrons and antihyperons for general projectile and target nucleus
        and for different XS parametrization as function of LAB frame kinetic variables.
        This cross section is integrated over all angles.
        
        This function should only be used for the parametrizations: Winkler, Korsmeier_I and Korsmeier_II.
        For diMauro apply a global factor 2.3 instead.
        
        Usage: The function is vectorized in both energies. dE_AA_pbar_LAB( [array], [arry] ) returns the correct corrisponding matrix.
        
        \param double Tn_proj_LAB      Kinetic energy per nucleus of the prjectile (in the LAB frame)
        \param doulbe T_pbar_LAB       Kinetic energy of the antiproton (in the LAB frame)
        \param int    A_projectile     Mass number of the projectile
        \param int    N_projectile     Number of neutrons in the projectile
        \param int    A_target         Mass number of the target
        \param int    N_target         Number of neutrons in the target
        \param string parametrization  Cross section parametrization [KORSMEIER_II (default), KORSMEIER_I, WINKLER, DI_MAURO_I, DI_MAURO_II]
        
        \return double XS              Cross section in mbarn/GeV
        """
    Tn_proj_LAB_v, T_pbar_LAB_v = np.meshgrid(Tn_proj_LAB, T_pbar_LAB, indexing='ij')
    ret = _dE_AA_pbar_LAB_incNbarAndHyperon_v(Tn_proj_LAB_v, T_pbar_LAB_v, A_projectile, N_projectile, A_target, N_target, parametrization)
    if len(ret[:,0])==1:
        ret=ret[0,:]
        if len(ret)==1:
            return ret[0]
        return ret
    if len(ret[0,:])==1:
        return ret[:,0]
    return ret

