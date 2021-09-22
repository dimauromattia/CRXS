import math
import numpy as np

fMass_proton = 0.9382720813


def tot_pp__diMauro( s ):
    if s<0:
        return 0
    Zpp  = 33.44
    Y1pp = 13.53
    Y2pp = 6.38
    n1   = 0.324
    n2   = 0.324
    hbar2= 0.38937966 # GeV^2 mbarn (PDG)
    M    = 2.06
    Bpp  = 3.1415926536 * hbar2/M/M
    sM   = np.power(2*fMass_proton+M, 2);
    sigmaPP = Zpp + Bpp*np.power( np.log(s/sM), 2) + Y1pp * np.power(sM/s, n1) - Y2pp * np.power(sM/s, n2)
    return sigmaPP;



def convert_LAB_to_CM__from__s_ppbar_eta(  s, p_pbar_LAB, eta_LAB, with_pL=False  ):
    global fMass_proton
    E_p_LAB         =   (s-2.*fMass_proton*fMass_proton)/(2.*fMass_proton)
    E_pbar_LAB      =   np.sqrt(p_pbar_LAB*p_pbar_LAB+fMass_proton*fMass_proton)
    
    beta            =   np.sqrt(E_p_LAB - fMass_proton)/np.sqrt(E_p_LAB + fMass_proton)
    gamma           =   1./np.sqrt(1 - beta*beta);
    gammabeta       =   gamma * beta
    
    E_pbar          =   gamma * E_pbar_LAB - gammabeta * p_pbar_LAB*np.tanh(eta_LAB)
    pL_pbar         = - gammabeta * E_pbar_LAB + gamma * p_pbar_LAB*np.tanh(eta_LAB)
    
    E_pbar_Max      =   ( s-8.*fMass_proton*fMass_proton )/2./np.sqrt( s )
    
    pT_pbar         =   p_pbar_LAB/np.cosh(eta_LAB)
    x_R             =   E_pbar/E_pbar_Max
    
    if with_pL:
        return s, pT_pbar, x_R, pL_pbar
    return s, pT_pbar, x_R

def convert_LAB_to_CM__from__s_ppbar_theta(  s, p_pbar_LAB, theta_LAB, with_pL=False  ):
    global fMass_proton
    E_p_LAB         =   (s-2.*fMass_proton*fMass_proton)/(2.*fMass_proton)
    E_pbar_LAB      =   np.sqrt(p_pbar_LAB*p_pbar_LAB+fMass_proton*fMass_proton)
    
    beta            =   np.sqrt(E_p_LAB - fMass_proton)/np.sqrt(E_p_LAB + fMass_proton)
    gamma           =   1./np.sqrt(1 - beta*beta);
    gammabeta       =   gamma * beta
    
    E_pbar          =   gamma * E_pbar_LAB - gammabeta * p_pbar_LAB*np.cos(theta_LAB)
    pL_pbar         = - gammabeta * E_pbar_LAB + gamma * p_pbar_LAB*np.cos(theta_LAB)

    E_pbar_Max      =   ( s-8.*fMass_proton*fMass_proton )/2./np.sqrt( s )
    
    pT_pbar         =   p_pbar_LAB*np.sin(theta_LAB)
    x_R             =   E_pbar/E_pbar_Max
    
    if with_pL:
        return s, pT_pbar, x_R, pL_pbar
    return s, pT_pbar, x_R


def convert_LAB_to_CM__from__Tp_Tpbar_eta(  T_p_LAB, T_par_LAB, eta_LAB, with_pL=False  ):
    global fMass_proton
    s               =   4*fMass_proton*fMass_proton + 2 * T_p_LAB * fMass_proton
    p_pbar_LAB      =   np.sqrt(  T_pbar_LAB*(T_pbar_LAB+2*fMass_proton)  )
    
    return convert_LAB_to_CM(  s, p_pbar_LAB, eta_LAB, with_pL=False  )

def theta_to_eta(theta):
    eta             = -np.log(np.tan(theta/2.))
    return eta

def eta_to_theta(eta):
    theta           = 2.*np.arctan(np.exp(-eta))
    return theta

def eta_to_cosTetha(eta):
    return np.arctanh(eta)

def eta_to_sinTetha(eta):
    return 1./np.arccosh(eta)

def y_to_eta(y, pT):
    global fMass_proton
    n = np.sqrt(   (pT*pT+fMass_proton*fMass_proton)*np.cosh(y)*np.cosh(y) - fMass_proton*fMass_proton) + np.sqrt(pT*pT+fMass_proton*fMass_proton) * np.sinh(y)
    z = pT
    eta = np.log(n/z)
    return eta


def eta_to_y(eta, pT):
    global fMass_proton
    n = np.sqrt( fMass_proton*fMass_proton + (pT*pT)*np.cosh(eta)*np.cosh(eta) ) + pT * np.sinh(eta)
    z = np.sqrt( fMass_proton*fMass_proton+pT*pT )
    y = np.log(n/z)
    return y

def x_R__from_s_ppbar(s, p_pbar):
    global fMass_proton
    E_pbar          =   np.sqrt(p_pbar*p_pbar+fMass_proton*fMass_proton)
    E_pbar_Max      =   ( s-8.*fMass_proton*fMass_proton )/2./np.sqrt( s )
    x_R             =   E_pbar/E_pbar_Max
    return x_R


def s__from_pLAB(pLAB):
    global fMass_proton
    E = np.sqrt(pLAB*pLAB + fMass_proton*fMass_proton)
    s = 2*fMass_proton*(E+fMass_proton)
    return s

def x_R__from_s_pT_xF(s, pT, xF):
    global fMass_proton
    pL = xF*np.sqrt(s)/2.
    E_pbar          =   np.sqrt(pT*pT+pL*pL+fMass_proton*fMass_proton)
    E_pbar_Max      =   ( s-8.*fMass_proton*fMass_proton )/2./np.sqrt( s )
    x_R             =   E_pbar/E_pbar_Max
    return x_R



