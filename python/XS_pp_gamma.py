#! /usr/bin/env python3

import sys
import math
import numpy as np

#
#   Taken from: arXiv 1406.7369v2
#       Parametrization of gamma-ray production cross-secteions for pp interactions in a broad proton energy range from the kinematic threshold to PeV energes
#       E. Kafexhiu, et al.
#

#   Units:
#
#       T_p__LAB:       GeV
#       E_gamma__LAB:   GeV
#       dE_sigma        mbar
#
#


m_p   = 0.938272       # GeV
m_pi  = 0.134976       # GeV

T_p__LAB_th = 2 * m_pi + m_pi**2/2/m_p

'''
    Eq. (8)
    '''
def dE_sigma__pp_gamma__LAB(T_p__LAB, E_gamma__LAB):
    if T_p__LAB<T_p__LAB_th:
        return 0.
    return A_max(T_p__LAB) * F(T_p__LAB, E_gamma__LAB)


'''
    Eq. (8)
    '''
def A_max(T_p__LAB):
    global T_p__LAB_th

    b_0 = 5.9                                                                                   # Eq. 14
    b_1 = 9.53                 # Tab. VII
    b_2 = 0.52                 # Tab. VII
    b_3 = 0.054                # Tab. VII
    if T_p__LAB >= 5 :
        b_1 = 9.13             # Tab. VII
        b_2 = 0.35             # Tab. VII
        b_3 = 9.7e-3           # Tab. VII
    
    s = 2 * m_p * ( T_p__LAB + 2*m_p )
    
    theta_p = T_p__LAB/m_p
    
    E_pi_max__CM = (s - 4 * m_p**2 + m_pi**2) / 2 / np.sqrt(s)                                  # Eq. 10
    gamma_CM = (T_p__LAB + 2 * m_p)/np.sqrt(s)                                                  # Eq. 10
    beta_CM  = np.sqrt(1-gamma_CM**-2)                                                          # Eq. 10
    E_pi_max__LAB = gamma_CM * ( E_pi_max__CM + np.sqrt(E_pi_max__CM**2-m_pi**2)*beta_CM )      # Eq. 10
    
    if T_p__LAB < T_p__LAB_th:
        return 0.
    if T_p__LAB < 1:
        return b_0 * sigma_pi(T_p__LAB)/E_pi_max__LAB
    else:
        return b_1 * theta_p**-b_2 * np.exp(b_3 * (np.log(theta_p))**2 ) * sigma_pi(T_p__LAB)/m_p

'''
    Eq. (7.5), Sec. II.4
    '''
def sigma_pi(T_p__LAB):
    global T_p__LAB_th
    if T_p__LAB < T_p__LAB_th:
        return 0.
    if T_p__LAB < 2.:
        return sigma_1pi(T_p__LAB) + sigma_2pi(T_p__LAB)
    else:
        return sigma_in(T_p__LAB) * n_pi(T_p__LAB)

'''
    Eq. (2)
    '''
def sigma_1pi(T_p__LAB):
    global T_p__LAB_th
    if T_p__LAB<T_p__LAB_th:
        return 0.
    sigma_0 = 7.66e-3
    s = 2 * m_p * ( T_p__LAB + 2*m_p )
    eta = np.sqrt((s-m_pi**2-4*m_p**2)**2-16*m_pi**2*m_p**2)/(2*m_pi*np.sqrt(s))                       # Eq. 3
    return sigma_0 * eta**1.95*(1+eta+eta**5) * (f_BW(np.sqrt(s)))**1.86

'''
    Eq. (4)
    '''
def f_BW(sqrtS):
    
    Gamma_res = 0.2264      # GeV
    M_res     = 1.1883      # GeV
    
    gamma = np.sqrt( M_res**2*(M_res**2 + Gamma_res**2) )
    K     = np.sqrt(8.) * M_res * Gamma_res * gamma / math.pi / np.sqrt(M_res**2 + gamma)
    
    return m_p * K / ( ((sqrtS-m_p)**2 - M_res**2)**2 + M_res**2 * Gamma_res**2 )

'''
    Eq. (5)
    '''
def sigma_2pi(T_p__LAB):
    if T_p__LAB<T_p__LAB_th:
        return 0.
    s0 = 5.7   # mbarn
    if T_p__LAB < 0.56:
        return 0.
    else:
        return s0/(1.+ np.exp(-9.3*(T_p__LAB-1.4)))


def E_pi_max__LAB(T_p__LAB):
    s = 2 * m_p * ( T_p__LAB + 2*m_p )
    E_pi_max__CM = (s - 4 * m_p**2 + m_pi**2) / 2 / np.sqrt(s)                                  # Eq. 10
    gamma_CM = (T_p__LAB + 2 * m_p)/np.sqrt(s)                                                  # Eq. 10
    beta_CM  = np.sqrt(1-gamma_CM**-2)                                                          # Eq. 10
    E_pi_max__LAB = gamma_CM * ( E_pi_max__CM + np.sqrt(E_pi_max__CM**2-m_pi**2)*beta_CM )      # Eq. 10
    return E_pi_max__LAB


'''
    Eq. (1)
    '''
def sigma_in(T_p__LAB):
    return ( 30.7 - 0.96 * np.log(T_p__LAB/T_p__LAB_th) + 0.18 * (np.log(T_p__LAB/T_p__LAB_th))**2 ) * ( 1- (T_p__LAB_th/T_p__LAB)**1.9)**3

'''
    Eq. (7)
    '''
def n_pi(T_p__LAB):
    if T_p__LAB < 5 :
        Q_p = (T_p__LAB-T_p__LAB_th)/m_p
        return -6e-3 + 0.237*Q_p -0.023*Q_p**2
    else:
        a_1 = 0.728                 # Tab. IV
        a_2 = 0.596                 # Tab. IV
        a_3 = 0.491                 # Tab. IV
        a_4 = 0.2503                # Tab. IV
        a_5 = 0.117                 # Tab. IV
        eps_p = (T_p__LAB-3)/m_p
        return a_1 * eps_p**a_4 * (1+np.exp(-a_2*eps_p**a_5)) * (1-np.exp(-a_3*eps_p**0.25))

'''
    Eq. (11)
    '''
def F(T_p__LAB, E_gamma__LAB):
    
    if T_p__LAB < T_p__LAB_th:
        return 0.
    
    lamb  = 3.0
    alpha = 1.0
    beta  = kappa(T_p__LAB)
    gamma = 0.
    if T_p__LAB > 1.:
        lamb  = 3.0
        alpha = 1.0
        beta  = mu(T_p__LAB) + 2.45
        gamma = mu(T_p__LAB) + 1.45
    if T_p__LAB > 4.:
        lamb  = 3.0
        alpha = 1.0
        beta  = 1.5*mu(T_p__LAB) + 4.95
        gamma =     mu(T_p__LAB) + 1.50
    if T_p__LAB > 20.:
        lamb  = 3.0
        alpha = 0.5
        beta  = 4.2
        gamma = 1.0
    if T_p__LAB > 100.:
        lamb  = 3.0
        alpha = 0.5
        beta  = 4.9
        gamma = 1.0
    
    s = 2 * m_p * ( T_p__LAB + 2*m_p )
    E_pi_max__CM = (s - 4 * m_p**2 + m_pi**2) / 2 / np.sqrt(s)                                  # Eq. 10
    if E_pi_max__CM <= m_pi:
        return 0.
    gamma_CM = (T_p__LAB + 2 * m_p)/np.sqrt(s)                                                  # Eq. 10
    beta_CM  = np.sqrt(1-gamma_CM**-2)                                                          # Eq. 10
    E_pi_max__LAB = gamma_CM * ( E_pi_max__CM + np.sqrt(E_pi_max__CM**2-m_pi**2)*beta_CM )      # Eq. 10
    gamma_pi__LAB = E_pi_max__LAB /m_pi
    beta_pi__LAB  = np.sqrt(1-gamma_pi__LAB**-2)
    #E_gamma_min = m_pi/2. * gamma_pi__LAB * ( 1 - beta_pi__LAB )
    E_gamma_max = m_pi/2. * gamma_pi__LAB * ( 1 + beta_pi__LAB )

    if E_gamma__LAB > E_gamma_max:
        return 0.
   
    
    E_gamma = E_gamma__LAB
    Y_gamma     = E_gamma     + m_pi**2/4./E_gamma
    Y_gamma_max = E_gamma_max + m_pi**2/4./E_gamma_max
    X_gamma     = ( Y_gamma - m_pi ) / ( Y_gamma_max - m_pi)
    
    C = lamb * m_pi/Y_gamma_max

    N  = (  1-X_gamma**alpha )**beta
    DN = (  1+X_gamma/C      )**gamma

    return N/DN

def kappa(T_p__LAB):
    theta_p = T_p__LAB/m_p
    return 3.29 - 0.2 * theta_p**-1.5

def mu(T_p__LAB):
    q = (T_p__LAB - 1)/m_p
    return 5./4 * q**(5./4) * np.exp( -5./4 * q )


####
####       TESTS
####

#
#import matplotlib                   as mpl
#mpl.use('Agg')
#import matplotlib.pyplot            as plt
#import MKastro.basic.PlotFunctions  as plf
#import matplotlib.colors            as colors
#
#sigma_pi_vec = np.vectorize(sigma_pi)
#sigma_1pi_vec = np.vectorize(sigma_1pi)
#sigma_2pi_vec = np.vectorize(sigma_2pi)
#
#
#
#T_p_2 = np.power( 10, np.linspace(-1,0.3, 20) )
#
#
#plot, fig = plf.new_plot(r'$T_p$ [GeV]', r'$\sigma_\pi$  [mbarn]', 'log', 'log')
#T_p   = np.power( 10, np.linspace(np.log10(T_p__LAB_th),1,  100) )
#T_p_2 = np.power( 10, np.linspace(np.log10(T_p__LAB_th),0.3, 20) )
#print(T_p_2)
#xs_pi = sigma_pi_vec(T_p )
#xs_1pi = sigma_1pi_vec( T_p_2 )
#xs_2pi = sigma_2pi_vec( T_p_2 )
#
#plot.plot(T_p,   xs_pi,    label=r'$\sigma_\pi$'   )
#plot.plot(T_p_2, xs_1pi,   label=r'$\sigma_{1\pi}$')
#plot.plot(T_p_2, xs_2pi,   label=r'$\sigma_{2\pi}$')
#
#plot.set_ylim( (0.01, 40) )
#
##plot.plot(T_p_2, 5.7/(1.+ np.exp(-9.3*(T_p_2-1.4))),   label=r'$\sigma_{2\pi}$ 2')
#plt.legend()
#plt.savefig( 'sigma_pi.png' )
#
#
#
#
#dE_sigma__pp_gamma__LAB_vec = np.vectorize(dE_sigma__pp_gamma__LAB)
#
#
#plot, fig = plf.new_plot(r'$T_p$ [GeV]', r'$E_\gamma$ [GeV]', 'log', 'log', label_size=15)
#plt.subplots_adjust(left=0.15, right=0.75, top=0.9, bottom=0.15)
#
#T_p   = np.power( 10, np.linspace(-1,5,  300) )
#E_g   = np.power( 10, np.linspace(-1,5,  300) )
#
#x, y = np.meshgrid( T_p, E_g, indexing='ij' )
#
#z = dE_sigma__pp_gamma__LAB_vec(x, y)
#
#print z
#
#
#cmap_str = 'magma_r'
#cmap    = plt.get_cmap(cmap_str)
#
#norm = colors.LogNorm(vmin=z.max()/1e3, vmax=z.max())
##norm = colors.Normalize(vmin=z.min(), vmax=z.max())
#
#colormesh = plot.pcolormesh(  x, y, z, norm=norm, cmap=cmap  )
#
#cbar_ax = fig.add_axes([0.8, 0.15, 0.02, 0.75])
#cbar = fig.colorbar(colormesh, cax=cbar_ax, cmap=cmap)
#cbar.ax.set_ylabel(r'$d\sigma/dE$ [mbar/GeV]' )
#
##plot.set_xticks(np.arange(1,n_bins+1))
##plot.set_yticks(np.arange(1,n_bins+1))
#
#plt.savefig('Michael.png' )
#
