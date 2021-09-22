#! /usr/bin/env python

import glob
import os
import math
import sys

import numpy as np
import scipy as sp
import scipy.special as sps
from   scipy.interpolate import interp1d

from   scipy.interpolate import interp2d

import matplotlib.pyplot    as plt
import matplotlib           as mpl
import matplotlib.colors    as colors
import matplotlib.patches   as mpatches
from   matplotlib.path  import Path
from   matplotlib       import rc






import argparse
parser = argparse.ArgumentParser(description='Propagation in analytic model.')
parser.add_argument('--species', help='Dbar or Hebar', action='store', dest='species', type=str, default='Dbar')
parser.add_argument('--step',        help='step',                           action='store', dest='step'   ,     type=int, default=1       )
args = parser.parse_args()
species=args.species
step   =args.step

suffix = ''

A_sp            = 2.
if species=='pbar':
    A_sp        = 1.
    Z_sp        = 1.
if species=='Dbar':
    A_sp        = 2.
    Z_sp        = 1.
if species=='Hebar':
    A_sp        = 3.
    Z_sp        = 2.

CRACS = os.environ['CRACS']

data                       = np.genfromtxt( species+'_SourceTerm.txt',           skip_header=1 )

data_cirelli_pbar          = np.genfromtxt( CRACS+'/data/project/AntiMatter/AtProduction_antiprotons.dat',      skip_header=1 )
#f_data_cirelli_pbar__bb    = interp2d(  data_cirelli_pbar[:,0], data_cirelli_pbar[:,1], data_cirelli_pbar[:,13], kind='cubic')


# XS:
data_dpbar_nar             = np.genfromtxt( CRACS+'/data/project/AntiMatter/table_dpbar_nar.txt',               skip_header=1 )
data_dpbar_tot             = np.genfromtxt( CRACS+'/data/project/AntiMatter/table_dpbar_tot.txt',               skip_header=1 )
data_ppbar_el              = np.genfromtxt( CRACS+'/data/project/AntiMatter/table_ppbar_el.txt',                skip_header=1 )
data_ppbar_tot             = np.genfromtxt( CRACS+'/data/project/AntiMatter/table_ppbar_tot.txt',               skip_header=1 )

data_Anderson              = np.genfromtxt( CRACS+'/data/CS_lab/dT_pp_p_LAB_Anderson.txt'                         )




c0='black'
c2='#04B404'
c4='#B40486'
c3='#0489B1'
c1='#0404B4'
c5='#B40404'

ff = 13

def getSpectrum(mDM):
    
    for i in range(1000):
        if data_cirelli_pbar[i*179,0]>mDM:
            break
    mDM_u = data_cirelli_pbar[ i   *179,0]
    mDM_l = data_cirelli_pbar[(i-1)*179,0]

    x       = data_cirelli_pbar[ i   *179:(i+1)*179, 1]
    dNdx_l  = data_cirelli_pbar[ i   *179:(i+1)*179, ff] + np.ones(179)*1e-90
    dNdx_u  = data_cirelli_pbar[(i+1)*179:(i+2)*179, ff] + np.ones(179)*1e-90

    log_dNdx_m  = np.log(dNdx_l)    +    (  np.log(dNdx_u)-np.log(dNdx_l) ) / (np.log(mDM_u)-np.log(mDM_l)) * (np.log(mDM)-np.log(mDM_l))

    return np.power(10,x)*mDM, np.exp(log_dNdx_m)/(  np.power(10,x)*mDM  )/np.log(10)


print_size=10
label_size=20

def plot_1D( xlabel, ylabel, xscale='log', yscale='log', sizex=1.3, sizey=1.0):
    
    global print_size, label_size
    plt.close('all')
    font_props = {"size":label_size}
    rc("font", **font_props)
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['mathtext.fontset']='stixsans'
    
    fig     = plt.figure(figsize=(print_size*sizex, print_size*sizey))
    plot    = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
    plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
    
    plot.set_xlabel ( xlabel, fontsize=label_size*1.4 )
    plot.set_ylabel ( ylabel, fontsize=label_size*1.4 )
    
    plot.set_xscale ( xscale )
    plot.set_yscale ( yscale )
    
    plot.tick_params('both', length=20, width=2, which='major', top=True, right=True, direction='in', pad=10 )
    plot.tick_params('both', length=10, width=1, which='minor', top=True, right=True, direction='in', pad=10 )
    
    plot.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
    plot.tick_params(axis='both', pad=10)
    
    return plot, fig



p_coal = 0.062

######################################
#   propagation parameters
######################################

h       =   0.1         # kpc
K_0     =   0.0112      # kpc^2/Myr          MED
V_c     =   12.         # km/s               MED
delta   =   0.7         #                    MED
L       =   4.          # kpc                MED

# transform
K_0     =   K_0 / (1e6 * 3.154e+7)           # in kpc^2/s
V_c     =   V_c * 3.24078e-17                # transform to kpc/s

######################################
#   DM parameters
######################################

sv              = 2.69e-26                  # cm^3/s        sigma v  (thermally averaged cross section
mDM             = 70.8                      # GeV           DM mass
rhoSun          = 0.43                      # GeV/cm^3      Local DM density
rSun            = 8.                        # kpc


######################################
#   grid parameters
######################################

fMass_He3   = 2.81          # GeV
fMass_d     = 1.875612928   # GeV
fMass_p     = 0.9382720813  # GeV
fMass_n     = 0.939565379   # GeV

R       = 20.           # kpc
Tn_min  = 1e-1          # GeV/n
Tn_max  = 1e3           # GeV/n

nE      =  40
nZ      =  60
nR      = 200
nI      =  80

nZ      =  20
nE      =  40
nR      =  80
nI      = 100

######################################
#   create grid
######################################
dR  = 1.*R/(nR-1)
dEf = 1.*np.log(1.*Tn_max/Tn_min)/(nE-1)

i_vec   = np.arange(0,nI)
r_vec   = np.arange(0, R+dR/2., dR)
Tn_vec  = np.arange( np.log(Tn_min), np.log(Tn_max)+ dEf/2., dEf)
Tn_vec  = np.exp(Tn_vec)


xi_i = sps.jn_zeros(0, nI+1)
    
def xi(i):
    global xi_i
    return xi_i[i]
vect_xi = np.vectorize(xi)

def delta_z(z):
    global dZ, h
    if z<1e-10:
        return 1.*h/dZ
    return 0
vect_delta_z = np.vectorize(delta_z)


def proper(case='DM', tertiary_loss=False, threshold_correction=False):

    global suffix
    global h, K_0, V_c, delta, L, sv, mDM, rhoSun
    global Z,  R, Tn_min, Tn_max, nE, nZ, nR, nI, dR, dZ, dEf, i_vec, r_vec, z_vec, Tn_vec, i_grid, Tn_grid, z_grid, r_grid
    global multiplicator_pp_XS
    
    global data_DM_energySpectrum, data, data_tert
    
    Z       = L             # kpc
    dZ  = 1.*Z/(nZ-1)
    z_vec   = np.arange(0, Z+dZ/2., dZ)
    i_grid, Tn_grid, z_grid, r_grid = np.meshgrid(i_vec, Tn_vec, z_vec, r_vec, indexing='ij')

    
    print 'Propergate: '+suffix
    #if mDM!=100:
    #    'WARNING:  adjust the file with dN/dT data to your DM mass!'

    num             = 1e-50
    rr_grid         = np.sqrt(r_grid*r_grid + z_grid*z_grid)+num
    
    ######################################
    #   create q_D (E, r, z)
    ######################################

    Rs              = 20.
    rho_0           = 1./ (rSun/Rs * np.power(1 + rSun/Rs, 2))
    rho_grid        = 1./ (rr_grid       /Rs * np.power(1 + rr_grid  /Rs, 2)) * rhoSun / rho_0

    if species=='pbar':
        p_vec       = np.sqrt(Tn_vec*Tn_vec+2*fMass_p*Tn_vec)
    if species=='Dbar':
        p_vec       = np.sqrt(A_sp*A_sp*Tn_vec*Tn_vec+2*fMass_d*A_sp*Tn_vec)
    if species=='Hebar':
        p_vec       = np.sqrt(A_sp*A_sp*Tn_vec*Tn_vec+2*fMass_He3*A_sp*Tn_vec)
    T_pbar_vec  = np.sqrt(p_vec*p_vec/A_sp/A_sp + fMass_p*fMass_p)-fMass_p
#    spec        =   f_data_cirelli_pbar__bb( mDM, np.log10(T_pbar_vec/mDM) )[:,0]
#    spec       *=   1./Tn_vec/np.log(10)
    T_cirelli, dNdT_cirelli = getSpectrum(mDM)
    spec        = np.interp(T_pbar_vec, T_cirelli, dNdT_cirelli)
    

    if species=='pbar':
        pref        =   1.
    if species=='Dbar':
            pref        =   4./3.*     pow(p_coal, 3)/p_vec        *    fMass_d/fMass_p/fMass_n
    if species=='Hebar':
        pref        =   3.   * pow(pow(p_coal, 3)/p_vec,2)     *    fMass_He3/fMass_p/fMass_n/fMass_n

    if threshold_correction and species=='Dbar':
        T_cirelli, dNdT_cirelli = getSpectrum(mDM)
        spec1        = np.interp(T_pbar_vec, T_cirelli, dNdT_cirelli)
        spec2 = []
        for i in range(len(T_pbar_vec)):
            T_cirelli, dNdT_cirelli = getSpectrum(mDM-2*T_pbar_vec[i]-2*fMass_p)
            spec2.append ( np.interp(T_pbar_vec[i], T_cirelli, dNdT_cirelli) )
        spec2 = np.array(spec2)
        dN_dTn_grid = np.interp( Tn_grid, Tn_vec, A_sp * pref * spec1*spec2  )
    else:
        dN_dTn_grid = np.interp( Tn_grid, Tn_vec, A_sp * pref * np.power(spec,A_sp)  )

    print 'Source at 600 MeV/n: ' +str(np.interp( 0.6, Tn_vec, A_sp * pref * np.power(spec,A_sp)  )*np.power( 0.080/0.062, 3*(A_sp-1)) )


    #dN_dTn_grid     = np.interp( Tn_grid, data_DM_energySpectrum[:,0], data_DM_energySpectrum[:,1]  )

    rho_q_D_grid  =  rho_grid*rho_grid/rhoSun/rhoSun


    ######################################
    #   create q_i (E, z)
    ######################################

    i_grid2     = i_grid [:,:,:,0]
    Tn_grid2    = Tn_grid[:,:,:,0]
    z_grid2     = z_grid [:,:,:,0]

    if case =='DM':
        integrand_q_i_grid  = r_grid  *  sps.j0( vect_xi(i_grid)*r_grid/R ) * rho_q_D_grid
    if case =='Secondary' or case =='Tertiary':
        integrand_q_i_grid  = r_grid  *  sps.j0( vect_xi(i_grid)*r_grid/R ) * vect_delta_z(z_grid)
    if case =='Secondary_Green' or case =='Tertiary_Green':
        integrand_q_i_grid  = r_grid  *  sps.j0( vect_xi(i_grid)*r_grid/R ) * vect_delta_z(z_grid)*np.power(r_grid/rSun, 1.09)*np.exp(-3.87*(r_grid-rSun)/rSun)

    pref_grid2          = 4./np.power(sps.j1( vect_xi(i_grid2) )*R, 2)
    q_i_grid2           = integrand_q_i_grid.sum(axis=3) * dR * pref_grid2


    ######################################
    #   calculate y_i (E)
    ######################################

    E_grid2     =   Tn_grid2*2 + fMass_d
    p_grid2     =   2.*np.sqrt(  Tn_grid2*(Tn_grid2+fMass_d)  )
    beta_grid2  =   p_grid2/E_grid2
    rig_grid2   =   p_grid2/Z_sp
    K_grid2     =   K_0 * beta_grid2 * np.power(rig_grid2, delta)

    S_i_grid2   =   np.sqrt(  np.power(  V_c/K_grid2, 2 ) + 4*np.power(vect_xi(i_grid2)/R, 2)  )            #   in 1/kpc

    nH  = 1.    *   1e6
    nHe = 0.1   *   1e6
    
    data_dpbar_in = (data_ppbar_tot-data_ppbar_el)*data_dpbar_tot/data_ppbar_tot  *  A_sp/2.
    if tertiary_loss:
        data_dpbar_in = ( (data_ppbar_tot[:,:2]-data_ppbar_el[:,:2])*data_dpbar_tot[:,:2]/data_ppbar_tot[:,:2] + data_dpbar_nar[:,:2])  *  A_sp/2.
    sigma_in_dbarp__grid2     = np.interp( Tn_grid2, data_dpbar_tot[:,0], data_dpbar_in[:,1]  ) * 1e-31 * np.power(A_sp/2., 2/3)
    Gamma_grid2               = (nH + np.power(4,2./3)*nHe)*sigma_in_dbarp__grid2*beta_grid2*2.99792e8
    
    A_i_grid2   =   V_c + 2*h*Gamma_grid2 + K_grid2*S_i_grid2  *  1./np.tanh(S_i_grid2*L/2.)

    integrand__y_i_grid2    = np.exp( V_c/2./K_grid2 * (L-z_grid2) )  * np.sinh(  S_i_grid2 * (L-z_grid2) /2.  ) * q_i_grid2

    y_i_grid3 = integrand__y_i_grid2.sum(axis=2) * dZ

    ######################################
    #   calculate R_D
    ######################################

    i_grid3     = i_grid [:,:,0,0]
    Tn_grid3    = Tn_grid[:,:,0,0]
    z_grid3     = z_grid [:,:,0,0]

    K_grid3     = K_grid2   [:,:,0]
    A_i_grid3   = A_i_grid2 [:,:,0]
    S_i_grid3   = S_i_grid2 [:,:,0]


    summand_R_D = sps.j0(vect_xi(i_grid3)*rSun/R) * np.exp(-V_c * L / 2. / K_grid3)* y_i_grid3/A_i_grid3/np.sinh( S_i_grid3 * L / 2. )
    R_D = summand_R_D.sum(axis=0)


    if case =='DM':
        n_D = 0.5* sv *rhoSun*rhoSun/mDM/mDM * dN_dTn_grid[0,:,0,0] * R_D
        n_D = 1e6*n_D                                                           # Transform all units to GeV, m, and s
    if case =='Secondary' or case=='Secondary_Green':
        n_D = np.interp( Tn_vec, data[:,0], data[:,1]  ) * R_D
    if case =='Tertiary'  or case=='Tertiary_Green':
        n_D = np.interp( Tn_vec, data_tert[:,0], data_tert[:,1]  ) * R_D


    ######################################
    #   calculate flux
    ######################################

    phi_D = n_D * beta_grid2[0,:,0] / 4. / math.pi * 2.99792e8

    phi=0.4
    if species=='Hebar':
        phi_D      *=   2    # Tritium + Helium3
        phi=0.6
    if species=='pbar':
        phi=0.8

    Tn_vec_mod = Tn_vec - phi * np.fabs(Z_sp) / A_sp
    phi_D_mod =  ( Tn_vec_mod*(Tn_vec_mod+2*fMass_p)  )/(  Tn_vec*(Tn_vec+2*fMass_p)  ) * phi_D

    ############################
    #
    #  Change p_coal from 124 to 160!!
    #
    ############################

    if species!='pbar':
        if case !='Tertiary':
            print np.power( 0.080/0.062, 3*(A_sp-1) )
            phi_D     *= np.power( 0.080/0.062, 3*(A_sp-1) )
            phi_D_mod *= np.power( 0.080/0.062, 3*(A_sp-1) )

    print 'Flux at 600 MeV/n: ' + str(np.interp(0.6, Tn_vec, phi_D))


    return (Tn_vec, phi_D), (Tn_vec_mod, phi_D_mod)

#################################
##         Tertiaries          ##
#################################
def heaviside(x):
    if x>0:
        return 1.
    else:
        return 0.
v_heaviside = np.vectorize(heaviside)


def getTertiaries( sec_Tn, sec_Sec):
    global Tn_min, Tn_max, Tn_vec
    nE_ter      = 10000
    dEf_ter     = 1.*np.log(1.*Tn_max/Tn_min)/(nE_ter-1)
    log_Tn_vec  = np.arange( np.log(Tn_min), np.log(Tn_max)+ dEf_ter/2., dEf_ter )

    Tn_vec_grid_tert, log_Tn_vec_grid_tert = np.meshgrid(Tn_vec, log_Tn_vec, indexing='ij')

    nH  = 1.    *   1e6
    nHe = 0.1   *   1e6

    dsigma_by_dT__grid_tert = np.power(A_sp/2., 2/3)*v_heaviside( np.exp(log_Tn_vec_grid_tert) - Tn_vec_grid_tert )/np.exp(log_Tn_vec_grid_tert)  *  np.interp( np.exp(log_Tn_vec_grid_tert), data_dpbar_nar[:,0], data_dpbar_nar[:,1]  *  A_sp/2.  ) * 1e-31
    integrand__grid_tert = dsigma_by_dT__grid_tert * np.exp(log_Tn_vec_grid_tert) * np.interp( np.exp(log_Tn_vec_grid_tert),  sec_Tn, sec_Sec )

    q_tert = integrand__grid_tert.sum(axis=1) * dEf_ter * 4 * 3.142 * (nH + np.power(4,2./3)*nHe)
    return q_tert

data_Anderson_Tprod = data_Anderson[0  , 1: ]
data_Anderson_Tproj = data_Anderson[1: , 0  ]
data_Anderson_XS    = data_Anderson[1: , 1: ]

dLog = np.log(data_Anderson_Tprod[1]/data_Anderson_Tprod[0])

v_Anderson_Tproj, v_Anderson_Tprod  = np.meshgrid(data_Anderson_Tproj, data_Anderson_Tprod, indexing='ij')

int_Tprod = ( (data_Anderson_XS*v_Anderson_Tprod).sum(axis=1) )*dLog

factor = np.interp( data_Anderson_Tproj, data_dpbar_nar[:,0], data_dpbar_nar[:,1]  )*1e-31/int_Tprod

data_Anderson_XS = data_Anderson_XS * np.interp( v_Anderson_Tproj, data_Anderson_Tproj, factor  )

int_anderson = interp2d( data_Anderson_Tprod, data_Anderson_Tproj, data_Anderson_XS, kind='cubic')


def getTertiaries_anderson( sec_Tn, sec_Sec):
    global Tn_min, Tn_max, Tn_vec, int_anders
    nE_ter      = 10000
    dEf_ter     = 1.*np.log(1.*Tn_max/Tn_min)/(nE_ter-1)
    log_Tn_vec  = np.arange( np.log(Tn_min), np.log(Tn_max)+ dEf_ter/2., dEf_ter )

    Tn_vec_grid_tert, log_Tn_vec_grid_tert = np.meshgrid(Tn_vec, log_Tn_vec, indexing='ij')

    nH  = 1.    *   1e6
    nHe = 0.1   *   1e6

    dsigma_by_dT__grid_tert = np.power(A_sp/2., 2/3)*np.transpose( int_anderson( Tn_vec, np.exp(log_Tn_vec)  ) )
    integrand__grid_tert = dsigma_by_dT__grid_tert * np.exp(log_Tn_vec_grid_tert) * np.interp( np.exp(log_Tn_vec_grid_tert),  sec_Tn, sec_Sec )

    q_tert = integrand__grid_tert.sum(axis=1) * dEf_ter * 4 * 3.142 * (nH + np.power(4,2./3)*nHe)
    return q_tert



#################################
##    Flux Transformations     ##
#################################


def T_to_R( T, Z, A ):
    m = A*0.938
    E = T + m
    p = np.sqrt(E*E-m*m)
    return p/np.fabs(1.*Z)

def R_to_T( R, Z, A ):
    m = A*0.938
    p = R*np.fabs(1.*Z)
    return np.sqrt(p*p+m*m)-m

def dT_by_dR( T, Z, A ):
    m = A*0.938
    E = T + m
    p = np.sqrt(E*E-m*m)
    return p*np.fabs(1.*Z)/E

def dR_by_dT( R,  Z, A):
    m = A*0.938
    p = R*np.fabs(1.*Z)
    E = np.sqrt(p*p+m*m)
    return E/(p*np.fabs(1.*Z))

def ConvertSpectrumFromTnToR( Tn, F, Z, A ):
    R = T_to_R(Tn*A, Z, A)
    Fnew = F*dT_by_dR(Tn*A,Z,A)/A;
    return R, Fnew

def ConvertSpectrumFromRToTn( R, F, Z, A ):
    Tn = R_to_T(R, Z, A)/A;
    Fnew = F*dR_by_dT(R, Z, A)*A;
    return Tn, Fnew

def SolarModulation( Tn, F, phi, Z, A):
    Tn_new    = Tn - np.fabs(1.*Z)/A*phi;
    F_new           = F*(Tn_new*Tn_new+2*Tn_new*0.938)/(Tn*Tn+2*Tn*0.938)
    return Tn_new, F_new

def SolarModulationR( R, Fr, phi, Z, A):
    Tn, F = ConvertSpectrumFromRToTn(R, Fr, Z, A)
    Tn_new    = Tn - np.fabs(1.*Z)/A*phi;
    F_new           = F*(Tn_new*Tn_new+2*Tn_new*0.938)/(Tn*Tn+2*Tn*0.938)
    R_new, Fr_new = ConvertSpectrumFromTnToR(Tn_new, F_new, Z, A)
    return R_new, Fr_new



######################################
#   PROPERGATE MED
######################################

h       =   0.1         # kpc
K_0     =   0.0112      # kpc^2/Myr          MED
V_c     =   12.         # km/s               MED
delta   =   0.7         #                    MED
L       =   4.          # kpc                MED

# transform
K_0     =   K_0 / (1e6 * 3.154e+7)           # in kpc^2/s
V_c     =   V_c * 3.24078e-17                # transform to kpc/s

suffix = '_MED'
med_is, med = proper()

suffix = '_Secondary_MED'
med_sec_is, med_sec = proper(case='Secondary')

if step==2:
    suffix = '_Secondary_MED_TertLoss'
    med_sec_is_loss, med_sec_loss = proper(case='Secondary', tertiary_loss=True)
    s_write='# T/n              sec MED              sec MED (+tert loss)'
    for i in range(len(med_sec[1])):
        s_write += '\n' + str( (med_sec[0])[i] ).rjust(20) + ' ' +str( (med_sec[1])[i] ).rjust(20) + ' ' + str( (med_sec_loss[1])[i] ).rjust(20) + ' '
    f = open( 'sec_loss.txt', 'w' )
    f.write( s_write )
    f.close()
    sys.exit(0)

if step==3:
    suffix = '_Secondary_MED_threshold_correction'
    med_is_thresh, med_thresh = proper(case='DM', threshold_correction=True)
    s_write='# T/n              DM MED              DM MED (+thresh corr)'
    for i in range(len(med_sec[1])):
        s_write += '\n' + str( (med[0])[i] ).rjust(20) + ' ' +str( (med[1])[i] ).rjust(20) + ' ' + str( (med_thresh[1])[i] ).rjust(20) + ' '
    f = open( 'DM_thres.txt', 'w' )
    f.write( s_write )
    f.close()
    sys.exit(0)


q_tert = getTertiaries_anderson(  med_sec_is[0], med_sec_is[1]  )
data_tert = np.transpose(np.array([Tn_vec, q_tert]))
med_tert_is, med_tert = proper(case='Tertiary')

for i in range(len(med[1])):
    if (med[1])[i]<=0:
        (med[1])[i]=1e-90
for i in range(len(med_sec[1])):
    if (med_sec[1])[i]<=0:
        (med_sec[1])[i]=1e-90
for i in range(len(med_tert[1])):
    if (med_tert[1])[i]<=0:
        (med_tert[1])[i]=1e-90


######################################
#   POPERGATE MAX
######################################

h       =   0.1         # kpc
K_0     =   0.0765      # kpc^2/Myr          MAX
V_c     =   5.          # km/s               MAX
delta   =   0.46        #                    MAX
L       =  15.          # kpc                MAX


# transform
K_0     =   K_0 / (1e6 * 3.154e+7)           # in kpc^2/s
V_c     =   V_c * 3.24078e-17                # transform to kpc/s

suffix = '_MAX'
max_is, max = proper()

suffix = '_Secondary_MAX'
max_sec_is, max_sec = proper(case='Secondary')

q_tert = getTertiaries_anderson(  max_sec_is[0], max_sec_is[1]  )
data_tert = np.transpose(np.array([Tn_vec, q_tert]))
max_tert_is, max_tert = proper(case='Tertiary')


######################################
#   POPERGATE KoCu
######################################

#h       =   0.1         # kpc
#D_0     =   9.84e28     # cm^2/s             KoCu
#V_c     =   45.3        # km/s               KoCu
#delta   =   0.245       #                    KoCu
#L       =   5.35        # kpc                KoCu
#
## transform
#D_0     =   D_0 * np.power(  (1e-3*3.24078e-17*1e-2),2  )       # in kpc^2/s
#K_0     =   D_0 / np.power( 4., delta )                         # from 4 to 1 GV
#print K_0*(1e6 * 3.154e+7)
#V_c     =   V_c * 3.24078e-17                                   # transform to kpc/s
#
#suffix = '_KoCu'
#KoCu_is, KoCu = proper()
#
#suffix = '_Secondary_KoCu'
#KoCu_sec_is, KoCu_sec = proper(case='Secondary')
#
#
#
##suffix = '_Secondary_KoCu_Green'
##KoCu_secR_is, KoCu_secR = proper(case='Secondary_Green')
#
#q_tert = getTertiaries_anderson(  KoCu_sec_is[0], KoCu_sec_is[1]  )
#data_tert = np.transpose(np.array([Tn_vec, q_tert]))
#KoCu_tert_is, KoCu_tert = proper(case='Tertiary')
#
##q_tert = getTertiaries_anderson(  KoCu_secR_is[0], KoCu_secR_is[1]  )
##data_tert = np.transpose(np.array([Tn_vec, q_tert]))
##KoCu_tertR_is, KoCu_tertR = proper(case='Tertiary')
#
#
##q_tertDM = getTertiaries_anderson(  KoCu_is[0], KoCu_is[1]  )
##data_tert  = np.transpose(np.array([Tn_vec, q_tertDM]))
##KoCu_tertDM_is, KoCu_tertDM = proper(case='Tertiary')


cmd = 'galprop -o . -g '+CRACS+'/data/project/AntiMatter -r standard  -f /usr/local/korsmeier/galprop/data/FITS'
os.system( cmd )


#   GALLPROP

if species=='Dbar':
    os.system('PROPERFIT_plotGalpropSpectrum nuclei_56_standard --A 2 --Z -1 --p 0 --t primary   --run false')
    os.system('PROPERFIT_plotGalpropSpectrum nuclei_56_standard --A 2 --Z -1 --p 0 --t secondary --run false')
    os.system('PROPERFIT_plotGalpropSpectrum nuclei_56_standard --A 2 --Z -1 --p 0 --t dm        --run false')
    galprop_dm   = np.genfromtxt('Spectrum_ekinpern_dm_Z_-1_A_2.txt')
    galprop_sec  = np.genfromtxt('Spectrum_ekinpern_primary_Z_-1_A_2.txt')
    galprop_tert = np.genfromtxt('Spectrum_ekinpern_secondary_Z_-1_A_2.txt')
    KoCu_is         = (galprop_dm[:,0]/1000.,galprop_dm[:,1]*1000.)
    KoCu            = SolarModulation( galprop_dm[:,0]/1000.,galprop_dm[:,1]*1000., 0.4, -1, 2)
    KoCu_sec_is     = (galprop_sec[:,0]/1000.,galprop_sec[:,1]*1000.)
    KoCu_sec        = SolarModulation( galprop_sec[:,0]/1000.,galprop_sec[:,1]*1000., 0.4, -1, 2)
    KoCu_tert_is    = (galprop_tert[:,0]/1000.,galprop_tert[:,1]*1000.)
    KoCu_tert       = SolarModulation( galprop_tert[:,0]/1000.,galprop_tert[:,1]*1000., 0.4, -1, 2)

if species=='Hebar':
    os.system('PROPERFIT_plotGalpropSpectrum nuclei_56_standard --A 3 --Z -2 --p 0 --t primary   --run false')
    os.system('PROPERFIT_plotGalpropSpectrum nuclei_56_standard --A 3 --Z -2 --p 0 --t secondary --run false')
    os.system('PROPERFIT_plotGalpropSpectrum nuclei_56_standard --A 3 --Z -2 --p 0 --t dm        --run false')
    galprop_dm   = np.genfromtxt('Spectrum_ekinpern_dm_Z_-2_A_3.txt')
    galprop_sec  = np.genfromtxt('Spectrum_ekinpern_primary_Z_-2_A_3.txt')
    galprop_tert = np.genfromtxt('Spectrum_ekinpern_secondary_Z_-2_A_3.txt')
    KoCu_is         = (galprop_dm[:,0]/1000.,galprop_dm[:,1]*1000.)
    KoCu            = SolarModulation( galprop_dm[:,0]/1000.,galprop_dm[:,1]*1000., 0.6, -2, 3)
    KoCu_sec_is     = (galprop_sec[:,0]/1000.,galprop_sec[:,1]*1000.)
    KoCu_sec        = SolarModulation( galprop_sec[:,0]/1000.,galprop_sec[:,1]*1000., 0.6, -2, 3)
    KoCu_tert_is    = (galprop_tert[:,0]/1000.,galprop_tert[:,1]*1000.)
    KoCu_tert       = SolarModulation( galprop_tert[:,0]/1000.,galprop_tert[:,1]*1000., 0.6, -2, 3)






######################################
#   PLOTTING EVERYTHING
######################################

R_27 = 0.
if species=='pbar':
    R_27 = 2.7

plot, fig = plot_1D(r'$\mathrm{T/n\,\,[GeV/n]}$', r'$\mathrm{\phi\,\,[(GeV/n)^{-1}m^{-2}s^{-1} sr^{-1}]}$')

plot.plot(KoCu       [0], np.power(KoCu       [0], R_27)*KoCu       [1],           color=c5, lw=3,                    zorder=10, label='DM CuKrKo'       )
#plot.plot(KoCu_tertDM[0], np.power(KoCu_tertDM[0], R_27)*KoCu_tertDM[1],           color=c5, lw=3,                    zorder=10, label=''                )
plot.plot(KoCu_is    [0], np.power(KoCu_is    [0], R_27)*KoCu_is    [1],           color=c5, lw=3, dashes=(10,3),     zorder= 9, label='DM CuKrKo LIS'   )
plot.fill_between(med[0], np.power(med        [0], R_27)*med[1], np.power(med       [0], R_27)*max[1],   color=c5, alpha=0.2, lw=0,         zorder= 5, label='DM MED-MAX'      )

plot.plot(KoCu_sec    [0], np.power(KoCu_sec    [0], R_27)*KoCu_sec    [1],               color=c1, lw=3,                        zorder=10, label= 'Secondary CuKrKo'                             )
plot.plot(KoCu_tert   [0], np.power(KoCu_tert   [0], R_27)*KoCu_tert   [1],               color=c2, lw=3,                        zorder=10, label= ''                             )
#plot.plot(KoCu_secR   [0], np.power(KoCu_secR   [0], R_27)*KoCu_secR   [1],               color=c1, lw=3, dashes=(3,3),          zorder= 8, label= 'Secondary CuKrKo (Green)'                     )
#plot.plot(KoCu_tertR  [0], np.power(KoCu_tertR  [0], R_27)*KoCu_tertR  [1],               color=c2, lw=3, dashes=(3,3),          zorder= 8, label= ''                     )
plot.plot(KoCu_sec_is [0], np.power(KoCu_sec_is [0], R_27)*KoCu_sec_is [1],               color=c1, lw=3, dashes=(10,3),         zorder= 9, label= 'Secondary CuKrKo LIS'                         )
plot.plot(KoCu_tert_is[0], np.power(KoCu_tert_is[0], R_27)*KoCu_tert_is[1],               color=c2, lw=3, dashes=(10,3),         zorder= 9, label= ''                         )

plot.fill_between(med_sec[0],  np.power(med_sec [0], R_27)*med_sec[1],  np.power(med_sec [0], R_27)*max_sec[1],    color=c1, alpha=0.2, lw=0,             zorder= 5, label= 'Secondary MED-MAX'                            )
plot.fill_between(med_tert[0], np.power(med_tert[0], R_27)*med_tert[1], np.power(med_tert[0], R_27)*max_tert[1],   color=c2, alpha=0.2, lw=0,             zorder= 5, label= ''                            )

plot.set_xlim( (1e-1,  1e2  ) )

handles, labels = plot.get_legend_handles_labels()

handles1 = handles[0:2]
labels1  = labels [0:2]
handles1.append(handles[-2])
labels1 .append(labels [-2])

handles2 = handles[2:5]
labels2  = labels [2:5]
handles2.append(handles[-1])
labels2 .append(labels [-1])



# print labels
# sort both labels and handles by labels

if species=='Dbar':
    
    plot.set_ylim( (1e-10, 1e-3 ) )


    leg1 = plot.legend( handles1, labels1,  loc='center left', bbox_to_anchor=(0.05, 0.45 ), ncol=1, frameon=False, fontsize=0.8*label_size)
    leg2 = plot.legend( handles2, labels2,  loc='upper right', bbox_to_anchor=(0.95, 0.70 ), ncol=1, frameon=False, fontsize=0.8*label_size)

    # 10.1103/PhysRevLett.95.081101
    bess = mpatches.Rectangle( (0.17, 1.9e-4), 1.15-0.17,  1e-3-1.9e-4, label='BESS limit (Fuke, et al. 2005)', linewidth=2, linestyle='dotted', zorder=0,  ec=c0, fc=(0,0,0,0.2) )
    bess2= mpatches.Rectangle( (0.17, 1.9e-4), 1.15-0.17,  1e-3-1.9e-4, label='BESS', linewidth=2, linestyle='dotted', zorder=20, ec=c0, fc=(0,0,0,0.0) )
    plot.add_patch(  bess  )
    plot.add_patch(  bess2 )

    #
    GAPS = mpatches.Rectangle( (0.1 , 2e-6),   0.25-0.1,   1e-3-2e-6,   label='GAPS expected sensitivity (Aramaki, et al. 2015)', linewidth=2, linestyle='dashed', zorder=1,  ec=c0, fc=(0,0,0,0.2) )
    GAPS2= mpatches.Rectangle( (0.1 , 2e-6),   0.25-0.1,   1e-3-2e-6,   label='GAPS', linewidth=2, linestyle='dashed', zorder=21, ec=c0, fc=(0,0,0,0.0) )
    plot.add_patch(  GAPS  )
    plot.add_patch(  GAPS2 )

    handles3=[bess, GAPS]
    plot.legend( handles=handles3,  loc='upper left', bbox_to_anchor=(0.17, 0.86 ), ncol=1, frameon=False, fontsize=0.8*label_size)

    plot.add_artist(leg1)
    plot.add_artist(leg2)

#   GALLPROP DBAR
    galprop_dm   = np.genfromtxt('Spectrum_ekinpern_dm_Z_-1_A_2.txt')
    Tn, F = SolarModulation( galprop_dm[:,0]/1000.,galprop_dm[:,1]*1000., 0.4, -1, 2)
    plot.plot(Tn, F , color=c0)
    galprop_sec  = np.genfromtxt('Spectrum_ekinpern_primary_Z_-1_A_2.txt')
    Tn, F = SolarModulation( galprop_sec[:,0]/1000.,galprop_sec[:,1]*1000., 0.4, -1, 2)
    plot.plot(Tn, F , color=c0)
    galprop_tert = np.genfromtxt('Spectrum_ekinpern_secondary_Z_-1_A_2.txt')
    Tn, F = SolarModulation( galprop_tert[:,0]/1000.,galprop_tert[:,1]*1000., 0.4, -1, 2)
    plot.plot(Tn, F , color=c0)



if species=='Hebar':
    
    plot.set_ylim( (1e-15, 1e-5 ) )


    leg1 = plot.legend( handles1, labels1,  loc='center left', bbox_to_anchor=(0.02, 0.65 ), ncol=1, frameon=False, fontsize=0.8*label_size)
    leg2 = plot.legend( handles2, labels2,  loc='center left', bbox_to_anchor=(0.02, 0.30 ), ncol=1, frameon=False, fontsize=0.8*label_size)
    
    He_flux         = np.genfromtxt(CRACS+'/data/project/AntiMatter/He_flux.txt', skip_header=1)
    
    phi = 0.6
    Tn_He_mod  =  He_flux[:,0] - phi / A_sp
    phi_He_mod =  ( Tn_He_mod*(Tn_He_mod+2*fMass_p)  )/(  He_flux[:,0]*(He_flux[:,0]+2*fMass_p)  ) * He_flux[:,1]
    for i in range(len(Tn_He_mod)):
        if Tn_He_mod[i]>0:
            break
    Tn_He_mod  = Tn_He_mod [i:]
    phi_He_mod = phi_He_mod[i:]

    R_He, F_He      = ConvertSpectrumFromTnToR( Tn_He_mod, phi_He_mod, 2, 4 )
    He_fluxR        = interp1d( R_He, F_He, 'cubic')

    R_bess          = np.arange( 0, 1.+1./60, 1./30)
    R_bess          = np.power( 10, R_bess);
    R_bess_limit    = He_fluxR(R_bess)*7e-8
    
    Tn_bess, Tn_bess_limit = ConvertSpectrumFromRToTn( R_bess, R_bess_limit, 2, 3 )

    verts = []
    codes = []
    
    verts.append( (Tn_bess[0], 1e0) )
    codes.append( Path.MOVETO )
    
    for i in range(len(Tn_bess)):
        verts.append( (Tn_bess[i], Tn_bess_limit[i]) )
        codes.append( Path.LINETO )

    verts.append( (Tn_bess[-1], 1e0) )
    codes.append( Path.LINETO )
    verts.append( (Tn_bess[0], 1e0) )
    codes.append( Path.CLOSEPOLY )

    path_bess = Path(verts, codes)


    R_ams          = np.arange( 0.28, 3.+1./600, 1./300)
    R_ams          = np.power( 10, R_ams)
    R_ams_limit    = He_fluxR(R_ams)
    for i in range(len(R_ams_limit)):
        if R_ams[i]<50:
            R_ams_limit[i] = R_ams_limit[i]*4e-10
        elif R_ams[i]<140:
            R_ams_limit[i] = R_ams_limit[i]*2e-9
        elif R_ams[i]<380:
            R_ams_limit[i] = R_ams_limit[i]*2e-6
        else:
            R_ams_limit[i] = R_ams_limit[i]*1e-4

    Tn_ams, Tn_ams_limit = ConvertSpectrumFromRToTn( R_ams, R_ams_limit, 2, 3 )
    Tn_ams_limit = Tn_ams_limit * 18/5  # scale 18 to 5 yeaers!

    verts = []
    codes = []

    verts.append( (Tn_ams[0], 1e0) )
    codes.append( Path.MOVETO )

    for i in range(len(Tn_ams)):
        verts.append( (Tn_ams[i], Tn_ams_limit[i]) )
        codes.append( Path.LINETO )

    verts.append( (Tn_ams[-1], 1e0) )
    codes.append( Path.LINETO )
    verts.append( (Tn_ams[0], 1e0) )
    codes.append( Path.CLOSEPOLY )

    path_AMS = Path(verts, codes)


    AMS = mpatches.PathPatch(path_AMS, label='AMS-02 expected sensitivity (5 yr, Kounine 2011)', linewidth=2, linestyle='dashed', zorder=1,  ec=c0, fc=(0,0,0,0.2) )
    plot.add_patch(  AMS  )

    bess = mpatches.PathPatch( path_bess, label='BESS limit (Fuke, et al. 2012)', linewidth=2, linestyle='dotted', zorder=0,  ec=c0, fc=(0,0,0,0.2) )
    plot.add_patch(  bess  )

    handles3=[AMS, bess]
    plot.legend( handles=handles3,  loc='upper left', bbox_to_anchor=(0.39, 0.88 ), ncol=1, frameon=False, fontsize=0.8*label_size)

    plot.add_artist(leg1)
    plot.add_artist(leg2)


if species=='pbar':
    plot.set_ylim( (1e-4, 1e1 ) )
    
    
    pbarAMS = np.genfromtxt(CRACS+'/data/AMS02/pbar_ams02.txt' )
    pbar_R      = np.sqrt( pbarAMS[:,0]*pbarAMS[:,1] )
    pbar_f      = pbarAMS[:,13]
    pbar_flux   = pbarAMS[:,10] * pbar_f
    pbar_err    = np.sqrt( np.power(pbarAMS[:,11], 2) +np.power(pbarAMS[:,12], 2) ) * pbar_f
    plot.errorbar( pbar_R, pbar_flux*np.power(pbar_R, 2.7), yerr=pbar_err*np.power(pbar_R, 2.7), fmt='o',  lw=1,  color='black', label=r'$\bar{p}$' )



plt.savefig('propagated_'+species+'.pdf')



######################################
#   PLOTTING PAPER
######################################



plot, fig = plot_1D(r'$\mathrm{T/n\,\,[GeV/n]}$', r'$\mathrm{\phi\,\,[(GeV/n)^{-1}m^{-2}s^{-1} sr^{-1}]}$')

plot.plot(KoCu       [0], KoCu       [1],   color=c5, lw=3,                    zorder=10, label='DM CuKrKo'       )
plot.fill_between(med[0], med[1], max[1],   color=(0.703,0.015,0.015, 0.2), edgecolor=c5, alpha=0.2, lw=1,         zorder= 5, label='MED-MAX'      )

plot.plot(KoCu_sec    [0], KoCu_sec    [1],               color=c1, lw=3,                        zorder=10, label= 'Secondary CuKrKo'                             )
plot.plot(KoCu_tert   [0], KoCu_tert   [1],               color=c2, lw=3,                        zorder=10, label= 'Tertiary  CuKrKo'                             )

plot.fill_between(med_sec[0],  med_sec[1],  max_sec[1],    color=(0.015,0.015,0.703, 0.2), edgecolor=c1, lw=1,             zorder= 5, label= 'MED-MAX'                            )
plot.fill_between(med_tert[0], med_tert[1], max_tert[1],   color=(0.015,0.703,0.015, 0.2), edgecolor=c2, lw=1,             zorder= 5, label= 'MED-MAX'                            )

plot.set_xlim( (1e-1,  1e2  ) )

handles, labels = plot.get_legend_handles_labels()

handles1 = handles[0:1]
labels1  = labels [0:1]
handles1.append(handles[-3])
labels1 .append(labels [-3])

handles2 = handles[1:2]
labels2  = labels [1:2]
handles2.append(handles[-2])
labels2 .append(labels [-2])

handles3 = handles[2:3]
labels3  = labels [2:3]
handles3.append(handles[-1])
labels3 .append(labels [-1])


# print labels
# sort both labels and handles by labels

if species=='Dbar':

    plot.set_ylim( (1e-10, 1e-3 ) )


    leg1 = plot.legend( handles1, labels1,  loc='center left', bbox_to_anchor=(0.05, 0.55 ), ncol=1, frameon=False, fontsize=0.8*label_size)
    leg2 = plot.legend( handles2, labels2,  loc='upper right', bbox_to_anchor=(0.95, 0.55 ), ncol=1, frameon=False, fontsize=0.8*label_size)
    leg3 = plot.legend( handles3, labels3,  loc='lower center',bbox_to_anchor=(0.30, 0.03 ), ncol=1, frameon=False, fontsize=0.8*label_size)


    plot.set_ylim( (1e-10, 1e-3 ) )

#    plot.text(3.0e-1,2.3e-04,   'Bess limit',                       color=c0, fontsize=0.8*label_size, rotation=0)
#    plot.text(1.3e-1,2.6e-06,   'GAPS expected\nsensitivity',       color=c0, fontsize=0.8*label_size, rotation='vertical', va='bottom', ha='left')
#    plot.text(2.7e+0,2.6e-06,   'AMS-02 expected\nsensitivity',     color=c0, fontsize=0.8*label_size, rotation='vertical', va='bottom', ha='left')
#    plot.text(3.0e-1,2.6e-06,   'AMS-02\nexpected\nsensitivity',    color=c0, fontsize=0.8*label_size )


    # 10.1103/PhysRevLett.95.081101
    bess = mpatches.Rectangle( (0.17, 1.9e-4), 1.15-0.17,  1e-3-1.9e-4, label='BESS limit', linewidth=2, linestyle='dashed', zorder=0,  ec=c0, fc=(0,0,0,0.4) )
    bess2= mpatches.Rectangle( (0.17, 1.9e-4), 1.15-0.17,  1e-3-1.9e-4, label='BESS', linewidth=2,                           linestyle='dashed', zorder=20, ec=c0, fc=(0,0,0,0.0) )
    plot.add_patch(  bess  )
    plot.add_patch(  bess2 )

    #
    GAPS = mpatches.Rectangle( (0.1 , 2e-6),   0.25-0.1,   1e-3-2e-6,   label='GAPS sensitivity', linewidth=2, linestyle='-', zorder=1,  ec=c0, fc=(0,0,0,0.3) )
    GAPS2= mpatches.Rectangle( (0.1 , 2e-6),   0.25-0.1,   1e-3-2e-6,   label='GAPS',                                             linewidth=2, linestyle='-', zorder=21, ec=c0, fc=(0,0,0,0.0) )
    plot.add_patch(  GAPS  )
    plot.add_patch(  GAPS2 )

    # AMS

    AMS = mpatches.Rectangle( (0.180342 , 2e-6),   0.722238-0.180342,   1e-3-2e-6,   label='AMS-02 sensitivity', linewidth=2, linestyle='dotted', zorder=1,  ec=c0, fc=(0,0,0,0.1) )
    AMS2= mpatches.Rectangle( (0.180342 , 2e-6),   0.722238-0.180342,   1e-3-2e-6,   label='AMS-02',                      linewidth=2, linestyle='dotted', zorder=21, ec=c0, fc=(0,0,0,0.0) )
    plot.add_patch(  AMS  )
    plot.add_patch(  AMS2 )


    verts = []
    codes = []

    verts.append( ( 2.40788 , 1.0e-3 ) )
    codes.append( Path.MOVETO )
    verts.append( ( 2.40788 , 1.0e-6 ) )
    codes.append( Path.LINETO )
    verts.append( ( 3.54403 , 1.0e-6 ) )
    codes.append( Path.LINETO )
    verts.append( ( 3.54403 , 1.7e-6 ) )
    codes.append( Path.LINETO )
    verts.append( ( 4.47348 , 1.7e-6 ) )
    codes.append( Path.LINETO )
    verts.append( ( 4.47348 , 1.0e-3 ) )
    codes.append( Path.LINETO )
    verts.append( ( 2.40788, 1.0e-3 ) )
    codes.append( Path.CLOSEPOLY )

    path_AMS = Path(verts, codes)

    AMS3 = mpatches.PathPatch(path_AMS, label='AMS-02 sensitivity', linewidth=2, linestyle='dotted', zorder=1,  ec=c0, fc=(0,0,0,0.1) )
    plot.add_patch(  AMS3  )
    
    handles3=[bess, GAPS, AMS]
    plot.legend( handles=handles3,  loc='upper left', bbox_to_anchor=(0.6, 0.90 ), ncol=1, frameon=False, fontsize=0.8*label_size)

#   GALLPROP DBAR
    galprop_dm   = np.genfromtxt('Spectrum_ekinpern_dm_Z_-1_A_2.txt')
    Tn, F = SolarModulation( galprop_dm[:,0]/1000.,galprop_dm[:,1]*1000., 0.4, -1, 2)
    plot.plot(Tn, F , color=c0)
    galprop_sec  = np.genfromtxt('Spectrum_ekinpern_primary_Z_-1_A_2.txt')
    Tn, F = SolarModulation( galprop_sec[:,0]/1000.,galprop_sec[:,1]*1000., 0.4, -1, 2)
    plot.plot(Tn, F , color=c0)
    galprop_tert = np.genfromtxt('Spectrum_ekinpern_secondary_Z_-1_A_2.txt')
    Tn, F = SolarModulation( galprop_tert[:,0]/1000.,galprop_tert[:,1]*1000., 0.4, -1, 2)
    plot.plot(Tn, F , color=c0)


    plot.add_artist(leg1)
    plot.add_artist(leg2)
    plot.add_artist(leg3)

if species=='Hebar':

    plot.set_ylim( (1e-15, 1e-5 ) )


    leg1 = plot.legend( handles1, labels1,  loc='center left', bbox_to_anchor=(0.05, 0.65 ), ncol=1, frameon=False, fontsize=0.8*label_size)
    leg2 = plot.legend( handles2, labels2,  loc='center left', bbox_to_anchor=(0.50, 0.45 ), ncol=1, frameon=False, fontsize=0.8*label_size)
    leg3 = plot.legend( handles3, labels3,  loc='center left', bbox_to_anchor=(0.02, 0.30 ), ncol=1, frameon=False, fontsize=0.8*label_size)

    He_flux         = np.genfromtxt(CRACS+'/data/project/AntiMatter/He_flux.txt', skip_header=1)
    phi = 0.6
    Tn_He_mod  =  He_flux[:,0] - phi / A_sp
    phi_He_mod =  ( Tn_He_mod*(Tn_He_mod+2*fMass_p)  )/(  He_flux[:,0]*(He_flux[:,0]+2*fMass_p)  ) * He_flux[:,1]
    for i in range(len(Tn_He_mod)):
        if Tn_He_mod[i]>0:
            break
    Tn_He_mod  = Tn_He_mod [i:]
    phi_He_mod = phi_He_mod[i:]

    R_He, F_He      = ConvertSpectrumFromTnToR( Tn_He_mod, phi_He_mod, 2, 4 )
    He_fluxR        = interp1d( R_He, F_He, 'cubic')

    R_bess          = np.arange( 0, 1.+1./60, 1./30)
    R_bess          = np.power( 10, R_bess);
    R_bess_limit    = He_fluxR(R_bess)*7e-8

    Tn_bess, Tn_bess_limit = ConvertSpectrumFromRToTn( R_bess, R_bess_limit, 2, 3 )

    verts = []
    codes = []

    verts.append( (Tn_bess[0], 1e0) )
    codes.append( Path.MOVETO )

    for i in range(len(Tn_bess)):
        verts.append( (Tn_bess[i], Tn_bess_limit[i]) )
        codes.append( Path.LINETO )

    verts.append( (Tn_bess[-1], 1e0) )
    codes.append( Path.LINETO )
    verts.append( (Tn_bess[0], 1e0) )
    codes.append( Path.CLOSEPOLY )

    path_bess = Path(verts, codes)


    R_ams          = np.arange( 0.28, 3.+1./600, 1./300)
    R_ams          = np.power( 10, R_ams)
    R_ams_limit    = He_fluxR(R_ams)
    for i in range(len(R_ams_limit)):
        if R_ams[i]<50:
            R_ams_limit[i] = R_ams_limit[i]*4e-10
        elif R_ams[i]<140:
            R_ams_limit[i] = R_ams_limit[i]*2e-9
        elif R_ams[i]<380:
            R_ams_limit[i] = R_ams_limit[i]*2e-6
        else:
            R_ams_limit[i] = R_ams_limit[i]*1e-4

    Tn_ams, Tn_ams_limit = ConvertSpectrumFromRToTn( R_ams, R_ams_limit, 2, 3 )
    Tn_ams_limit = Tn_ams_limit * 18/5  # scale 18 to 5 yeaers!

    verts = []
    codes = []

    verts.append( (Tn_ams[0], 1e0) )
    codes.append( Path.MOVETO )

    for i in range(len(Tn_ams)):
        verts.append( (Tn_ams[i], Tn_ams_limit[i]) )
        codes.append( Path.LINETO )

    verts.append( (Tn_ams[-1], 1e0) )
    codes.append( Path.LINETO )
    verts.append( (Tn_ams[0], 1e0) )
    codes.append( Path.CLOSEPOLY )

    path_AMS = Path(verts, codes)


    AMS = mpatches.PathPatch(path_AMS, label='AMS-02 expected sensitivity (5 yr, Kounine 2011)', linewidth=2, linestyle='dashed', zorder=1,  ec=c0, fc=(0,0,0,0.2) )
    plot.add_patch(  AMS  )

    bess = mpatches.PathPatch( path_bess, label='BESS limit (Fuke, et al. 2012)', linewidth=2, linestyle='dotted', zorder=0,  ec=c0, fc=(0,0,0,0.2) )
    plot.add_patch(  bess  )

    handles3=[bess, AMS]
    plot.legend( handles=handles3,  loc='upper left', bbox_to_anchor=(0.39, 0.88 ), ncol=1, frameon=False, fontsize=0.8*label_size)

    plot.add_artist(leg1)
    plot.add_artist(leg2)
    plot.add_artist(leg3)

if species=='pbar':
    pbarAMS = np.genfromtxt(CRACS+'/data/AMS02/pbar_ams02.txt' )
    
    pbar_R      = np.sqrt( pbarAMS[:,0]*pbarAMS[:,1] )
    pbar_f      = pbarAMS[:,13]
    pbar_flux   = pbarAMS[:,10] * pbar_f
    pbar_err    = np.sqrt( np.power(pbarAMS[:,11], 2) +np.power(pbarAMS[:,12], 2) ) * pbar_f
    plot.errorbar( pbar_R, pbar_flux, yerr=pbar_err, fmt='o',  lw=1,  color='black', label=r'$\bar{p}$' )



plt.savefig('propagated_'+species+'_paper.pdf')



######################################
#   PLOTTING PAPER     p_coal = 248
######################################

f= np.power(0.124/0.08, (A_sp-1)*3 )

plot, fig = plot_1D(r'$\mathrm{T/n\,\,[GeV/n]}$', r'$\mathrm{\phi\,\,[(GeV/n)^{-1}m^{-2}s^{-1} sr^{-1}]}$')

plot.plot(KoCu       [0], KoCu         [1]*f,   color=c5, lw=3,                    zorder=10, label='DM CuKrKo'       )
plot.fill_between(med[0], med[1]*f, max[1]*f,   color=(0.703,0.015,0.015, 0.2), edgecolor=c5, alpha=0.2, lw=1,         zorder= 5, label='MED-MAX'      )

plot.plot(KoCu_sec    [0], KoCu_sec    [1]*f,               color=c1, lw=3,                        zorder=10, label= 'Secondary CuKrKo'                             )
plot.plot(KoCu_tert   [0], KoCu_tert   [1]*f,               color=c2, lw=3,                        zorder=10, label= 'Tertiary  CuKrKo'                             )

plot.fill_between(med_sec[0],  med_sec[1]*f,  max_sec[1]*f,    color=(0.015,0.015,0.703, 0.2), edgecolor=c1, lw=1,             zorder= 5, label= 'MED-MAX'                            )
plot.fill_between(med_tert[0], med_tert[1]*f, max_tert[1]*f,   color=(0.015,0.703,0.015, 0.2), edgecolor=c2, lw=1,             zorder= 5, label= 'MED-MAX'                            )

plot.set_xlim( (1e-1,  1e2  ) )

handles, labels = plot.get_legend_handles_labels()

handles1 = handles[0:1]
labels1  = labels [0:1]
handles1.append(handles[-3])
labels1 .append(labels [-3])

handles2 = handles[1:2]
labels2  = labels [1:2]
handles2.append(handles[-2])
labels2 .append(labels [-2])

handles3 = handles[2:3]
labels3  = labels [2:3]
handles3.append(handles[-1])
labels3 .append(labels [-1])


# print labels
# sort both labels and handles by labels

if species=='Dbar':

    plot.set_ylim( (1e-10, 1e-3 ) )


    leg1 = plot.legend( handles1, labels1,  loc='center left', bbox_to_anchor=(0.10, 0.55 ), ncol=1, frameon=False, fontsize=0.8*label_size)
    leg2 = plot.legend( handles2, labels2,  loc='upper right', bbox_to_anchor=(0.95, 0.65 ), ncol=1, frameon=False, fontsize=0.8*label_size)
    #leg3 = plot.legend( handles3, labels3,  loc='center left', bbox_to_anchor=(0.02, 0.40 ), ncol=1, frameon=False, fontsize=0.8*label_size)
    leg3 = plot.legend( handles3, labels3,  loc='lower center',bbox_to_anchor=(0.30, 0.03 ), ncol=1, frameon=False, fontsize=0.8*label_size)


#    plot.text(3.0e-1,2.3e-04,   'Bess limit',                       color=c0, fontsize=0.8*label_size, rotation=0)
#    plot.text(1.3e-1,1.0e-05,   'GAPS expected\nsensitivity',       color=c0, fontsize=0.8*label_size, rotation='vertical', va='bottom', ha='left')
#    plot.text(3.0e-1,1.0e-05,   'AMS-02\nexpected\nsensitivity',    color=c0, fontsize=0.8*label_size )
#    plot.text(2.7e+0,2.6e-06,   'AMS-02 expected\nsensitivity',     color=c0, fontsize=0.8*label_size, rotation='vertical', va='bottom', ha='left')


    # 10.1103/PhysRevLett.95.081101
    bess = mpatches.Rectangle( (0.17, 1.9e-4), 1.15-0.17,  1e-3-1.9e-4, label='BESS limit', linewidth=2, linestyle='dashed', zorder=0,  ec=c0, fc=(0,0,0,0.4) )
    bess2= mpatches.Rectangle( (0.17, 1.9e-4), 1.15-0.17,  1e-3-1.9e-4, label='BESS', linewidth=2,                           linestyle='dashed', zorder=20, ec=c0, fc=(0,0,0,0.0) )
    plot.add_patch(  bess  )
    plot.add_patch(  bess2 )

    #
    GAPS = mpatches.Rectangle( (0.1 , 2e-6),   0.25-0.1,   1e-3-2e-6,   label='GAPS sensitivity', linewidth=2, linestyle='-', zorder=1,  ec=c0, fc=(0,0,0,0.3) )
    GAPS2= mpatches.Rectangle( (0.1 , 2e-6),   0.25-0.1,   1e-3-2e-6,   label='GAPS',                                             linewidth=2, linestyle='-', zorder=21, ec=c0, fc=(0,0,0,0.0) )
    plot.add_patch(  GAPS  )
    plot.add_patch(  GAPS2 )

    # AMS

    AMS = mpatches.Rectangle( (0.180342 , 2e-6),   0.722238-0.180342,   1e-3-2e-6,   label='AMS-02 sensitivity', linewidth=2, linestyle='dotted', zorder=1,  ec=c0, fc=(0,0,0,0.1) )
    AMS2= mpatches.Rectangle( (0.180342 , 2e-6),   0.722238-0.180342,   1e-3-2e-6,   label='AMS',                      linewidth=2, linestyle='dotted', zorder=21, ec=c0, fc=(0,0,0,0.0) )
    plot.add_patch(  AMS  )
    plot.add_patch(  AMS2 )


    verts = []
    codes = []

    verts.append( ( 2.40788 , 1.0e-3 ) )
    codes.append( Path.MOVETO )
    verts.append( ( 2.40788 , 1.0e-6 ) )
    codes.append( Path.LINETO )
    verts.append( ( 3.54403 , 1.0e-6 ) )
    codes.append( Path.LINETO )
    verts.append( ( 3.54403 , 1.7e-6 ) )
    codes.append( Path.LINETO )
    verts.append( ( 4.47348 , 1.7e-6 ) )
    codes.append( Path.LINETO )
    verts.append( ( 4.47348 , 1.0e-3 ) )
    codes.append( Path.LINETO )
    verts.append( ( 2.40788, 1.0e-3 ) )
    codes.append( Path.CLOSEPOLY )

    path_AMS = Path(verts, codes)

    AMS3 = mpatches.PathPatch(path_AMS, label='AMS expected sensitivity', linewidth=2, linestyle='dotted', zorder=1,  ec=c0, fc=(0,0,0,0.1) )
    plot.add_patch(  AMS3  )
    
    handles3=[bess, GAPS, AMS]
    plot.legend( handles=handles3,  loc='upper left', bbox_to_anchor=(0.6, 0.90 ), ncol=1, frameon=False, fontsize=0.8*label_size)


    plot.add_artist(leg1)
    plot.add_artist(leg2)
    plot.add_artist(leg3)

if species=='Hebar':

    plot.set_ylim( (1e-15, 1e-5 ) )


    leg1 = plot.legend( handles1, labels1,  loc='center left', bbox_to_anchor=(0.05, 0.60 ), ncol=1, frameon=False, fontsize=0.8*label_size)
    leg2 = plot.legend( handles2, labels2,  loc='center left', bbox_to_anchor=(0.60, 0.35 ), ncol=1, frameon=False, fontsize=0.8*label_size)
    leg3 = plot.legend( handles3, labels3,  loc='center left', bbox_to_anchor=(0.02, 0.20 ), ncol=1, frameon=False, fontsize=0.8*label_size)

    He_flux         = np.genfromtxt(CRACS+'/data/project/AntiMatter/He_flux.txt', skip_header=1)
    phi = 0.6
    Tn_He_mod  =  He_flux[:,0] - phi / A_sp
    phi_He_mod =  ( Tn_He_mod*(Tn_He_mod+2*fMass_p)  )/(  He_flux[:,0]*(He_flux[:,0]+2*fMass_p)  ) * He_flux[:,1]
    for i in range(len(Tn_He_mod)):
        if Tn_He_mod[i]>0:
            break
    Tn_He_mod  = Tn_He_mod [i:]
    phi_He_mod = phi_He_mod[i:]

    R_He, F_He      = ConvertSpectrumFromTnToR( Tn_He_mod, phi_He_mod, 2, 4 )
    He_fluxR        = interp1d( R_He, F_He, 'cubic')

    R_bess          = np.arange( 0, 1.+1./60, 1./30)
    R_bess          = np.power( 10, R_bess);
    R_bess_limit    = He_fluxR(R_bess)*7e-8

    Tn_bess, Tn_bess_limit = ConvertSpectrumFromRToTn( R_bess, R_bess_limit, 2, 3 )

    verts = []
    codes = []

    verts.append( (Tn_bess[0], 1e0) )
    codes.append( Path.MOVETO )

    for i in range(len(Tn_bess)):
        verts.append( (Tn_bess[i], Tn_bess_limit[i]) )
        codes.append( Path.LINETO )

    verts.append( (Tn_bess[-1], 1e0) )
    codes.append( Path.LINETO )
    verts.append( (Tn_bess[0], 1e0) )
    codes.append( Path.CLOSEPOLY )

    path_bess = Path(verts, codes)


    R_ams          = np.arange( 0.28, 3.+1./600, 1./300)
    R_ams          = np.power( 10, R_ams)
    R_ams_limit    = He_fluxR(R_ams)
    for i in range(len(R_ams_limit)):
        if R_ams[i]<50:
            R_ams_limit[i] = R_ams_limit[i]*4e-10
        elif R_ams[i]<140:
            R_ams_limit[i] = R_ams_limit[i]*2e-9
        elif R_ams[i]<380:
            R_ams_limit[i] = R_ams_limit[i]*2e-6
        else:
            R_ams_limit[i] = R_ams_limit[i]*1e-4

    Tn_ams, Tn_ams_limit = ConvertSpectrumFromRToTn( R_ams, R_ams_limit, 2, 3 )
    Tn_ams_limit = Tn_ams_limit * 18/5  # scale 18 to 5 yeaers!

    verts = []
    codes = []

    verts.append( (Tn_ams[0], 1e0) )
    codes.append( Path.MOVETO )

    for i in range(len(Tn_ams)):
        verts.append( (Tn_ams[i], Tn_ams_limit[i]) )
        codes.append( Path.LINETO )

    verts.append( (Tn_ams[-1], 1e0) )
    codes.append( Path.LINETO )
    verts.append( (Tn_ams[0], 1e0) )
    codes.append( Path.CLOSEPOLY )

    path_AMS = Path(verts, codes)


    AMS = mpatches.PathPatch(path_AMS, label='AMS-02 expected sensitivity (5 yr, Kounine 2011)', linewidth=2, linestyle='dashed', zorder=1,  ec=c0, fc=(0,0,0,0.2) )
    plot.add_patch(  AMS  )

    bess = mpatches.PathPatch( path_bess, label='BESS limit (Fuke, et al. 2012)', linewidth=2, linestyle='dotted', zorder=0,  ec=c0, fc=(0,0,0,0.2) )
    plot.add_patch(  bess  )

    handles3=[bess, AMS]
    plot.legend( handles=handles3,  loc='upper left', bbox_to_anchor=(0.39, 0.88 ), ncol=1, frameon=False, fontsize=0.8*label_size)

    plot.add_artist(leg1)
    plot.add_artist(leg2)
    plot.add_artist(leg3)


plt.savefig('propagated_'+species+'_paper_largePcoal.pdf')




#################################################
#   PLOTTING PAPER p_coal uncertainty    CuKrKo
#################################################



plot, fig = plot_1D(r'$\mathrm{T/n\,\,[GeV/n]}$', r'$\mathrm{\phi\,\,[(GeV/n)^{-1}m^{-2}s^{-1} sr^{-1}]}$')

factor = pow( 0.124/0.080,3*(A_sp-1) )
plot.fill_between(KoCu[0], KoCu     [1]+1e-90, KoCu     [1]*factor+1e-90,  color=(0.703,0.015,0.015, 0.2), edgecolor=c5, lw=1,         zorder= 1,  label='DM CuKrKo'       )
plot.fill_between(KoCu[0], KoCu_sec [1]+1e-90, KoCu_sec [1]*factor+1e-90,  color=(0.015,0.015,0.703, 0.2), edgecolor=c1, lw=1,         zorder= 4,  label='DM CuKrKo'       )
plot.fill_between(KoCu[0], KoCu_tert[1]+1e-90, KoCu_tert[1]*factor+1e-90,  color=(0.015,0.703,0.015, 0.2), edgecolor=c2, lw=1,         zorder= 5,  label='DM CuKrKo'       )

plot.set_xlim( (1e-1,  1e2  ) )



# print labels
# sort both labels and handles by labels

if species=='Dbar':

    plot.set_ylim( (1e-10, 1e-3 ) )

    plot.text(8.0e-1,1.9e-06, 'DM',            color=c5, fontsize=1.2*label_size, rotation=-15)
    plot.text(1.1e00,3.2e-09, 'Tertiary',      color=c2, fontsize=1.2*label_size, rotation=-18)
    plot.text(1.2e+1,7e-08,   'Secondary',     color=c1, fontsize=1.2*label_size, rotation=-38)

    plot.text(3.0e-1,2.3e-04,   'Bess limit',                       color=c0, fontsize=0.8*label_size, rotation=0)
    plot.text(1.2e-1,2.6e-06,   'GAPS expected\nsensitivity',       color=c0, fontsize=0.8*label_size, rotation='vertical', va='bottom', ha='left')
    plot.text(2.9e+0,2.6e-06,   'AMS-02 expected\nsensitivity',     color=c0, fontsize=0.8*label_size, rotation='vertical', va='bottom', ha='left')

    plot.text(3.0e-1,2.6e-06,   'AMS-02\nexpected\nsensitivity',    color=c0, fontsize=0.8*label_size )


    # 10.1103/PhysRevLett.95.081101
    bess = mpatches.Rectangle( (0.17, 1.9e-4), 1.15-0.17,  1e-3-1.9e-4, label='BESS limit (Fuke, et al. 2005)', linewidth=2, linestyle='dashed', zorder=0,  ec=c0, fc=(0,0,0,0.2) )
    bess2= mpatches.Rectangle( (0.17, 1.9e-4), 1.15-0.17,  1e-3-1.9e-4, label='BESS', linewidth=2,                           linestyle='dashed', zorder=20, ec=c0, fc=(0,0,0,0.0) )
    plot.add_patch(  bess  )
    plot.add_patch(  bess2 )

    #
    GAPS = mpatches.Rectangle( (0.1 , 2e-6),   0.25-0.1,   1e-3-2e-6,   label='GAPS expected sensitivity (Aramaki, et al. 2015)', linewidth=2, linestyle='-', zorder=1,  ec=c0, fc=(0,0,0,0.5) )
    GAPS2= mpatches.Rectangle( (0.1 , 2e-6),   0.25-0.1,   1e-3-2e-6,   label='GAPS',                                             linewidth=2, linestyle='-', zorder=21, ec=c0, fc=(0,0,0,0.0) )
    plot.add_patch(  GAPS  )
    plot.add_patch(  GAPS2 )

    # AMS

    AMS = mpatches.Rectangle( (0.180342 , 2e-6),   0.722238-0.180342,   1e-3-2e-6,   label='AMS expected sensitivity', linewidth=2, linestyle='dotted', zorder=1,  ec=c0, fc=(0,0,0,0.1) )
    AMS2= mpatches.Rectangle( (0.180342 , 2e-6),   0.722238-0.180342,   1e-3-2e-6,   label='AMS',                      linewidth=2, linestyle='dotted', zorder=21, ec=c0, fc=(0,0,0,0.0) )
    plot.add_patch(  AMS  )
    plot.add_patch(  AMS2 )


    verts = []
    codes = []

    verts.append( ( 2.40788 , 1.0e-3 ) )
    codes.append( Path.MOVETO )
    verts.append( ( 2.40788 , 1.0e-6 ) )
    codes.append( Path.LINETO )
    verts.append( ( 3.54403 , 1.0e-6 ) )
    codes.append( Path.LINETO )
    verts.append( ( 3.54403 , 1.7e-6 ) )
    codes.append( Path.LINETO )
    verts.append( ( 4.47348 , 1.7e-6 ) )
    codes.append( Path.LINETO )
    verts.append( ( 4.47348 , 1.0e-3 ) )
    codes.append( Path.LINETO )
    verts.append( ( 2.40788, 1.0e-3 ) )
    codes.append( Path.CLOSEPOLY )

    path_AMS = Path(verts, codes)

    AMS3 = mpatches.PathPatch(path_AMS, label='AMS expected sensitivity', linewidth=2, linestyle='dotted', zorder=1,  ec=c0, fc=(0,0,0,0.1) )
    plot.add_patch(  AMS3  )


#    handles3=[bess, GAPS, AMS]
#    plot.legend( handles=handles3,  loc='upper left', bbox_to_anchor=(0.17, 0.86 ), ncol=1, frameon=False, fontsize=0.8*label_size)


if species=='Hebar':

    plot.set_ylim( (1e-15, 1e-5 ) )

    plot.text(2.0e-1,7e-09, 'DM',            color=c5, fontsize=1.2*label_size)
    plot.text(1.5e+0,4e-13, 'Tertiary',      color=c2, fontsize=1.2*label_size, rotation= -5)
    plot.text(1.2e+1,3e-11, 'Secondary',     color=c1, fontsize=1.2*label_size, rotation=-23)

#    plot.text(7.0e-1,6.2e-08, 'AMS-02 expected sensitivity (5 yr, Kounine 2011)',    color=c0, fontsize=0.8*label_size, rotation=-22.5)
#    plot.text(1.2e+0,3e-06, 'BESS limit (Fuke, et al. 2012)',             color=c0, fontsize=0.8*label_size, rotation=-22)
    plot.text(3.5e+0,4.2e-08, 'AMS-02 expected sensitivity',    color=c0, fontsize=0.8*label_size, rotation=-27)
    plot.text(2.4e+0,4.0e-06, 'BESS limit',                     color=c0, fontsize=0.8*label_size, rotation=-22)
    plot.text(6.0e-1,2.5e-07,  '5-yr',     color=c0, fontsize=0.8*label_size, rotation=0, horizontalalignment='right')
    plot.text(6.0e-1,8.5e-08, '13-yr',     color=c0, fontsize=0.8*label_size, rotation=0, horizontalalignment='right')
    

    He_flux         = np.genfromtxt(CRACS+'/data/project/AntiMatter/He_flux.txt', skip_header=1)
    phi = 0.6
    Tn_He_mod  =  He_flux[:,0] - phi / A_sp
    phi_He_mod =  ( Tn_He_mod*(Tn_He_mod+2*fMass_p)  )/(  He_flux[:,0]*(He_flux[:,0]+2*fMass_p)  ) * He_flux[:,1]
    for i in range(len(Tn_He_mod)):
        if Tn_He_mod[i]>0:
            break
    Tn_He_mod  = Tn_He_mod [i:]
    phi_He_mod = phi_He_mod[i:]

    R_He, F_He      = ConvertSpectrumFromTnToR( Tn_He_mod, phi_He_mod, 2, 4 )
    He_fluxR        = interp1d( R_He, F_He, 'cubic')

    R_bess          = np.arange( 0, 1.+1./60, 1./30)
    R_bess          = np.power( 10, R_bess);
    R_bess_limit    = He_fluxR(R_bess)*7e-8

    Tn_bess, Tn_bess_limit = ConvertSpectrumFromRToTn( R_bess, R_bess_limit, 2, 3 )

    verts = []
    codes = []

    verts.append( (Tn_bess[0], 1e0) )
    codes.append( Path.MOVETO )

    for i in range(len(Tn_bess)):
        verts.append( (Tn_bess[i], Tn_bess_limit[i]) )
        codes.append( Path.LINETO )

    verts.append( (Tn_bess[-1], 1e0) )
    codes.append( Path.LINETO )
    verts.append( (Tn_bess[0], 1e0) )
    codes.append( Path.CLOSEPOLY )

    path_bess = Path(verts, codes)


    R_ams          = np.arange( 0.28, 3.+1./600, 1./300)
    R_ams          = np.power( 10, R_ams)
    R_ams_limit    = He_fluxR(R_ams)
    for i in range(len(R_ams_limit)):
        if R_ams[i]<50:
            R_ams_limit[i] = R_ams_limit[i]*4e-10
        elif R_ams[i]<140:
            R_ams_limit[i] = R_ams_limit[i]*2e-9
        elif R_ams[i]<380:
            R_ams_limit[i] = R_ams_limit[i]*2e-6
        else:
            R_ams_limit[i] = R_ams_limit[i]*1e-4

    Tn_ams, Tn_ams_limit = ConvertSpectrumFromRToTn( R_ams, R_ams_limit, 2, 3 )
    Tn_ams_limit = Tn_ams_limit * 18/5  # scale 18 to 5 yeaers!

    verts = []
    codes = []

    verts.append( (Tn_ams[0], 1e0) )
    codes.append( Path.MOVETO )

    for i in range(len(Tn_ams)):
        verts.append( (Tn_ams[i], Tn_ams_limit[i]) )
        codes.append( Path.LINETO )

    verts.append( (Tn_ams[-1], 1e0) )
    codes.append( Path.LINETO )
    verts.append( (Tn_ams[0], 1e0) )
    codes.append( Path.CLOSEPOLY )

    path_AMS = Path(verts, codes)

    plot.plot(Tn_ams[:], Tn_ams_limit[:]*5./13., color=c0, linewidth=2, linestyle='dashed')


    AMS = mpatches.PathPatch(path_AMS, label='AMS-02 expected sensitivity (5 yr, Kounine 2011)', linewidth=2, linestyle='dashed', zorder=1,  ec=c0, fc=(0,0,0,0.2) )
    plot.add_patch(  AMS  )

    bess = mpatches.PathPatch( path_bess, label='BESS limit (Fuke, et al. 2012)', linewidth=2, linestyle='dotted', zorder=0,  ec=c0, fc=(0,0,0,0.2) )
    plot.add_patch(  bess  )

#    handles3=[bess, AMS]
#    plot.legend( handles=handles3,  loc='upper left', bbox_to_anchor=(0.39, 0.88 ), ncol=1, frameon=False, fontsize=0.8*label_size)


plt.savefig('propagated_'+species+'_large_pcoal.pdf')



#################################################
#   PLOTTING PAPER p_coal uncertainty    MED
#################################################



plot, fig = plot_1D(r'$\mathrm{T/n\,\,[GeV/n]}$', r'$\mathrm{\phi\,\,[(GeV/n)^{-1}m^{-2}s^{-1} sr^{-1}]}$')

factor = pow( 0.124/0.080,3*(A_sp-1) )
plot.fill_between(med[0], med     [1]+1e-90, med     [1]*factor+1e-90,  color=(0.703,0.015,0.015, 0.2), edgecolor=c5, lw=1,         zorder= 1,  label='DM CuKrKo'       )
plot.fill_between(med[0], med_sec [1]+1e-90, med_sec [1]*factor+1e-90,  color=(0.015,0.015,0.703, 0.2), edgecolor=c1, lw=1,         zorder= 4,  label='DM CuKrKo'       )
plot.fill_between(med[0], med_tert[1]+1e-90, med_tert[1]*factor+1e-90,  color=(0.015,0.703,0.015, 0.2), edgecolor=c2, lw=1,         zorder= 5,  label='DM CuKrKo'       )


plot.set_xlim( (1e-1,  1e2  ) )



# print labels
# sort both labels and handles by labels

if species=='Dbar':

    plot.set_ylim( (1e-10, 1e-3 ) )

    plot.text(8.0e-1,1.9e-06, 'DM',            color=c5, fontsize=1.2*label_size, rotation=-15)
    plot.text(1.1e00,3.2e-09, 'Tertiary',      color=c2, fontsize=1.2*label_size, rotation=-18)
    plot.text(1.2e+1,7e-08,   'Secondary',     color=c1, fontsize=1.2*label_size, rotation=-38)

    plot.text(3.0e-1,2.3e-04,   'Bess limit',                       color=c0, fontsize=0.8*label_size, rotation=0)
    plot.text(1.2e-1,2.6e-06,   'GAPS expected\nsensitivity',       color=c0, fontsize=0.8*label_size, rotation='vertical', va='bottom', ha='left')
    plot.text(2.9e+0,2.6e-06,   'AMS-02 expected\nsensitivity',     color=c0, fontsize=0.8*label_size, rotation='vertical', va='bottom', ha='left')

    plot.text(3.0e-1,2.6e-06,   'AMS-02\nexpected\nsensitivity',    color=c0, fontsize=0.8*label_size )


    # 10.1103/PhysRevLett.95.081101
    bess = mpatches.Rectangle( (0.17, 1.9e-4), 1.15-0.17,  1e-3-1.9e-4, label='BESS limit (Fuke, et al. 2005)', linewidth=2, linestyle='dashed', zorder=0,  ec=c0, fc=(0,0,0,0.2) )
    bess2= mpatches.Rectangle( (0.17, 1.9e-4), 1.15-0.17,  1e-3-1.9e-4, label='BESS', linewidth=2,                           linestyle='dashed', zorder=20, ec=c0, fc=(0,0,0,0.0) )
    plot.add_patch(  bess  )
    plot.add_patch(  bess2 )

    #
    GAPS = mpatches.Rectangle( (0.1 , 2e-6),   0.25-0.1,   1e-3-2e-6,   label='GAPS expected sensitivity (Aramaki, et al. 2015)', linewidth=2, linestyle='-', zorder=1,  ec=c0, fc=(0,0,0,0.5) )
    GAPS2= mpatches.Rectangle( (0.1 , 2e-6),   0.25-0.1,   1e-3-2e-6,   label='GAPS',                                             linewidth=2, linestyle='-', zorder=21, ec=c0, fc=(0,0,0,0.0) )
    plot.add_patch(  GAPS  )
    plot.add_patch(  GAPS2 )
    
    # AMS
    
    AMS = mpatches.Rectangle( (0.180342 , 2e-6),   0.722238-0.180342,   1e-3-2e-6,   label='AMS expected sensitivity', linewidth=2, linestyle='dotted', zorder=1,  ec=c0, fc=(0,0,0,0.1) )
    AMS2= mpatches.Rectangle( (0.180342 , 2e-6),   0.722238-0.180342,   1e-3-2e-6,   label='AMS',                      linewidth=2, linestyle='dotted', zorder=21, ec=c0, fc=(0,0,0,0.0) )
    plot.add_patch(  AMS  )
    plot.add_patch(  AMS2 )
    
    
    verts = []
    codes = []

    verts.append( ( 2.40788 , 1.0e-3 ) )
    codes.append( Path.MOVETO )
    verts.append( ( 2.40788 , 1.0e-6 ) )
    codes.append( Path.LINETO )
    verts.append( ( 3.54403 , 1.0e-6 ) )
    codes.append( Path.LINETO )
    verts.append( ( 3.54403 , 1.7e-6 ) )
    codes.append( Path.LINETO )
    verts.append( ( 4.47348 , 1.7e-6 ) )
    codes.append( Path.LINETO )
    verts.append( ( 4.47348 , 1.0e-3 ) )
    codes.append( Path.LINETO )
    verts.append( ( 2.40788, 1.0e-3 ) )
    codes.append( Path.CLOSEPOLY )
    
    path_AMS = Path(verts, codes)

    AMS3 = mpatches.PathPatch(path_AMS, label='AMS expected sensitivity', linewidth=2, linestyle='dotted', zorder=1,  ec=c0, fc=(0,0,0,0.1) )
    plot.add_patch(  AMS3  )

    handles3=[bess, GAPS]
    #plot.legend( handles=handles3,  loc='upper left', bbox_to_anchor=(0.17, 0.86 ), ncol=1, frameon=False, fontsize=0.8*label_size)


if species=='Hebar':

    plot.set_ylim( (1e-15, 1e-5 ) )

    plot.text(2.0e-1,4.0e-10, 'DM',            color=c5, fontsize=1.2*label_size)
    plot.text(2.0e-1,1.0e-12, 'Tertiary',      color=c2, fontsize=1.2*label_size, rotation= +5)
    plot.text(1.2e+1,3.0e-11, 'Secondary',     color=c1, fontsize=1.2*label_size, rotation=-23)

    #plot.text(7.0e-1,6.2e-08, 'AMS-02 expected sensitivity (5 yr, Kounine 2011)',    color=c0, fontsize=0.8*label_size, rotation=-22.5)
    #plot.text(1.2e+0,3.0e-06, 'BESS limit (Fuke, et al. 2012)',             color=c0, fontsize=0.8*label_size, rotation=-22)
    plot.text(3.5e+0,4.2e-08, 'AMS-02 expected sensitivity',    color=c0, fontsize=0.8*label_size, rotation=-27)
    plot.text(2.4e+0,4.0e-06, 'BESS limit',                     color=c0, fontsize=0.8*label_size, rotation=-22)
    plot.text(6.0e-1,2.5e-07,  '5-yr',     color=c0, fontsize=0.8*label_size, rotation=0, horizontalalignment='right')
    plot.text(6.0e-1,8.5e-08, '13-yr',     color=c0, fontsize=0.8*label_size, rotation=0, horizontalalignment='right')
    

    He_flux         = np.genfromtxt(CRACS+'/data/project/AntiMatter/He_flux.txt', skip_header=1)
    phi = 0.6
    Tn_He_mod  =  He_flux[:,0] - phi / A_sp
    phi_He_mod =  ( Tn_He_mod*(Tn_He_mod+2*fMass_p)  )/(  He_flux[:,0]*(He_flux[:,0]+2*fMass_p)  ) * He_flux[:,1]
    for i in range(len(Tn_He_mod)):
        if Tn_He_mod[i]>0:
            break
    Tn_He_mod  = Tn_He_mod [i:]
    phi_He_mod = phi_He_mod[i:]

    R_He, F_He      = ConvertSpectrumFromTnToR( Tn_He_mod, phi_He_mod, 2, 4 )
    He_fluxR        = interp1d( R_He, F_He, 'cubic')

    R_bess          = np.arange( 0, 1.+1./60, 1./30)
    R_bess          = np.power( 10, R_bess);
    R_bess_limit    = He_fluxR(R_bess)*7e-8

    Tn_bess, Tn_bess_limit = ConvertSpectrumFromRToTn( R_bess, R_bess_limit, 2, 3 )

    verts = []
    codes = []

    verts.append( (Tn_bess[0], 1e0) )
    codes.append( Path.MOVETO )

    for i in range(len(Tn_bess)):
        verts.append( (Tn_bess[i], Tn_bess_limit[i]) )
        codes.append( Path.LINETO )

    verts.append( (Tn_bess[-1], 1e0) )
    codes.append( Path.LINETO )
    verts.append( (Tn_bess[0], 1e0) )
    codes.append( Path.CLOSEPOLY )

    path_bess = Path(verts, codes)


    R_ams          = np.arange( 0.28, 3.+1./600, 1./300)
    R_ams          = np.power( 10, R_ams)
    R_ams_limit    = He_fluxR(R_ams)
    for i in range(len(R_ams_limit)):
        if R_ams[i]<50:
            R_ams_limit[i] = R_ams_limit[i]*4e-10
        elif R_ams[i]<140:
            R_ams_limit[i] = R_ams_limit[i]*2e-9
        elif R_ams[i]<380:
            R_ams_limit[i] = R_ams_limit[i]*2e-6
        else:
            R_ams_limit[i] = R_ams_limit[i]*1e-4

    Tn_ams, Tn_ams_limit = ConvertSpectrumFromRToTn( R_ams, R_ams_limit, 2, 3 )

    Tn_ams_limit = Tn_ams_limit * 18/5  # scale 18 to 5 yeaers!

    verts = []
    codes = []

    verts.append( (Tn_ams[0], 1e0) )
    codes.append( Path.MOVETO )

    for i in range(len(Tn_ams)):
        verts.append( (Tn_ams[i], Tn_ams_limit[i]) )
        codes.append( Path.LINETO )

    verts.append( (Tn_ams[-1], 1e0) )
    codes.append( Path.LINETO )
    verts.append( (Tn_ams[0], 1e0) )
    codes.append( Path.CLOSEPOLY )

    path_AMS = Path(verts, codes)

    plot.plot(Tn_ams[:], Tn_ams_limit[:]*5./13., color=c0, linewidth=2, linestyle='dashed')


    AMS = mpatches.PathPatch(path_AMS, label='AMS-02 expected sensitivity (5 yr, Kounine 2011)', linewidth=2, linestyle='dashed', zorder=1,  ec=c0, fc=(0,0,0,0.2) )
    plot.add_patch(  AMS  )

    bess = mpatches.PathPatch( path_bess, label='BESS limit (Fuke, et al. 2012)', linewidth=2, linestyle='dotted', zorder=0,  ec=c0, fc=(0,0,0,0.2) )
    plot.add_patch(  bess  )

#    handles3=[bess, AMS]
#    plot.legend( handles=handles3,  loc='upper left', bbox_to_anchor=(0.39, 0.88 ), ncol=1, frameon=False, fontsize=0.8*label_size)


plt.savefig('propagated_'+species+'_MED_large_pcoal.pdf')




