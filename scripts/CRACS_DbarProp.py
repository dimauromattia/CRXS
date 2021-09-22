#! /usr/bin/env python

import glob
import math
import sys

import numpy as np
import scipy as sp
import scipy.special as sps

import matplotlib.pyplot    as plt
import matplotlib           as mpl
import matplotlib.colors    as colors
import matplotlib.patches   as mpatches
from   matplotlib import rc

from scipy.optimize import curve_fit


c0='black'
c2='#04B404'
c4='#B40486'
c3='#0489B1'
c1='#0404B4'
c5='#B40404'


suffix=''

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
    
    plot.tick_params('both', length=20, width=2, which='major')
    plot.tick_params('both', length=10, width=1, which='minor')
    
    plot.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
    plot.tick_params(axis='both', pad=10)
    
    return plot, fig


data_Dbar                       = np.genfromtxt( 'Dbar_SourceTerm.txt',           skip_header=1 )
data_Dbar_DM                    = np.genfromtxt( 'Dbar_SourceTerm_DM.txt',        skip_header=1 )
data_Dbar_DM_energySpectrum     = np.genfromtxt( 'Dbar_dN_by_dTn.txt',            skip_header=1 )
data_XS_pp                      = np.genfromtxt( 'XS_tot_pp.txt',                 skip_header=1 )
data_pbar_Cirelli               = np.genfromtxt( 'pbar_dN_by_dT.txt',             skip_header=1 )

data_pbar_vittino               = np.genfromtxt( 'antiproton_vittino.dat',        skip_header=1 )

#col 0 :            sqrts (CoM energy)
#col 1 :            x = kinetic energy/sqrts
#col 2 , ..., 12 :  dN/dx for each of the 11 channels, as:
#                   ee  mumu  tautau  uu  dd  ss  cc  bb  toptop  WW  ZZ
#                   2    3    4       5   6   7   8   9   10      11  12
#find range of sqrtS = 200
i_from = 0
i_to   = len(data_pbar_vittino[:,0])
for i in range(i_from, i_to ):
    if  data_pbar_vittino[i,0] < 200:
        i_from = i
    if  data_pbar_vittino[i,0] > 200:
        i_to   = i
        break

T_p_vittino  = data_pbar_vittino[i_from:i_to,1]*data_pbar_vittino[i_from:i_to,0]
dN_dT_vittino= data_pbar_vittino[i_from:i_to,9]/data_pbar_vittino[i_from:i_to,0]


plot, fig = plot_1D('T [GeV]', 'dN/dE  [1/GeV]')
plot.plot(data_pbar_Cirelli[:,0],data_pbar_Cirelli[:,1], color=c1, lw=3, label='Cirelli')
plot.plot(T_p_vittino, dN_dT_vittino, color=c2, dashes=(7,7), lw=3, label='Vittino')
plot.set_ylim( (1e-3, 1e0) )
plot.set_xlim( (1e-2, 1e2) )
plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), ncol=1, frameon=False, fontsize=0.7*label_size)
plt.savefig('pbar_dN_by_dT'+suffix+'.png')


plot, fig = plot_1D('T/n [GeV/n]', 'dN/(dT/n)  [n/GeV]')
plot.plot(data_Dbar_DM_energySpectrum[:,0],data_Dbar_DM_energySpectrum[:,1], color=c1, lw=3, label='cpp')

T_p     = data_pbar_Cirelli[:,0]        #  = T_D/nucleon
dNdT    = data_pbar_Cirelli[:,1]
p_D     = 2. * np.sqrt(  T_p*T_p +  T_p*0.938*2  )
dNDdT   = 2. * np.power(0.160,3)  / 6. / p_D * 0.938*2/0.938/0.938 * dNdT*dNdT
plot.plot( T_p, dNDdT, color=c2, lw=3, dashes=(7,7), label='python')


plot.set_ylim( (1e-7, 1e-3) )
plot.set_xlim( (1e-2, 1e2 ) )
plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), ncol=1, frameon=False, fontsize=0.7*label_size)
plt.savefig('Dbar_dN_by_dT'+suffix+'.png')



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

fMass_d = 1.875612928   # GeV

Z       = L             # kpc
R       = 20.           # kpc
Tn_min  = 0.1           # GeV/n
Tn_max  = 100.          # GeV/n

nE      = 50
nZ      = 50
nR      = 500
nI      = 100

nZ      =  20
nE      =  30
nR      =  80
nI      =  40

######################################
#   create grid
######################################
dR  = 1.*R/(nR-1)
dZ  = 1.*Z/(nZ-1)
dEf = 1.*np.log(1.*Tn_max/Tn_min)/(nE-1)

i_vec   = np.arange(0,nI)
r_vec   = np.arange(0, R+dR/2., dR)
z_vec   = np.arange(0, Z+dZ/2., dZ)
Tn_vec  = np.arange( np.log(Tn_min), np.log(Tn_max)+ dEf/2., dEf)
Tn_vec  = np.exp(Tn_vec)

i_grid, Tn_grid, z_grid, r_grid = np.meshgrid(i_vec, Tn_vec, z_vec, r_vec, indexing='ij')

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


def proper(case='DM'):
    
    global suffix
    global h, K_0, V_c, delta, L, sv, mDM, rhoSun
    global Z,  R, Tn_min, Tn_max, nE, nZ, nR, nI, dR, dZ, dEf, i_vec, r_vec, z_vec, Tn_vec, i_grid, Tn_grid, z_grid, r_grid
    
    global data_Dbar_DM_energySpectrum, data_Dbar
    
    print 'Propergate: '+suffix
    if mDM!=100:
        'WARNING:  adjust the file with dN/dT data to your DM mass!'

    num             = 1e-50
    rr_grid         = np.sqrt(r_grid*r_grid + z_grid*z_grid)+num
    
    ######################################
    #   create q_D (E, r, z)
    ######################################

    Rs              = 20.
    rho_0           = 1./ (rSun/Rs * np.power(1 + rSun/Rs, 2))
    rho_grid        = 1./ (rr_grid       /Rs * np.power(1 + rr_grid  /Rs, 2)) * rhoSun / rho_0

    dN_dTn_grid     = np.interp( Tn_grid, data_Dbar_DM_energySpectrum[:,0], data_Dbar_DM_energySpectrum[:,1]  )

    rho_q_D_grid  =  rho_grid*rho_grid/rhoSun/rhoSun

    plot, fig = plot_1D(r'$\mathrm{r\,\,[kpc]}$', r'$\mathrm{\rho_{NFW}}/\rho_\odot$')
    plot.plot(r_vec, rho_grid[0,0,0,:], color=c1, lw=3, label='MK')
    plot.set_xlim( (1e-1, 1e2 ) )
    plot.set_ylim( (1e-1, 1e3 ) )
    plt.savefig('rho_NFW'+suffix+'.png')

    ######################################
    #   create q_i (E, z)
    ######################################

    i_grid2     = i_grid [:,:,:,0]
    Tn_grid2    = Tn_grid[:,:,:,0]
    z_grid2     = z_grid [:,:,:,0]

    if case =='DM':
        integrand_q_i_grid  = r_grid  *  sps.j0( vect_xi(i_grid)*r_grid/R ) * rho_q_D_grid
    if case =='Secondary':
        integrand_q_i_grid  = r_grid  *  sps.j0( vect_xi(i_grid)*r_grid/R ) * vect_delta_z(z_grid)
    if case =='Secondary_Green':
        integrand_q_i_grid  = r_grid  *  sps.j0( vect_xi(i_grid)*r_grid/R ) * vect_delta_z(z_grid)*np.power(r_grid/rSun, 1.09)*np.exp(-3.87*(r_grid-rSun)/rSun)

    pref_grid2          = 4./np.power(sps.j1( vect_xi(i_grid2) )*R, 2)
    q_i_grid2           = integrand_q_i_grid.sum(axis=3) * dR * pref_grid2




    ######################################
    #   calculate y_i (E)
    ######################################

    E_grid2     =   Tn_grid2*2 + fMass_d
    p_grid2     =   2.*np.sqrt(  Tn_grid2*(Tn_grid2+fMass_d)  )
    beta_grid2  =   p_grid2/E_grid2
    rig_grid2   =   p_grid2
    K_grid2     =   K_0 * beta_grid2 * np.power(rig_grid2, delta)

    S_i_grid2   =   np.sqrt(  np.power(  V_c/K_grid2, 2 ) + 4*np.power(vect_xi(i_grid2)/R, 2)  )            #   in 1/kpc


    plot, fig = plot_1D(r'$\mathrm{T/n\,\,[GeV/n]}$', r'$K$')
    plot.plot(Tn_vec, K_grid2[0,:,0], color=c1, lw=3, label='MK')
    plt.savefig('K'+suffix+'.png')

    plot, fig = plot_1D(r'$\mathrm{T/n\,\,[GeV/n]}$', r'$\beta$')
    plot.plot(Tn_vec, beta_grid2[0,:,0], color=c1, lw=3, label='MK')
    plt.savefig('beta'+suffix+'.png')


    #########  TO DO      ########
    nH  = 1.    *   1e6
    nHe = 0.1   *   1e6
    sigma_in_pp__grid2     = 2.* np.interp( Tn_grid2, data_XS_pp[:,0], data_XS_pp[:,1]  ) * 1e-31
    Gamma_grid2            = (nH + np.power(4,2./3)*nHe)*sigma_in_pp__grid2*beta_grid2*2.99792e8

    #########  END TO DO  ########


    A_i_grid2   =   V_c + 2*h*Gamma_grid2 + K_grid2*S_i_grid2  *  1./np.tanh(S_i_grid2*L/2.)


    plot, fig = plot_1D(r'$\mathrm{T/n\,\,[GeV/n]}$', r'$A_i$')
    plot.plot(Tn_vec, A_i_grid2[0,:,0], color=c1, lw=3, label='MK')
    plt.savefig('A_i'+suffix+'.png')


#    print "q_i"
#    print q_i_grid2

    integrand__y_i_grid2    = np.exp( V_c/2./K_grid2 * (L-z_grid2) )  * np.sinh(  S_i_grid2 * (L-z_grid2) /2.  ) * q_i_grid2

    y_i_grid3 = integrand__y_i_grid2.sum(axis=2) * dZ

    plot, fig = plot_1D(r'$\mathrm{T/n\,\,[GeV/n]}$', r'$y_i$')
    plot.plot(Tn_vec, y_i_grid3[0,:], color=c1, lw=3, label='i=0')
    plot.plot(Tn_vec, y_i_grid3[1,:], color=c2, lw=3, label='i=1')
    plot.plot(Tn_vec, y_i_grid3[2,:], color=c3, lw=3, label='i=2')
    plot.plot(Tn_vec, y_i_grid3[10,:],color=c4, lw=3, label='i=10')
    plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), ncol=1, frameon=False, fontsize=0.7*label_size)
    plt.savefig('y_i'+suffix+'.png')

    ######################################
    #   calculate R_D
    ######################################


    i_grid3     = i_grid [:,:,0,0]
    Tn_grid3    = Tn_grid[:,:,0,0]
    z_grid3     = z_grid [:,:,0,0]

    K_grid3     = K_grid2   [:,:,0]
    A_i_grid3   = A_i_grid2 [:,:,0]
    S_i_grid3   = S_i_grid2 [:,:,0]


    plot, fig = plot_1D(r'$\mathrm{T/n\,\,[GeV/n]}$', r'$\mathrm{S_i\,\,[]}$')
    plot.plot(Tn_vec, S_i_grid3[0,:], color=c1, lw=3, label='MK')
    plt.savefig('S_i'+suffix+'.png')


    summand_R_D = sps.j0(vect_xi(i_grid3)*rSun/R) * np.exp(-V_c * L / 2. / K_grid3)* y_i_grid3/A_i_grid3/np.sinh( S_i_grid3 * L / 2. )
    R_D = summand_R_D.sum(axis=0)

    #    print "y"
    #    print y_i_grid3
    #
    #    print "A"
    #    print A_i_grid3
    #
    #
    #    print "S"
    #    print S_i_grid3
    #
    #    print "sum R_D"
    #    print summand_R_D
    #
    #    print "R_D"
    #    print R_D
    plot, fig = plot_1D(r'$\mathrm{T/n\,\,[GeV/n]}$', r'$\mathrm{R_D\,\,[yr]}$')
    plot.plot(Tn_vec, R_D*3.17098e-8, color=c1, lw=3, label='MK')
    plt.savefig('R_D'+suffix+'.png')


    if case =='DM':
        n_D = 0.5* sv *rhoSun*rhoSun/mDM/mDM * dN_dTn_grid[0,:,0,0] * R_D
        n_D = 1e6*n_D                                                           # Transform all units to GeV, m, and s
    if case =='Secondary' or case=='Secondary_Green':
        n_D = np.interp( Tn_vec, data_Dbar[:,0], data_Dbar[:,1]  ) * R_D

    plot, fig = plot_1D(r'$\mathrm{T/n\,\,[GeV/n]}$', r'$\mathrm{n_D\,\,[]}$')
    plot.plot(Tn_vec, n_D, color=c1, lw=3, label='MK')
    plt.savefig('n_D'+suffix+'.png')



    ######################################
    #   calculate flux
    ######################################

    phi_D = n_D * beta_grid2[0,:,0] / 4. / math.pi * 2.99792e8

    plot, fig = plot_1D(r'$\mathrm{T/n\,\,[GeV/n]}$', r'$\mathrm{\phi_D\,\,[(GeV/n)^{-1}m^{-2}s^{-1} sr^{-1}]}$')
    plot.plot(Tn_vec, phi_D, color=c1, lw=3, label='MK')
    plot.set_xlim( (1e-1, 1e2 ) )
    plot.set_ylim( (1e-10, 1e-4 ) )
    plt.savefig('phi_D_IS'+suffix+'.png')

    phi=0.4

    Tn_vec_mod = Tn_vec - phi / 2.

    phi_D_mod =  np.sqrt(  Tn_vec_mod*(Tn_vec_mod+fMass_d)  )/np.sqrt(  Tn_vec*(Tn_vec_mod+fMass_d)  ) * phi_D

    plot, fig = plot_1D(r'$\mathrm{T/n\,\,[GeV/n]}$', r'$\mathrm{\phi_D\,\,[(GeV/n)^{-1}m^{-2}s^{-1} sr^{-1}]}$')
    plot.plot(Tn_vec_mod, phi_D_mod, color=c1, lw=3, label='MK')
    plot.set_xlim( (1e-1, 1e2 ) )
    plot.set_ylim( (1e-10, 1e-4 ) )
    plt.savefig('phi_D_SMod'+suffix+'.png')

    return (Tn_vec, phi_D), (Tn_vec_mod, phi_D_mod)




#######################################
##   PROPERGATE MIN                              FIXME: Here is some numerical instability...
#######################################
#
#h       =   0.1         # kpc
#K_0     =   0.0016      # kpc^2/Myr          MIN
#V_c     =   13.5        # km/s               MIN
#delta   =   0.85        #                    MIN
#L       =   1.          # kpc                MIN
#
## transform
#K_0     =   K_0 / (1e6 * 3.154e+7)           # in kpc^2/s
#V_c     =   V_c * 3.24078e-17                # transform to kpc/s
#
#suffix = '_MIN'
#proper()

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

######################################
#   POPERGATE KoCu
######################################

h       =   0.2         # kpc
D_0     =   9.84e28     # cm^2/s             KoCu
V_c     =   45.3        # km/s               KoCu
delta   =   0.245       #                    KoCu
L       =   5.35        # kpc                KoCu

# transform
D_0     =   D_0 * np.power(  (1e-3*3.24078e-17*1e-2),2  )      # in kpc^2/s
K_0     =   D_0 / np.power( 4., delta )              # from 4 to 1 GV
V_c     =   V_c * 3.24078e-17                        # transform to kpc/s

suffix = '_KoCu'
KoCu_is, KoCu = proper()

suffix = '_Secondary_KoCu'
KoCu_sec_is, KoCu_sec = proper(case='Secondary')

suffix = '_Secondary_KoCu_Green'
KoCu_secR_is, KoCu_secR = proper(case='Secondary_Green')



plot, fig = plot_1D(r'$\mathrm{T/n\,\,[GeV/n]}$', r'$\mathrm{\phi_D\,\,[(GeV/n)^{-1}m^{-2}s^{-1} sr^{-1}]}$')

plot.plot(KoCu   [0], KoCu   [1], color=c1, lw=3,                    zorder=10, label='DM KoCu'         )
plot.plot(KoCu_is[0], KoCu_is[1], color=c1, lw=3, dashes=(7,7),      zorder= 9, label='DM KoCu (IS)'    )
plot.fill_between(med[0], med[1], max[1], color=c0, alpha=0.2, lw=0, zorder= 0, label='MED - MAX'       )

plot.set_xlim( (1e-1, 1e2 ) )
plot.set_ylim( (1e-10, 1e-4 ) )
plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), ncol=1, frameon=False, fontsize=1.0*label_size)
plt.savefig('propagated_DM_flux.png')



plot, fig = plot_1D(r'$\mathrm{T/n\,\,[GeV/n]}$', r'$\mathrm{\phi_D\,\,[(GeV/n)^{-1}m^{-2}s^{-1} sr^{-1}]}$')

plot.plot(KoCu_sec   [0], KoCu_sec   [1], color=c1, lw=3,                    zorder=10, label='Secondary CuKrKo'         )
plot.plot(KoCu_sec_is[0], KoCu_sec_is[1], color=c1, lw=3, dashes=(7,3),      zorder= 9, label='Secondary CuKrKo LIS'    )
plot.plot(KoCu_secR  [0], KoCu_secR  [1], color=c1, lw=3, dashes=(3,3),      zorder= 8, label='Secondary CuKrKo (Green)' )

plot.set_xlim( (1e-1, 1e2 ) )
plot.set_ylim( (1e-10, 1e-4 ) )
plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), ncol=1, frameon=False, fontsize=1.0*label_size)
plt.savefig('propagated_secondary_flux.png')



plot, fig = plot_1D(r'$\mathrm{T/n\,\,[GeV/n]}$', r'$\mathrm{\phi_D\,\,[(GeV/n)^{-1}m^{-2}s^{-1} sr^{-1}]}$')

plot.plot(KoCu   [0], KoCu   [1],           color=c5, lw=3,                    zorder=10, label='DM CuKrKo'       )
#plot.plot(np.NaN, np.NaN, '-',              color='none',                                 label=' ')
plot.plot(KoCu_is[0], KoCu_is[1],           color=c5, lw=3, dashes=(10,3),     zorder= 9, label='DM CuKrKo LIS'   )
plot.fill_between(med[0], med[1], max[1],   color=c5, alpha=0.2, lw=0,         zorder= 5, label='DM MED-MAX'      )


plot.plot(KoCu_sec   [0], KoCu_sec   [1],               color=c1, lw=3,                        zorder=10, label= 'Secondary CuKrKo'                             )
plot.plot(KoCu_secR  [0], KoCu_secR  [1],               color=c1, lw=3, dashes=(3,3),          zorder= 8, label=r'Secondary CuKrKo (Green)'                     )
plot.plot(KoCu_sec_is[0], KoCu_sec_is[1],               color=c1, lw=3, dashes=(10,3),         zorder= 9, label= 'Secondary CuKrKo LIS'                         )
plot.fill_between(med_sec[0], med_sec[1], max_sec[1],   color=c1, alpha=0.2, lw=0,             zorder= 5, label= 'Secondary MED-MAX'                            )

#plot.plot(KoCu_sec   [0], KoCu_sec[1]+KoCu[1],          color=c0, lw=3,                        zorder= 6, label=  ''                             )


plot.set_xlim( (1e-1,  1e2  ) )
plot.set_ylim( (1e-10, 1e-3 ) )


handles, labels = plot.get_legend_handles_labels()

handles1 = handles[0:2]
labels1  = labels [0:2]
handles1.append(handles[-2])
labels1 .append(labels [-2])

handles2 = handles[2:5]
labels2  = labels [2:5]
handles2.append(handles[-1])
labels2 .append(labels [-1])



#print labels
# sort both labels and handles by labels

leg1 = plot.legend( handles1, labels1,  loc='center left', bbox_to_anchor=(0.05, 0.45 ), ncol=1, frameon=False, fontsize=0.8*label_size)
leg2 = plot.legend( handles2, labels2,  loc='upper right', bbox_to_anchor=(0.95, 0.70 ), ncol=1, frameon=False, fontsize=0.8*label_size)


# 10.1103/PhysRevLett.95.081101
bess = mpatches.Rectangle( (0.17, 1.9e-4), 1.15-0.17,  1e-3-1.9e-4, label='BESS limit (Fuke, et al. 2005)', linewidth=2, linestyle='dotted', zorder=0,  ec=c0, fc=(0,0,0,0.2) )
bess2= mpatches.Rectangle( (0.17, 1.9e-4), 1.15-0.17,  1e-3-1.9e-4, label='BESS', linewidth=2, linestyle='dotted', zorder=20, ec=c0, fc=(0,0,0,0.0) )
plot.add_patch(  bess  )
plot.add_patch(  bess2 )

#
GAPS = mpatches.Rectangle( (0.1 , 2e-6),   0.25-0.1,   1e-3-2e-6,   label='GAPS expected sensitivit (Aramaki, et al. 2015)', linewidth=2, linestyle='dashed', zorder=1,  ec=c0, fc=(0,0,0,0.2) )
GAPS2= mpatches.Rectangle( (0.1 , 2e-6),   0.25-0.1,   1e-3-2e-6,   label='GAPS', linewidth=2, linestyle='dashed', zorder=21, ec=c0, fc=(0,0,0,0.0) )
plot.add_patch(  GAPS  )
plot.add_patch(  GAPS2 )

handles3=[bess, GAPS]
plot.legend( handles=handles3,  loc='upper left', bbox_to_anchor=(0.17, 0.86 ), ncol=1, frameon=False, fontsize=0.8*label_size)


plot.add_artist(leg1)
plot.add_artist(leg2)


plt.savefig('propagated.pdf')




