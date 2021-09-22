#! /usr/bin/env python

import glob
import math
import sys
import os

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
parser.add_argument('--step'   , help='step',          action='store', dest='step'   , type=int, default=1     )
args = parser.parse_args()
species=args.species
step   =args.step

suffix = ''

A_sp            = 2.
if species=='Dbar':
    A_sp        = 2.
if species=='Hebar':
    A_sp        = 3.


CRACS = os.environ['CRACS']

data                       = np.genfromtxt( species+'_SourceTerm.txt',           skip_header=1 )
#data_DM_energySpectrum     = np.genfromtxt( species+'_dN_by_dTn.txt',            skip_header=1 )

# XS:
# XS:
data_dpbar_nar             = np.genfromtxt( CRACS+'/data/project/AntiMatter/table_dpbar_nar.txt',               skip_header=1 )
data_dpbar_tot             = np.genfromtxt( CRACS+'/data/project/AntiMatter/table_dpbar_tot.txt',               skip_header=1 )
data_ppbar_el              = np.genfromtxt( CRACS+'/data/project/AntiMatter/table_ppbar_el.txt',                skip_header=1 )
data_ppbar_tot             = np.genfromtxt( CRACS+'/data/project/AntiMatter/table_ppbar_tot.txt',               skip_header=1 )

data_Anderson              = np.genfromtxt( CRACS+'/data/CS_lab/dT_pp_p_LAB_Anderson.txt'                         )


data_cirelli_pbar          = np.genfromtxt( CRACS+'/data/project/AntiMatter/AtProduction_antiprotons.dat',      skip_header=1 )
f_data_cirelli_pbar__bb    = interp2d(  data_cirelli_pbar[:,0], data_cirelli_pbar[:,1], data_cirelli_pbar[:,13], kind='cubic')

#data_vittino_dbar          = np.genfromtxt( 'Dbar_dN_dx__bb_100_vittino.txt' )
#data_vittino_dbar_multip   = np.genfromtxt( 'Dbar_multiplicity_bb_vittino.txt' )
#f_data_vittino_dbar_multip = interp1d(  data_vittino_dbar_multip[:,0], data_vittino_dbar_multip[:,1], kind='cubic'  )

c0='black'
c2='#04B404'
c4='#B40486'
c3='#0489B1'
c1='#0404B4'
c5='#B40404'




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



######################################
#   grid parameters
######################################

fMass_d = 1.875612928   # GeV
fMass_p = 0.9382720813  # GeV
fMass_n = 0.939565379   # GeV

R       = 20.           # kpc
Tn_min  = 1e-1          # GeV/n
Tn_max  = 1e3           # GeV/n

nE      = 50
nZ      = 50
nR      = 500
nI      = 100

nZ      =  20
nE      =  40
nR      =  80
nI      =  40

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

def delta_z(z, h):
    global dZ
    if z<1e-10:
        return 1.*h/dZ
    return 0
vect_delta_z = np.vectorize(delta_z)



def proper(ret='spectra', propagation='MED', dark_matter='', p_coal=0.078, injection='analytic'):

    global data, data_tert, f_data_cirelli_pbar__bb
    
    if propagation=='galprop':
        cmd = 'galprop -o . -g '+CRACS+'/data/project/AntiMatter -r standard  -f /usr/local/korsmeier/galprop/data/FITS'
        os.system( cmd )
        os.system('PROPERFIT_plotGalpropSpectrum nuclei_56_standard --A 2 --Z -1 --p 0 --t dm        --run false')
        galprop_dm   = np.genfromtxt('Spectrum_ekinpern_dm_Z_-1_A_2.txt')

        phi_D = np.interp(Tn_vec, galprop_dm[:,0]/1000., galprop_dm[:,1]*1000.)

    else:
        ######################################
        #   propagation parameters   (default MED)
        ######################################

        h       =   0.1             # kpc
        K_0     =   0.0112          # kpc^2/Myr          MED
        V_c     =   12.             # km/s               MED
        delta   =   0.7             #                    MED
        L       =   4.              # kpc                MED
        
        if propagation=='MAX':
            h       =   0.1         # kpc
            K_0     =   0.0765      # kpc^2/Myr          MAX
            V_c     =   5.          # km/s               MAX
            delta   =   0.46        #                    MAX
            L       =  15.          # kpc                MAX

        if propagation!='MED' and propagation!='MAX':
            h       =   propagation[0]       # kpc
            K_0     =   propagation[1]       # kpc^2/Myr
            V_c     =   propagation[2]       # km/s
            delta   =   propagation[3]       #
            L       =   propagation[4]       # kpc


        Z       = L             # kpc
        dZ  = 1.*Z/(nZ-1)
        z_vec   = np.arange(0, Z+dZ/2., dZ)

        i_grid, Tn_grid, z_grid, r_grid = np.meshgrid(i_vec, Tn_vec, z_vec, r_vec, indexing='ij')




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

        if dark_matter!='':
            mDM     = dark_matter[0]
            sv      = dark_matter[1]
            rhoSun  = dark_matter[2]


        #print 'Propergate: '+suffix
        #    if mDM!=100:
        #        'WARNING:  adjust the file with dN/dT data to your DM mass!'

        num             = 1e-50
        rr_grid         = np.sqrt(r_grid*r_grid + z_grid*z_grid)+num
        
        ######################################
        #   create q_D (E, r, z)
        ######################################


        Rs              = 20.
        rho_0           = 1./ (rSun/Rs * np.power(1 + rSun/Rs, 2))
        rho_grid        = 1./ (rr_grid       /Rs * np.power(1 + rr_grid  /Rs, 2)) * rhoSun / rho_0



        spec = 1.
        pref = 1.
        if injection=='vittino':
            dN_dTn_grid         = np.interp( Tn_grid, data_vittino_dbar[:,0] * mDM / 2., data_vittino_dbar[:,1] / mDM *2. )  * f_data_vittino_dbar_multip(mDM)/f_data_vittino_dbar_multip(100.)
        


        if injection=='analytic':
            spec        =   f_data_cirelli_pbar__bb( mDM, np.log10(Tn_vec/mDM) )[:,0]
            spec       *=   1./Tn_vec/np.log(10)
            pref        =   2*4./3.*pow(p_coal, 3)/np.sqrt(2*2*Tn_vec*Tn_vec+2*fMass_d*2*Tn_vec)*fMass_d/fMass_p/fMass_n
            dN_dTn_grid = np.interp( Tn_grid, Tn_vec, pref*spec*spec  )
        

        rho_q_D_grid    =  rho_grid*rho_grid/rhoSun/rhoSun


        ######################################
        #   create q_i (E, z)
        ######################################

        i_grid2     = i_grid [:,:,:,0]
        Tn_grid2    = Tn_grid[:,:,:,0]
        z_grid2     = z_grid [:,:,:,0]

        integrand_q_i_grid  = r_grid  *  sps.j0( vect_xi(i_grid)*r_grid/R ) * rho_q_D_grid

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

        nH  = 1.    *   1e6
        nHe = 0.1   *   1e6
        
        data_dpbar_in           = (data_ppbar_tot-data_ppbar_el)*data_dpbar_tot/data_ppbar_tot
        sigma_in_dbarp__grid2   = np.interp( Tn_grid2, data_dpbar_tot[:,0], data_dpbar_in[:,1]  ) * 1e-31 * np.power(A_sp/2., 2/3)
        Gamma_grid2             = (nH + np.power(4,2./3)*nHe)*sigma_in_dbarp__grid2*beta_grid2*2.99792e8

        A_i_grid2               =   V_c + 2*h*Gamma_grid2 + K_grid2*S_i_grid2  *  1./np.tanh(S_i_grid2*L/2.)

        integrand__y_i_grid2    = np.exp( V_c/2./K_grid2 * (L-z_grid2) )  * np.sinh(  S_i_grid2 * (L-z_grid2) /2.  ) * q_i_grid2

        y_i_grid3               = integrand__y_i_grid2.sum(axis=2) * dZ

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


        n_D = 0.5* sv *rhoSun*rhoSun/mDM/mDM * dN_dTn_grid[0,:,0,0] * R_D
        n_D = 1e6*n_D                                                           # Transform all units to GeV, m, and s

        ######################################
        #   calculate flux
        ######################################

        phi_D = n_D * beta_grid2[0,:,0] / 4. / math.pi * 2.99792e8


    n_integral = []

    phi=0.0
    Tn_vec_mod = Tn_vec - phi / A_sp
    phi_D_mod =  ( Tn_vec_mod*(Tn_vec_mod+2*fMass_p)  )/(  Tn_vec*(Tn_vec+2*fMass_p)  ) * phi_D
    T_int     = np.arange(0.1,0.25,0.0001)
    phi_int   = np.interp( T_int, Tn_vec_mod, phi_D_mod )
    n_Dbar    = 2*phi_int.sum() * 0.0001/(2e-6*0.15)
    n_integral.append(n_Dbar)

    phi=0.1
    Tn_vec_mod = Tn_vec - phi / A_sp
    phi_D_mod =  ( Tn_vec_mod*(Tn_vec_mod+2*fMass_p)  )/(  Tn_vec*(Tn_vec+2*fMass_p)  ) * phi_D
    T_int     = np.arange(0.1,0.25,0.0001)
    phi_int   = np.interp( T_int, Tn_vec_mod, phi_D_mod )
    n_Dbar    = 2*phi_int.sum() * 0.0001/(2e-6*0.15)
    n_integral.append(n_Dbar)

    phi=0.2
    Tn_vec_mod = Tn_vec - phi / A_sp
    phi_D_mod =  ( Tn_vec_mod*(Tn_vec_mod+2*fMass_p)  )/(  Tn_vec*(Tn_vec+2*fMass_p)  ) * phi_D
    T_int     = np.arange(0.1,0.25,0.0001)
    phi_int   = np.interp( T_int, Tn_vec_mod, phi_D_mod )
    n_Dbar    = 2*phi_int.sum() * 0.0001/(2e-6*0.15)
    n_integral.append(n_Dbar)

    phi=0.3
    Tn_vec_mod = Tn_vec - phi / A_sp
    phi_D_mod =  ( Tn_vec_mod*(Tn_vec_mod+2*fMass_p)  )/(  Tn_vec*(Tn_vec+2*fMass_p)  ) * phi_D
    T_int     = np.arange(0.1,0.25,0.0001)
    phi_int   = np.interp( T_int, Tn_vec_mod, phi_D_mod )
    n_Dbar    = 2*phi_int.sum() * 0.0001/(2e-6*0.15)
    n_integral.append(n_Dbar)

    phi=0.4
    Tn_vec_mod = Tn_vec - phi / A_sp
    phi_D_mod =  ( Tn_vec_mod*(Tn_vec_mod+2*fMass_p)  )/(  Tn_vec*(Tn_vec+2*fMass_p)  ) * phi_D
    T_int     = np.arange(0.1,0.25,0.0001)
    phi_int   = np.interp( T_int, Tn_vec_mod, phi_D_mod )
    n_Dbar    = 2*phi_int.sum() * 0.0001/(2e-6*0.15)
    n_integral.append(n_Dbar)

    phi=0.5
    Tn_vec_mod = Tn_vec - phi / A_sp
    phi_D_mod =  ( Tn_vec_mod*(Tn_vec_mod+2*fMass_p)  )/(  Tn_vec*(Tn_vec+2*fMass_p)  ) * phi_D
    T_int     = np.arange(0.1,0.25,0.0001)
    phi_int   = np.interp( T_int, Tn_vec_mod, phi_D_mod )
    n_Dbar    = 2*phi_int.sum() * 0.0001/(2e-6*0.15)
    n_integral.append(n_Dbar)


    if ret=='spectra':
        return (Tn_vec, phi_D), (Tn_vec_mod, phi_D_mod)
    if ret=='integral':
        return n_integral




data_bbbar_dm = np.genfromtxt('bbbar_DM/bbbar.txt')
x = np.argsort(data_bbbar_dm[:, 1])

data_bbbar_dm = data_bbbar_dm[x,::]

chiSq_bb    =   data_bbbar_dm[:, 1]
K_0_bb      =   data_bbbar_dm[:, 8]
V_c_bb      =   data_bbbar_dm[:,11]
delta_bb    =   data_bbbar_dm[:, 9]
L_bb        =   data_bbbar_dm[:,12]

m_DM_bb     =   np.power( 10, data_bbbar_dm[:,13]-3 )
sv_DM_bb    =   np.power( 10, data_bbbar_dm[:,14]   )



K_0_bb     *=   np.power(0.25,delta_bb)* np.power(3.24e-22,2)/(1e-6*3.171e-8)      # Transform from (galprop) D0_xx to K_0

m_DM        = []
n_particles = []

injection = 'analytic'
if step==4:
    injection = 'vittino'
if step==1 or step==4:
    not1sigma = True
    not2sigma = True
    not3sigma = True
    for i in range( 0, len(chiSq_bb) ):
        m_DM.append(m_DM_bb[i])
        propagation = ( 0.1,  K_0_bb[i],  V_c_bb[i],  delta_bb[i],  L_bb[i] )
        dm          = ( m_DM_bb[i], sv_DM_bb[i], 0.43 )
        n =  proper(ret='integral', propagation=propagation, dark_matter=dm, p_coal=0.078, injection=injection)
        n_particles.append(n)
        print str(i) + '   ' + str(chiSq_bb[i]) + '/' + str(chiSq_bb[0])
        print n_particles

        #case='CuKoKr'
        case='galprop'
        if step == 4:
            case = 'CuKoKr_vittino'

        if chiSq_bb[i]-chiSq_bb[0] > 2.3 and not1sigma:
            f = open(case+'_npar_1sigma.txt', 'w')
            s = '#' + str('m_DM/GeV').ljust(19) + ' ' + str('sm: phi = 0.0 GV').ljust(20) + ' ' + str('sm: phi = 0.1 GV').ljust(20) + ' ' + str('sm: phi = 0.2 GV').ljust(20) + ' ' + str('sm: phi = 0.3 GV').ljust(20) + ' ' + str('sm: phi = 0.4 GV').ljust(20) + ' ' + str('sm: phi = 0.5 GV').ljust(20)
            for j in range(len(m_DM)):
                s += '\n' + str(m_DM[j]).ljust(19)
                for k in range(len(n_particles[j])):
                    s += ' ' + str( (n_particles[j])[k]).ljust(20)
            f.write(s)
            f.close()
            m_DM        = []
            n_particles = []
            not1sigma = False
        if chiSq_bb[i]-chiSq_bb[0] > 6.18 and not2sigma:
            f = open(case+'_npar_2sigma.txt', 'w')
            s = '#' + str('m_DM/GeV').ljust(19) + ' ' + str('#').ljust(20)
            for j in range(len(m_DM)):
                s += '\n' + str(m_DM[j]).ljust(19)
                for k in range(len(n_particles[j])):
                    s += ' ' + str( (n_particles[j])[k]).ljust(20)
            f.write(s)
            f.close()
            m_DM        = []
            n_particles = []
            not2sigma = False
        if chiSq_bb[i]-chiSq_bb[0] > 11.83 and not3sigma:
            f = open(case+'_npar_3sigma.txt', 'w')
            s = '#' + str('m_DM/GeV').ljust(19) + ' ' + str('#').ljust(20)
            for j in range(len(m_DM)):
                s += '\n' + str(m_DM[j]).ljust(19)
                for k in range(len(n_particles[j])):
                    s += ' ' + str( (n_particles[j])[k]).ljust(20)
            f.write(s)
            f.close()
            not3sigma = False
            break
if step==2:
    not1sigma = True
    not2sigma = True
    not3sigma = True
    for i in range( 0, len(chiSq_bb) ):
        m_DM.append(m_DM_bb[i])
        #propagation = ( 0.2,  K_0_bb[i],  V_c_bb[i],  delta_bb[i],  L_bb[i] )
        dm          = ( m_DM_bb[i], sv_DM_bb[i], 0.43 )
        n =  proper(ret='integral', propagation='MED', dark_matter=dm, p_coal=0.078)
        n_particles.append(n)
        print str(i) + '   ' + str(chiSq_bb[i]) + '/' + str(chiSq_bb[0])
        case='MED'

        if chiSq_bb[i]-chiSq_bb[0] > 2.3 and not1sigma:
            f = open(case+'_npar_1sigma.txt', 'w')
            s = '#' + str('m_DM/GeV').ljust(19) + ' ' + str('#').ljust(20)
            for j in range(len(m_DM)):
                s += '\n' + str(m_DM[j]).ljust(19)
                for k in range(len(n_particles[j])):
                    s += ' ' + str( (n_particles[j])[k]).ljust(20)
            f.write(s)
            f.close()
            m_DM        = []
            n_particles = []
            not1sigma = False
        if chiSq_bb[i]-chiSq_bb[0] > 6.18 and not2sigma:
            f = open(case+'_npar_2sigma.txt', 'w')
            s = '#' + str('m_DM/GeV').ljust(19) + ' ' + str('#').ljust(20)
            for j in range(len(m_DM)):
                s += '\n' + str(m_DM[j]).ljust(19)
                for k in range(len(n_particles[j])):
                    s += ' ' + str( (n_particles[j])[k]).ljust(20)
            f.write(s)
            f.close()
            m_DM        = []
            n_particles = []
            not2sigma = False
        if chiSq_bb[i]-chiSq_bb[0] > 11.83 and not3sigma:
            f = open(case+'_npar_3sigma.txt', 'w')
            s = '#' + str('m_DM/GeV').ljust(19) + ' ' + str('#').ljust(20)
            for j in range(len(m_DM)):
                s += '\n' + str(m_DM[j]).ljust(19)
                for k in range(len(n_particles[j])):
                    s += ' ' + str( (n_particles[j])[k]).ljust(20)
            f.write(s)
            f.close()
            not3sigma = False
            break
if step==3:
    not1sigma = True
    not2sigma = True
    not3sigma = True
    for i in range( 0, len(chiSq_bb) ):
        m_DM.append(m_DM_bb[i])
        propagation = ( 0.1,  K_0_bb[i],  V_c_bb[i],  delta_bb[i],  L_bb[i] )
        dm          = ( m_DM_bb[i], sv_DM_bb[i], 0.43 )
        n =  proper(ret='integral', propagation='MAX', dark_matter=dm, p_coal=0.078)
        n_particles.append(n)
        print str(i) + '   ' + str(chiSq_bb[i]) + '/' + str(chiSq_bb[0])

        case='MAX'

        if chiSq_bb[i]-chiSq_bb[0] > 2.3 and not1sigma:
            f = open(case+'_npar_1sigma.txt', 'w')
            s = '#' + str('m_DM/GeV').ljust(19) + ' ' + str('#').ljust(20)
            for j in range(len(m_DM)):
                s += '\n' + str(m_DM[j]).ljust(19)
                for k in range(len(n_particles[j])):
                    s += ' ' + str( (n_particles[j])[k]).ljust(20)
            f.write(s)
            f.close()
            m_DM        = []
            n_particles = []
            not1sigma = False
        if chiSq_bb[i]-chiSq_bb[0] > 6.18 and not2sigma:
            f = open(case+'_npar_2sigma.txt', 'w')
            s = '#' + str('m_DM/GeV').ljust(19) + ' ' + str('#').ljust(20)
            for j in range(len(m_DM)):
                s += '\n' + str(m_DM[j]).ljust(19)
                for k in range(len(n_particles[j])):
                    s += ' ' + str( (n_particles[j])[k]).ljust(20)
            f.write(s)
            f.close()
            m_DM        = []
            n_particles = []
            not2sigma = False
        if chiSq_bb[i]-chiSq_bb[0] > 11.83 and not3sigma:
            f = open(case+'_npar_3sigma.txt', 'w')
            s = '#' + str('m_DM/GeV').ljust(19) + ' ' + str('#').ljust(20)
            for j in range(len(m_DM)):
                s += '\n' + str(m_DM[j]).ljust(19) + ' ' + str(n_particles[j]).ljust(20)
            f.write(s)
            f.close()
            not3sigma = False
            break

from matplotlib.ticker import FormatStrFormatter

def get_contour(x, y, n=10, t='cubic'):
    min = np.amin(x);
    max = np.amax(x);
    
    s    = (max-min)/n
    bins = np.arange(min, max+s+s/2., s )

    lower = np.ones(len(bins))*(+1e90)
    upper = np.ones(len(bins))*(-1e90)

    y_min = 1
    y_max = 1
    for i in range(len(lower)-2):
        for j in range(len(x)):
            if x[j] == min:
                y_min = y[j]
            if x[j] == max:
                y_max = y[j]
#            if y[j]>4:
#                print y[j]
            if x[j]>bins[i] and x[j]<bins[i+1]:
                if y[j] < lower[i+1]:
                    lower[i+1] = y[j]
                if y[j] > upper[i+1]:
                    upper[i+1] = y[j]


    x_r = bins-s/2
    x_r[ 0] = min
    x_r[-1] = max

    lower[ 0] = 0.333*(lower[ 1]+upper[ 1] + y_min)
    upper[ 0] = 0.333*(lower[ 1]+upper[ 1] + y_min)
    lower[-1] = 0.333*(lower[-2]+upper[-2] + y_max)
    upper[-1] = 0.333*(lower[-2]+upper[-2] + y_max)


    if t=='cubic':
        for i in range(3):
            for j in range(1, len(lower)-1 ):
                if lower[j]>(lower[j-1] + lower[j+1] )*0.5:
                    lower[j] = (lower[j-1] + lower[j+1] ) * 0.5
                if upper[j] < (upper[j-1] + upper[j+1] ) * 0.5 :
                    upper[j] = (upper[j-1] + upper[j+1] ) * 0.5
        x_r[ 0] *= 0.999
        x_r[-1] *= 1.001
        f_l = interp1d( x_r, np.log(lower), 'cubic' )
        f_u = interp1d( x_r, np.log(upper), 'cubic' )
        s = s/10
        x_r = np.arange(min, max+s/2, s)
        
#        print upper[0]
#        print lower[0]
#        print np.exp(f_l(x_r[0]))
#        print np.exp(f_u(x_r[0]))
#        print '****'

        return x_r, np.exp(f_l(x_r)), np.exp(f_u(x_r))

    return x_r, lower, upper



if step==5:
    not1sigma = True
    not2sigma = True
    not3sigma = True
    for i in range( 0, len(chiSq_bb) ):
        if chiSq_bb[i]-chiSq_bb[0] > 2.3 and not1sigma:
            print '1 sigma: i = ' + str(i)
            not1sigma = False
        if chiSq_bb[i]-chiSq_bb[0] > 6.18 and not2sigma:
            print '1 sigma: i = ' + str(i)
            not2sigma = False
        if chiSq_bb[i]-chiSq_bb[0] > 11.83 and not3sigma:
            print '1 sigma: i = ' + str(i)
            not3sigma = False


    print 'plot'
    MED_1       = np.genfromtxt('MED_npar_1sigma.txt')
    MAX_1       = np.genfromtxt('MAX_npar_1sigma.txt')
    CuKoKr_1    = np.genfromtxt('CuKoKr_npar_1sigma.txt')
    Vittino_1   = np.genfromtxt('CuKoKr_vittino_npar_1sigma.txt')

    MED_2       = np.genfromtxt('MED_npar_2sigma.txt')
    MAX_2       = np.genfromtxt('MAX_npar_2sigma.txt')
    CuKoKr_2    = np.genfromtxt('CuKoKr_npar_2sigma.txt')
    Vittino_2   = np.genfromtxt('CuKoKr_vittino_npar_2sigma.txt')

#    MED_3 = np.genfromtxt('MED_npar_3sigma.txt')
#    MAX_3 = np.genfromtxt('MAX_npar_3sigma.txt')
#    CuKoKr_3 = np.genfromtxt('CuKoKr_npar_3sigma.txt')

    plot, fig = plot_1D(r'$m_{DM}$ [GeV]', r'# $\bar{D}$', 'linear', 'log')

#    plot.scatter( MED_2   [:,0], MED_2   [:,1], c='blue' )
    m, l_n, u_n = get_contour(  MED_2[:,0],  MED_2[:,5]  )
    plot.fill_between(m, l_n, u_n, color='blue', alpha=0.2, label='MED')
    m, l_n, u_n = get_contour(  MAX_2[:,0],  MAX_2[:,5]  )
    plot.fill_between(m, l_n, u_n, color='green', alpha=0.2, label='MAX')
    m, l_n, u_n = get_contour(  CuKoKr_2[:,0],  CuKoKr_2[:,5]  )
    plot.fill_between(m, l_n, u_n, color='red', alpha=0.2, label='CuKoKr')


    plot.legend     (  loc='lower left', numpoints=1, bbox_to_anchor=(0.05, 0.05), frameon=False, fontsize=0.8*label_size  )
    plot.axhline( y=1, lw=2, c='black' )
    plot.grid(b=None, axis='y')
    plot.set_xlim( ( 50, 90) )
    plot.set_ylim( (0.2, 80) )
    plt.savefig('number_prop.pdf')


    plot, fig = plot_1D(r'$m_{DM}$ [GeV]', r'# $\bar{D}$', 'linear', 'log')

    m, l_n, u_n = get_contour(  CuKoKr_2[:,0],  CuKoKr_2[:,1]  )
    plot.fill_between(m, l_n, u_n, color='blue',  alpha=0.2, label=r'Solar Modulation, $\phi = 0\quad \mathrm{MV}$')
    m, l_n, u_n = get_contour(  CuKoKr_2[:,0],  CuKoKr_2[:,3]  )
    plot.fill_between(m, l_n, u_n, color='red',   alpha=0.2, label=r'Solar Modulation, $\phi = 200\quad \mathrm{MV}$')
    m, l_n, u_n = get_contour(  CuKoKr_2[:,0],  CuKoKr_2[:,5]  )
    plot.fill_between(m, l_n, u_n, color='green', alpha=0.2, label=r'Solar Modulation, $\phi = 400\quad \mathrm{MV}$')

    plot.legend     (  loc='upper right', numpoints=1, bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size  )
    plot.axhline( y=1, lw=2, c='black' )
    plot.grid(b=None, axis='y')
    plot.set_xlim( ( 50, 90) )
    plot.set_ylim( (0.2, 80) )
    plt.savefig('number_solarmod.pdf')


    plot, fig = plot_1D(r'$m_{DM}$ [GeV]', r'# $\bar{D}$', 'linear', 'log')
    
    m, l_n, u_n = get_contour(  CuKoKr_2[:,0],  CuKoKr_2[:,5]*np.power(0.3/0.43, 2)  )
    plot.fill_between(m, l_n, u_n, color='red',   alpha=0.2, label=r'Local DM density, $\rho_\odot = 0.30 \mathrm{GeV/cm^3}$')
    m, l_n, u_n = get_contour(  CuKoKr_2[:,0],  CuKoKr_2[:,5]*np.power(0.43/0.43, 2) )
    plot.fill_between(m, l_n, u_n, color='green', alpha=0.2, label=r'Local DM density, $\rho_\odot = 0.43 \mathrm{GeV/cm^3}$')
    m, l_n, u_n = get_contour(  CuKoKr_2[:,0],  CuKoKr_2[:,5]*np.power(0.6/0.43, 2)  )
    plot.fill_between(m, l_n, u_n, color='blue',  alpha=0.2, label=r'Local DM density, $\rho_\odot = 0.60\quad \mathrm{GeV/cm^3}$')
    
    plot.legend     (  loc='upper right', numpoints=1, bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size  )
    plot.axhline( y=1, lw=2, c='black' )
    plot.grid(b=None, axis='y')
    plot.set_xlim( ( 50, 90) )
    plot.set_ylim( (0.2, 80) )
    plt.savefig('number_DM_density.pdf')


    plot, fig = plot_1D(r'$m_{DM}$ [GeV]', r'# $\bar{D}$', 'linear', 'log')

    m, l_n, u_n = get_contour(  CuKoKr_2[:,0],  CuKoKr_2[:,5]*np.power(70./78, 3)  )
    plot.fill_between(m, l_n, u_n, color='red',   alpha=0.2, label=r'Coalescence momentum, $p_{\mathrm{coal}} = 70 \mathrm{MeV}$')
    m, l_n, u_n = get_contour(  CuKoKr_2[:,0],  CuKoKr_2[:,5]*np.power(78./78, 3) )
    plot.fill_between(m, l_n, u_n, color='green', alpha=0.2, label=r'Coalescence momentum, $p_{\mathrm{coal}} = 78 \mathrm{MeV}$')
    m, l_n, u_n = get_contour(  CuKoKr_2[:,0],  CuKoKr_2[:,5]*np.power(90./78, 3)  )
    plot.fill_between(m, l_n, u_n, color='blue',  alpha=0.2, label=r'Coalescence momentum, $p_{\mathrm{coal}} = 90 \mathrm{MeV}$')
    
    plot.legend     (  loc='upper right', numpoints=1, bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size  )
    plot.axhline( y=1, lw=2, c='black' )
    plot.grid(b=None, axis='y')
    plot.set_xlim( ( 50, 90) )
    plot.set_ylim( (0.2, 80) )
    plt.savefig('number_coalescence_momentum.pdf')



    plot, fig = plot_1D(r'$m_{DM}$ [GeV]', r'# $\bar{D}$', 'linear', 'log')
    
    m, l_n, u_n = get_contour(  CuKoKr_2[:,0],  CuKoKr_2[:,5] )
    plot.fill_between(m, l_n, u_n, color='red',   alpha=0.2, label=r'Analytic coalescence')
    m, l_n, u_n = get_contour(  Vittino_2[:,0],  Vittino_2[:,5]  )
    plot.fill_between(m, l_n, u_n, color='green', alpha=0.2, label=r'Monte Carlo coalescence')

    plot.legend     (  loc='upper right', numpoints=1, bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size  )
    plot.axhline( y=1, lw=2, c='black' )
    plot.grid(b=None, axis='y')
    plot.set_xlim( ( 50, 90) )
    plot.set_ylim( (0.2, 80) )
    plt.savefig('number_coalescence_model.pdf')








