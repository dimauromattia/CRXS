#! /usr/bin/env python

import glob
import math
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib import rc

c0='black'
c1='#0404B4'
c2='#B40404'
c3='#04B404'
c4='#0484B4'
c5='#B40484'
c6='#B48404'

c1_a = (0.015,0.015,0.703, 0.2)
c2_a = (0.703,0.015,0.015, 0.2)
c3_a = (0.015,0.703,0.015, 0.2)
c4_a = (0.015,0.535,0.703, 0.2)
c5_a = (0.703,0.015,0.535, 0.2)
c6_a = (0.703,0.703,0.015, 0.2)

c2_aa = (0.015,0.703,0.015, 0.1)

print_size=10
label_size=20

species = { 'proton':0, 'helium':1, 'carbon':2, 'nitrogen':3, 'oxygen':4, 'lithium':5, 'beryllium':6, 'boron':7, 'pbar':8}


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



param       = np.genfromtxt( 'power_law_param.txt',         skip_header=1 )
brokenPowerLaw_param = param;

T = np.arange(-1,3.3+0.01,0.01)
T = np.power(10,T);

T_long = np.arange(-1,5.2+0.01,0.01)
T_long = np.power(10,T_long);




def brokenPowerLaw(Tn, species):
    global brokenPowerLaw_param
    p = brokenPowerLaw_param[species, :]
    if Tn!=Tn:
        return 0
    if Tn<=0:
        return 0
    if Tn<=p[5]:
        delta21 = p[2]-p[1]
        ret =  p[0]  *  np.power(  Tn, -p[1] )   /   np.power( p[4], -p[1] ) * np.power(  np.power(Tn, 1./p[6])+np.power(p[4], 1./p[6]), -p[6]*delta21 )   / np.power( pow(2, p[6])*p[4], -delta21 )
        return  ret
    return     brokenPowerLaw(p[5], species) * pow( Tn, -p[3] ) / pow( p[5], -p[3] )



SM_potential = 0.6

# assume fluxes in T/n
def modulated_flux(Tn, species):
    
    mp = 0.938
    global SM_potential
    
    Z_over_A = 0.5
    if species==0 or species==4:
        Z_over_A = 1.
    Tn_LIS    = Tn + Z_over_A*SM_potential
    phi_LIS   = brokenPowerLaw(Tn_LIS, species)*np.power(Tn_LIS, -2.7)
    factor    = (Tn*Tn+2*Tn*mp)/(Tn_LIS*Tn_LIS+2*Tn_LIS*mp)
    
    return phi_LIS*factor*np.power(Tn, 2.7)



def proton_FIT(Tn):
    return brokenPowerLaw(Tn, species['proton'])
proton_FIT_v = np.vectorize(proton_FIT)
def proton_SM(Tn):
    return modulated_flux(Tn, species['proton'])
proton_SM_v = np.vectorize(proton_SM)


def helium_FIT(Tn):
    return brokenPowerLaw(Tn, species['helium'])
helium_FIT_v = np.vectorize(helium_FIT)
def helium_SM(Tn):
    return modulated_flux(Tn, species['helium'])
helium_SM_v = np.vectorize(helium_SM)



def carbon_FIT(Tn):
    return brokenPowerLaw(Tn, species['carbon'])
carbon_FIT_v = np.vectorize(carbon_FIT)
def carbon_SM(Tn):
    return modulated_flux(Tn, species['carbon'])
carbon_SM_v = np.vectorize(carbon_SM)



def oxygen_FIT(Tn):
    return brokenPowerLaw(Tn, species['oxygen'])
oxygen_FIT_v = np.vectorize(oxygen_FIT)
def oxygen_SM(Tn):
    return modulated_flux(Tn, species['oxygen'])
oxygen_SM_v = np.vectorize(oxygen_SM)



def antiproton_FIT(Tn):
    return brokenPowerLaw(Tn, species['pbar'])
antiproton_FIT_v = np.vectorize(antiproton_FIT)
def antiproton_SM(Tn):
    return modulated_flux(Tn, species['pbar'])
antiproton_SM_v = np.vectorize(antiproton_SM)


proton          = np.genfromtxt( 'proton.txt',          skip_header=0 )
helium          = np.genfromtxt( 'helium.txt',          skip_header=0 )
carbon          = np.genfromtxt( 'carbon.txt',          skip_header=0 )
oxygen          = np.genfromtxt( 'oxygen.txt',          skip_header=0 )

nitrogen        = np.genfromtxt( 'nitrogen.txt',        skip_header=0 )
lithium         = np.genfromtxt( 'lithium.txt',         skip_header=0 )
beryllium       = np.genfromtxt( 'beryllium.txt',       skip_header=0 )
boron           = np.genfromtxt( 'boron.txt',           skip_header=0 )

antiproton      = np.genfromtxt( 'antiproton.txt',      skip_header=0 )

proton_LIS      = np.genfromtxt( 'proton_LIS.txt',      skip_header=0 )
helium_LIS      = np.genfromtxt( 'helium_LIS.txt',      skip_header=0 )
carbon_LIS      = np.genfromtxt( 'carbon_LIS.txt',      skip_header=0 )
oxygen_LIS      = np.genfromtxt( 'oxygen_LIS.txt',      skip_header=0 )

antiproton_LIS  = np.genfromtxt( 'antiproton_LIS.txt',  skip_header=0 )


#proton_pamela          = np.genfromtxt( 'proton_pamela.txt',          skip_header=0 )
#helium_pamela          = np.genfromtxt( 'helium_pamela.txt',          skip_header=0 )
#antiproton_pamela      = np.genfromtxt( 'antiproton_pamela.txt',      skip_header=0 )
#proton_cream           = np.genfromtxt( 'proton_cream.txt',           skip_header=0 )
#helium_cream           = np.genfromtxt( 'helium_cream.txt',           skip_header=0 )

proton_fit      = np.genfromtxt( 'proton_fit.txt',          skip_header=0 )
helium_fit      = np.genfromtxt( 'helium_fit.txt',          skip_header=0 )
carbon_fit      = np.genfromtxt( 'carbon_fit.txt',          skip_header=0 )
nitrogen_fit    = np.genfromtxt( 'nitrogen_fit.txt',        skip_header=0 )
oxygen_fit      = np.genfromtxt( 'oxygen_fit.txt',          skip_header=0 )

lithium_fit     = np.genfromtxt( 'lithium_fit.txt',         skip_header=0 )
beryllium_fit   = np.genfromtxt( 'beryllium_fit.txt',       skip_header=0 )
boron_fit       = np.genfromtxt( 'boron_fit.txt',           skip_header=0 )


antiproton_fit  = np.genfromtxt( 'antiproton_fit.txt',      skip_header=0 )


############################################################
############################################################
############################################################


#plot, fig = plot_1D( r'$\mathrm{T/n \quad [GeV/n]}$', r'$\mathrm{\phi \cdot (T/n)^{2.7} \quad[(GeV/n)^{1.7} m^{-2} s^{-1} sr^{-1}] }$' )
#
#
#
#plot.errorbar( proton_LIS    [:,0], proton_LIS       [:,1], yerr=proton_LIS        [:,3], zorder=10,  fmt='o',  lw=1,  color='black'  , label=r'AMS-02'        )
#plot.plot(proton_fit[:,0],proton_fit[:,1] )
#
#plot.errorbar( helium_LIS    [:,0], helium_LIS       [:,1], yerr=helium_LIS        [:,3], zorder=10,  fmt='o',  lw=1,  color='black'  , label=r'AMS-02'        )
#plot.plot(helium_fit[:,0],helium_fit[:,1] )
#
#plot.errorbar( carbon_LIS    [:,0], carbon_LIS       [:,1], yerr=carbon_LIS        [:,3], zorder=10,  fmt='o',  lw=1,  color='black'  , label=r'AMS-02'        )
#plot.plot(carbon_fit[:,0],carbon_fit[:,1] )
#
#plot.errorbar( oxygen_LIS    [:,0], oxygen_LIS       [:,1], yerr=oxygen_LIS        [:,3], zorder=10,  fmt='o',  lw=1,  color='black'  , label=r'AMS-02'        )
#plot.plot(oxygen_fit[:,0],oxygen_fit[:,1] )
#
#
#plot.errorbar( antiproton_LIS    [:,0], antiproton_LIS       [:,1], yerr=antiproton_LIS        [:,3], zorder=10,  fmt='o',  lw=1,  color='black'  , label=r'AMS-02'        )
#plot.plot(antiproton_fit[:,0],antiproton_fit[:,1] )
#
#plot.set_xlim(5e-1, 5e4)
#plot.set_ylim(2e-1, 3e4)
#
#plt.savefig('proton.pdf')



################

plot, fig = plot_1D( r'$\mathrm{T/n \quad [GeV/n]}$', r'$\mathrm{\phi \cdot (T/n)^{2.7} \quad[(GeV/n)^{1.7} m^{-2} s^{-1} sr^{-1}] }$' )


plot.errorbar( proton        [:,0], proton           [:,1], yerr=proton        [:,3], zorder=10,  fmt='o',  lw=1,  color=c1, label=r'p'         )
plot.errorbar( helium        [:,0], helium           [:,1], yerr=helium        [:,3], zorder=10,  fmt='o',  lw=1,  color=c2, label=r'He'        )
plot.errorbar( carbon        [:,0], carbon           [:,1], yerr=carbon        [:,3], zorder=10,  fmt='o',  lw=1,  color=c3, label=r'C'         )
plot.errorbar( nitrogen      [:,0], nitrogen         [:,1], yerr=nitrogen      [:,3], zorder=10,  fmt='o',  lw=1,  color=c4, label=r'N'         )
plot.errorbar( oxygen        [:,0], oxygen           [:,1], yerr=oxygen        [:,3], zorder=10,  fmt='o',  lw=1,  color=c5, label=r'O'         )

plot.errorbar( lithium       [:,0], lithium          [:,1], yerr=lithium       [:,3], zorder=10,  fmt='o',  lw=1,  color=c3, label=r'Li'        )
plot.errorbar( beryllium     [:,0], beryllium        [:,1], yerr=beryllium     [:,3], zorder=10,  fmt='o',  lw=1,  color=c4, label=r'Be'        )
plot.errorbar( boron         [:,0], boron            [:,1], yerr=boron         [:,3], zorder=10,  fmt='o',  lw=1,  color=c5, label=r'B'         )

plot.errorbar( antiproton    [:,0], antiproton       [:,1], yerr=antiproton    [:,3], zorder=10,  fmt='o',  lw=1,  color=c0, label=r'$\bar{p}$' )

plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), ncol=1, frameon=False, fontsize=0.7*label_size)

plot.set_xlim(5e-1, 5e4)
#plot.set_ylim(2e-1, 3e4)

plt.savefig('fluxes.pdf')


plot.plot(proton_fit[:,0],      proton_fit[:,1],           c=c1)
plot.plot(helium_fit[:,0],      helium_fit[:,1],           c=c2 )
plot.plot(carbon_fit[:,0],      carbon_fit[:,1],           c=c3 )
plot.plot(nitrogen_fit[:,0],    nitrogen_fit[:,1],         c=c4 )
plot.plot(oxygen_fit[:,0],      oxygen_fit[:,1],           c=c5 )


plot.plot(lithium_fit[:,0],     lithium_fit[:,1],           c=c3 )
plot.plot(beryllium_fit[:,0],   beryllium_fit[:,1],         c=c4 )
plot.plot(boron_fit[:,0],       boron_fit[:,1],             c=c5 )


plot.plot(antiproton_fit[:,0],  antiproton_fit[:,1],       c=c0 )

plt.savefig('fluxes_fit.pdf')


#################

plot, fig = plot_1D( r'$\mathrm{T/n \quad [GeV/n]}$', r'$\mathrm{\phi \cdot (T/n)^{2.7} \quad[(GeV/n)^{1.7} m^{-2} s^{-1} sr^{-1}] }$' )


plot.errorbar( proton        [:,0], proton           [:,1], yerr=proton        [:,3], zorder=10,  fmt='o',  lw=1,  color=c1, label=r'p'         )
plot.errorbar( helium        [:,0], helium           [:,1], yerr=helium        [:,3], zorder=10,  fmt='o',  lw=1,  color=c2, label=r'He'        )

plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), ncol=1, frameon=False, fontsize=0.7*label_size)

plot.set_xlim(5e-1, 5e4)
plot.set_ylim(2e+1, 2e4)

plt.savefig('pHe.pdf')


plot.plot(proton_fit[:,0],      proton_fit[:,1],           c=c1)
plot.plot(helium_fit[:,0],      helium_fit[:,1],           c=c2 )

plt.savefig('pHe_fit.pdf')

################

plot, fig = plot_1D( r'$\mathrm{T/n \quad [GeV/n]}$', r'$\mathrm{\phi \cdot (T/n)^{2.7} \quad[(GeV/n)^{1.7} m^{-2} s^{-1} sr^{-1}] }$' )


plot.errorbar( carbon        [:,0], carbon           [:,1], yerr=carbon        [:,3], zorder=10,  fmt='o',  lw=1,  color=c3, label=r'C'         )
plot.errorbar( nitrogen      [:,0], nitrogen         [:,1], yerr=nitrogen      [:,3], zorder=10,  fmt='o',  lw=1,  color=c4, label=r'N'         )
plot.errorbar( oxygen        [:,0], oxygen           [:,1], yerr=oxygen        [:,3], zorder=10,  fmt='o',  lw=1,  color=c5, label=r'O'         )


plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), ncol=1, frameon=False, fontsize=0.7*label_size)

plot.set_xlim(5e-1, 5e4)
plot.set_ylim(2e-1, 5e1)

plt.savefig('CNO.pdf')


plot.plot(carbon_fit[:,0],      carbon_fit[:,1],           c=c3 )
plot.plot(nitrogen_fit[:,0],    nitrogen_fit[:,1],         c=c4 )
plot.plot(oxygen_fit[:,0],      oxygen_fit[:,1],           c=c5 )

plt.savefig('CNO_fit.pdf')

################

plot, fig = plot_1D( r'$\mathrm{T/n \quad [GeV/n]}$', r'$\mathrm{\phi \cdot (T/n)^{2.7} \quad[(GeV/n)^{1.7} m^{-2} s^{-1} sr^{-1}] }$' )

plot.errorbar( lithium       [:,0], lithium          [:,1], yerr=lithium       [:,3], zorder=10,  fmt='o',  lw=1,  color=c3, label=r'Li'        )
plot.errorbar( beryllium     [:,0], beryllium        [:,1], yerr=beryllium     [:,3], zorder=10,  fmt='o',  lw=1,  color=c4, label=r'Be'        )
plot.errorbar( boron         [:,0], boron            [:,1], yerr=boron         [:,3], zorder=10,  fmt='o',  lw=1,  color=c5, label=r'B'         )

plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), ncol=1, frameon=False, fontsize=0.7*label_size)

plot.set_xlim(5e-1, 5e4)
plot.set_ylim(2e-1, 1e1)

plt.savefig('LiBeB.pdf')


plot.plot(lithium_fit[:,0],     lithium_fit[:,1],           c=c3 )
plot.plot(beryllium_fit[:,0],   beryllium_fit[:,1],         c=c4 )
plot.plot(boron_fit[:,0],       boron_fit[:,1],             c=c5 )


plt.savefig('LiBeB_fit.pdf')


sys.exit(0)


############################################################
############################################################
############################################################

plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))
upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

upperPlt.errorbar( T, (proton_FIT_v(T)-proton_SM_v(T))/proton_FIT_v(T),  lw=2,  color='black'  ,   dashes=(     ), label=r'Proton'    )

upperPlt.errorbar( T, (helium_FIT_v(T)-helium_SM_v(T))/helium_FIT_v(T),  lw=2,  color='black'  ,   dashes=(15,5 ), label=r'Helium'    )

upperPlt.errorbar( T, (antiproton_FIT_v(T)-antiproton_SM_v(T))/antiproton_FIT_v(T),  lw=2,  color='black'  ,   dashes=(3,3  ), label=r'Antiproton'    )


upperPlt.set_xscale('log')
upperPlt.set_yscale('linear')

upperPlt.set_xlim(5e-1, 1e3)

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')

upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)

upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T/n \quad [GeV/n]}$')
upperPlt.set_ylabel(r'$\mathrm{(LIS-SM)/LIS}$')

upperPlt.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), frameon=False, fontsize=0.8*label_size)

plt.savefig('relative_SM.pdf')

upperPlt.set_ylim(0, 0.1)
plt.savefig('relative_SM_zoom.pdf')

############################################################

plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*0.9, print_size))
upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.15)

upperPlt.errorbar( antiproton[:,0], antiproton[:,1], yerr=antiproton[:,3], fmt='o',  lw=1,  color='black'  , label=r'AMS-02'                      )
upperPlt.set_xscale('log')
upperPlt.set_yscale('log')
upperPlt.set_xlim(5e-1, 1e3)
upperPlt.set_ylim(1e-3, 1e1)

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')

upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)


upperPlt.set_xlabel(r'$\mathrm{T \quad [GeV]}$')
upperPlt.set_ylabel(r'$\mathrm{\phi \cdot T^{2.7} \quad[(GeV)^{1.7} m^{-2} s^{-1} sr^{-1}]}$')

plt.savefig('pbar_flux.pdf')

upperPlt.set_ylim(5e-1, 1e1)
plt.savefig('pbar_flux_zoom.pdf')

plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))
upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

upperPlt.errorbar( antiproton_LIS[:,0], antiproton_LIS[:,3]/antiproton_LIS[:,1], fmt='o',  lw=1,  color='black'  , label=r''                      )

upperPlt.set_xscale('log')
upperPlt.set_yscale('linear')
upperPlt.set_xlim(5e-1, 1e3)
upperPlt.set_ylim(0, 0.3)

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')

upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)

upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T \quad [GeV]}$')
upperPlt.set_ylabel(r'$\mathrm{\sigma_{\phi}/\phi}$')

plt.savefig('relativ_pbar_uncertainty.pdf')



plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))

upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)

plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)


upperPlt.errorbar( proton        [:,0], proton    [:,1], yerr=proton        [:,3], fmt='o',  lw=1,  color='black'  , label=r'AMS-02 Measurement'                      )
#upperPlt.errorbar( proton_LIS    [:,0], proton_LIS[:,1], yerr=proton_LIS    [:,3], fmt='s',  lw=1,  color='black'  , label=r'LIS spectrum AMS-02 ($\phi=500\,\mathrm{MV}$)'   )
brokenPowerLaw_param = param_upper
protonUpper = proton_FIT_v(T)
brokenPowerLaw_param = param_lower
protonLower = proton_FIT_v(T)
upperPlt.fill_between( T, protonUpper, protonLower,   lw=2,  color='black'  , alpha=0.25, label=r''   )
brokenPowerLaw_param = param
upperPlt.errorbar( T, proton_FIT_v(T),   lw=2,  color='black'  ,   dashes=(     ), label=r'LIS Fux ($\phi=600_{-200}^{+100}\,\mathrm{MV}$)'   )
SM_potential = 0.6
upperPlt.errorbar( T, proton_SM_v(T),    lw=2,  color='black'  ,   dashes=(15,5), label=r''   )
upperPlt.text    ( 3e3,    proton_FIT_v(T)[-1]*0.93 , 'p',   color='black' )

upperPlt.errorbar( helium        [:,0], helium        [:,1], yerr=helium        [:,3], fmt='o',  lw=1,  color='black'  )
#upperPlt.errorbar( helium_LIS    [:,0], helium_LIS    [:,1], yerr=helium_LIS    [:,3], fmt='s',  lw=1,  color='black'  )
brokenPowerLaw_param = param_upper
Upper = helium_FIT_v(T)
brokenPowerLaw_param = param_lower
Lower = helium_FIT_v(T)
upperPlt.fill_between( T, Upper, Lower,       lw=2,  color='black'  , alpha=0.25, label=r''   )
brokenPowerLaw_param = param
upperPlt.errorbar( T, helium_FIT_v(T),   lw=2,  color='black'  ,   dashes=(     ), label=r''   )
SM_potential = 0.6
upperPlt.errorbar( T, helium_SM_v(T),    lw=2,  color='black'  ,   dashes=(15,5), label=r''   )
upperPlt.text    ( 3e3,    helium_FIT_v(T)[-1]*0.93 , 'He',   color='black' )



upperPlt.set_xscale('log')
upperPlt.set_yscale('log')
upperPlt.set_xlim(5e-1, 1e4)
upperPlt.set_ylim(5e1, 3e4)

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')

upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)

upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T/n \quad [GeV/n]}$')
upperPlt.set_ylabel(r'$\mathrm{\phi \cdot (T/n)^{2.7} \quad[(GeV/n)^{1.7} m^{-2} s^{-1} sr^{-1}] }$')

upperPlt.legend(loc='lower right', bbox_to_anchor=(0.95, 0.05), frameon=False, fontsize=0.8*label_size)

plt.savefig('lis_fluxes_pHe.pdf')




upperPlt.errorbar( antiproton    [:,0], antiproton    [:,1]*100, yerr=antiproton    [:,3]*100, fmt='o',  lw=1,  color='black'  )
#upperPlt.errorbar( antiproton_LIS[:,0], antiproton_LIS[:,1]*100, yerr=antiproton_LIS[:,3]*100, fmt='s',  lw=1,  color='magenta'  )
brokenPowerLaw_param = param_upper
Upper = antiproton_FIT_v(T)
brokenPowerLaw_param = param_lower
Lower = antiproton_FIT_v(T)
upperPlt.fill_between( T, Upper*100, Lower*100,   lw=2,  color='black'  , alpha=0.25, label=r''   )
brokenPowerLaw_param = param
upperPlt.errorbar( T, antiproton_FIT_v(T)*100,   lw=2,  color='black'  ,   dashes=(     ), label=r''   )
SM_potential = 0.6
upperPlt.errorbar( T, antiproton_SM_v(T)*100,    lw=2,  color='black'  ,   dashes=(15,5), label=r''   )
upperPlt.text    ( 2.2e3, antiproton_FIT_v(T)[-1]*100*0.93 , r'$\mathrm{100 \times \bar{p}}$',   color='black' )



upperPlt.set_xlim(5e-1, 1e4)
upperPlt.set_ylim(5e-0, 3e4)

plt.savefig('lis_fluxes_pHePbar.pdf')


plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))

upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)

plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)


upperPlt.errorbar( proton        [:,0], proton           [:,1], yerr=proton        [:,3], zorder=10,  fmt='o',  lw=1,  color='red'  , label=r'AMS-02'        )
upperPlt.errorbar( proton_pamela [:,0], proton_pamela    [:,1], yerr=proton_pamela [:,3], zorder= 9,  fmt='v',  lw=1,  color='blue'   , label=r'PAMELA'        )
upperPlt.errorbar( proton_cream  [:,0], proton_cream     [:,1], yerr=proton_cream  [:,3], zorder= 8,  fmt='s',  lw=1,  color='green'  , label=r'CREAM'         )
#upperPlt.errorbar( proton_LIS    [:,0], proton_LIS[:,1], yerr=proton_LIS    [:,3], fmt='s',  lw=1,  color='black'  , label=r'LIS spectrum AMS-02 ($\phi=500\,\mathrm{MV}$)'   )
brokenPowerLaw_param = param_upper
protonUpper = proton_FIT_v(T     )
brokenPowerLaw_param = param_lower
protonLower = proton_FIT_v(T     )
upperPlt.fill_between( T     , protonUpper, protonLower,   lw=2,  color='black'  , alpha=0.25, zorder= 0,   label=r''   )
brokenPowerLaw_param = param
upperPlt.errorbar( T_long, proton_FIT_v(T_long),   lw=2,       color='black'  ,   dashes=(3,3  ), zorder=1   )
upperPlt.errorbar( T     , proton_FIT_v(T     ),   lw=2,       color='black'  ,   dashes=(     ), zorder=1, label=r'LIS Fux ($\phi=600_{-200}^{+100}\,\mathrm{MV}$)' )
SM_potential = 0.6
upperPlt.errorbar( T     , proton_SM_v(T     ),    lw=2,       color='black'  ,   dashes=(15,5), zorder=2,  label=r''   )
upperPlt.text    ( 3e5, proton_FIT_v(T_long)[-1]*0.93 , 'p',   color='black'    )

upperPlt.errorbar( helium        [:,0], helium        [:,1], yerr=helium        [:,3], zorder=10,  fmt='o',  lw=1,  color='red'  )
upperPlt.errorbar( helium_pamela [:,0], helium_pamela [:,1], yerr=helium_pamela [:,3], zorder= 9,  fmt='v',  lw=1,  color='blue'   )
upperPlt.errorbar( helium_cream  [:,0], helium_cream  [:,1], yerr=helium_cream  [:,3], zorder= 8,  fmt='s',  lw=1,  color='green'  )
#upperPlt.errorbar( helium_LIS    [:,0], helium_LIS    [:,1], yerr=helium_LIS    [:,3], fmt='s',  lw=1,  color='black'  )
brokenPowerLaw_param = param_upper
Upper = helium_FIT_v(T     )
brokenPowerLaw_param = param_lower
Lower = helium_FIT_v(T     )
upperPlt.fill_between( T     , Upper, Lower,   lw=2,  color='black'  , alpha=0.25, zorder=1,   label=r''   )
brokenPowerLaw_param = param
upperPlt.errorbar( T_long, helium_FIT_v(T_long),   lw=2,  color='black'  ,   zorder=2,  dashes=(3,3  ), label=r''   )
upperPlt.errorbar( T     , helium_FIT_v(T     ),   lw=2,  color='black'  ,   zorder=2,  dashes=(     ), )
SM_potential = 0.6
upperPlt.errorbar( T     , helium_SM_v(T     ),    lw=2,  color='black'  ,   zorder=3,  dashes=(15,5), label=r''   )
upperPlt.text    ( 3e5, helium_FIT_v(T_long)[-1]*0.93 , 'He',   color='black' )



upperPlt.set_xscale('log')
upperPlt.set_yscale('log')

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')

upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)

upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T/n \quad [GeV/n]}$', fontsize=label_size*1.4 )
upperPlt.set_ylabel(r'$\mathrm{\phi \cdot (T/n)^{2.7} \quad[(GeV/n)^{1.7} m^{-2} s^{-1} sr^{-1}] }$', fontsize=label_size*1.4 )

upperPlt.legend(loc='lower right', numpoints=1, bbox_to_anchor=(0.75, 0.05), frameon=False, fontsize=0.8*label_size)



upperPlt.errorbar( antiproton        [:  ,0], antiproton       [:  ,1]*100, yerr=antiproton        [:  ,3]*100,zorder=10, fmt='o',  lw=1,  color='red'  )
upperPlt.errorbar( antiproton_pamela [:-1,0], antiproton_pamela[:-1,1]*100, yerr=antiproton_pamela [:-1,3]*100,zorder= 9, fmt='o',  lw=1,  color='blue'  )
brokenPowerLaw_param = param_upper
Upper = antiproton_FIT_v(T)
brokenPowerLaw_param = param_lower
Lower = antiproton_FIT_v(T)
upperPlt.fill_between( T, Upper*100, Lower*100,   lw=2,  color='black'  ,zorder= 1, alpha=0.25, label=r''   )
brokenPowerLaw_param = param
upperPlt.errorbar( T, antiproton_FIT_v(T)*100,   lw=2,  color='black'  ,zorder= 2,   dashes=(     ), label=r''   )
SM_potential = 0.6
upperPlt.errorbar( T, antiproton_SM_v(T)*100,    lw=2,  color='black'  ,zorder= 3,   dashes=(15,5), label=r''   )
upperPlt.text    ( 2.2e3, antiproton_FIT_v(T)[-1]*100*0.93 , r'$\mathrm{100 \times \bar{p}}$',   color='black' )



upperPlt.set_xlim(5e-1, 5e6)
upperPlt.set_ylim(2e-0, 3e4)

plt.savefig('lis_fluxes_pHePbar_AMS_CREAM_PAMELA.pdf')



plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))

upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)

plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)


upperPlt.errorbar( proton        [:,0], proton           [:,1], yerr=proton        [:,3], zorder=10,  fmt='o',  lw=1,  color=c1  , label=r'AMS-02'        )
upperPlt.errorbar( proton_cream  [:,0], proton_cream     [:,1], yerr=proton_cream  [:,3], zorder= 8,  fmt='s',  lw=1,  color=c1  , label=r'CREAM'         )
brokenPowerLaw_param = param_upper
protonUpper = proton_FIT_v(T     )
brokenPowerLaw_param = param_lower
protonLower = proton_FIT_v(T     )
upperPlt.fill_between( T     , protonUpper, protonLower,   lw=2,  color=c1  , alpha=0.25, zorder= 0,   label=r''   )
brokenPowerLaw_param = param
upperPlt.errorbar( T_long, proton_FIT_v(T_long),   lw=2,       color=c1  ,   dashes=(3,3  ), zorder=1   )
upperPlt.errorbar( T     , proton_FIT_v(T     ),   lw=2,       color=c1  ,   dashes=(     ), zorder=1, label=r'LIS Fux ($\phi=600_{-200}^{+100}\,\mathrm{MV}$)' )
SM_potential = 0.6
upperPlt.errorbar( T     , proton_SM_v(T     ),    lw=2,       color=c1  ,   dashes=(15,5), zorder=2,  label=r''   )
upperPlt.text    ( 3e5, proton_FIT_v(T_long)[-1]*0.93 , 'p',   color=c1    )

upperPlt.errorbar( helium        [:,0], helium        [:,1], yerr=helium        [:,3], zorder=10,  fmt='o',  lw=1,  color=c1  )
upperPlt.errorbar( helium_cream  [:,0], helium_cream  [:,1], yerr=helium_cream  [:,3], zorder= 8,  fmt='s',  lw=1,  color=c1  )
#upperPlt.errorbar( helium_LIS    [:,0], helium_LIS    [:,1], yerr=helium_LIS    [:,3], fmt='s',  lw=1,  color='black'  )
brokenPowerLaw_param = param_upper
Upper = helium_FIT_v(T     )
brokenPowerLaw_param = param_lower
Lower = helium_FIT_v(T     )
upperPlt.fill_between( T     , Upper, Lower,   lw=2,  color=c1  , alpha=0.25, zorder=1,   label=r''   )
brokenPowerLaw_param = param
upperPlt.errorbar( T_long, helium_FIT_v(T_long),   lw=2,  color=c1  ,   zorder=2,  dashes=(3,3  ), label=r''   )
upperPlt.errorbar( T     , helium_FIT_v(T     ),   lw=2,  color=c1  ,   zorder=2,  dashes=(     ), )
SM_potential = 0.6
upperPlt.errorbar( T     , helium_SM_v(T     ),    lw=2,  color=c1  ,   zorder=3,  dashes=(15,5), label=r''   )
upperPlt.text    ( 3e5, helium_FIT_v(T_long)[-1]*0.93 , 'He',   color=c1 )



upperPlt.set_xscale('log')
upperPlt.set_yscale('log')

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')

upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)

upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T/n \quad [GeV/n]}$')
upperPlt.set_ylabel(r'$\mathrm{\phi \cdot (T/n)^{2.7} \quad[(GeV/n)^{1.7} m^{-2} s^{-1} sr^{-1}] }$')

upperPlt.legend(loc='lower right', bbox_to_anchor=(0.65, 0.05), frameon=False, fontsize=0.8*label_size)


upperPlt.set_xlim(5e-1, 5e6)
upperPlt.set_ylim(3e+1, 3e4)

plt.savefig('lis_fluxes_pHe_TALK__AMS_CREAM.pdf')





######################


plt.close('all')
font_props = {"size":label_size}
rc("font", **font_props)
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['mathtext.fontset']='stixsans'
fig = plt.figure(figsize=(print_size*1.5, print_size))

upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)

plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)


upperPlt.errorbar( proton        [:,0], proton    [:,1], yerr=proton        [:,3], fmt='o',  lw=1,  color='black'  , label=r'AMS-02 Measurement'                      )
#upperPlt.errorbar( proton_LIS    [:,0], proton_LIS[:,1], yerr=proton_LIS    [:,3], fmt='s',  lw=1,  color='black'  , label=r'LIS spectrum AMS-02 ($\phi=500\,\mathrm{MV}$)'   )
brokenPowerLaw_param = param_upper
protonUpper = proton_FIT_v(T)
brokenPowerLaw_param = param_lower
protonLower = proton_FIT_v(T)
upperPlt.fill_between( T, protonUpper, protonLower,   lw=2,  color='black'  , alpha=0.25, label=r''   )
brokenPowerLaw_param = param
upperPlt.errorbar( T, proton_FIT_v(T),   lw=2,  color='black'  ,   dashes=(     ), label=r'LIS Fux ($\phi=600_{-200}^{+100}\,\mathrm{MV}$)'   )
SM_potential = 0.6
upperPlt.errorbar( T, proton_SM_v(T),    lw=2,  color='black'  ,   dashes=(15,5), label=r''   )
upperPlt.text    ( 3e3, proton_FIT_v(T)[-1]*0.93 , 'p',   color='black' )

upperPlt.errorbar( helium        [:,0], helium        [:,1], yerr=helium        [:,3], fmt='o',  lw=1,  color='black'  )
#upperPlt.errorbar( helium_LIS    [:,0], helium_LIS    [:,1], yerr=helium_LIS    [:,3], fmt='s',  lw=1,  color='black'  )
brokenPowerLaw_param = param_upper
Upper = helium_FIT_v(T)
brokenPowerLaw_param = param_lower
Lower = helium_FIT_v(T)
upperPlt.fill_between( T, Upper, Lower,   lw=2,  color='black'  , alpha=0.25, label=r''   )
brokenPowerLaw_param = param
upperPlt.errorbar( T, helium_FIT_v(T),   lw=2,  color='black'  ,   dashes=(     ), label=r''   )
SM_potential = 0.6
upperPlt.errorbar( T, helium_SM_v(T),    lw=2,  color='black'  ,   dashes=(15,5), label=r''   )
upperPlt.text    ( 3e3, helium_FIT_v(T)[-1]*0.93 , 'He',   color='black' )


upperPlt.errorbar( antiproton    [:,0], antiproton    [:,1], yerr=antiproton    [:,3], fmt='o',  lw=1,  color='black'  )
#upperPlt.errorbar( antiproton_LIS[:,0], antiproton_LIS[:,1]*100, yerr=antiproton_LIS[:,3]*100, fmt='s',  lw=1,  color='magenta'  )
brokenPowerLaw_param = param_upper
Upper = antiproton_FIT_v(T)
brokenPowerLaw_param = param_lower
Lower = antiproton_FIT_v(T)
upperPlt.fill_between( T, Upper, Lower,   lw=2,  color='black'  , alpha=0.25, label=r''   )
brokenPowerLaw_param = param
upperPlt.errorbar( T, antiproton_FIT_v(T),   lw=2,  color='black'  ,   dashes=(     ), label=r''   )
SM_potential = 0.6
upperPlt.errorbar( T, antiproton_SM_v(T),    lw=2,  color='black'  ,   dashes=(15,5), label=r''   )
upperPlt.text    ( 3e3, antiproton_FIT_v(T)[-1]*0.93 , r'$\mathrm{\bar{p}}$',   color='black' )



upperPlt.errorbar( carbon        [:,0], carbon        [:,1]*5, yerr=carbon        [:,3]*5, fmt='o',  lw=1,  color='black'  )
#upperPlt.errorbar( carbon_LIS    [:,0], carbon_LIS    [:,1]*5, yerr=carbon_LIS    [:,3]*5, fmt='s',  lw=1,  color='black'  )
brokenPowerLaw_param = param_upper
Upper = carbon_FIT_v(T)
brokenPowerLaw_param = param_lower
Lower = carbon_FIT_v(T)
upperPlt.fill_between( T, Upper*5, Lower*5,   lw=2,  color='black'  , alpha=0.25, label=r''   )
brokenPowerLaw_param = param
upperPlt.errorbar( T, carbon_FIT_v(T)*5,   lw=2,  color='black'  ,   dashes=(     ), label=r''   )
SM_potential = 0.6
upperPlt.errorbar( T, carbon_SM_v(T)*5,    lw=2,  color='black'  ,   dashes=(15,5),  label=r''   )
upperPlt.text    ( 2.2e3, carbon_FIT_v(T)[-1]*5  , r'$\times 5$ C',   color='black' )


upperPlt.errorbar( oxygen        [:,0], oxygen        [:,1]*0.5, yerr=oxygen        [:,3]*0.5, fmt='o',  lw=1,  color='black'  )
#upperPlt.errorbar( oxygen_LIS    [:,0], oxygen_LIS    [:,1]*0.5, yerr=oxygen_LIS    [:,3]*0.5, fmt='s',  lw=1,  color='green'  )
brokenPowerLaw_param = param_upper
Upper = oxygen_FIT_v(T)
brokenPowerLaw_param = param_lower
Lower = oxygen_FIT_v(T)
upperPlt.fill_between( T, Upper*0.5, Lower*0.5,   lw=2,  color='black'  , alpha=0.25, label=r''   )
brokenPowerLaw_param = param
upperPlt.errorbar( T, oxygen_FIT_v(T)*0.5,   lw=2,  color='black'  ,   dashes=(     ), label=r''   )
SM_potential = 0.6
upperPlt.errorbar( T, oxygen_SM_v(T)*0.5,    lw=2,  color='black'  ,   dashes=(15,5),  label=r''   )

upperPlt.text    ( 2.2e3, oxygen_FIT_v(T)[-1]*0.5 , r'$\times 0.5$ O',   color='black' )



upperPlt.set_xscale('log')
upperPlt.set_yscale('log')
upperPlt.set_xlim(5e-1, 1e4)
upperPlt.set_ylim(5e-1, 2e5)

upperPlt.tick_params('both', length=20, width=2, which='major')
upperPlt.tick_params('both', length=10, width=1, which='minor')

upperPlt.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)

upperPlt.tick_params(axis='both', pad=10)

upperPlt.set_xlabel(r'$\mathrm{T/n \quad [GeV/n]}$')
upperPlt.set_ylabel(r'$\mathrm{\phi \cdot (T/n)^{2.7} \quad[(GeV/n)^{1.7} m^{-2} s^{-1} sr^{-1}] }$')

upperPlt.legend(loc='upper left', numpoints=1,  bbox_to_anchor=(0.05, 0.95), frameon=False, fontsize=0.8*label_size, ncol=2)

plt.savefig('lis_fluxes_all.pdf')




#
#
#plt.close('all')
#font_props = {"size":label_size}
#rc("font", **font_props)
#mpl.rcParams['axes.linewidth'] = 2
#mpl.rcParams['mathtext.fontset']='stixsans'
#fig = plt.figure(figsize=(print_size*1.5, print_size))
#
#upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
#
#plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
#
#upperPlt.errorbar( proton        [:,0], proton    [:,1]*1000, yerr=proton        [:,3], fmt='o',  lw=1,  color='grey'  , label=r'Spectrum AMS-02'                        )
#upperPlt.errorbar( proton_LIS    [:,0], proton_LIS[:,1]*1000, yerr=proton_LIS    [:,3], fmt='s',  lw=1,  color='grey'  , label=r'LIS spectrum AMS-02 ($\phi=500\,\mathrm{MV}$)'   )
#upperPlt.errorbar( T, proton_FIT_v(T)*100,   lw=2,  color='grey'  , alpha=0.75,   dashes=(     ), label=r'Fit'   )
#
#upperPlt.errorbar( proton        [:,0], proton        [:,1], yerr=proton        [:,3], fmt='o',  lw=1,  color='blue'  )
#upperPlt.errorbar( proton_LIS    [:,0], proton_LIS    [:,1], yerr=proton_LIS    [:,3], fmt='s',  lw=1,  color='blue'  )
#upperPlt.errorbar( T,   proton_FIT_v(T)     , lw=2,  color='blue'  , alpha=0.75,   dashes=(     ) )
#upperPlt.text    ( 3e3, proton_FIT_v(T)[-1] , 'p',   color='blue' )
#
#upperPlt.errorbar( helium        [:,0], helium        [:,1], yerr=helium        [:,3], fmt='o',  lw=1,  color='red'  )
#upperPlt.errorbar( helium_LIS    [:,0], helium_LIS    [:,1], yerr=helium_LIS    [:,3], fmt='s',  lw=1,  color='red'  )
#upperPlt.errorbar( T,   helium_FIT_v(T)     , lw=2,  color='red'  , alpha=0.75,   dashes=(     ) )
#upperPlt.text    ( 3e3, helium_FIT_v(T)[-1] , 'He',  color='red' )
#
#upperPlt.errorbar( antiproton    [:,0], antiproton    [:,1]*10, yerr=antiproton    [:,3]*10, fmt='o',  lw=1,  color='magenta'  )
#upperPlt.errorbar( antiproton_LIS[:,0], antiproton_LIS[:,1]*10, yerr=antiproton_LIS[:,3]*10, fmt='s',  lw=1,  color='magenta'  )
#upperPlt.errorbar( T,   antiproton_FIT_v(T)*10     , lw=2,  color='magenta'  , alpha=0.75,   dashes=(     ) )
#upperPlt.text    ( 3e3, antiproton_FIT_v(T)[-1]*10 , r'$\mathrm{\bar{p}} \quad \times 10$',   color='magenta' )
#
#
#upperPlt.set_xscale('log')
#upperPlt.set_yscale('log')
#upperPlt.set_xlim(5e-1, 1e4)
#upperPlt.set_ylim(5e-1, 3e4)
#
#upperPlt.tick_params('both', length=20, width=2, which='major')
#upperPlt.tick_params('both', length=10, width=1, which='minor')
#
#upperPlt.grid(b=True, which='major', alpha=0.5, linestyle='-', linewidth=1)
#
#upperPlt.tick_params(axis='both', pad=10)
#
#upperPlt.set_xlabel(r'$\mathrm{T/n [GeV/n]}$')
#upperPlt.set_ylabel(r'$\mathrm{\phi (T/n)^{2.7} [(GeV/n)^{1.7} m^{-2} s^{-1} sr^{-1}] }$')
#
#upperPlt.legend(loc='lower right', fontsize=20)
#
#plt.savefig('lis_fluxes_pHePbar.pdf')
#
#
#
#plt.close('all')
#font_props = {"size":label_size}
#rc("font", **font_props)
#mpl.rcParams['axes.linewidth'] = 2
#mpl.rcParams['mathtext.fontset']='stixsans'
#fig = plt.figure(figsize=(print_size*1.5, print_size))
#
#upperPlt = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
#
#plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
#
#upperPlt.errorbar( proton        [:,0], proton    [:,1]*1000, yerr=proton        [:,3], fmt='o',  lw=1,  color='grey'  , label=r'Spectrum AMS-02'                        )
#upperPlt.errorbar( proton_LIS    [:,0], proton_LIS[:,1]*1000, yerr=proton_LIS    [:,3], fmt='s',  lw=1,  color='grey'  , label=r'LIS spectrum AMS-02 ($\phi=500\,\mathrm{MV}$)'   )
#upperPlt.errorbar( T, proton_FIT_v(T)*100,   lw=2,  color='grey'  , alpha=0.75,   dashes=(     ), label=r'Fit'   )
#
#
#upperPlt.errorbar( proton        [:,0], proton        [:,1], yerr=proton        [:,3], fmt='o',  lw=1,  color='blue'  )
#upperPlt.errorbar( proton_LIS    [:,0], proton_LIS    [:,1], yerr=proton_LIS    [:,3], fmt='s',  lw=1,  color='blue'  )
#upperPlt.errorbar( T,   proton_FIT_v(T)     , lw=2,  color='blue'  , alpha=0.75,   dashes=(     ) )
#upperPlt.text    ( 3e3, proton_FIT_v(T)[-1] , 'p',   color='blue' )
#
#
#
#upperPlt.errorbar( helium        [:,0], helium        [:,1], yerr=helium        [:,3], fmt='o',  lw=1,  color='red'  )
#upperPlt.errorbar( helium_LIS    [:,0], helium_LIS    [:,1], yerr=helium_LIS    [:,3], fmt='s',  lw=1,  color='red'  )
#upperPlt.errorbar( T,   helium_FIT_v(T)     , lw=2,  color='red'  , alpha=0.75,   dashes=(     ) )
#upperPlt.text    ( 3e3, helium_FIT_v(T)[-1] , 'He',  color='red' )
#
#
#upperPlt.errorbar( carbon        [:,0], carbon        [:,1]*10, yerr=carbon        [:,3]*10, fmt='o',  lw=1,  color='black'  )
#upperPlt.errorbar( carbon_LIS    [:,0], carbon_LIS    [:,1]*10, yerr=carbon_LIS    [:,3]*10, fmt='s',  lw=1,  color='black'  )
#upperPlt.errorbar( T,   carbon_FIT_v(T)*10     , lw=2,  color='black'  , alpha=0.75,   dashes=(     ) )
#upperPlt.text    ( 3e3, carbon_FIT_v(T)[-1]*10  , r'C $\times 10$',   color='black' )
#
#
#upperPlt.errorbar( oxygen        [:,0], oxygen        [:,1]*3, yerr=oxygen        [:,3]*3, fmt='o',  lw=1,  color='green'  )
#upperPlt.errorbar( oxygen_LIS    [:,0], oxygen_LIS    [:,1]*3, yerr=oxygen_LIS    [:,3]*3, fmt='s',  lw=1,  color='green'  )
#upperPlt.errorbar( T,   oxygen_FIT_v(T)*3     , lw=2,  color='green'  , alpha=0.75,   dashes=(     ) )
#upperPlt.text    ( 3e3, oxygen_FIT_v(T)[-1]*3 , r'O $\times 3$',   color='green' )
#
#
#upperPlt.errorbar( antiproton    [:,0], antiproton    [:,1]*10, yerr=antiproton    [:,3]*10, fmt='o',  lw=1,  color='magenta'  )
#upperPlt.errorbar( antiproton_LIS[:,0], antiproton_LIS[:,1]*10, yerr=antiproton_LIS[:,3]*10, fmt='s',  lw=1,  color='magenta'  )
#upperPlt.errorbar( T,   antiproton_FIT_v(T)*10     , lw=2,  color='magenta'  , alpha=0.75,   dashes=(     ) )
#upperPlt.text    ( 3e3, antiproton_FIT_v(T)[-1]*10 , r'$\mathrm{\bar{p}} \quad \times 10$',   color='magenta' )
#
#
#upperPlt.set_xscale('log')
#upperPlt.set_yscale('log')
#upperPlt.set_xlim(5e-1, 1e4)
#upperPlt.set_ylim(5e-1, 3e4)
#
#upperPlt.tick_params('both', length=20, width=2, which='major')
#upperPlt.tick_params('both', length=10, width=1, which='minor')
#
#upperPlt.grid(b=True, which='major', alpha=0.5, linestyle='-', linewidth=1)
#
#upperPlt.tick_params(axis='both', pad=10)
#
#upperPlt.set_xlabel(r'$\mathrm{T/n [GeV/n]}$')
#upperPlt.set_ylabel(r'$\mathrm{\phi (T/n)^{2.7} [(GeV/n)^{1.7} m^{-2} s^{-1} sr^{-1}] }$')
#
#upperPlt.legend(loc='lower right', fontsize=20)
#
#plt.savefig('lis_fluxes.pdf')
#
#
#
#
#
