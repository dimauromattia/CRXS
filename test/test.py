#! /usr/bin/env python3

import numpy                         as np

import CRXS.XS_wrapper               as XS
import KDD18.interpolate_Param_II_B  as kdd18


#
#   Use matplotlib to produce plots
#
import matplotlib.pyplot    as       plt
import matplotlib           as       mpl
import matplotlib.colors    as       colors
from   matplotlib           import   rc

#
#   Function to set some nice layout
#
print_size=10
label_size=35
def plot_1D( xlabel, ylabel, xscale='log', yscale='log', sizex=1.3, sizey=1.0):
    global print_size, label_size
    plt.close('all')
    font_props = {"size":label_size}
    rc("font", **font_props)
    mpl.rcParams['axes.linewidth'] = 2
    fig     = plt.figure(figsize=(print_size*sizex, print_size*sizey))
    plot    = plt.subplot2grid((1, 1), (0, 0), rowspan=1)
    plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
    plot.set_xlabel ( xlabel )
    plot.set_ylabel ( ylabel )
    plot.set_xscale ( xscale )
    plot.set_yscale ( yscale )
    plot.tick_params('both', length=20, width=2, which='major', top=True, right=True, direction='in', pad=10)
    plot.tick_params('both', length=10, width=1, which='minor', top=True, right=True, direction='in', pad=10)
    plot.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
    return plot, fig


dat_pp_winkler = np.genfromtxt('pp_cross.dat.txt')

dat_pHe_winkler = np.genfromtxt('pHe_cross.dat.txt')


#print(dat_pp_winkler[9789:9789+251, 0])
#print(dat_pp_winkler[9789:9789+251, 1])
#sys.exit(0)


##########################################################################
#   Plot a profile of the cross sections d sigma/d T (T=200 GeV, Tpbar)
##########################################################################

#
#   p + p -> pbar + X
#

dE_AA_pbar_LAB_incNbarAndHyperon = np.vectorize(XS.dE_AA_pbar_LAB_incNbarAndHyperon)
dE_AA_pbar_LAB                   = np.vectorize(XS.dE_AA_pbar_LAB                  )


list_T  = [19.95, 104.7, 199.5]

index_w = [3521-7, 8039-7, 9796-7 ]

#T      = 200
name = [ 'pp',  'pHe' ]
nuc  = [ [1,0], [4,2] ]


m_proton = 0.938

def deltaHyperon(s):
    
    C1 = 0.31
    C2 = 0.30
    C3 = 21316.
    C4 = 0.9

    factor   = 0.81

    hyperon = C1 + C2/(1+pow(C3/s,C4))
    hyperon       *= factor
    return hyperon

    
def deltaIsospin(s):

    C14 = 0.114
    C15 = 20736.
    C16 = 0.51

    return C14/(1+pow(s/C15, C16));



for i, n in enumerate(name):
    A_t = nuc[i][0]
    N_t = nuc[i][1]
    print(n)
    print(A_t)
    print(N_t)

    for i, T in enumerate(list_T):
        
        s = 4*m_proton**2 + 2 * T * m_proton;
        
        fac = 1 # 2 + 2 * deltaHyperon(s) + deltaIsospin(s)
        
        
        Tpbar  = np.power(10, np.arange(-1,2.01,0.05) )

        if A_t == 1:
            xs_kdd18    = kdd18.dT_pp_pbar                    (T, Tpbar)[:,0]
            dat_winkler = dat_pp_winkler
        if A_t == 4:
            xs_kdd18    = kdd18.dT_pHe_pbar                   (T, Tpbar)[:,0]
            dat_winkler = dat_pHe_winkler

        xs_wrapper        = dE_AA_pbar_LAB_incNbarAndHyperon   (T, Tpbar, A_projectile=1, N_projectile=0, A_target=A_t, N_target=N_t, parametrization='KORSMEIER_II') * 1e-31   # factor 1e-31, conversion from mbarn to m^2
        xs_wrapper_winkler= dE_AA_pbar_LAB_incNbarAndHyperon   (T, Tpbar, A_projectile=1, N_projectile=0, A_target=A_t, N_target=N_t, parametrization='WINKLER'     ) * 1e-31   # factor 1e-31, conversion from mbarn to m^2
        xs_wrapper_noHyp  = dE_AA_pbar_LAB                     (T, Tpbar, A_projectile=1, N_projectile=0, A_target=A_t, N_target=N_t, parametrization='KORSMEIER_II') * 1e-31   # factor 1e-31, conversion from mbarn to m^2

        plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'$d\sigma/dT_{\bar{p}} \;[\mathrm{m^2GeV^{-1}}]$', 'log', 'log' )
        plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)
        plot.plot(Tpbar, xs_kdd18,             c='#aa0000', lw=3,               label='KDD18 (table)' )
        plot.plot(Tpbar, xs_wrapper,           c='#0000aa', lw=3, dashes=(5,5), label='CRXS         ' )
        plot.plot(Tpbar, xs_wrapper_noHyp*fac, c='#00aa00', lw=3, dashes=(2,2), label='CRXS, no Iso, no Hyp         ' )

        plot.plot(Tpbar, xs_wrapper_winkler, c='#00aaaa', lw=3, dashes=(10,2), label='CRXS, Winkler         ' )

        i_w = index_w[i]
        plot.plot(dat_winkler[i_w:i_w+251, 1], dat_winkler[i_w:i_w+251, 2]*1e-31, c='#aaaa00', lw=3, dashes=(5,2,2,2), label='Winkler (table)         ' )

        
        plot.legend( loc='upper left', bbox_to_anchor=(0.03, 0.97), frameon=False, ncol=1, fontsize=0.5*label_size )
        plot.text(0.7, 2e-34, r'$T_p=%.2f\;\mathrm{GeV}$'%T)
        plot.set_xlim( ( 1e-1,  1e2  ) )
        plot.set_ylim( ( 5e-35, 8e-32) )
        plt.savefig('profile_'+n+'_T_%i.png'%T)


        plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'$d\sigma/dT_{\bar{p}} \;[\mathrm{m^2GeV^{-1}}]$', 'log', 'linear' )
        plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)
        plot.plot(Tpbar, xs_kdd18/xs_wrapper,         c='#aa0000', lw=3,               label='KDD18 (table/CRXS)' )
        
        #plot.plot(Tpbar, xs_wrapper_winkler, c='#00aaaa', lw=3, dashes=(10,2), label='CRXS, Winkler         ' )
        
        i_w = index_w[i]
        
        xs_wrapper_winkler_interp = np.interp(dat_winkler[i_w:i_w+251, 1], Tpbar, xs_wrapper_winkler)
        plot.plot(dat_winkler[i_w:i_w+251, 1], dat_winkler[i_w:i_w+251, 2]*1e-31/xs_wrapper_winkler_interp, c='#aaaa00', lw=3, dashes=(5,2,2,2), label='Winkler (table/CRXS) ' )
        
        
        plot.legend( loc='upper left', bbox_to_anchor=(0.03, 0.97), frameon=False, ncol=1, fontsize=0.5*label_size )
        plot.text(0.7, 2e-34, r'$T_p=%.2f\;\mathrm{GeV}$'%T)
        plot.set_xlim( ( 1e-1,  1e2  ) )
        plot.set_ylim( ( 0.8, 1.2) )
        plt.savefig('ratio_'+n+'_T_%i.png'%T)


