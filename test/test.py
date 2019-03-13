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

for i, T in enumerate(list_T):
    Tpbar  = np.power(10, np.arange(-1,2.01,0.05) )

    xs_kdd18    = kdd18.dT_pp_pbar                   (T, Tpbar)[:,0]

    xs_wrapper        = dE_AA_pbar_LAB_incNbarAndHyperon   (T, Tpbar, A_projectile=1, N_projectile=0, A_target=1, N_target=0, parametrization='KORSMEIER_II') * 1e-31   # factor 1e-31, conversion from mbarn to m^2
    xs_wrapper_winkler= dE_AA_pbar_LAB_incNbarAndHyperon   (T, Tpbar, A_projectile=1, N_projectile=0, A_target=1, N_target=0, parametrization='WINKLER'     ) * 1e-31   # factor 1e-31, conversion from mbarn to m^2
    xs_wrapper_noHyp  = dE_AA_pbar_LAB                     (T, Tpbar, A_projectile=1, N_projectile=0, A_target=1, N_target=0, parametrization='KORSMEIER_II') * 1e-31   # factor 1e-31, conversion from mbarn to m^2

    plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'$d\sigma/dT_{\bar{p}} \;[\mathrm{m^2GeV^{-1}}]$', 'log', 'log' )
    plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)
    plot.plot(Tpbar, xs_kdd18,         c='#aa0000', lw=3,               label='KDD18 (table)' )
    plot.plot(Tpbar, xs_wrapper,       c='#0000aa', lw=3, dashes=(5,5), label='CRXS         ' )
    plot.plot(Tpbar, xs_wrapper_noHyp, c='#00aa00', lw=3, dashes=(2,2), label='CRXS, no Iso, no Hyp         ' )

    plot.plot(Tpbar, xs_wrapper_winkler, c='#00aaaa', lw=3, dashes=(10,2), label='CRXS, Winkler         ' )

    i_w = index_w[i]
    plot.plot(dat_pp_winkler[i_w:i_w+251, 1], dat_pp_winkler[i_w:i_w+251, 2]*1e-31, c='#aaaa00', lw=3, dashes=(5,2,2,2), label='Winkler (table)         ' )

    
    plot.legend( loc='upper left', bbox_to_anchor=(0.03, 0.97), frameon=False, ncol=1, fontsize=0.5*label_size )
    plot.text(0.7, 2e-34, r'$T_p=%.2f\;\mathrm{GeV}$'%T)
    plot.set_xlim( ( 1e-1,  1e2  ) )
    plot.set_ylim( ( 5e-35, 8e-32) )
    plt.savefig('profile_pp_T_%i.png'%T)

