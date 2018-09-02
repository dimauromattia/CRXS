#! /usr/bin/env python3

import numpy  as np
import scipy.integrate   as     integrate
import math

import python.XS_tools as XS

#
#   Some examples how to call XS_tools:
#
#   dE_AA_pbar_LAB_incNbarAndHyperon( Tn_proj_lab, T_pbar_lab, type )
#       - the function is vectorized
#       - if you use vectors for either T_p_lab or T_pbar_lab it returns a vector
#       - if you use vectors for both T_p_lab and T_pbar_lab it returns a matrix



print('Evaluate: XS.dE_AA_pbar_LAB_incNbarAndHyperon(  100        , 3" )  automaticaly assumes pp scattering and Korsmeier_II parametrization')
print(           XS.dE_AA_pbar_LAB_incNbarAndHyperon(  100        , 3  )  )
print('Evaluate: XS.dE_AA_pbar_LAB_incNbarAndHyperon( [10, 100]   , 3  )' )
print(           XS.dE_AA_pbar_LAB_incNbarAndHyperon( [10, 100]   , 3  )  )
print('Evaluate: XS.dE_AA_pbar_LAB_incNbarAndHyperon(  100        , [5,10] )' )
print(           XS.dE_AA_pbar_LAB_incNbarAndHyperon(  100        , [5,10] )  )
print('Evaluate: XS.dE_AA_pbar_LAB_incNbarAndHyperon( [100,50]    , [5,10] )' )
print(           XS.dE_AA_pbar_LAB_incNbarAndHyperon( [100,50]    , [5,10] )  )





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


##########################################################################
#   Plot a profile of the cross sections d sigma/d T (T=4000 GeV, Tpbar)
##########################################################################

T      = 4000
Tpbar  = np.power(10, np.arange(-0.9,2.91,0.05) )


#
#   p + p -> pbar + X
#

xs_tools_KorsmeierII    = XS.dE_AA_pbar_LAB_incNbarAndHyperon   (T, Tpbar, parametrization='Korsmeier_II') * 1e-31   # factor 1e-31, conversion from mbarn to m^2

plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'$d\sigma/dT_{\bar{p}} \;[\mathrm{m^2GeV^{-1}}]$', 'log', 'log' )
plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)
plot.plot(Tpbar, xs_tools_KorsmeierII, c='#0000aa', lw=3, dashes=(5,5), label='XS Tools     ' )
plot.legend( loc='upper right', bbox_to_anchor=(0.97, 0.97), frameon=False, ncol=1, fontsize=0.65*label_size )
plot.text(0.7, 2e-34, r'$T_p='+str(T)+r'\;\mathrm{GeV}$')
plot.set_xlim( ( 1e-1,  1e3  ) )
plot.set_ylim( ( 5e-35, 8e-32) )
plt.savefig('profile_pp.png')

xs_tools_Winkler       = 1e-31 * XS.dE_AA_pbar_LAB_incNbarAndHyperon   (T, Tpbar, parametrization='Winkler'     )
xs_tools_KorsmeierI    = 1e-31 * XS.dE_AA_pbar_LAB_incNbarAndHyperon   (T, Tpbar, parametrization='Korsmeier_I' )
xs_tools_diMauroI      = 1e-31 * XS.dE_AA_pbar_LAB                     (T, Tpbar, parametrization='diMauro_I'   ) * 2.3
xs_tools_diMauroII     = 1e-31 * XS.dE_AA_pbar_LAB                     (T, Tpbar, parametrization='diMauro_II'  ) * 2.3

plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'ratio', 'log', 'linear' )
plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)
plot.plot(Tpbar, xs_tools_KorsmeierI/xs_tools_KorsmeierII,  c='#aa0000', lw=3, label='KorsmeierI/KorsmeierII' )
plot.plot(Tpbar, xs_tools_Winkler   /xs_tools_KorsmeierII,  c='#00aa00', lw=3, label='Winkler/KorsmeierII'    )
plot.plot(Tpbar, xs_tools_diMauroI  /xs_tools_KorsmeierII,  c='#0000aa', lw=3, label='diMauroI/KormseierII'   )
plot.plot(Tpbar, xs_tools_diMauroII /xs_tools_KorsmeierII,  c='#990099', lw=3, label='diMauroII/KormseierII'  )
plot.legend( loc='upper right', bbox_to_anchor=(0.97, 0.97), frameon=False, ncol=1, fontsize=0.65*label_size  )
plot.set_xlim( ( 1e-1,  1e3  ) )
plot.set_ylim( ( 0.5, 1.5) )
plot.text(0.7, 0.7, r'$T_p='+str(T)+r'\;\mathrm{GeV}$')
plt.savefig('ratio_pp.png')


#
#   p + He -> pbar + X
#

xs_pHe_tools_KorsmeierII    = XS.dE_AA_pbar_LAB_incNbarAndHyperon   (T, Tpbar, A_projectile=1, N_projectile=0, A_target=4, N_target=2, parametrization='Korsmeier_II') * 1e-31

plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'$d\sigma/dT_{\bar{p}} \;[\mathrm{m^2GeV^{-1}}]$', 'log', 'log' )
plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)
plot.plot(Tpbar, xs_pHe_tools_KorsmeierII, c='#0000aa', lw=3, dashes=(5,5), label='XS Tools     ' )
plot.legend( loc='upper right', bbox_to_anchor=(0.97, 0.97), frameon=False, ncol=1, fontsize=0.65*label_size )
plot.text(0.7, 2e-33, r'$T_p='+str(T)+r'\;\mathrm{GeV}$')
plot.set_xlim( ( 1e-1,  1e3  ) )
plot.set_ylim( ( 5e-34, 8e-31) )
plt.savefig('profile_pHe.png')


#
#   He + p -> pbar + X
#

xs_Hep_tools_KorsmeierII    = XS.dE_AA_pbar_LAB_incNbarAndHyperon   (T, Tpbar, A_projectile=4, N_projectile=2, A_target=1, N_target=0, parametrization='Korsmeier_II') * 1e-31

plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'$d\sigma/dT_{\bar{p}} \;[\mathrm{m^2GeV^{-1}}]$', 'log', 'log' )
plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)
plot.plot(Tpbar, xs_Hep_tools_KorsmeierII, c='#0000aa', lw=3, dashes=(5,5), label='XS Tools     ' )
plot.legend( loc='upper right', bbox_to_anchor=(0.97, 0.97), frameon=False, ncol=1, fontsize=0.65*label_size )
plot.text(0.7, 2e-33, r'$T_{He}='+str(T)+r'\;\mathrm{GeV}$')
plot.set_xlim( ( 1e-1,  1e3  ) )
plot.set_ylim( ( 5e-34, 8e-31) )
plt.savefig('profile_Hep.png')
