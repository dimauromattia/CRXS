#! /usr/bin/env python3

import numpy  as np
import scipy.integrate   as     integrate
import math

# Include the python wrapper for the cpp functions
import python.XS_wrapper as XS

#
#   Some examples how to call functions from the wrapper:
#

print('Evaluate: "XS.dEn_AA_Dbar_LAB(  100        , 3           )"  automaticaly assumes pp scattering and KORSMEIER_II parametrization')
print(            XS.dEn_AA_Dbar_LAB(  100        , 3           ) )

# You may also vectorize the function:
_dE_AA_pbar_LAB_incNbarAndHyperon = np.vectorize(XS.dE_AA_pbar_LAB_incNbarAndHyperon)


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
#   Plot a profile of the cross sections d sigma/d T (T=200 GeV, Tn_Dbar)
##########################################################################

T          = 6000.
Tn_Dbar_V  = np.power(10, np.arange(-0.9,2.91,0.05) )

_dEn_AA_Dbar_LAB = np.vectorize(XS.dEn_AA_Dbar_LAB)
#
#print(_dEn_AA_Dbar_LAB( [100., 10.]   , [30., 3.] ) )

#
#   p + p -> pbar + X
#

xs_KorsmeierII__coal_fix   = _dEn_AA_Dbar_LAB (T, Tn_Dbar_V, coalescence='FIXED_P0'                  ) * 1e-31   # factor 1e-31, conversion from mbarn to m^2
xs_diMauroI__coal_fix      = _dEn_AA_Dbar_LAB (T, Tn_Dbar_V, coalescence='FIXED_P0', parametrization='DI_MAURO_I') * 1e-31   # factor 1e-31, conversion from mbarn to m^2
xs_KorsmeierII__coal_E_dep = _dEn_AA_Dbar_LAB (T, Tn_Dbar_V, coalescence='ENERGY_DEP__VAN_DOETINCHEM') * 1e-31   # factor 1e-31, conversion from mbarn to m^2

plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'$d\sigma/dT_{\bar{p}} \;[\mathrm{m^2GeV^{-1}}]$', 'log', 'log' )
plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)
plot.plot(Tn_Dbar_V, xs_KorsmeierII__coal_fix,   c='#0000aa', lw=3, dashes=(   ), label=r'$p_c=80$ MeV fixed' )
plot.plot(Tn_Dbar_V, xs_diMauroI__coal_fix,      c='#00aa00', lw=3, dashes=(   ), label=r'$p_c=80$ MeV fixed, di Mauro' )
plot.plot(Tn_Dbar_V, xs_KorsmeierII__coal_E_dep, c='#aa0000', lw=3, dashes=(5,5), label=r'$p_c$ energy dependent' )
plot.legend( loc='upper right', bbox_to_anchor=(0.97, 0.97), frameon=False, ncol=1, fontsize=0.65*label_size )
#plot.text(0.7, 2e-34, r'$T_p='+str(T)+r'\;\mathrm{GeV}$')
plot.set_xlim( ( 5e-1,  5e2  ) )
plot.set_ylim( ( 1e-38, 1e-35) )
plt.savefig('profile_pp_dbar.png')

