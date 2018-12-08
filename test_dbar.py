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

_dEn_AA_Dbar_LAB = np.vectorize(XS.dEn_AA_Dbar_LAB)


for d in np.arange(1,3,0.1):
    E = np.power(10, d)+0.938;
    print( _dEn_AA_Dbar_LAB (E, 50+0.938, coalescence='FIXED_P0', parametrization='DI_MAURO_I') * 1e-31 )

T_list      = [10., 20., 40., 80., 160., 320., 640., 1280., 2540., 5120.]

Tn_Dbar_V   = np.power(10, np.arange(-0.9,2.01,0.05) )

plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'$d\sigma/dT_{\bar{p}} \;[\mathrm{m^2GeV^{-1}}]$', 'log', 'log' )
plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)

for T in T_list:
    #
    #   p + p -> pbar + X
    #

    xs_KorsmeierII__coal_fix   = _dEn_AA_Dbar_LAB (T, Tn_Dbar_V, coalescence='FIXED_P0'                  ) * 1e-31   # factor 1e-31, conversion from mbarn to m^2
    xs_diMauroI__coal_fix      = _dEn_AA_Dbar_LAB (T, Tn_Dbar_V, coalescence='FIXED_P0', parametrization='DI_MAURO_I') * 1e-31   # factor 1e-31, conversion from mbarn to m^2
    xs_KorsmeierII__coal_E_dep = _dEn_AA_Dbar_LAB (T, Tn_Dbar_V, coalescence='ENERGY_DEP__VAN_DOETINCHEM') * 1e-31   # factor 1e-31, conversion from mbarn to m^2

    plot.plot(Tn_Dbar_V, xs_KorsmeierII__coal_fix, lw=3, dashes=(   ), label=r'$T_p = %5.1f$ GeV' % T )
    #plot.plot(Tn_Dbar_V, xs_KorsmeierII__coal_fix,   c='#0000aa', lw=3, dashes=(   ), label=r'$p_c=80$ MeV fixed' )
    #plot.plot(Tn_Dbar_V, xs_diMauroI__coal_fix,      c='#00aa00', lw=3, dashes=(   ), label=r'$p_c=80$ MeV fixed, di Mauro' )
    #plot.plot(Tn_Dbar_V, xs_KorsmeierII__coal_E_dep, c='#aa0000', lw=3, dashes=(5,5), label=r'$p_c$ energy dependent' )
    #plot.legend( loc='upper right', bbox_to_anchor=(0.97, 0.97), frameon=False, ncol=1, fontsize=0.65*label_size )
    #plot.text(0.7, 2e-34, r'$T_p='+str(T)+r'\;\mathrm{GeV}$')
    plot.set_xlim( ( 5e-1,  5e2  ) )
    plot.set_ylim( ( 1e-38, 1e-35) )
plot.legend( loc='upper right', bbox_to_anchor=(0.97, 0.97), frameon=False, ncol=2, fontsize=0.45*label_size )
plt.savefig('profile_pp_dbar.png')


_dE_AA_pbar_LAB = np.vectorize(XS.dE_AA_pbar_LAB)

plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'$d\sigma/dT_{\bar{p}} \;[\mathrm{m^2GeV^{-1}}]$', 'log', 'log' )
plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)

for T in T_list:
    #
    #   p + p -> pbar + X
    #
    
    xs_KorsmeierII   = _dE_AA_pbar_LAB (T, Tn_Dbar_V ) * 1e-31   # factor 1e-31, conversion from mbarn to m^2
   
    plot.plot(Tn_Dbar_V, xs_KorsmeierII, lw=3, dashes=(   ), label=r'$T_p = %5.1f$ GeV' % T )
    #plot.plot(Tn_Dbar_V, xs_KorsmeierII__coal_fix,   c='#0000aa', lw=3, dashes=(   ), label=r'$p_c=80$ MeV fixed' )
    #plot.plot(Tn_Dbar_V, xs_diMauroI__coal_fix,      c='#00aa00', lw=3, dashes=(   ), label=r'$p_c=80$ MeV fixed, di Mauro' )
    #plot.plot(Tn_Dbar_V, xs_KorsmeierII__coal_E_dep, c='#aa0000', lw=3, dashes=(5,5), label=r'$p_c$ energy dependent' )
    #plot.legend( loc='upper right', bbox_to_anchor=(0.97, 0.97), frameon=False, ncol=1, fontsize=0.65*label_size )
    #plot.text(0.7, 2e-34, r'$T_p='+str(T)+r'\;\mathrm{GeV}$')
plot.set_xlim( ( 5e-1,  5e2  ) )
plot.set_ylim( ( 1e-34, 1e-31) )
plot.legend( loc='upper right', bbox_to_anchor=(0.97, 0.97), frameon=False, ncol=2, fontsize=0.45*label_size )
plt.savefig('profile_pp_pbar.png')





ams02_p = np.genfromtxt('proton_ams02.txt')
ams02_p__R  = np.sqrt(ams02_p[:,0]*ams02_p[:,1])
ams02_p__F  =         (ams02_p[:,2]*ams02_p[:,9])
ams02_p__Fe = np.sqrt(ams02_p[:,3]*ams02_p[:,3]+ams02_p[:,8]*ams02_p[:,8])*ams02_p[:,9]

from scipy.interpolate import interp1d

int_fun_ams_p = interp1d(ams02_p__R, ams02_p__R**2.8 * ams02_p__F, bounds_error=False, fill_value=ams02_p__R[-1]**2.8 * ams02_p__F[-1])

def fun_ams_p(R):
    return int_fun_ams_p(R)*np.power(R,-2.8)


plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'$f$', 'log', 'log' )
plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)
R = np.power(10,np.arange(0, 5, 0.1))
plot.plot(R, fun_ams_p(R))
plot.plot(ams02_p__R, ams02_p__F)
plt.savefig('flux_p.png')

Tn_Dbar= np.power(10,np.arange(0.4, 3.01, 0.1))
pref =  4*math.pi* 0.1e6
q = []
for Tn_d in Tn_Dbar:
    q0 = integrate.quad(lambda T: _dEn_AA_Dbar_LAB (T, Tn_d, coalescence='FIXED_P0', parametrization='DI_MAURO_I') * 1e-31 * fun_ams_p(T), 2, 1e5)
    q.append(q0[0])

q1 = []
for Tn_d in Tn_Dbar:
    q0 = integrate.quad(lambda logT: _dEn_AA_Dbar_LAB (np.exp(logT), Tn_d, coalescence='FIXED_P0', parametrization='DI_MAURO_I') * 1e-31 * fun_ams_p(np.exp(logT))*np.exp(logT), np.log(2), np.log(1e6))
    q1.append(q0[0])


plot, fig = plot_1D ( r'$T_{\bar{p}}\;[\mathrm{GeV}]$', r'$f$', 'log', 'log' )
plt.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.15)
plot.plot(Tn_Dbar, q*Tn_Dbar)
plot.plot(Tn_Dbar, q1*Tn_Dbar)
plt.savefig('q_dbar.png')


