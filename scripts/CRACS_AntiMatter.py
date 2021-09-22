#! /usr/bin/env python

import glob
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib import rc



import argparse
from PlotFunctions import plot2D


parser = argparse.ArgumentParser(description='Michael Korsmeier - Plotting script for source terms.'        )
parser.add_argument('--file',  help='Filename of table.', action='store', dest='file', type=str, default='' )
parser.add_argument('--label', help='Label.', action='store', dest='label', type=str, default=''            )
args = parser.parse_args()
file   	    = args.file
label 	    = args.label


c0='black'
c2='#04B404'
c4='#B40486'
c3='#0489B1'
c1='#0404B4'
c5='#B40404'

print_size=10
label_size=30

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

def profile(file, Tn, type='prod', Tn_min=1e-1, Tn_max=1e5, y_fac=1. ):
    T_prod = Tn
    data = np.genfromtxt( file )
    if not type=='prod':
        data = np.transpose(data)
    min  = 1e90;
    col  = 1
    for i in range(1,len(data[0,:])):
        diff = (data[0,i]-T_prod)*(data[0,i]-T_prod)
        if diff<min:
            min = diff
            col = i
    x_i = 1
    for i in range(1,len(data[:,0])):
        if data[i,0]<Tn_max:
            x_i = i
    x = data[1:x_i,0 ]
    y = data[1:x_i,col]*y_fac
    return x,y




data_Dbar    = np.genfromtxt( 'Dbar_SourceTerm.txt',           skip_header=1 )
data_Dbar_DM = np.genfromtxt( 'Dbar_SourceTerm_DM.txt',        skip_header=1 )

data_pbar    = np.genfromtxt( 'pbar_SourceTerm_diMauro.txt',   skip_header=1 )
data_pbar_DM = np.genfromtxt( 'pbar_SourceTerm_DM.txt',        skip_header=1 )

data_Hebar    = np.genfromtxt( 'Hebar_SourceTerm.txt',           skip_header=1 )
data_Hebar_DM = np.genfromtxt( 'Hebar_SourceTerm_DM.txt',        skip_header=1 )


T_power = 0

##############################
##############################
#   Anti Deuteron
##############################
##############################


plot, fig = plot_1D( r'$\mathrm{T_{\bar{D}}/n [GeV/n]}$', r'$\mathrm{q^{(\bar{D})}  \, [(GeV/n)^{-1}m^{-3}s^{-1}]}}$' )

plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.15)


plot.errorbar( data_Dbar[:,0], data_Dbar[:, 1]*np.power(data_Dbar[:,0],T_power),   lw=4,  color=c0   , dashes=(        ), label=r'Secondary'   )
plot.errorbar( data_Dbar[:,0], data_Dbar[:, 2]*np.power(data_Dbar[:,0],T_power),   lw=2,  color=c1   , dashes=(        ), label=r'p H'   )
plot.errorbar( data_Dbar[:,0], data_Dbar[:, 3]*np.power(data_Dbar[:,0],T_power),   lw=2,  color=c1   , dashes=(3, 3    ), label=r'p He'  )
#plot.errorbar( data_Dbar[:,0], data_Dbar[:, 6]*np.power(data_Dbar[:,0],T_power),   lw=2,  color=c1   , dashes=(14,3    ), label=r'p C'   )
#plot.errorbar( data_Dbar[:,0], data_Dbar[:, 7]*np.power(data_Dbar[:,0],T_power),   lw=2,  color=c1   , dashes=(20,5,5,5), label=r'p O'   )
plot.errorbar( data_Dbar[:,0], data_Dbar[:, 4]*np.power(data_Dbar[:,0],T_power),   lw=2,  color=c2   , dashes=(10, 3   ), label=r'He H'  )
plot.errorbar( data_Dbar[:,0], data_Dbar[:, 5]*np.power(data_Dbar[:,0],T_power),   lw=2,  color=c2   , dashes=(10,3,3,3), label=r'He He' )

plot.errorbar( data_Dbar[:,0], data_Dbar[:, 8]*np.power(data_Dbar[:,0],T_power),   lw=2,  color=c3   , dashes=(6, 6     ), label=r'$\mathrm{\bar{p}}$ H'   )
plot.errorbar( data_Dbar[:,0], data_Dbar[:, 9]*np.power(data_Dbar[:,0],T_power),   lw=2,  color=c3   , dashes=(14,3,10,5), label=r'$\mathrm{\bar{p}}$ He'  )

plot.set_xlim(1e-1, 2e3)
plot.set_ylim(1e-33, 2e-26)


plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), ncol=1, frameon=False, fontsize=0.7*label_size)

plt.savefig('Dbar_source.pdf')

plot, fig = plot_1D( r'$\mathrm{T_{\bar{D}}/n [GeV/n]}$', r'$\mathrm{q^{(\bar{D})}  \, [(GeV/n)^{-1}m^{-3}s^{-1}]}}$' )

plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.15)


plot.errorbar( data_Dbar[:,0], data_Dbar[:, 1]*np.power(data_Dbar[:,0],T_power),       lw=4,  color=c0   , dashes=(        ), label=r'Secondary'   )
plot.errorbar( data_Dbar_DM[:,0], data_Dbar_DM[:, 1]*np.power(data_Dbar[:,0],T_power), lw=4,  color=c5   , dashes=(4,4     ), label=r'DM (CuKrKo)'  )

plot.errorbar( data_Dbar[:,0], data_Dbar[:, 2]*np.power(data_Dbar[:,0],T_power),   lw=2,  color=c1   , dashes=(        ), label=r'p H'   )
plot.errorbar( data_Dbar[:,0], data_Dbar[:, 3]*np.power(data_Dbar[:,0],T_power),   lw=2,  color=c1   , dashes=(3, 3    ), label=r'p He'  )
#plot.errorbar( data_Dbar[:,0], data_Dbar[:, 6]*np.power(data_Dbar[:,0],T_power),   lw=2,  color=c1   , dashes=(14,3    ), label=r'p C'   )
#plot.errorbar( data_Dbar[:,0], data_Dbar[:, 7]*np.power(data_Dbar[:,0],T_power),   lw=2,  color=c1   , dashes=(20,5,5,5), label=r'p O'   )
plot.errorbar( data_Dbar[:,0], data_Dbar[:, 4]*np.power(data_Dbar[:,0],T_power),   lw=2,  color=c2   , dashes=(10, 3   ), label=r'He H'  )
plot.errorbar( data_Dbar[:,0], data_Dbar[:, 5]*np.power(data_Dbar[:,0],T_power),   lw=2,  color=c2   , dashes=(10,3,3,3), label=r'He He' )

plot.errorbar( data_Dbar[:,0], data_Dbar[:, 8]*np.power(data_Dbar[:,0],T_power),   lw=2,  color=c3   , dashes=(6, 6     ), label=r'$\mathrm{\bar{p}}$ H'   )
plot.errorbar( data_Dbar[:,0], data_Dbar[:, 9]*np.power(data_Dbar[:,0],T_power),   lw=2,  color=c3   , dashes=(14,3,10,5), label=r'$\mathrm{\bar{p}}$ He'  )

plot.set_xlim(1e-1, 2e3)
plot.set_ylim(1e-33, 2e-26)


plot.legend(loc='upper center', bbox_to_anchor=(0.5, 0.95), ncol=4, frameon=False, fontsize=0.6*label_size)

plt.savefig('Dbar_source_DM.pdf')


#plot2D('dT_pp_Dbar_LAB.txt',          resfile='XS_pp.pdf'       , x_label=r'$\mathrm{T_{p      } \quad [GeV]}$', y_label=r'$\mathrm{T_{\bar{D}}/n \quad [GeV/n]}$', Tmin_proj=1e1, Tmax_proj=1e5, Tmin_prod=1e-1, Tmax_prod=1e3, Zlabel=r'$\mathrm{d\sigma/(dT_{\bar{D}}/n)\quad[m^2/(GeV/n)]}$'  )
#plot2D('dT_ppbar_Dbar_LAB.txt',       resfile='XS_pbarp.pdf'    , x_label=r'$\mathrm{T_{\bar{p}} \quad [GeV]}$', y_label=r'$\mathrm{T_{\bar{D}}/n \quad [GeV/n]}$', Tmin_proj=1e1, Tmax_proj=1e5, Tmin_prod=1e-1, Tmax_prod=1e3, Zlabel=r'$\mathrm{d\sigma/(dT_{\bar{D}}/n)\quad[m^2/(GeV/n)]}$'  )


# Profiles

plot, fig = plot_1D (r'$\mathrm{T_{\bar{D}}/n \quad [GeV/n]}$',      r'$\mathrm{d\sigma/(dT_{\bar{D}}/n)\quad[m^2/(GeV/n)]}$', 'log', 'log' )

plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.15)

x1,   y1   = profile('dT_pp_Dbar_LAB.txt',                             30, type='proj', y_fac=2.3     )
x2,   y2   = profile('dT_pp_Dbar_LAB.txt',                             50, type='proj', y_fac=2.3     )
x3,   y3   = profile('dT_pp_Dbar_LAB.txt',                            100, type='proj', y_fac=2.3     )
x4,   y4   = profile('dT_pp_Dbar_LAB.txt',                            300, type='proj', y_fac=2.3     )
x5,   y5   = profile('dT_pp_Dbar_LAB.txt',                           1000, type='proj', y_fac=2.3     )

plot.plot(x1,y1, color=c1, lw=3)
plot.plot(x2,y2, color=c2, lw=3)
plot.plot(x3,y3, color=c3, lw=3)
plot.plot(x4,y4, color=c4, lw=3)
plot.plot(x5,y5, color=c5, lw=3)

plot.set_xlim( (5e-1, 5e2)          )
plot.set_ylim( (1e-38, 1e-35)       )

plot.text    ( 1.5e0,  5e-38 , r'$\mathrm{T_{p}=                   30\,GeV}$',   color=c1,     rotation= 0  )
plot.text    (   8e0, 15e-38 , r'$\mathrm{T_{p}=                   50\,GeV}$',   color=c2,     rotation=-75 )
plot.text    (  17e0, 2e-37  , r'$\mathrm{T_{p}=                  100\,GeV}$',   color=c3,     rotation=-75 )
plot.text    (   4e1, 2.7e-37, r'$\mathrm{T_{p}=                  300\,GeV}$',   color=c4,     rotation=-75 )
plot.text    (  10e1, 4e-37  , r'$\mathrm{T_{p}=                 1000\,GeV}$',   color=c5,     rotation=-75 )


plt.savefig('XS_comparison_pp_proj.pdf')

plot, fig = plot_1D (r'$\mathrm{T_{\bar{D}}/n \quad [GeV/n]}$',      r'$\mathrm{d\sigma/(dT_{\bar{D}}/n)\quad[m^2/(GeV/n)]}$', 'log', 'log' )

plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.15)

x1,   y1   = profile('dT_ppbar_Dbar_LAB.txt',                             10, type='proj', y_fac=2.3     )
x2,   y2   = profile('dT_ppbar_Dbar_LAB.txt',                             20, type='proj', y_fac=2.3     )
x3,   y3   = profile('dT_ppbar_Dbar_LAB.txt',                             30, type='proj', y_fac=2.3     )
x4,   y4   = profile('dT_ppbar_Dbar_LAB.txt',                             50, type='proj', y_fac=2.3     )
x5,   y5   = profile('dT_ppbar_Dbar_LAB.txt',                            100, type='proj', y_fac=2.3     )

plot.plot(x1,y1, color=c1, lw=3)
plot.plot(x2,y2, color=c2, lw=3)
plot.plot(x3,y3, color=c3, lw=3)
plot.plot(x4,y4, color=c4, lw=3)
plot.plot(x5,y5, color=c5, lw=3)

plot.set_xlim( (1e-1, 1e2)          )
plot.set_ylim( (2e-38, 1e-34)       )

plot.text    (   1e0,  3e-37 , r'$\mathrm{T_{\bar{p}}=                   10\,GeV}$',   color=c1,     rotation=0 )
plot.text    ( 3.4e0,2.3e-36 , r'$\mathrm{T_{\bar{p}}=                   20\,GeV}$',   color=c2,     rotation=-60 )
plot.text    ( 5.5e0,  4e-36 , r'$\mathrm{T_{\bar{p}}=                   30\,GeV}$',   color=c3,     rotation=-60 )
plot.text    (   9e0,  6e-36 , r'$\mathrm{T_{\bar{p}}=                   50\,GeV}$',   color=c4,     rotation=-60 )
plot.text    (  15e0,  1e-35 , r'$\mathrm{T_{\bar{p}}=                  100\,GeV}$',   color=c5,     rotation=-60 )



plt.savefig('XS_comparison_pbarp_proj.pdf')



##############################
##############################
#   Anti Helium
##############################
##############################

T_power = 0.


plot, fig = plot_1D( r'$\mathrm{T_{\bar{He}}/n [GeV/n]}$', r'$\mathrm{q^{(\bar{He})}  \, [(GeV/n)^{-1}m^{-3}s^{-1}]}}$' )

plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.15)


plot.errorbar( data_Hebar[:,0], data_Hebar[:, 1]*np.power(data_Hebar[:,0],T_power),   lw=4,  color=c0   , dashes=(        ), label=r'Secondary'   )
plot.errorbar( data_Hebar[:,0], data_Hebar[:, 2]*np.power(data_Hebar[:,0],T_power),   lw=2,  color=c1   , dashes=(        ), label=r'p H'   )
plot.errorbar( data_Hebar[:,0], data_Hebar[:, 3]*np.power(data_Hebar[:,0],T_power),   lw=2,  color=c1   , dashes=(3, 3    ), label=r'p He'  )
#plot.errorbar( data_Hebar[:,0], data_Hebar[:, 6]*np.power(data_Hebar[:,0],T_power),   lw=2,  color=c1   , dashes=(14,3    ), label=r'p C'   )
#plot.errorbar( data_Hebar[:,0], data_Hebar[:, 7]*np.power(data_Hebar[:,0],T_power),   lw=2,  color=c1   , dashes=(20,5,5,5), label=r'p O'   )
plot.errorbar( data_Hebar[:,0], data_Hebar[:, 4]*np.power(data_Hebar[:,0],T_power),   lw=2,  color=c2   , dashes=(10, 3   ), label=r'He H'  )
plot.errorbar( data_Hebar[:,0], data_Hebar[:, 5]*np.power(data_Hebar[:,0],T_power),   lw=2,  color=c2   , dashes=(10,3,3,3), label=r'He He' )

plot.errorbar( data_Hebar[:,0], data_Hebar[:, 8]*np.power(data_Hebar[:,0],T_power),   lw=2,  color=c3   , dashes=(6, 6     ), label=r'$\mathrm{\bar{p}}$ H'   )
plot.errorbar( data_Hebar[:,0], data_Hebar[:, 9]*np.power(data_Hebar[:,0],T_power),   lw=2,  color=c3   , dashes=(14,3,10,5), label=r'$\mathrm{\bar{p}}$ He'  )

plot.set_xlim(1e-1, 2e3)
plot.set_ylim(1e-38, 2e-30)


plot.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95), ncol=1, frameon=False, fontsize=0.7*label_size)

plt.savefig('Hebar_source.pdf')

plot, fig = plot_1D( r'$\mathrm{T_{\bar{He}}/n [GeV/n]}$', r'$\mathrm{q^{(\bar{He})}  \, [(GeV/n)^{-1}m^{-3}s^{-1}]}}$' )

plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.15)


plot.errorbar( data_Hebar[:,0], data_Hebar[:, 1]*np.power(data_Hebar[:,0],T_power),       lw=4,  color=c0   , dashes=(        ), label=r'Secondary'   )
plot.errorbar( data_Hebar_DM[:,0], data_Hebar_DM[:, 1]*np.power(data_Hebar[:,0],T_power), lw=4,  color=c5   , dashes=(4,4     ), label=r'DM (CuKrKo)'  )

plot.errorbar( data_Hebar[:,0], data_Hebar[:, 2]*np.power(data_Hebar[:,0],T_power),   lw=2,  color=c1   , dashes=(        ), label=r'p H'   )
plot.errorbar( data_Hebar[:,0], data_Hebar[:, 3]*np.power(data_Hebar[:,0],T_power),   lw=2,  color=c1   , dashes=(3, 3    ), label=r'p He'  )
#plot.errorbar( data_Hebar[:,0], data_Hebar[:, 6]*np.power(data_Hebar[:,0],T_power),   lw=2,  color=c1   , dashes=(14,3    ), label=r'p C'   )
#plot.errorbar( data_Hebar[:,0], data_Hebar[:, 7]*np.power(data_Hebar[:,0],T_power),   lw=2,  color=c1   , dashes=(20,5,5,5), label=r'p O'   )
plot.errorbar( data_Hebar[:,0], data_Hebar[:, 4]*np.power(data_Hebar[:,0],T_power),   lw=2,  color=c2   , dashes=(10, 3   ), label=r'He H'  )
plot.errorbar( data_Hebar[:,0], data_Hebar[:, 5]*np.power(data_Hebar[:,0],T_power),   lw=2,  color=c2   , dashes=(10,3,3,3), label=r'He He' )

plot.errorbar( data_Hebar[:,0], data_Hebar[:, 8]*np.power(data_Hebar[:,0],T_power),   lw=2,  color=c3   , dashes=(6, 6     ), label=r'$\mathrm{\bar{p}}$ H'   )
plot.errorbar( data_Hebar[:,0], data_Hebar[:, 9]*np.power(data_Hebar[:,0],T_power),   lw=2,  color=c3   , dashes=(14,3,10,5), label=r'$\mathrm{\bar{p}}$ He'  )

plot.set_xlim(1e-1, 2e3)
plot.set_ylim(1e-38, 2e-30)


plot.legend(loc='upper center', bbox_to_anchor=(0.5, 0.95), ncol=4, frameon=False, fontsize=0.6*label_size)

plt.savefig('Hebar_source_DM.pdf')


plot2D('dT_pp_Hebar_LAB.txt',          resfile='XS_pp_Hebar.pdf'       , x_label=r'$\mathrm{T_{p      } \quad [GeV]}$', y_label=r'$\mathrm{T_{\bar{He}}/n \quad [GeV/n]}$', Tmin_proj=1e1, Tmax_proj=1e5, Tmin_prod=1e-1, Tmax_prod=1e3, Zlabel=r'$\mathrm{d\sigma/(dT_{\bar{He}}/n)\quad[m^2/(GeV/n)]}$'  )
#plot2D('dT_ppbar_Hebar_LAB.txt',       resfile='XS_pbarp_Hebar.pdf'    , x_label=r'$\mathrm{T_{\bar{p}} \quad [GeV]}$', y_label=r'$\mathrm{T_{\bar{He}}/n \quad [GeV/n]}$', Tmin_proj=1e1, Tmax_proj=1e5, Tmin_prod=1e-1, Tmax_prod=1e3, Zlabel=r'$\mathrm{d\sigma/(dT_{\bar{He}}/n)\quad[m^2/(GeV/n)]}$'  )



##############################
##############################
#   Anti Proton
##############################
##############################

T_power = 2.7

plot, fig = plot_1D( r'$\mathrm{T_{\bar{p}} [GeV]}$', r'$\mathrm{q^{(\bar{p})} \cdot (T_{\bar{p}})^{2.7}  \, [(GeV)^{1.7}m^{-3}s^{-1}]}}$' )

plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.15)


plot.errorbar( data_pbar[:,0], data_pbar[:, 1]*np.power(data_pbar[:,0],T_power),   lw=4,  color=c0   , dashes=(        ), label=r'Secondary' )
plot.errorbar( data_pbar[:,0], data_pbar[:, 2]*np.power(data_pbar[:,0],T_power),   lw=2,  color=c1   , dashes=(        ), label=r'p H'   )
plot.errorbar( data_pbar[:,0], data_pbar[:, 3]*np.power(data_pbar[:,0],T_power),   lw=2,  color=c1   , dashes=(3, 3    ), label=r'p He'  )
#plot.errorbar( data_pbar[:,0], data_pbar[:, 6]*np.power(data_pbar[:,0],T_power),   lw=2,  color=c1   , dashes=(14,3    ), label=r'p C'   )
#plot.errorbar( data_pbar[:,0], data_pbar[:, 7]*np.power(data_pbar[:,0],T_power),   lw=2,  color=c1   , dashes=(20,5,5,5), label=r'p O'   )
plot.errorbar( data_pbar[:,0], data_pbar[:, 4]*np.power(data_pbar[:,0],T_power),   lw=2,  color=c2   , dashes=(10, 3   ), label=r'He H'  )
plot.errorbar( data_pbar[:,0], data_pbar[:, 5]*np.power(data_pbar[:,0],T_power),   lw=2,  color=c2   , dashes=(10,3,3,3), label=r'He He' )


plot.set_xlim(1e-1, 2e3)
plot.set_ylim(5e-26, 5e-20)


plot.legend(loc='upper center', bbox_to_anchor=(0.5, 0.95), ncol=3, frameon=False, fontsize=0.7*label_size)

plt.savefig('pbar_source.pdf')

plot, fig = plot_1D( r'$\mathrm{T_{\bar{p}} [GeV]}$', r'$\mathrm{q^{(\bar{p})} \cdot (T_{\bar{p}})^{2.7}  \, [(GeV)^{1.7}m^{-3}s^{-1}]}}$' )

plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.15)


plot.errorbar( data_pbar[:,0], data_pbar[:, 1]*np.power(data_pbar[:,0],T_power),         lw=4,  color=c0   , dashes=(        ), label=r'Secondary'    )
plot.errorbar( data_pbar_DM[:,0], data_pbar_DM[:, 1]*np.power(data_pbar[:,0],T_power),   lw=4,  color=c5   , dashes=( 4,4    ), label=r'DM (CuKrKo)'  )

plot.errorbar( data_pbar[:,0], data_pbar[:, 2]*np.power(data_pbar[:,0],T_power),   lw=2,  color=c1   , dashes=(        ), label=r'p H'   )
plot.errorbar( data_pbar[:,0], data_pbar[:, 3]*np.power(data_pbar[:,0],T_power),   lw=2,  color=c1   , dashes=(3, 3    ), label=r'p He'  )
#plot.errorbar( data_pbar[:,0], data_pbar[:, 6]*np.power(data_pbar[:,0],T_power),   lw=2,  color=c1   , dashes=(14,3    ), label=r'p C'   )
#plot.errorbar( data_pbar[:,0], data_pbar[:, 7]*np.power(data_pbar[:,0],T_power),   lw=2,  color=c1   , dashes=(20,5,5,5), label=r'p O'   )
plot.errorbar( data_pbar[:,0], data_pbar[:, 4]*np.power(data_pbar[:,0],T_power),   lw=2,  color=c2   , dashes=(10, 3   ), label=r'He H'  )
plot.errorbar( data_pbar[:,0], data_pbar[:, 5]*np.power(data_pbar[:,0],T_power),   lw=2,  color=c2   , dashes=(10,3,3,3), label=r'He He' )


plot.set_xlim(1e-1, 2e3)
plot.set_ylim(5e-26, 5e-20)



plot.legend(loc='upper center', bbox_to_anchor=(0.5, 0.95), ncol=4, frameon=False, fontsize=0.6*label_size)

plt.savefig('pbar_source_DM.pdf')


