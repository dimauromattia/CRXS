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
parser.add_argument('--species',     help='Dbar',                           action='store', dest='species',     type=str, default='Dbar'  )
parser.add_argument('--final_state', help='bbbar',                          action='store', dest='final_state', type=str, default='bbbar' )
parser.add_argument('--step',        help='step',                           action='store', dest='step'   ,     type=int, default=1       )
parser.add_argument('--e',           help='galprop_src',                    action='store', dest='galprop_src', type=str, default='galprop'     )

args        = parser.parse_args()
species     = args.species
final_state = args.final_state
step        = args.step
galprop_src = args.galprop_src
suffix      = ''

################################################################################
##                                                                            ##
##                                                                            ##
##          Careful: Use this only for antideuteron, not for Hebar            ##
##                                                                            ##
##                                                                            ##
################################################################################


c0='black'
c1='#0404B4'
c2='#04B404'
c3='#0484B4'
c4='#B40484'
c5='#B40404'
c6='#B48404'

c1_a = (0.015,0.015,0.703, 0.2)
c2_a = (0.015,0.703,0.015, 0.2)
c3_a = (0.015,0.535,0.703, 0.2)
c4_a = (0.703,0.015,0.535, 0.2)
c5_a = (0.703,0.015,0.015, 0.2)
c6_a = (0.703,0.703,0.015, 0.2)

c2_aa = (0.015,0.703,0.015, 0.1)

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
    
    return plot, fig

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
    
    no_value = []
    for i in range(len(lower)-2):
        for j in range(len(x)):
            if x[j] == min:
                y_min = y[j]
            if x[j] == max:
                y_max = y[j]
            if x[j]>bins[i] and x[j]<bins[i+1]:
                if y[j] < lower[i+1]:
                    lower[i+1] = y[j]
                if y[j] > upper[i+1]:
                    upper[i+1] = y[j]
        if lower[i+1]>1e80:
            no_value.append( i+1 )
    
    x_r = bins-s/2
    x_r[ 0] = min
    x_r[-1] = max

    for i in no_value[::-1]:
        x_r   = np.delete(x_r,  i)
        lower = np.delete(lower,i)
        upper = np.delete(upper,i)

    lower[ 0] = 0.333*(lower[ 1]+upper[ 1] + y_min)
    upper[ 0] = 0.333*(lower[ 1]+upper[ 1] + y_min)
    lower[-1] = 0.333*(lower[-2]+upper[-2] + y_max)
    upper[-1] = 0.333*(lower[-2]+upper[-2] + y_max)

    if t=='cubic':
        for i in range(5):
            for j in range(1, len(lower)-1 ):
                if lower[j]>(lower[j-1] + lower[j+1] )*0.5:
                    lower[j] = (lower[j-1] + lower[j+1] ) * 0.5
                if upper[j] < (upper[j-1] + upper[j+1] ) * 0.5 :
                    upper[j] = (upper[j-1] + upper[j+1] ) * 0.5
        x_r[ 0] -= np.fabs(0.001*x_r[ 0])
        x_r[-1] += np.fabs(0.001*x_r[-1])
        f_l = interp1d( x_r, lower, 'cubic' )
        f_u = interp1d( x_r, upper, 'cubic' )
        s = s/10
        x_r = np.arange(min, max+s/2, s)
        return x_r, f_l(x_r), f_u(x_r)
    
    return x_r, lower, upper

def get_contour_logy(x, y, n=10, t='cubic'):
    x_r, lower, upper = get_contour(x, np.log(y), n, t)
    return x_r, np.exp(lower), np.exp(upper)


if step==5:
    
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

    plot, fig = plot_1D(r'$m_\mathrm{DM}$ [GeV]', r'Average flux over GAPS sensitivity', 'linear', 'log')

    m, l_n, u_n = get_contour_logy(  MED_2[:,0],  MED_2[:,5]  )
    plot.fill_between(m, l_n, u_n, color=c5_a, edgecolor=c5, lw=2, label='MED')
    m, l_n, u_n = get_contour_logy(  MAX_2[:,0],  MAX_2[:,5]  )
    plot.fill_between(m, l_n, u_n, color=c2_aa, edgecolor=c2, lw=2, linestyle=':', label='MAX')
    m, l_n, u_n = get_contour_logy(  CuKoKr_2[:,0],  CuKoKr_2[:,5]  )
    plot.fill_between(m, l_n, u_n, color=c1_a, edgecolor=c1, lw=2, label='CuKoKr')

    plot.text( 58, 19.0,   r'MAX'   ,    color=c2,  fontsize=label_size, horizontalalignment='right'  )
    plot.text( 59,  1.8,   r'MED'   ,    color=c5,  fontsize=label_size, horizontalalignment='right'  )
    plot.text( 57,  6.0,   r'CuKrKo',    color=c1,  fontsize=label_size, horizontalalignment='right'  )

    plot.axhline( y=1, lw=2, c='black' )
    plot.grid(b=None, axis='y')
    plot.set_xlim( ( 50, 90) )
    plot.set_ylim( (0.5,100) )
    plt.savefig('number_prop.pdf')


    plot, fig = plot_1D(r'$m_\mathrm{DM}$ [GeV]', r'Average flux over GAPS sensitivity', 'linear', 'log')

    m, l_n, u_n = get_contour_logy(  CuKoKr_2[:,0],  CuKoKr_2[:,4]  )
    plot.fill_between(m, l_n, u_n, color=c5_a, edgecolor=c5, lw=2, label=r'Solar Modulation, $\phi = 200\quad \mathrm{MV}$')
    m, l_n, u_n = get_contour_logy(  CuKoKr_2[:,0],  CuKoKr_2[:,5]  )
    plot.fill_between(m, l_n, u_n, color=c2_a, edgecolor=c2, lw=2, label=r'Solar Modulation, $\phi = 300\quad \mathrm{MV}$')
    m, l_n, u_n = get_contour_logy(  CuKoKr_2[:,0],  CuKoKr_2[:,6]  )
    plot.fill_between(m, l_n, u_n, color=c1_a, edgecolor=c1, lw=2, label=r'Solar Modulation, $\phi = 400\quad \mathrm{MV}$')

    plot.text( 73,  9.9,  r'$\phi = 300  \,\, \mathrm{MV}$'  , color=c5,  fontsize=label_size  )
    plot.text( 52,  6.0,  r'$\phi = 400  \,\, \mathrm{MV}$'  , color=c2,  fontsize=label_size  )
    plot.text( 67,  2.4,  r'$\phi = 500  \,\, \mathrm{MV}$'  , color=c1,  fontsize=label_size  )

    plot.axhline( y=1, lw=2, c='black' )
    plot.grid(b=None, axis='y')
    plot.set_xlim( ( 50, 90) )
    plot.set_ylim( (0.5,100) )
    plt.savefig('number_solarmod.pdf')

    plot, fig = plot_1D(r'$m_\mathrm{DM}$ [GeV]', r'Average flux over GAPS sensitivity', 'linear', 'log')

    m, l_n, u_n = get_contour_logy(  CuKoKr_2[:,0],  CuKoKr_2[:,5]*np.power( 80./80.,3)  )
    plot.fill_between(m, l_n, u_n, color=c1_a, edgecolor=c1, lw=2 )
    m, l_n, u_n = get_contour_logy(  CuKoKr_2[:,0],  CuKoKr_2[:,5]*np.power(124./80, 3)  )
    plot.fill_between(m, l_n, u_n, color=c5_a, edgecolor=c5, lw=2,)
    
    plot.text( 64,  2.50,  r'$p_{\mathrm{coal}} = 160  \,\, \mathrm{MeV}$'  , color=c1,  fontsize=label_size  )
    plot.text( 70,  35.0,  r'$p_{\mathrm{coal}} = 248  \,\, \mathrm{MeV}$'  , color=c5,  fontsize=label_size  )

    plot.axhline( y=1, lw=2, c='black' )
    plot.grid(b=None, axis='y')
    plot.set_xlim( ( 50, 90) )
    plot.set_ylim( (0.5,100) )
    plt.savefig('number_coalescence_momentum.pdf')

    plot, fig = plot_1D(r'$m_\mathrm{DM}$ [GeV]', r'Average flux over GAPS sensitivity', 'linear', 'log')

    m, l_n, u_n = get_contour_logy(  CuKoKr_2[:,0],  CuKoKr_2[:,5] )
    plot.fill_between(m, l_n, u_n, color=c1_a, edgecolor=c1, lw=2) #, label=r'Analytic coalescence')
    m, l_n, u_n = get_contour_logy(  Vittino_2[:,0],  Vittino_2[:,5]  )
    plot.fill_between(m, l_n, u_n, color=c5_a, edgecolor=c5, lw=2) #, label=r'Monte Carlo coalescence')

    plot.text( 73, 9  ,  'analytic\ncoalescence'                    , color=c1,  fontsize=label_size  )
    plot.text( 52, 1.5,  'Monte Carlo\ncoalescence\n(Pythia)'       , color=c5,  fontsize=label_size  )

    plot.axhline( y=1, lw=2, c='black' )
    plot.grid(b=None, axis='y')
    plot.set_xlim( ( 50, 90) )
    plot.set_ylim( (0.5,100) )
    plt.savefig('number_coalescence_model.pdf')


if step==7:
    print 'plot'

    CuKoKr_2_bb    = np.genfromtxt('CuKoKr_npar_2sigma.txt')
    CuKoKr_2_gg    = np.genfromtxt('CuKoKr_gg_npar_2sigma.txt')
    CuKoKr_2_ww    = np.genfromtxt('CuKoKr_ww_npar_2sigma.txt')
    CuKoKr_2_zz    = np.genfromtxt('CuKoKr_zz_npar_2sigma.txt')
    CuKoKr_2_hh    = np.genfromtxt('CuKoKr_hh_npar_2sigma.txt')
    CuKoKr_2_tt    = np.genfromtxt('CuKoKr_ttbar_npar_2sigma.txt')


    plot, fig = plot_1D(r'$m_\mathrm{DM}$ [GeV]', r'Average flux over GAPS sensitivity', 'log', 'log')

    m, l_n, u_n = get_contour_logy(  CuKoKr_2_gg[:,0],  CuKoKr_2_gg[:,5]  )
    plot.fill_between(m, l_n, u_n, color=c5_a, edgecolor=c5, lw=2, label='gg')

    #m, l_n, u_n = get_contour_logy(  CuKoKr_2_ww[:,0],  CuKoKr_2_ww[:,5]  )
    #plot.fill_between(m, l_n, u_n, color=c5_a, edgecolor=c5, lw=2, label='ww')
    #plot.scatter(CuKoKr_2_ww[:,0],  CuKoKr_2_ww[:,5], color=c5)

    m, l_n, u_n = get_contour_logy(  CuKoKr_2_zz[:,0],  CuKoKr_2_zz[:,5]  )
    plot.fill_between(m, l_n, u_n, color=c2_a, edgecolor=c2, lw=2, label='zz')
    #plot.scatter(CuKoKr_2_zz[:,0],  CuKoKr_2_zz[:,5], color=c2)


    m, l_n, u_n = get_contour_logy(  CuKoKr_2_bb[:,0],  CuKoKr_2_bb[:,5]  )
    plot.fill_between(m, l_n, u_n, color=c1_a, edgecolor=c1, lw=2, label='bb')
    

    y, l_x, u_x = get_contour(  np.log(CuKoKr_2_hh[:,5]),  np.log(CuKoKr_2_hh[:,0])  )

    x_vec = np.ones(2*len(y))
    y_vec = np.ones(2*len(y))
    x_vec[:len(y)] = np.exp(l_x      )
    x_vec[len(y):] = np.exp(u_x[::-1])
    y_vec[:len(y)] = np.exp(y      )
    y_vec[len(y):] = np.exp(y[::-1])
    verts = []
    codes = []
    verts.append( (x_vec[0], y_vec[0]) )
    codes.append( Path.MOVETO )
    for i in range(len(x_vec)):
        verts.append( (x_vec[i], y_vec[i]) )
        codes.append( Path.LINETO )
    verts.append( (x_vec[0], y_vec[0]) )
    codes.append( Path.CLOSEPOLY )
    path_HH = Path(verts, codes)
    HH = mpatches.PathPatch(path_HH, label='hh', linewidth=2,  ec=c3, fc=c3_a )
    plot.add_patch(  HH  )

    y, l_x, u_x = get_contour(  np.log(CuKoKr_2_tt[:,5]),  np.log(CuKoKr_2_tt[:,0])  )
    
    x_vec = np.ones(2*len(y))
    y_vec = np.ones(2*len(y))
    x_vec[:len(y)] = np.exp(l_x      )
    x_vec[len(y):] = np.exp(u_x[::-1])
    y_vec[:len(y)] = np.exp(y      )
    y_vec[len(y):] = np.exp(y[::-1])
    verts = []
    codes = []
    verts.append( (x_vec[0], y_vec[0]) )
    codes.append( Path.MOVETO )
    for i in range(len(x_vec)):
        verts.append( (x_vec[i], y_vec[i]) )
        codes.append( Path.LINETO )
    verts.append( (x_vec[0], y_vec[0]) )
    codes.append( Path.CLOSEPOLY )
    path_TT = Path(verts, codes)
    TT = mpatches.PathPatch(path_TT, label='hh', linewidth=2,  ec=c4, fc=c4_a )
    plot.add_patch(  TT  )



    plot.text( 3.0e+1, 9.6e+0,   r'$gg$'   ,         color=c5,  fontsize=1.2*label_size, ha='center', va='center'  )
    plot.text( 1.6e+2, 6.0e+0,   r'$hh$'   ,         color=c3,  fontsize=1.2*label_size, ha='center', va='center'  )
    plot.text( 2.0e+2, 1.2e+0,   r'$t\bar{t}$'  ,    color=c4,  fontsize=1.2*label_size, ha='center', va='center'  )
    plot.text( 8.0e+1, 8.0e+0,   r'$b\bar{b}$'  ,    color=c1,  fontsize=1.2*label_size, ha='center', va='center'  )
    plot.text( 5.4e+1, 2.0e+0,   r'$ZZ^*$'  ,        color=c2,  fontsize=1.2*label_size, ha='center', va='center'  )
    #plot.text( 4.5e+1, 2.1e+0,   r'$WW^*$'  ,        color=c5,  fontsize=1.2*label_size, ha='center', va='center'  )

    plot.axhline( y=1, lw=2, c='black' )
    plot.grid(b=None, axis='y')
    plot.set_xlim( ( 10, 1000) )
    plot.set_ylim( (0.5, 20) )
    plt.savefig('number_finalstate.pdf')

    sys.exit(0)



A_sp            = 2.
if species=='Dbar':
    A_sp        = 2.
if species=='Hebar':
    A_sp        = 3.


CRACS = os.environ['CRACS']

if not step==0:
    data_cirelli_pbar          = np.genfromtxt( CRACS+'/data/project/AntiMatter/AtProduction_antiprotons_with_ZZstar_WWstar.dat',      skip_header=1 )


ff = 13 # bbbar

if final_state=='gg':
    ff=21
if final_state=='bbbar':
    ff=13
if final_state=='ww':
    ff=17
if final_state=='zz':
    ff=20
if final_state=='hh':
    ff=23
if final_state=='ttbar':
    ff=14




    
if ff==14:  ## ttbar
    intM    = 20
    intX    = 179
    i_from   =  intM     * intX
    print 'ttbar'
    print data_cirelli_pbar[i_from,0]
    i_to     = (intM+1)  * intX
    data_cirelli_pbar[i_from:i_to,ff] = data_cirelli_pbar[i_from+intX:i_to+intX,ff]
if ff==23:  ## hh
    intM    = 16
    intX    = 179
    i_from   =  intM     * intX
    i_to     = (intM+1)  * intX
    print 'hh'
    print data_cirelli_pbar[i_from,0]
    data_cirelli_pbar[i_from:i_to,ff] = data_cirelli_pbar[i_from+intX:i_to+intX,ff]

#intM    = 16
#intX    = 179
#i_from   =  intM     * intX
#i_to     = (intM+1)  * intX
#print data_cirelli_pbar[i_from:i_to,0]
#print data_cirelli_pbar[i_from:i_to,1]


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

#ff=17
#
#plot, fig = plot_1D(r'x', 'dN/dx', 'log', 'log')
#x, N = getSpectrum(40)
#plot.plot(x, N, label='40')
#x, N = getSpectrum(50)
#plot.plot(x, N, label='50')
#x, N = getSpectrum(60)
#plot.plot(x, N, label='60')
#x, N = getSpectrum(70)
#plot.plot(x, N, label='70')
#x, N = getSpectrum(80)
#plot.plot(x, N, label='80')
#x, N = getSpectrum(90)
#plot.plot(x, N, label='90')
#
#plot.set_ylim( (1e-3, 1e+1) )
#
#plot.legend()
#plt.savefig('test_ww.png')
#
#sys.exit(0)





#f_data_cirelli_pbar__ff    = interp2d(  np.log(data_cirelli_pbar[:,0]), data_cirelli_pbar[:,1], data_cirelli_pbar[:,ff], kind='linear')


# XS:
data_dpbar_nar             = np.genfromtxt( CRACS+'/data/project/AntiMatter/table_dpbar_nar.txt',               skip_header=1 )
data_dpbar_tot             = np.genfromtxt( CRACS+'/data/project/AntiMatter/table_dpbar_tot.txt',               skip_header=1 )
data_ppbar_el              = np.genfromtxt( CRACS+'/data/project/AntiMatter/table_ppbar_el.txt',                skip_header=1 )
data_ppbar_tot             = np.genfromtxt( CRACS+'/data/project/AntiMatter/table_ppbar_tot.txt',               skip_header=1 )

data_Anderson              = np.genfromtxt( CRACS+'/data/CS_lab/dT_pp_p_LAB_Anderson.txt'                           )

data_vittino_dbar          = np.genfromtxt( CRACS+'/data/project/AntiMatter/Dbar_dN_dx__bb_100_vittino.txt'         )
data_vittino_dbar_multip   = np.genfromtxt( CRACS+'/data/project/AntiMatter/Dbar_multiplicity_bb_vittino.txt'       )

f_data_vittino_dbar_multip = interp1d(  data_vittino_dbar_multip[:,0], data_vittino_dbar_multip[:,1], kind='cubic'  )

if step==100:
    for dm in np.arange(1,3.01, 0.1 ):
        for x in np.arange(-8.9, 0.01, 0.05):
            m = np.power(10,dm)
            dNdx = np.interp( np.power(10,x), data_vittino_dbar[:,0]/2., data_vittino_dbar[:,1]*2., 0, 0 )*f_data_vittino_dbar_multip(m)/f_data_vittino_dbar_multip(100)
            cir  =  np.log(10)*np.power(10,x)*dNdx
            print str(m).ljust(20) + ' ' + str(x).ljust(20) + ' ' + str(cir).ljust(20)

    plot, fig = plot_1D(r'$x', 'dN/dx', 'log', 'log')

    plot.plot(data_vittino_dbar[:,0],data_vittino_dbar[:,1])

    new_bb = np.genfromtxt('bbbar_vittino.txt')
    plot.plot(new_bb[:,0]/100.*2,new_bb[:,1]*100./2., dashes=(5,5))

#    for i in range(len(new_bb[:,0])):
#        print str(new_bb[len(new_bb[:,0])-1-i,0]/100.*2).ljust(20) + ' ' + str(new_bb[len(new_bb[:,0])-1-i,1]*100./2.).ljust(20)


    plt.savefig('dNdx_vittino.png')


    sys.exit(0)




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

def delta_z(z, h):
    global dZ
    if z<1e-10:
        return 1.*h/dZ
    return 0
vect_delta_z = np.vectorize(delta_z)

if step==10:
    
    p_coal = 0.080
    
    Tn_min = 0.5 * ( np.sqrt( 0.62*0.62 + fMass_d*fMass_d ) - fMass_d )
    Tn_max = 0.5 * ( np.sqrt( 1.03*1.03 + fMass_d*fMass_d ) - fMass_d )
    
    dT = 1e-5
    Tn_vec = np.arange(Tn_min, Tn_max, dT )
    p_vec       = np.sqrt(2*2*Tn_vec*Tn_vec  +  2*fMass_d*2*Tn_vec)
    T_pbar_vec  = np.sqrt(p_vec*p_vec/A_sp/A_sp + fMass_p*fMass_p)-fMass_p
    
    pref =  4./3.*     pow(p_coal, 3)/(p_vec)   *   fMass_d/fMass_p/fMass_n

#    ff = 20
#    T_cirelli, dNdT_cirelli = getSpectrum(91)
#    spec =  np.interp( T_pbar_vec, T_cirelli, dNdT_cirelli ) * 0.5  # factor 0.5 for to correct ZZ -> Z
#    dN_dTn_grid = 2 * pref * np.power(spec,2)                       # factor 2 to make dN_dT -> dN_dTn
#
#    dN_dTn_grid /= 0.69    # hadronic only

    ff = 13
    T_cirelli, dNdT_cirelli = getSpectrum(45.5)
    spec_bb        =  np.interp(T_pbar_vec, T_cirelli, dNdT_cirelli)
    dN_dTn_grid_bb = 2 * pref * np.power(spec_bb,2) # factor 2 to make dN_dT -> dN_dTn
    ff = 12
    T_cirelli, dNdT_cirelli = getSpectrum(45.5)
    spec_cc        =  np.interp(T_pbar_vec, T_cirelli, dNdT_cirelli)
    dN_dTn_grid_cc = 2 * pref * np.power(spec_cc,2) # factor 2 to make dN_dT -> dN_dTn
    ff = 11
    T_cirelli, dNdT_cirelli = getSpectrum(45.5)
    spec_qq        =  np.interp(T_pbar_vec, T_cirelli, dNdT_cirelli)
    dN_dTn_grid_qq = 2 * pref * np.power(spec_qq,2) # factor 2 to make dN_dT -> dN_dTn

    # add up accrding to branching ratios
    dN_dTn_grid    = dN_dTn_grid_qq * (0.6991-0.1512-0.1203)/0.6991 + dN_dTn_grid_cc * 0.1203/0.6991 + dN_dTn_grid_bb * 0.1512/0.6991

    Rd       = dN_dTn_grid.sum() * dT * 0.95
    Rd_aleph = 5.9e-6
    
    p_coal = p_coal * np.power( Rd_aleph/Rd , 1./3. )

    print "Rd:"
    print Rd
    
    print "Rd:"
    print Rd_aleph
    
    print "Pcoal in MeV (new):"
    print p_coal*1e3*2
    
    sys.exit(0)

with open(CRACS+'/data/project/AntiMatter/galdef_56_onlyDM') as f:
    f_galdef = f.readlines()

def proper(ret='spectra', propagation='MED', dark_matter='', p_coal=0.080, injection='analytic', n=0, galprop_param=np.zeros(13), DM_antideuteron=1 ):

    global data_tert, f_data_cirelli_pbar__ff
    
    if propagation=='galprop':
        os.system('mkdir galprop')
        
        s_galdef=''
        for l in f_galdef:
            if 'nuc_g_0_01_001   ' in l:
#                print l
#                print galprop_param[ 0]
                s_galdef += 'nuc_g_0_01_001       = '+str( galprop_param[ 0] ) + '\n'
            elif 'nuc_g_0          ' in l:
#                print l
#                print galprop_param[ 1]
                s_galdef += 'nuc_g_0              = '+str( galprop_param[ 1] ) + '\n'
            elif 'nuc_g_1_01_001   ' in l:
#                print l
#                print galprop_param[ 2]
                s_galdef += 'nuc_g_1_01_001       = '+str( galprop_param[ 2] ) + '\n'
            elif 'nuc_g_1   ' in l:
#                print l
#                print galprop_param[ 3]
                s_galdef += 'nuc_g_1              = '+str( galprop_param[ 3] ) + '\n'
            elif 'nuc_rigid_br0    ' in l:
#                print l
#                print galprop_param[ 4]
                s_galdef += 'nuc_rigid_br0        = '+str( galprop_param[ 4] ) + '\n'
            elif 'nuc_rigid_brS    ' in l:
#                print l
#                print galprop_param[ 5]
                s_galdef += 'nuc_rigid_brS        = '+str( galprop_param[ 5] ) + '\n'
            elif 'D0_xx            ' in l:
#                print l
#                print galprop_param[ 6]
                s_galdef += 'D0_xx                = '+str( galprop_param[ 6] ) + '\n'
            elif 'D_g_1            ' in l:
#                print l
#                print galprop_param[ 7]
                s_galdef += 'D_g_1                = '+str( galprop_param[ 7] ) + '\n'
            elif 'D_g_2            ' in l:
#                print l
#                print galprop_param[ 7]
                s_galdef += 'D_g_2                = '+str( galprop_param[ 7] ) + '\n'
            elif 'v_Alfven         ' in l:
#                print l
#                print galprop_param[ 8]
                s_galdef += 'v_Alfven             = '+str( galprop_param[ 8] ) + '\n'
            elif 'v0_conv          ' in l:
#                print l
#                print galprop_param[ 9]
                s_galdef += 'v0_conv              = '+str( galprop_param[ 9] ) + '\n'
            elif 'z_min            ' in l:
#                print l
#                print galprop_param[10]
                s_galdef += 'z_min                = -'+str(galprop_param[10] ) + '\n'
            elif 'z_max            ' in l:
#                print l
#                print galprop_param[10]
                s_galdef += 'z_max                = +'+str(galprop_param[10] ) + '\n'
            elif 'DM_double0' in l:
#                print l
#                print galprop_param[11]
                s_galdef += 'DM_double0           = '+str( galprop_param[11] ) + '\n'
            elif 'DM_double1' in l:
#                print l
#                print galprop_param[12]
                s_galdef += 'DM_double1           = '+str( galprop_param[12] ) + '\n'
            elif 'DM_antideuteron' in l:
                #                print l
                #                print galprop_param[12]
                s_galdef += 'DM_antideuteron      = '+str( DM_antideuteron ) + '\n'
            else:
                s_galdef += l


        n_galdef = open( 'galprop/galdef_56_'+str(n).rjust(5,'0'), 'w' )
        n_galdef.write(s_galdef)
        n_galdef.close()

#sys.exit(0)
#os.system('cp '+CRACS+'/data/project/AntiMatter/galdef_56_onlyDM galprop/galdef_56_'+str(n).rjust(5,'0') )
        
        cmd = galprop_src + ' -o galprop -g galprop -r '+str(n).rjust(5,'0')+'  -f /usr/local/korsmeier/galprop/data/FITS'
        print cmd
        os.system( cmd )
        os.system('PROPERFIT_plotGalpropSpectrum galprop/nuclei_56_'+str(n).rjust(5,'0')+' --A 2 --Z -1 --p 0 --t dm        --run false')
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
            p_vec       = np.sqrt(A_sp*A_sp*Tn_vec*Tn_vec+2*fMass_d*A_sp*Tn_vec)
            T_pbar_vec  = np.sqrt(p_vec*p_vec/A_sp/A_sp + fMass_p*fMass_p)-fMass_p
    #        spec        =   f_data_cirelli_pbar__ff( np.log(mDM), np.log10(T_pbar_vec/mDM) )[:,0]
    #        spec       *=   1./Tn_vec/np.log(10)
            T_cirelli, dNdT_cirelli = getSpectrum(mDM)
            spec        = np.interp(T_pbar_vec, T_cirelli, dNdT_cirelli)
            if species=='Dbar':
                print p_coal
                pref        =   4./3.*     pow(p_coal, 3)/p_vec        *    fMass_d/fMass_p/fMass_n
            if species=='Hebar':
                pref        =   3.   * pow(pow(p_coal, 3)/p_vec,2)     *    fMass_He3/fMass_p/fMass_n/fMass_n
            dN_dTn_grid = np.interp( Tn_grid, Tn_vec, A_sp * pref * np.power(spec,A_sp)  )
            print 'Source at 600 MeV/n: ' +str(np.interp( 0.6, Tn_vec, A_sp * pref * np.power(spec,A_sp)  ))



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
        
        data_dpbar_in           = (data_ppbar_tot-data_ppbar_el)*data_dpbar_tot/data_ppbar_tot * A_sp/2.
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

    print mDM
    print sv
    print 'Flux at 600 MeV/n: ' + str(np.interp(0.6, Tn_vec, phi_D))

    n_integral = []

    phi=0.0
    Tn_vec_mod = Tn_vec - phi / A_sp
    phi_D_mod =  ( Tn_vec_mod*(Tn_vec_mod+2*fMass_p)  )/(  Tn_vec*(Tn_vec+2*fMass_p)  ) * phi_D
    T_int     = np.arange(0.1,0.25,0.0001)
    phi_int   = np.interp( T_int, Tn_vec_mod, phi_D_mod )
    n_Dbar    = phi_int.sum() * 0.0001/(2e-6*0.15)
    n_integral.append(n_Dbar)

    phi=0.1
    Tn_vec_mod = Tn_vec - phi / A_sp
    phi_D_mod =  ( Tn_vec_mod*(Tn_vec_mod+2*fMass_p)  )/(  Tn_vec*(Tn_vec+2*fMass_p)  ) * phi_D
    T_int     = np.arange(0.1,0.25,0.0001)
    phi_int   = np.interp( T_int, Tn_vec_mod, phi_D_mod )
    n_Dbar    = phi_int.sum() * 0.0001/(2e-6*0.15)
    n_integral.append(n_Dbar)

    phi=0.2
    Tn_vec_mod = Tn_vec - phi / A_sp
    phi_D_mod =  ( Tn_vec_mod*(Tn_vec_mod+2*fMass_p)  )/(  Tn_vec*(Tn_vec+2*fMass_p)  ) * phi_D
    T_int     = np.arange(0.1,0.25,0.0001)
    phi_int   = np.interp( T_int, Tn_vec_mod, phi_D_mod )
    n_Dbar    = phi_int.sum() * 0.0001/(2e-6*0.15)
    n_integral.append(n_Dbar)

    phi=0.3
    Tn_vec_mod = Tn_vec - phi / A_sp
    phi_D_mod =  ( Tn_vec_mod*(Tn_vec_mod+2*fMass_p)  )/(  Tn_vec*(Tn_vec+2*fMass_p)  ) * phi_D
    T_int     = np.arange(0.1,0.25,0.0001)
    phi_int   = np.interp( T_int, Tn_vec_mod, phi_D_mod )
    n_Dbar    = phi_int.sum() * 0.0001/(2e-6*0.15)
    n_integral.append(n_Dbar)

    phi=0.4
    Tn_vec_mod = Tn_vec - phi / A_sp
    phi_D_mod =  ( Tn_vec_mod*(Tn_vec_mod+2*fMass_p)  )/(  Tn_vec*(Tn_vec+2*fMass_p)  ) * phi_D
    T_int     = np.arange(0.1,0.25,0.0001)
    phi_int   = np.interp( T_int, Tn_vec_mod, phi_D_mod )
    n_Dbar    = phi_int.sum() * 0.0001/(2e-6*0.15)
    n_integral.append(n_Dbar)

    phi=0.5
    Tn_vec_mod = Tn_vec - phi / A_sp
    phi_D_mod =  ( Tn_vec_mod*(Tn_vec_mod+2*fMass_p)  )/(  Tn_vec*(Tn_vec+2*fMass_p)  ) * phi_D
    T_int     = np.arange(0.1,0.25,0.0001)
    phi_int   = np.interp( T_int, Tn_vec_mod, phi_D_mod )
    n_Dbar    = phi_int.sum() * 0.0001/(2e-6*0.15)
    n_integral.append(n_Dbar)


    if ret=='spectra':
        return (Tn_vec, phi_D), (Tn_vec_mod, phi_D_mod)
    if ret=='integral':
        return n_integral




data_gg_dm    = np.genfromtxt( CRACS+'/data/project/AntiMatter/gg/gg.txt'       )
data_bbbar_dm = np.genfromtxt( CRACS+'/data/project/AntiMatter/bbbar/bbbar.txt' )
data_ww_dm    = np.genfromtxt( CRACS+'/data/project/AntiMatter/ww/ww.txt'       )
data_zz_dm    = np.genfromtxt( CRACS+'/data/project/AntiMatter/zz/zz.txt'       )
data_hh_dm    = np.genfromtxt( CRACS+'/data/project/AntiMatter/hh/hh.txt'       )
data_ttbar_dm = np.genfromtxt( CRACS+'/data/project/AntiMatter/ttbar/ttbar.txt' )

x = np.argsort( data_gg_dm[:, 1]    )
data_gg_dm    = data_gg_dm[x,::]
x = np.argsort( data_bbbar_dm[:, 1] )
data_bbbar_dm = data_bbbar_dm[x,::]
x = np.argsort( data_ww_dm[:, 1]    )
data_ww_dm    = data_ww_dm[x,::]
x = np.argsort( data_zz_dm[:, 1]    )
data_zz_dm    = data_zz_dm[x,::]
x = np.argsort( data_hh_dm[:, 1]    )
data_hh_dm    = data_hh_dm[x,::]
x = np.argsort( data_ttbar_dm[:, 1] )
data_ttbar_dm = data_ttbar_dm[x,::]


chiSq_bb    =   data_bbbar_dm[:, 1]
K_0_bb      =   data_bbbar_dm[:, 8]
V_c_bb      =   data_bbbar_dm[:,11]
delta_bb    =   data_bbbar_dm[:, 9]
L_bb        =   data_bbbar_dm[:,12]

m_DM_bb     =   np.power( 10, data_bbbar_dm[:,13]-3 )
sv_DM_bb    =   np.power( 10, data_bbbar_dm[:,14]   )



m_DM        = []
n_particles = []

injection = 'analytic'
if step==4:
    injection = 'vittino'
if step==1 or step==4 or step==0:
    if not step==0 and not step==4:
        K_0_bb     *=   np.power(0.25,delta_bb)* np.power(3.24e-22,2)/(1e-6*3.171e-8)      # Transform from (galprop) D0_xx to K_0

    not1sigma = True
    not2sigma = True
    not3sigma = True
    for i in range( 0, len(chiSq_bb) ):
        m_DM.append(m_DM_bb[i])
        propagation = ( 0.1,  K_0_bb[i],  V_c_bb[i],  delta_bb[i],  L_bb[i] )
        if step==0 or step==4:
            propagation='galprop'
        if step==0:
            DM_antideuteron=1
        if step==4:
            DM_antideuteron=2
        dm          = ( m_DM_bb[i], sv_DM_bb[i], 0.43 )
        n =  proper(ret='integral', propagation=propagation, dark_matter=dm, p_coal=0.080, injection=injection, galprop_param=data_bbbar_dm[i,2:], n=i, DM_antideuteron=DM_antideuteron )
        print n
        n_particles.append(n)
        print str(i) + '   ' + str(chiSq_bb[i]) + '/' + str(chiSq_bb[0])
        case='CuKoKr'
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
        n =  proper(ret='integral', propagation='MED', dark_matter=dm, p_coal=0.080)
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
        n =  proper(ret='integral', propagation='MAX', dark_matter=dm, p_coal=0.080)
        n_particles.append(n)
        print n
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




if step==6:
    
    
    chiSq_ff    =   data_bbbar_dm[:, 1]
    K_0_ff      =   data_bbbar_dm[:, 8]
    V_c_ff      =   data_bbbar_dm[:,11]
    delta_ff    =   data_bbbar_dm[:, 9]
    L_ff        =   data_bbbar_dm[:,12]
    m_DM_ff     =   np.power( 10, data_bbbar_dm[:,13]-3 )
    sv_DM_ff    =   np.power( 10, data_bbbar_dm[:,14]   )

    data_ff_dm  =   data_bbbar_dm
    print final_state
    if final_state=='gg':
        print "Hallo Welt"
        chiSq_ff    =   data_gg_dm[:, 1]
        K_0_ff      =   data_gg_dm[:, 8]
        V_c_ff      =   data_gg_dm[:,11]
        delta_ff    =   data_gg_dm[:, 9]
        L_ff        =   data_gg_dm[:,12]
        m_DM_ff     =   np.power( 10, data_gg_dm[:,13]-3 )
        sv_DM_ff    =   np.power( 10, data_gg_dm[:,14]   )
        data_ff_dm  =   data_gg_dm
    if final_state=='ww':
        chiSq_ff    =   data_ww_dm[:, 1]
        K_0_ff      =   data_ww_dm[:, 8]
        V_c_ff      =   data_ww_dm[:,11]
        delta_ff    =   data_ww_dm[:, 9]
        L_ff        =   data_ww_dm[:,12]
        m_DM_ff     =   np.power( 10, data_ww_dm[:,13]-3 )
        sv_DM_ff    =   np.power( 10, data_ww_dm[:,14]   )
        data_ff_dm  =   data_ww_dm
    if final_state=='zz':
        chiSq_ff    =   data_zz_dm[:, 1]
        K_0_ff      =   data_zz_dm[:, 8]
        V_c_ff      =   data_zz_dm[:,11]
        delta_ff    =   data_zz_dm[:, 9]
        L_ff        =   data_zz_dm[:,12]
        m_DM_ff     =   np.power( 10, data_zz_dm[:,13]-3 )
        sv_DM_ff    =   np.power( 10, data_zz_dm[:,14]   )
        data_ff_dm  =   data_zz_dm
    if final_state=='hh':
        chiSq_ff    =   data_hh_dm[:, 1]
        K_0_ff      =   data_hh_dm[:, 8]
        V_c_ff      =   data_hh_dm[:,11]
        delta_ff    =   data_hh_dm[:, 9]
        L_ff        =   data_hh_dm[:,12]
        m_DM_ff     =   np.power( 10, data_hh_dm[:,13]-3 )
        sv_DM_ff    =   np.power( 10, data_hh_dm[:,14]   )
        data_ff_dm  =   data_hh_dm
    if final_state=='ttbar':
        chiSq_ff    =   data_ttbar_dm[:, 1]
        K_0_ff      =   data_ttbar_dm[:, 8]
        V_c_ff      =   data_ttbar_dm[:,11]
        delta_ff    =   data_ttbar_dm[:, 9]
        L_ff        =   data_ttbar_dm[:,12]
        m_DM_ff     =   np.power( 10, data_ttbar_dm[:,13]-3 )
        sv_DM_ff    =   np.power( 10, data_ttbar_dm[:,14]   )
        data_ff_dm  =   data_ttbar_dm


    #K_0_ff     *=   np.power(0.25,delta_ff)* np.power(3.24e-22,2)/(1e-6*3.171e-8)      # Transform from (galprop) D0_xx to K_0


    not1sigma = True
    not2sigma = True
    not3sigma = True
    for i in range( 0, len(chiSq_ff) ):
        print K_0_ff[i]
        print data_ff_dm[i, 8]
        m_DM.append(m_DM_ff[i])
        propagation = ( 0.1,  K_0_ff[i],  V_c_ff[i],  delta_ff[i],  L_ff[i] )
        dm          = ( m_DM_ff[i], sv_DM_ff[i], 0.43 )
        n =  proper(ret='integral', propagation='galprop', dark_matter=dm, p_coal=0.080, injection=injection, galprop_param=data_ff_dm[i,2:], n=i)
        n_particles.append(n)
        print str(i) + '   ' + str(chiSq_ff[i]) + '/' + str(chiSq_ff[0])
        print n
        if i==0:
            print 'm_DM: ' + str(m_DM_ff[i])
            print 'm_DM: ' + str(sv_DM_ff[i])

        case='CuKoKr'
       
        if chiSq_ff[i]-chiSq_ff[0] > 2.3 and not1sigma:
            f = open(case+'_'+final_state+'_npar_1sigma.txt', 'w')
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
        if chiSq_ff[i]-chiSq_ff[0] > 6.18 and not2sigma:
            f = open(case+'_'+final_state+'_npar_2sigma.txt', 'w')
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
        if chiSq_ff[i]-chiSq_ff[0] > 11.83 and not3sigma:
            f = open(case+'_'+final_state+'_npar_3sigma.txt', 'w')
            s = '#' + str('m_DM/GeV').ljust(19) + ' ' + str('#').ljust(20)
            for j in range(len(m_DM)):
                s += '\n' + str(m_DM[j]).ljust(19)
                for k in range(len(n_particles[j])):
                    s += ' ' + str( (n_particles[j])[k]).ljust(20)
            f.write(s)
            f.close()
            not3sigma = False
            break






