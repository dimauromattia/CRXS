#! /usr/bin/env python

import glob
import math
import sys
import os
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.patches as mpatches
from matplotlib import rc
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import ScalarFormatter

from PlotFunctions import plot2D
from PlotFunctions import plot2D


#   Argument Handling
#######################################################
parser = argparse.ArgumentParser(description='CRACS - Plotting the fit results of the Cross Section Scan by MultiNest.')
parser.add_argument('--step',     help='Step .',                        action='store', dest='step',         type=int, default=1)
parser.add_argument('--dir',      help='Directory ',                    action='store', dest='dir',          type=str, default='.')
parser.add_argument('--sdir',     help='Directory ',                    action='store', dest='sdir',         type=str, default='.')
parser.add_argument('--exp',      help='Experiment',                    action='store', dest='exp',          type=str, default='lhcb')
parser.add_argument('--data',     help='Data',                          action='store', dest='data',         type=str, default='na49,na61,dekkers')
parser.add_argument('--data_pHe', help='Data pHe',                      action='store', dest='data_pHe',     type=str, default='')
parser.add_argument('--data_pC',  help='Data pC',                       action='store', dest='data_pC',      type=str, default='na49')
parser.add_argument('--cs ',      help='Cross section parametrization', action='store', dest='cs',           type=str, default='0')
parser.add_argument('--A1',       help='A_1 .',                         action='store', dest='A1',          type=int, default=1)
parser.add_argument('--A2',       help='A_2 .',                         action='store', dest='A2',          type=int, default=1)

parser.add_argument('--cluster',  help='write only script .',           action='store_true',    dest='cluster' )


args        = parser.parse_args()
step   	    = args.step
exp         = args.exp
dir         = args.dir
sdir        = args.sdir
data        = args.data
data_pHe    = args.data_pHe
data_pC     = args.data_pC
cs          = args.cs
A1          = args.A1
A2          = args.A2

cluster     = args.cluster

#######################################################

if sdir=='':
    sdir = os.environ('CRACS')


c0='black'
c1='#0404B4'
c2='#B40404'
c3='#04B404'
c4='#0484B4'
c5='#B40484'
c6='#B48404'
c7='#444444'
c8='#F48404'
c9='#8404B4'

c1_a = (0.015,0.015,0.703, 0.2)
c2_a = (0.703,0.015,0.015, 0.2)
c3_a = (0.015,0.703,0.015, 0.2)
c4_a = (0.015,0.535,0.703, 0.2)
c5_a = (0.703,0.015,0.535, 0.2)
c6_a = (0.703,0.703,0.015, 0.2)


c=[c0,c4,c2,c3,c1,c5]

dashes_array= [(),(3,3),(15,5),(15,5,5,5),(10,3,3,3,3,3)]

def var_to_label(var):
    if 'pT'         in var: return r'p^\mathrm{T}'
    if 'xR'         in var: return r'x_R'
    if 'xf'         in var: return r'x_f'
    if 'y'          in var: return r'y'
    if 'theta_LAB'  in var: return r'\theta_\mathrm{LAB}'
    if 'p_pbar_LAB' in var: return r'p_{\bar{p}, \mathrm{LAB}}'
    return var
def var_to_unit (var):
    if 'pT'         in var: return 'GeV'
    if 'xR'         in var: return ''
    if 'xf'         in var: return ''
    if 'y'          in var: return ''
    if 'theta_LAB'  in var: return 'rad'
    if 'p_pbar_LAB' in var: return 'GeV'
    return ''

def var_to_scale (var):
    if 'pT'         in var: return 'log'
    if 'xR'         in var: return 'linear'
    if 'xf'         in var: return 'linear'
    if 'y'          in var: return 'linear'
    if 'theta_LAB'  in var: return 'linear'
    if 'p_pbar_LAB' in var: return 'log'
    return 'log'


cb='#EFEFFF'

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
    
    plot.set_xlabel ( xlabel )
    plot.set_ylabel ( ylabel )
    
    plot.set_xscale ( xscale )
    plot.set_yscale ( yscale )
    
    plot.tick_params('both', length=20, width=2, which='major', top=True, right=True, direction='in', pad=10)
    plot.tick_params('both', length=10, width=1, which='minor', top=True, right=True, direction='in', pad=10)
    
    plot.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
    
    return plot, fig


experiment  = ['NA49', 'NA61', 'Allaby', 'Antreasyan', 'Guettler', 'Capiluppi', 'Johnson', 'Dekkers', 'BRAHMS', 'Phenix', 'LHCb', 'Abramov', 'Amann', 'Barton', 'Sugaya', 'Veronin']
name        = ['na49', 'na61', 'allaby', 'antreasyan', 'guettler', 'capiluppi', 'johnson', 'dekkers', 'brahms', 'phenix', 'lhcb', 'abramov', 'amann', 'barton', 'sugaya', 'veronin']


def draw( file, prefix, fmt='o', new_plot=True, sqrtS='' ):
    global experiment, name, plot, fig

    lines = [line.rstrip('\n') for line in open(file)]

    for l in lines:
        if l.startswith( '*  ' ):
            var_list = l.split(' ')
            while '' in var_list:
                var_list.remove('')
            var_list.remove('*')
            break

    d = np.genfromtxt( file , skip_header=1)
    e = ''
    A = ''
    for i in range(len(experiment)):
        if name[i] in file:
            e = experiment[i]
            break

    if '_4_' in file:
        A = 'He'
    if '_9_' in file:
        A = 'Be'
    if '_12_' in file:
        A = 'C'
    if '_minuit' in file:
        A += '_minuit_'
    if '_scan' in file:
        A += '_scan_'


    ind = np.lexsort((d[:,4],d[:,3],d[:,0]))
    d = d[ind]

    s   = [d[0,0]]
    s_i = [0     ]
    for j in range(len(d[:,0])):
        if d[j,0]!=s[-1]:
            s.  append(d[j,0])
            s_i.append(j)
    if s_i[-1]!=len(d[:,0])-1:
        s_i.append(len(d[:,0]))

    print e

    for k in range(len(s)):

        print s[k]

        if sqrtS!='':
            if s[k]> 1.1*sqrtS or s[k]<0.9*sqrtS:
                continue

        start = s_i[k]
        stop  = s_i[k+1]

        d_var_fix = d[start:stop,:]

        ind = np.lexsort((d_var_fix[:,4],d_var_fix[:,3],d_var_fix[:,0]))
        d_var_fix= d_var_fix[ind]


        var_fix   = [d_var_fix[0,3]]
        var_fix_u = []
        var_fix_i = [0     ]

        xR = d_var_fix[0,4]
        for j in range(len(d_var_fix[:,3])):
            if xR>d_var_fix[j,4]:
                var_fix.  append(d_var_fix[j,3])
                var_fix_u.append(d_var_fix[j-1,3])
                var_fix_i.append(j)
            xR = d_var_fix[j,4]
        var_fix_i.append(len(d_var_fix[:,3]))
        var_fix_u.append(d_var_fix[-1,3])

        var_x = '$'+var_to_label(var_list[4])
        if var_to_unit(var_list[4])!='':
            var_x += r'\, \mathrm{['+var_to_unit(var_list[4])+']}'
        var_x += '$'
        if new_plot:
            plot, fig = plot_1D( var_x, r'$\sigma_{\mathrm{inv}}\quad \mathrm{[mb/GeV^2]}$', xscale=var_to_scale(var_list[4]) )

        plt.subplots_adjust(left=0.15, right=0.5, top=0.9, bottom=0.15)

        yscale = 1.
        ys = [1,5,2]
        for j in range(len(var_fix)):

            yscale = ys[j%3]*pow(10,-int((j+2)/3))
            start = var_fix_i[j]
            stop  = var_fix_i[j+1]
            col = c[j%6]
            dashes = dashes_array[j/6]

            l = ''
            var_fix_min = var_fix  [j]
            var_fix_max = var_fix_u[j]

            if np.fabs(1-var_fix_min/var_fix_max)<0.01 or var_fix_max<0.001:
                l = '$'+var_to_label(var_list[3])+r'='+format(var_fix_min, '.2f')+'\, \mathrm{'+var_to_unit(var_list[3])+'}$'
            else:
                l = '$'+var_to_label(var_list[3])+r'=['+format(var_fix_min, '.2f')+'-'+format(var_fix_max, '.2f')+']\, \mathrm{'+var_to_unit(var_list[3])+'}$'

            if j>0:
                l += r'$\qquad(\times {'+str(ys[j%3])+'\cdot 10^{-'+str(int((j+2)/3))+'}})$'
            plot.errorbar   ( d_var_fix[start:stop,4], d_var_fix[start:stop,5]*d_var_fix[start:stop,15]*yscale, yerr=d_var_fix[start:stop,6]*d_var_fix[start:stop,15]*yscale, fmt=fmt,  lw=1,  color=col   )

            plot.errorbar       ( d_var_fix[start:stop,4], d_var_fix[start:stop, 10]*yscale,                              lw=2, alpha=1.0,   color=col, dashes=dashes  , label=l   )

            print l + '   ' + var_x + '   ' + str(d_var_fix[start,4]) + ' to ' + str(d_var_fix[stop-1,4])

            plot.fill_between   ( d_var_fix[start:stop,4], d_var_fix[start:stop,13]*yscale, d_var_fix[start:stop,14]*yscale, lw=0, alpha=0.1,  color=col,   )
            plot.fill_between   ( d_var_fix[start:stop,4], d_var_fix[start:stop,11]*yscale, d_var_fix[start:stop,12]*yscale, lw=0, alpha=0.3,  color=col,   )

            if len(d_var_fix[start:stop,4])==1:
                x = [0.9*d_var_fix[start,4],1.1*d_var_fix[start,2]]
                y = [d_var_fix[start,10]*yscale, d_var_fix[start,10]*yscale]
                plot.errorbar   ( x, y,  lw=3, alpha=0.6,  color=col, dashes=dashes  )

        x_min = np.amin(d_var_fix[:,4])
        x_max = np.amax(d_var_fix[:,4])
        fac = 1.1


        if np.log10(x_max/x_min)<1:
            ax = plot.get_xaxis()
            plot.set_xticks([0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500])
            ax.set_major_formatter(ScalarFormatter())

        fac1 = fac
        fac2 = fac
        if x_min<0:
            fac1 = 1./fac
        if x_max<0:
            fac2 = 1./fac
        plot.set_xlim((x_min/fac1, x_max*fac2))

        legend = plot.legend     ( loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False, fontsize=0.6*label_size, title=(e+r'$\qquad \sqrt{s}='+format(s[k], '.1f')+'\,\mathrm{GeV}$') )
        plt.setp        (legend.get_title(),fontsize=0.6*label_size)
        plt.savefig     (prefix + e+A+'_s_'+str(s[k])+'.png')


def draw_lhcb( ):
    global experiment, name, plot, fig

    if not os.path.isfile('bestFit_4_lhcb.txt'):
        return

    d = np.genfromtxt( 'bestFit_4_lhcb.txt' , skip_header=1)
    
    ind = np.lexsort((d[:,4],d[:,3],d[:,0]))
    d = d[ind]

    s = ''
    for i in range(1,len(d[:,0])):
        if d[i,3]!=d[i-1,3]:
#            print str(i) + '  ' + str(d[i-1,3])
            s += str(i)+', '
#    print len(d[:,0])
    s += str(len(d[:,0]))+', '
#    print s

    sep = [0,4, 15, 27, 41, 56, 72, 90, 108, 126, 136]

    plot, fig = plot_1D ( '$x_f$', r'$\sigma_{\mathrm{inv}}\quad \mathrm{[mb/GeV^2]}$', 'linear', 'log' )

    for i in range(1,len(sep)):
        start = sep[i-1]
        stop  = sep[i]
        yscale = np.power(0.6,+i-1)
        plot.errorbar     ( d[start:stop,4], d[start:stop, 5]*d[start:stop,15]*yscale, yerr=d[start:stop,6]*d[start:stop,15]*yscale, fmt='o',  lw=1,  color='black'   )
        plot.errorbar     ( d[start:stop,4], d[start:stop,10]*yscale, lw=2, alpha=1.0,   color='black' )
        #plot.fill_between ( d[start:stop,4], d[start:stop,13]*yscale, d[start:stop,14]*yscale, lw=0, alpha=0.2,  color='black'   )
        plot.fill_between ( d[start:stop,4], d[start:stop,11]*yscale, d[start:stop,12]*yscale, lw=0, alpha=0.3,  color='black'   )

    plot.text( -0.080, 25e-1, r'$p_T=0.47\,\,\mathrm{GeV}$', rotation=+4, fontsize=0.7*label_size)
    plot.text( -0.110, 80e-2, r'$p_T=0.62\,\,\mathrm{GeV}$', rotation=+4, fontsize=0.7*label_size)
    plot.text( -0.118, 30e-2, r'$p_T=0.75\,\,\mathrm{GeV}$', rotation=+4, fontsize=0.7*label_size)
    plot.text( -0.125, 10e-2, r'$p_T=0.85\,\,\mathrm{GeV}$', rotation=+4, fontsize=0.7*label_size)
    plot.text( -0.135, 35e-3, r'$p_T=0.97\,\,\mathrm{GeV}$', rotation=+4, fontsize=0.7*label_size)
    plot.text( -0.147, 10e-3, r'$p_T=1.12\,\,\mathrm{GeV}$', rotation=+4, fontsize=0.7*label_size)
    plot.text( -0.165, 27e-4, r'$p_T=1.32\,\,\mathrm{GeV}$', rotation=+4, fontsize=0.7*label_size)
    plot.text( -0.240, 18e-5, r'$p_T=[1.67-1.69]\,\,\mathrm{GeV}$', rotation=+6, fontsize=0.7*label_size)
    plot.text( -0.240, 43e-6, r'$p_T=[2.21-2.25]\,\,\mathrm{GeV}$', rotation=+6, fontsize=0.7*label_size)
    plot.text( -0.230, 10e-7, r'$p_T=[3.06-3.12]\,\,\mathrm{GeV}$', rotation=+6, fontsize=0.7*label_size)

    plot.text(  0.000, 10e-1, r'$\times\,\, 0.6^1$', rotation=+0, fontsize=0.7*label_size)
    plot.text(  0.002, 40e-2, r'$\times\,\, 0.6^2$', rotation=+0, fontsize=0.7*label_size)
    plot.text(  0.004, 15e-2, r'$\times\,\, 0.6^3$', rotation=+0, fontsize=0.7*label_size)
    plot.text(  0.006, 55e-3, r'$\times\,\, 0.6^4$', rotation=+0, fontsize=0.7*label_size)
    plot.text(  0.008, 17e-3, r'$\times\,\, 0.6^5$', rotation=+0, fontsize=0.7*label_size)
    plot.text(  0.010, 50e-4, r'$\times\,\, 0.6^6$', rotation=+0, fontsize=0.7*label_size)
    plot.text(  0.003, 80e-5, r'$\times\,\, 0.6^7$', rotation=+0, fontsize=0.7*label_size)
    plot.text( -0.005, 70e-6, r'$\times\,\, 0.6^8$', rotation=+0, fontsize=0.7*label_size)
    plot.text( -0.030, 25e-7, r'$\times\,\, 0.6^9$', rotation=+0, fontsize=0.7*label_size)

    plot.set_xlim( (-0.25, 0.05) )
    plot.set_ylim( (1e-7, 1e1) )
              
              

    plt.savefig('lhcb.pdf')

def plot_XS_comparison(prefix, dir_Winkler, dir_diMauro, label_Winkler='Param. II', label_diMauro='Param. I', type='1-1'):

    n_max=1e8

    plot, fig = plot_1D ( r'$\mathrm{T_{\bar{p}} \quad [GeV]}$', r'$\mathrm{d\sigma/dT_{\bar{p}}\quad[m^2/GeV]}$', 'log', 'log' )

    dir=dir_Winkler
    
    x,y = profile (dir+'/sourceTerm_'+type+'/dT_'+type+'_pbar.txt',   20 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
    y_min = np.copy(y)
    y_max = np.copy(y)
    list_1s = glob.glob(dir+'/sourceTerm_'+type+'/dT_'+type+'_pbar_1s*.txt')
    j = 0
    for f in list_1s:
        j+=1
        if j>n_max:
            break
        xx, yy = profile (f,   20 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
        for i in range(len(xx)):
            if yy[i] < y_min[i]:
                y_min[i] = yy[i]
            if yy[i] > y_max[i]:
                y_max[i] = yy[i]
    y_max +=1e-90
    y_min +=1e-90
    plot.plot        ( x, y,              lw=3, zorder=1, dashes=( ), color=c1, label=label_Winkler)
    plot.fill_between( x, y_min, y_max,   lw=0, zorder=1,             color=c1, alpha=0.2 )
    plot.plot        ( x, y_min,          lw=1, zorder=1, dashes=( ), color=c1            )
    plot.plot        ( x, y_max,          lw=1, zorder=1, dashes=( ), color=c1            )

    save = np.zeros( (len(x),1+2*3*3) )
    
    save[:,0]=x
    save[:,1]=y + 1e-90
    save[:,2]=y_min
    save[:,3]=y_max


    x,y = profile (dir+'/sourceTerm_'+type+'/dT_'+type+'_pbar.txt',  450 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
    y_min = np.copy(y)
    y_max = np.copy(y)
    list_1s = glob.glob(dir+'/sourceTerm_'+type+'/dT_'+type+'_pbar_1s*.txt')
    j = 0
    for f in list_1s:
        j+=1
        if j>n_max:
            break
        xx, yy = profile (f,  450 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
        for i in range(len(xx)):
            if yy[i] < y_min[i]:
                y_min[i] = yy[i]
            if yy[i] > y_max[i]:
                y_max[i] = yy[i]
    y_max +=1e-90
    y_min +=1e-90
    plot.plot        ( x, y,              lw=3, zorder=1, dashes=( ), color=c1            )
    plot.fill_between( x, y_min, y_max,   lw=0, zorder=1,             color=c1, alpha=0.2 )
    plot.plot        ( x, y_min,          lw=1, zorder=1, dashes=( ), color=c1            )
    plot.plot        ( x, y_max,          lw=1, zorder=1, dashes=( ), color=c1            )
    
    save[:,4]=y + 1e-90
    save[:,5]=y_min
    save[:,6]=y_max
    
    x,y = profile (dir+'/sourceTerm_'+type+'/dT_'+type+'_pbar.txt', 6500 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
    y_min = np.copy(y)
    y_max = np.copy(y)
    list_1s = glob.glob(dir+'/sourceTerm_'+type+'/dT_'+type+'_pbar_1s*.txt')
    j = 0
    for f in list_1s:
        j+=1
        if j>n_max:
            break
        xx, yy = profile (f, 6500 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
        for i in range(len(xx)):
            if yy[i] < y_min[i]:
                y_min[i] = yy[i]
            if yy[i] > y_max[i]:
                y_max[i] = yy[i]
    y    *=5
    y_max*=5
    y_min*=5
    y_max +=1e-90
    y_min +=1e-90
    plot.plot        ( x, y,              lw=3, zorder=1, dashes=( ), color=c1            )
    plot.fill_between( x, y_min, y_max,   lw=0, zorder=1,             color=c1, alpha=0.2 )
    plot.plot        ( x, y_min,          lw=1, zorder=1, dashes=( ), color=c1            )
    plot.plot        ( x, y_max,          lw=1, zorder=1, dashes=( ), color=c1            )


    save[:,7]=y + 1e-90
    save[:,8]=y_min
    save[:,9]=y_max

    dir=dir_diMauro

    x,y = profile (dir+'/sourceTerm_'+type+'/dT_'+type+'_pbar.txt',   20 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
    y_min = np.copy(y)
    y_max = np.copy(y)
    list_1s = glob.glob(dir+'/sourceTerm_'+type+'/dT_'+type+'_pbar_1s*.txt')
    j = 0
    for f in list_1s:
        j+=1
        if j>n_max:
            break
        xx, yy = profile (f,   20 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
        for i in range(len(xx)):
            if yy[i] < y_min[i]:
                y_min[i] = yy[i]
            if yy[i] > y_max[i]:
                y_max[i] = yy[i]
    y_max +=1e-90
    y_min +=1e-90
    plot.plot        ( x, y,              lw=3, zorder=2, dashes=(7,4), color=c2, label=label_diMauro)
    plot.fill_between( x, y_min, y_max,   lw=0, zorder=2,             color=c2, alpha=0.2 )
    plot.plot        ( x, y_min,          lw=1, zorder=2, dashes=(21,12), color=c2            )
    plot.plot        ( x, y_max,          lw=1, zorder=2, dashes=(21,12), color=c2            )


    save[:,10]=y
    save[:,11]=y_min
    save[:,12]=y_max

    x,y = profile (dir+'/sourceTerm_'+type+'/dT_'+type+'_pbar.txt',  450 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
    y_min = np.copy(y)
    y_max = np.copy(y)
    list_1s = glob.glob(dir+'/sourceTerm_'+type+'/dT_'+type+'_pbar_1s*.txt')
    j = 0
    for f in list_1s:
        j+=1
        if j>n_max:
            break
        xx, yy = profile (f,  450 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
        for i in range(len(xx)):
            if yy[i] < y_min[i]:
                y_min[i] = yy[i]
            if yy[i] > y_max[i]:
                y_max[i] = yy[i]
    y_max +=1e-90
    y_min +=1e-90
    plot.plot        ( x, y,              lw=3, zorder=2, dashes=(7,4), color=c2            )
    plot.fill_between( x, y_min, y_max,   lw=0, zorder=2,             color=c2, alpha=0.2 )
    plot.plot        ( x, y_min,          lw=1, zorder=2, dashes=(21,12), color=c2            )
    plot.plot        ( x, y_max,          lw=1, zorder=2, dashes=(21,12), color=c2            )


    save[:,13]=y
    save[:,14]=y_min
    save[:,15]=y_max

    x,y = profile (dir+'/sourceTerm_'+type+'/dT_'+type+'_pbar.txt', 6500 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
    y_min = np.copy(y)
    y_max = np.copy(y)
    list_1s = glob.glob(dir+'/sourceTerm_'+type+'/dT_'+type+'_pbar_1s*.txt')
    j = 0
    for f in list_1s:
        j+=1
        if j>n_max:
            break
        xx, yy = profile (f, 6500 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
        for i in range(len(xx)):
            if yy[i] < y_min[i]:
                y_min[i] = yy[i]
            if yy[i] > y_max[i]:
                y_max[i] = yy[i]
    y    *=5
    y_max*=5
    y_min*=5
    y_max +=1e-90
    y_min +=1e-90
    plot.plot        ( x, y,              lw=3, zorder=1, dashes=(7,4), color=c2            )
    plot.fill_between( x, y_min, y_max,   lw=0, zorder=1,             color=c2, alpha=0.2 )
    plot.plot        ( x, y_min,          lw=1, zorder=1, dashes=(21,12), color=c2            )
    plot.plot        ( x, y_max,          lw=1, zorder=1, dashes=(21,12), color=c2            )


    save[:,16]=y
    save[:,17]=y_min
    save[:,18]=y_max


    CRACS = os.getenv('CRACS')

    lim_y = 1.

    if type=='1-1' or type=='4-1' or type=='1-4':
        
        AA='pp'
        proj = 'T_{p}'
        if type=='4-1':
            AA='Hep'
            lim_y=4.
            proj = r'T_{\rm He}/n'
        if type=='1-4':
            AA='pHe'
            lim_y=4.
        
        if type=='1-1':
            x,y = profile (CRACS+'/data/CS_lab/dT_'+AA+'_pbar_LAB_diMauro12.txt',   20 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
            plot.plot( x, y, lw=2, zorder=0, dashes=(4,4           ), color=c0, label='di Mauro')
            x,y = profile (CRACS+'/data/CS_lab/dT_'+AA+'_pbar_LAB_diMauro12.txt',  450 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
            plot.plot( x, y, lw=2, zorder=0, dashes=(4,4           ), color=c0  )
            x,y = profile (CRACS+'/data/CS_lab/dT_'+AA+'_pbar_LAB_diMauro12.txt', 6500 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=5. )
            plot.plot( x, y, lw=2, zorder=0, dashes=(4,4           ), color=c0 )
        
        x,y = profile (CRACS+'/data/CS_lab/dT_'+AA+'_pbar_LAB_Winkler.txt',   20 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
        plot.plot( x, y, lw=2, zorder=0, dashes=(17,4,4,4       ), color=c0, label='Winkler')
        x,y = profile (CRACS+'/data/CS_lab/dT_'+AA+'_pbar_LAB_Winkler.txt',  450 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
        plot.plot( x, y, lw=2, zorder=0, dashes=(17,4,4,4       ), color=c0 )
        x,y = profile (CRACS+'/data/CS_lab/dT_'+AA+'_pbar_LAB_Winkler.txt', 6500 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=5. )
        plot.plot( x, y, lw=2, zorder=0, dashes=(17,4,4,4       ), color=c0 )

        x,y = profile (CRACS+'/data/CS_lab/dT_'+AA+'_pbar_LAB_KR_PPFRAG.txt',   20 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
        plot.plot( x, y, lw=2, zorder=0, dashes=(12,3,3,3,3,3   ), color=c0, label='KMO')
        x,y = profile (CRACS+'/data/CS_lab/dT_'+AA+'_pbar_LAB_KR_PPFRAG.txt',  450 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
        plot.plot( x, y, lw=2, zorder=0, dashes=(12,3,3,3,3,3   ), color=c0 )
        x,y = profile (CRACS+'/data/CS_lab/dT_'+AA+'_pbar_LAB_KR_PPFRAG.txt', 6500 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=5. )
        plot.plot( x, y, lw=2, zorder=0, dashes=(12,3,3,3,3,3   ), color=c0 )
        

    plot.text    ( 0.8e1, lim_y*2e-34 , r'$\mathrm{'+proj+r'=                   20\,GeV}$',   color=c0,     rotation=-80 )
    plot.text    (  15e1, lim_y*2e-34 , r'$\mathrm{'+proj+r'=                  450\,GeV}$',   color=c0,     rotation=-80 )
    plot.text    (   1e3, lim_y*1e-33 , r'$\mathrm{'+proj+r'=  6.5\,TeV\,\,(\times\,5)   }$', color=c0,     rotation=-70 )

    plot.set_xlim( (5e-1,  1e+4 ) )
    plot.set_ylim( (4e-36*lim_y, 4e-31*lim_y) )

    handles,   labels   = plot  .get_legend_handles_labels()
    leg1 = plot  .legend( ((handles)[:2])[::-1]+(handles)[2:], ((labels)[:2])[::-1]+(labels)[2:], loc='upper right', bbox_to_anchor=(0.97, 0.97), ncol=1, frameon=False, fontsize=0.9*label_size)

    plt.savefig('XS_comparison_'+prefix+'_'+type+'.pdf')


    plot, fig = plot_1D ( r'$\mathrm{T_{\bar{p}} \quad [GeV]}$', r'$\mathrm{d\sigma/dT_{\bar{p}}\quad[m^2/GeV]}$', 'log', 'linear', 1.3, 0.5 )

    plot.plot        ( save[:,0], save[:,1]/save[:,1],                       lw=3, zorder=1, dashes=( ), color=c1, label=label_Winkler)
    #plot.fill_between( save[:,0], save[:,2]/save[:,1], save[:,3]/save[:,1],  lw=0, zorder=1,             color=c1, alpha=0.2 )
    plot.plot        ( save[:,0], save[:,2]/save[:,1],                       lw=1, zorder=1, dashes=( ), color=c1            )
    plot.plot        ( save[:,0], save[:,3]/save[:,1],                       lw=1, zorder=1, dashes=( ), color=c1            )

    plot.plot        ( save[:,0], save[:,10]/save[:,1],                       lw=3, zorder=2, dashes=( 7,4 ), color=c2, label=label_diMauro)
    #plot.fill_between( save[:,0], save[:,11]/save[:,1], save[:,11]/save[:,1], lw=0, zorder=2,                 color=c2, alpha=0.2 )
    plot.plot        ( save[:,0], save[:,11]/save[:,1],                       lw=1, zorder=2, dashes=(21,12), color=c2            )
    plot.plot        ( save[:,0], save[:,12]/save[:,1],                       lw=1, zorder=2, dashes=(21,12), color=c2            )

    if type=='1-1' or type=='4-1' or type=='1-4':

        AA='pp'
        proj = 'T_{p}'
        if type=='4-1':
            AA='Hep'
            lim_y=4.
            proj = r'T_{\rm He}/n'
        if type=='1-4':
            AA='pHe'
            lim_y=4.

        if type=='1-1':
            x,y = profile (CRACS+'/data/CS_lab/dT_'+AA+'_pbar_LAB_diMauro12.txt',   20 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
            plot.plot( x, y/save[:,1], lw=2, zorder=0, dashes=(4,4           ), color=c0, label='di Mauro')

        x,y = profile (CRACS+'/data/CS_lab/dT_'+AA+'_pbar_LAB_Winkler.txt',   20 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
        plot.plot( x, y/save[:,1], lw=2, zorder=0, dashes=(17,4,4,4       ), color=c0, label='Winkler')

        x,y = profile (CRACS+'/data/CS_lab/dT_'+AA+'_pbar_LAB_KR_PPFRAG.txt',   20 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
        plot.plot( x, y/save[:,1], lw=2, zorder=0, dashes=(12,3,3,3,3,3   ), color=c0, label='KMO')

    plot.text    ( 8e-1, 0.6 , r'$\mathrm{'+proj+r' = 20\,GeV}$',   color=c0,     rotation=0 )

    plot.set_xlim( (5e-1,  1e+4 ) )
    plot.set_ylim( ( 0.5, 1.5) )

    plt.savefig('XS_comparison_ratio_20_'+prefix+'_'+type+'.pdf')

    plot, fig = plot_1D ( r'$\mathrm{T_{\bar{p}} \quad [GeV]}$', r'$\mathrm{d\sigma/dT_{\bar{p}}\quad[m^2/GeV]}$', 'log', 'linear', 1.3, 0.5 )

    plot.plot        ( save[:,0], save[:,4]/save[:,4],                       lw=3, zorder=1, dashes=( ), color=c1, label=label_Winkler)
    #plot.fill_between( save[:,0], save[:,5]/save[:,4], save[:,6]/save[:,1],  lw=0, zorder=1,             color=c1, alpha=0.2 )
    plot.plot        ( save[:,0], save[:,5]/save[:,4],                       lw=1, zorder=1, dashes=( ), color=c1            )
    plot.plot        ( save[:,0], save[:,6]/save[:,4],                       lw=1, zorder=1, dashes=( ), color=c1            )

    plot.plot        ( save[:,0], save[:,13]/save[:,4],                       lw=3, zorder=2, dashes=( 7,4 ), color=c2, label=label_diMauro)
    #plot.fill_between( save[:,0], save[:,14]/save[:,4], save[:,15]/save[:,1], lw=0, zorder=2,                 color=c2, alpha=0.2 )
    plot.plot        ( save[:,0], save[:,14]/save[:,4],                       lw=1, zorder=2, dashes=(21,12), color=c2            )
    plot.plot        ( save[:,0], save[:,15]/save[:,4],                       lw=1, zorder=2, dashes=(21,12), color=c2            )

    if type=='1-1' or type=='4-1' or type=='1-4':

        AA='pp'
        proj = 'T_{p}'
        if type=='4-1':
            AA='Hep'
            lim_y=4.
            proj = r'T_{\rm He}/n'
        if type=='1-4':
            AA='pHe'
            lim_y=4.

        if type=='1-1':
            x,y = profile (CRACS+'/data/CS_lab/dT_'+AA+'_pbar_LAB_diMauro12.txt',  450 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
            plot.plot( x, y/save[:,1], lw=2, zorder=0, dashes=(4,4           ), color=c0, label='di Mauro')

        x,y = profile (CRACS+'/data/CS_lab/dT_'+AA+'_pbar_LAB_Winkler.txt',  450 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
        plot.plot( x, y/save[:,1], lw=2, zorder=0, dashes=(17,4,4,4       ), color=c0, label='Winkler')

        x,y = profile (CRACS+'/data/CS_lab/dT_'+AA+'_pbar_LAB_KR_PPFRAG.txt',  450 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=1. )
        plot.plot( x, y/save[:,1], lw=2, zorder=0, dashes=(12,3,3,3,3,3   ), color=c0, label='KMO')

    plot.text    ( 8e-1, 0.6 , r'$\mathrm{'+proj+r' = 450\,GeV}$',   color=c0,     rotation=0 )

    plot.set_xlim( (5e-1,  1e+4 ) )
    plot.set_ylim( ( 0.5, 1.5) )

    plt.savefig('XS_comparison_ratio_450_'+prefix+'_'+type+'.pdf')

    plot, fig = plot_1D ( r'$\mathrm{T_{\bar{p}} \quad [GeV]}$', r'$\mathrm{d\sigma/dT_{\bar{p}}\quad[m^2/GeV]}$', 'log', 'linear', 1.3, 0.5 )

    plot.plot        ( save[:,0], save[:,7]/save[:,7],                       lw=3, zorder=1, dashes=( ), color=c1, label=label_Winkler)
    #plot.fill_between( save[:,0], save[:,8]/save[:,7], save[:,9]/save[:,7],  lw=0, zorder=1,             color=c1, alpha=0.2 )
    plot.plot        ( save[:,0], save[:,8]/save[:,7],                       lw=1, zorder=1, dashes=( ), color=c1            )
    plot.plot        ( save[:,0], save[:,9]/save[:,7],                       lw=1, zorder=1, dashes=( ), color=c1            )

    plot.plot        ( save[:,0], save[:,16]/save[:,7],                       lw=3, zorder=2, dashes=( 7,4 ), color=c2, label=label_diMauro)
    #plot.fill_between( save[:,0], save[:,17]/save[:,7], save[:,18]/save[:,7], lw=0, zorder=2,                 color=c2, alpha=0.2 )
    plot.plot        ( save[:,0], save[:,17]/save[:,7],                       lw=1, zorder=2, dashes=(21,12), color=c2            )
    plot.plot        ( save[:,0], save[:,18]/save[:,7],                       lw=1, zorder=2, dashes=(21,12), color=c2            )

    if type=='1-1' or type=='4-1' or type=='1-4':

        AA='pp'
        proj = 'T_{p}'
        if type=='4-1':
            AA='Hep'
            lim_y=4.
            proj = r'T_{\rm He}/n'
        if type=='1-4':
            AA='pHe'
            lim_y=4.

        if type=='1-1':
            x,y = profile (CRACS+'/data/CS_lab/dT_'+AA+'_pbar_LAB_diMauro12.txt', 6500 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=5. )
            plot.plot( x, y/save[:,1], lw=2, zorder=0, dashes=(4,4           ), color=c0, label='di Mauro')

        x,y = profile (CRACS+'/data/CS_lab/dT_'+AA+'_pbar_LAB_Winkler.txt', 6500 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=5. )
        plot.plot( x, y/save[:,1], lw=2, zorder=0, dashes=(17,4,4,4       ), color=c0, label='Winkler')

        x,y = profile (CRACS+'/data/CS_lab/dT_'+AA+'_pbar_LAB_KR_PPFRAG.txt', 6500 , type='proj', Tn_min=1e-1, Tn_max=1e5, y_fac=5. )
        plot.plot( x, y/save[:,1], lw=2, zorder=0, dashes=(12,3,3,3,3,3   ), color=c0, label='KMO')

    plot.text    ( 8e-1, 0.6 , r'$\mathrm{'+proj+r' = 6.5\,TeV}$',   color=c0,     rotation=0 )

    plot.set_xlim( (5e-1,  1e+4 ) )
    plot.set_ylim( ( 0.5, 1.5) )

    plt.savefig('XS_comparison_ratio_6500_'+prefix+'_'+type+'.pdf')


def plot_sourceterm_comparison(prefix, dir_Winkler, dir_diMauro, label_Winkler='Param. Winkler', label_diMauro='Param. di Mauro', type='1-1'):

    sizex=1.3
    sizey=1.3

    plt.close('all')
    font_props = {"size":label_size}
    rc("font", **font_props)
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['mathtext.fontset']='stixsans'

    fig     = plt.figure(figsize=(print_size*sizex, print_size*sizey))
    plot    = plt.subplot2grid((5, 1), (0, 0), rowspan=3)
    plot_u  = plt.subplot2grid((5, 1), (3, 0), rowspan=2)
    plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

    plot.set_ylabel ( r'$\mathrm{q^{(\bar{p})} \, T_{\bar{p}}^{2.7}  \,  \, [GeV^{1.7}m^{-3}s^{-1}]}}$'     )
    plot.tick_params(labelbottom='off')

    plot.set_xscale ( 'log' )
    plot.set_yscale ( 'log' )

    plot.tick_params('both', length=20, width=2, which='major', top=True, right=True, direction='in', pad=10)
    plot.tick_params('both', length=10, width=1, which='minor', top=True, right=True, direction='in', pad=10)

    plot.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
    
    plot_u.set_xlabel ( r'$\mathrm{T_{\bar{p}} [GeV]}$'                                                       )
    plot_u.set_ylabel ( r'$k$'     )

    plot_u.set_xscale ( 'log' )
    plot_u.set_yscale ( 'linear' )

    plot_u.tick_params('both', length=20, width=2, which='major', top=True, right=True, direction='in', pad=10)
    plot_u.tick_params('both', length=10, width=1, which='minor', top=True, right=True, direction='in', pad=10)

    plot_u.grid(b=True, which='major', axis='x', alpha=0.1, linestyle='-', linewidth=2)
    
    dir=dir_Winkler

    pp_source = np.genfromtxt(dir+'/sourceTerm_'+type+'/sourceTerm_'+type+'_pbar.txt', skip_header=1)

    y1_0 =pp_source[:, 1]*np.power(pp_source[:,0],2.7)

    if type=='1-1' or type=='4-1' or type=='1-4':
        y1_1 =pp_source[:, 2]*np.power(pp_source[:,0],2.7)
        y1_2 =pp_source[:, 4]*np.power(pp_source[:,0],2.7)
        y1_3 =pp_source[:, 3]*np.power(pp_source[:,0],2.7)


    pp_source_min = np.copy( pp_source[:, 1] )
    pp_source_max = np.copy( pp_source[:, 1] )

    j = 0
    list_1s = glob.glob(dir+'/sourceTerm_'+type+'/sourceTerm_'+type+'_pbar_1s*.txt')
    for f in list_1s:
        source = np.genfromtxt(f, skip_header=1)
#        plot    .plot( source[:,0], source[:, 1]*np.power(source[:,0],2.7) , color=c3, lw=1, alpha=0.5)
#        if j<5:
#            j += 1
#        plot_u  .plot( source[:,0], source[:, 1]/pp_source[:,1] ,            color=c3, lw=1, alpha=0.5)
        for i in range(len(source[:,1])):
            if source[i,1] < pp_source_min[i]:
                pp_source_min[i] = source[i,1]
            if source[i,1] > pp_source_max[i]:
                pp_source_max[i] = source[i,1]


    plot.   plot        ( pp_source[:,0], pp_source[:, 1]   *np.power(pp_source[:,0],2.7) ,                                        lw=3, zorder=1, dashes=( ), color=c1, label=label_Winkler)

    plot.   fill_between( pp_source[:,0], pp_source_min     *np.power(source[:,0],2.7), pp_source_max*np.power(source[:,0],2.7),   lw=0, zorder=1,             color=c1, alpha=0.2 )
    plot.   plot        ( pp_source[:,0], pp_source_min     *np.power(source[:,0],2.7) ,                                           lw=1, zorder=1, dashes=( ), color=c1            )
    plot.   plot        ( pp_source[:,0],                                               pp_source_max*np.power(source[:,0],2.7),   lw=1, zorder=1, dashes=( ), color=c1            )

    plot_u. plot        ( pp_source[:,0], pp_source[:, 1]   /pp_source[:,1],                                                       lw=3, zorder=1, dashes=( ), color=c1)
    plot_u. fill_between( pp_source[:,0], pp_source_min     /pp_source[:,1],            pp_source_max/pp_source[:,1],              lw=0, zorder=1, color=c1, alpha=0.2 )
    plot_u. plot        ( pp_source[:,0], pp_source_min     /pp_source[:,1],                                                       lw=1, zorder=1, dashes=( ), color=c1            )
    plot_u. plot        ( pp_source[:,0],                                               pp_source_max/pp_source[:,1],              lw=1, zorder=1, dashes=( ), color=c1            )

    pp_source_min2 = np.copy( pp_source_min )
    pp_source_max2 = np.copy( pp_source_max )

    norm = np.copy( pp_source[:,1] )


    dir=dir_diMauro

    pp_source = np.genfromtxt(dir+'/sourceTerm_'+type+'/sourceTerm_'+type+'_pbar.txt', skip_header=1)

    pp_source_min = np.copy( pp_source[:, 1] )
    pp_source_max = np.copy( pp_source[:, 1] )

    j = 0
    list_1s = glob.glob(dir+'/sourceTerm_'+type+'/sourceTerm_'+type+'_pbar_1s*.txt')
    for f in list_1s:
        source = np.genfromtxt(f, skip_header=1)
        #plot    .plot( source[:,0], source[:, 1]*np.power(source[:,0],2.7) , color=c3, lw=1, alpha=0.5)
        #if j<5:
        #    j += 1
        #    plot_u  .plot( source[:,0], source[:, 1]/norm ,            color=c3, lw=1, alpha=0.5)
        for i in range(len(source[:,1])):
            if source[i,1] < pp_source_min[i]:
                pp_source_min[i] = source[i,1]
            if source[i,1] > pp_source_max[i]:
                pp_source_max[i] = source[i,1]


    plot.   plot        ( pp_source[:,0], pp_source[:, 1]   *np.power(pp_source[:,0],2.7) ,                                        lw=3, zorder=2, dashes=(7,4), color=c2, label=label_diMauro)
    plot_u. plot        ( pp_source[:,0], pp_source[:, 1]   /norm                         ,                                        lw=3, zorder=2, dashes=(7,4), color=c2, label='')
    plot.   fill_between( pp_source[:,0], pp_source_min     *np.power(source[:,0],2.7), pp_source_max*np.power(source[:,0],2.7),   lw=0, zorder=2, color=c2, alpha=0.2 )
    plot.   plot        ( pp_source[:,0], pp_source_min     *np.power(source[:,0],2.7) ,                                           lw=1, zorder=2, dashes=(21,8), color=c2            )
    plot.   plot        ( pp_source[:,0],                                               pp_source_max*np.power(source[:,0],2.7),   lw=1, zorder=2, dashes=(21,8), color=c2            )

    plot_u. fill_between( pp_source[:,0], pp_source_min     /norm,                      pp_source_max/norm,                        lw=0, zorder=2, color=c2, alpha=0.2 )
    plot_u. plot        ( pp_source[:,0], pp_source_min     /norm,                                                                 lw=1, zorder=2, dashes=(21,8), color=c2            )
    plot_u. plot        ( pp_source[:,0],                                               pp_source_max/norm,                        lw=1, zorder=2, dashes=(21,8), color=c2            )

    pp_source_min2 = np.copy( pp_source_min )
    pp_source_max2 = np.copy( pp_source_max )


    plot.set_xlim( (5e-1, 5e3)         )
    plot.set_ylim( (1e-23, 1e-21)   )

    if type=='1-1' or type=='4-1' or type=='1-4':
        if type=='1-1':
            plot.plot(pp_source[:,0], y1_1, lw=2, color=c0, zorder=0, dashes=(4,4))
        plot.plot(pp_source[:,0], y1_2, lw=2, color=c0, zorder=0, dashes=(17,4,4,4))
        plot.plot(pp_source[:,0], y1_3, lw=2, color=c0, zorder=0, dashes=(12,3,3,3,3,3))

        if type=='1-1':
            plot_u.plot(pp_source[:,0], y1_1/y1_0, lw=2, label='di Mauro', color=c0, zorder=3, dashes=(4,4)  )
        plot_u.plot(pp_source[:,0], y1_2/y1_0, lw=2, label='Winkler' , color=c0, zorder=3, dashes=(17,4,4,4))
        plot_u.plot(pp_source[:,0], y1_3/y1_0, lw=2, label='KMO'     , color=c0, zorder=3, dashes=(12,4,4,4,4,4))

    #plot_u.set_yticks( [0.5, 1.0, 1.5] )

    plot_u.set_xlim( (5e-1, 5e3)         )
    plot_u.set_ylim( ( 0.7, 1.3)         )

    handles,   labels   = plot  .get_legend_handles_labels()

    leg1 = plot  .legend( (handles)[::-1], (labels)[::-1], loc='upper left', bbox_to_anchor=(0.03, 0.95), ncol=1, frameon=False, fontsize=0.9*label_size)

    if type=='1-1' or type=='4-1' or type=='1-4':
        leg2 = plot_u.legend(                                                     loc='lower right', bbox_to_anchor=(0.97 , 1.20), ncol=1, frameon=False, fontsize=0.9*label_size)


    plt.savefig('sourceterm_comparison_'+prefix+'_'+type+'.pdf')





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


def get_min_array( list ):
    y_min = np.zeros(len(list[0]))
    new = 0
    for l in list:
        for j in range(len(l)):
            if (l[j]<y_min[j] and l[j]>1e-90) or y_min[j]<1e-90:
                y_min[j] = l[j]
    return y_min

def get_max_array( list ):
    y_max = np.zeros(len(list[0]))
    new = 0
    for l in list:
        for j in range(len(l)):
            if ( l[j]>y_max[j] ):
                y_max[j] = l[j]
    return y_max


from numpy import linalg as LA
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator


if step==1:
    bestFit = glob.glob(dir+'/bestFit_*.txt')
    for f in bestFit:
        draw(f, 'fit_')
    comparison = glob.glob(dir+'/comparison_*.txt')
    for f in comparison:
        draw(f, 'comparison_')

if step==10:
    bV = np.genfromtxt(dir+'/covariance_matrix_pp.txt')
    b  = bV[:,0]
    V  = bV[:,1:]
    iV = LA.inv(V)
    
    e_val, RT = LA.eigh(iV)
    R  = np.transpose(RT)
    dim = len(e_val)
    
    fSigma_1_region = [1.000043426, 2.295815161, 3.526822180, 4.719568761, 5.887700474, 7.038515492, 8.176359741, 9.304044023, 10.42350189, 11.53612748, 12.64296378, 13.74481437, 14.84231367, 15.93597259, 17.02620974, 18.11337308, 19.19775552, 20.27960639, 21.35913992, 22.43654182];
    fSigma_2_region = [4.000009776, 6.180085905, 8.024894670, 9.715641142, 11.31387084, 12.84885057, 14.33712678, 15.78911024, 17.21184692, 18.61036514, 19.98840090, 21.34881933, 22.69387460, 24.02537791, 25.34481051, 26.65340209, 27.95218680, 29.24204418, 30.52372968, 31.79789790];

    count_1s = 0;
    count_2s = 0;
    
    f_all = open( 'parameters_pp.txt', 'w' )
    for i in range(dim):
        diff = RT[:,i]
        chiSq = np.vdot(   diff, np.dot(iV, diff)  )
    for n in range(1000):
        diff = np.zeros(dim)
        rand = np.random.normal(0, 1, dim)
        for i in range(dim):
            diff = diff + 1./np.sqrt(e_val[i])*rand[i]*RT[:,i]
        chiSq = np.vdot(   diff, np.dot(iV, diff)  )

        val = b + diff
        write = ''
        for v in val:
            write += str(v).ljust(20) + ' '
        f_all.write(write+'\n')
    f_all.close()


if step==2:
    bestFit = glob.glob('bestFit_*.txt')
    for f in bestFit:
        draw(f, 'fit_')
    comparison = glob.glob('comparison_*.txt')
    for f in comparison:
        draw(f, 'comparison_')
    draw_lhcb()


if step==20:
    bV_pp = np.genfromtxt(dir+'/covariance_matrix_pp.txt')
    bV_pA = np.genfromtxt(dir+'/covariance_matrix_pA.txt')
    
    dim_pp = len(bV_pp[:,0])
    dim_pA = len(bV_pA[:,0])
    dim    = dim_pA + dim_pp
    bV    = np.zeros( (dim,dim+1) )
    
    bV[:dim_pp,:dim_pp+1] = bV_pp
    bV[dim_pp:,dim_pp+1:] = bV_pA[:,1:]
    bV[dim_pp:,0] = bV_pA[:,0]
    
    b  = bV[:,0]
    V  = bV[:,1:]
    iV = LA.inv(V)
    e_val, RT = LA.eigh(iV)
    R  = np.transpose(RT)
    
    dim = len(e_val)
    fSigma_1_region = [1.000043426, 2.295815161, 3.526822180, 4.719568761, 5.887700474, 7.038515492, 8.176359741, 9.304044023, 10.42350189, 11.53612748, 12.64296378, 13.74481437, 14.84231367, 15.93597259, 17.02620974, 18.11337308, 19.19775552, 20.27960639, 21.35913992, 22.43654182];
    fSigma_2_region = [4.000009776, 6.180085905, 8.024894670, 9.715641142, 11.31387084, 12.84885057, 14.33712678, 15.78911024, 17.21184692, 18.61036514, 19.98840090, 21.34881933, 22.69387460, 24.02537791, 25.34481051, 26.65340209, 27.95218680, 29.24204418, 30.52372968, 31.79789790];
    
    count_1s = 0;
    count_2s = 0;
    
    f_all = open( 'parameters_pA.txt', 'w' )
    for i in range(dim):
        diff = RT[:,i]
        chiSq = np.vdot(   diff, np.dot(iV, diff)  )
    for n in range(1000):
        diff = np.zeros(dim)
        rand = np.random.normal(0, 1, dim)
        for i in range(dim):
            diff = diff + 1./np.sqrt(e_val[i])*rand[i]*RT[:,i]

        chiSq = np.vdot(   diff, np.dot(iV, diff)  )

        val = b + diff
        write = ''
        for v in val:
            write += str(v).ljust(20) + ' '
        f_all.write(write+'\n')
    f_all.close()


if step==30:
    d = dir
    if cluster:
        d = '.'
    if A1==1 and A2==1:
        f = np.genfromtxt(d+'/parameters_1_sigma.txt')
    else:
        f = np.genfromtxt(d+'/parameters_pA_1_sigma.txt')
    l = len(f[:,0])
    if not ( (A1==1 and A2==1) or (A1==1 and A2==4) or (A1==4 and A2==1) ):
        l = np.minimum(l,100)
    write='#!/bin/bash \n\n'
    for i in range(0,l):
        s = sdir+'/bin/CRACS_crossSectionFitter --step 30 --npoint '+str(i)+' --cs ' + str(cs)+' --data ' + str(data)+' --data_pHe ' + str(data_pHe) +' --data_pC ' + str(data_pC) +' --A1 ' + str(A1) +' --A2 ' + str(A2) + ' --dir ' + dir + ' --softwarePath ' + sdir
        print s
        if cluster:
            write += s+'\n'
        else:
            os.system(s)
    if cluster:
        os.system( 'mkdir cluster' )
        file = open( 'cluster/run_sourceTerm_'+str(A1)+'-'+str(A2)+'.sh', 'w' )
        file.write(write)
        file.close()
        os.system( 'echo '+dir+'/cluster/run_sourceTerm_'+str(A1)+'-'+str(A2)+'.sh  >> runlist.txt' )



if step==3:

    sizex=1.3
    sizey=1.3

    plt.close('all')
    font_props = {"size":label_size}
    rc("font", **font_props)
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['mathtext.fontset']='stixsans'

    fig     = plt.figure(figsize=(print_size*sizex, print_size*sizey))
    plot    = plt.subplot2grid((5, 1), (0, 0), rowspan=3)
    plot_u  = plt.subplot2grid((5, 1), (3, 0), rowspan=2)
    plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)

    #plot.set_xlabel ( r'$\mathrm{T_{\bar{p}} [GeV]}$'                                                       )
    plot.set_ylabel ( r'$\mathrm{q^{(\bar{p})} \, T_{\bar{p}}^{2.7}  \,  \, [GeV^{1.7}m^{-3}s^{-1}]}}$'     )
    plot.tick_params(labelbottom='off')

    plot.set_xscale ( 'log' )
    plot.set_yscale ( 'log' )

    plot.tick_params('both', length=20, width=2, which='major', top=True, right=True, direction='in', pad=10)
    plot.tick_params('both', length=10, width=1, which='minor', top=True, right=True, direction='in', pad=10)

    plot.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
    
    plot_u.set_xlabel ( r'$\mathrm{T_{\bar{p}} [GeV]}$'                                                       )
    plot_u.set_ylabel ( r'$k$'     )

    plot_u.set_xscale ( 'log' )
    plot_u.set_yscale ( 'linear' )

    plot_u.tick_params('both', length=20, width=2, which='major', top=True, right=True, direction='in', pad=10)
    plot_u.tick_params('both', length=10, width=1, which='minor', top=True, right=True, direction='in', pad=10)

    plot_u.grid(b=True, which='major', axis='x', alpha=0.1, linestyle='-', linewidth=2)
    plot_u.axhline(y=1.0, color='black', alpha=1, lw=2)
    
    AA = str(A1) + '-' + str(A2)

    pp_source = np.genfromtxt(dir+'/sourceTerm_'+AA+'/sourceTerm_'+AA+'_pbar.txt', skip_header=1)

    y1_0 =pp_source[:, 1]*np.power(pp_source[:,0],2.7)
    if A1==1 and A2==1 or A1==1 and A2==4 or A1==4 and A2==1:
        y1_1 =pp_source[:, 2]*np.power(pp_source[:,0],2.7)
        y1_2 =pp_source[:, 4]*np.power(pp_source[:,0],2.7)
        y1_3 =pp_source[:, 3]*np.power(pp_source[:,0],2.7)


    pp_source_min = np.copy( pp_source[:, 1] )
    pp_source_max = np.copy( pp_source[:, 1] )

    j = 0
    list_1s = glob.glob('sourceTerm_'+AA+'/sourceTerm_'+AA+'_pbar_1s*.txt')
    for f in list_1s:
        source = np.genfromtxt(f, skip_header=1)
        #plot    .plot( source[:,0], source[:, 1]*np.power(source[:,0],2.7) , color=c3, lw=1, alpha=0.5)
        if j<10:
            j += 1
            plot_u  .plot( source[:,0], source[:, 1]/pp_source[:,1] ,            color=c0, lw=1, alpha=0.5)
        for i in range(len(source[:,1])):
            if source[i,1] < pp_source_min[i]:
                pp_source_min[i] = source[i,1]
            if source[i,1] > pp_source_max[i]:
                pp_source_max[i] = source[i,1]
    cs = int(cs)
    if cs==0:
        label = 'Param. I'
    else:
        label = 'Param. II'
    print label
    plot.   plot        ( pp_source[:,0], pp_source[:, 1]   *np.power(pp_source[:,0],2.7) ,                                        lw=3, color=c0, label='$pp$ '+label)
    plot.   fill_between( pp_source[:,0], pp_source_min     *np.power(source[:,0],2.7), pp_source_max*np.power(source[:,0],2.7),   lw=0, color=c0, alpha=0.1 )
    plot_u. fill_between( pp_source[:,0], pp_source_min     /pp_source[:,1],            pp_source_max/pp_source[:,1],              lw=0, color=c0, alpha=0.1 )

    plot.set_xlim( (5e-1, 5e3)         )
    plot.set_ylim( (1e-23, 1e-21)   )

    if A1==1 and A2==1 or A1==1 and A2==4 or A1==4 and A2==1:
        plot.plot(pp_source[:,0], y1_1, lw=2, color=c1, dashes=dashes_array[1])
        plot.plot(pp_source[:,0], y1_2, lw=2, color=c2, dashes=dashes_array[2])
        plot.plot(pp_source[:,0], y1_3, lw=2, color=c3, dashes=dashes_array[3])

        plot_u.plot(pp_source[:,0], y1_1/y1_0, lw=2, label='di Mauro', color=c1, dashes=dashes_array[1])
        plot_u.plot(pp_source[:,0], y1_2/y1_0, lw=2, label='Winkler' , color=c2, dashes=dashes_array[2])
        plot_u.plot(pp_source[:,0], y1_3/y1_0, lw=2, label='KMO'     , color=c3, dashes=dashes_array[3])

    #plot_u.set_yticks( [0.5, 1.0, 1.5] )

    plot_u.set_xlim( (5e-1, 5e3)         )
    plot_u.set_ylim( ( 0.5, 1.5)         )

    handles_u, labels_u = plot_u.get_legend_handles_labels()
    handles,   labels   = plot  .get_legend_handles_labels()

    plot.legend( handles+handles_u, labels+labels_u, loc='lower right', bbox_to_anchor=(0.95, 0.05), ncol=1, frameon=False, fontsize=0.9*label_size)


    plt.savefig('sourceterm_'+AA+'.pdf')


def get_source_term( AA='1-1', type='_tot' ):

    s_tot  = np.genfromtxt(dir+'/sourceTerm_'+AA+'/sourceTerm_'+AA+'_pbar_tot.txt',                skip_header=1)
#    s_isoL = np.genfromtxt(dir+'/sourceTerm_'+AA+'/sourceTerm_'+AA+'_pbar_tot_lowerIsoSpin.txt',   skip_header=1)
#    s_isoU = np.genfromtxt(dir+'/sourceTerm_'+AA+'/sourceTerm_'+AA+'_pbar_tot_upperIsoSpin.txt',   skip_header=1)
    s_L    = np.genfromtxt(dir+'/sourceTerm_'+AA+'/sourceTerm_'+AA+'_pbar_tot_lower.txt',          skip_header=1)
    s_U    = np.genfromtxt(dir+'/sourceTerm_'+AA+'/sourceTerm_'+AA+'_pbar_tot_upper.txt',          skip_header=1)
    
    s       = np.genfromtxt(dir+'/sourceTerm_'+AA+'/sourceTerm_'+AA+'_pbar.txt',              skip_header=1)
    s_min   = np.copy( s[:, 1] )
    s_max   = np.copy( s[:, 1] )

    list_1s = glob.glob('sourceTerm_'+AA+'/sourceTerm_'+AA+'_pbar_1s*.txt')
    for f in list_1s:
        s_1s = np.genfromtxt(f, skip_header=1)
        for i in range(len(s_1s[:,1])):
            if s_1s[i,1] < s_min[i]:
                s_min[i] = s_1s[i,1]
            if s_1s[i,1] > s_max[i]:
                s_max[i] = s_1s[i,1]

    ret = np.zeros( (len(s_tot[:,0]), 8) )
    ret[:,:2] = s_tot
    ret[:, 2] = s_tot[:,1]*s_min/s[:,1]
    ret[:, 3] = s_tot[:,1]*s_max/s[:,1]

    ret[:, 4] = s_L     [:,1]
    ret[:, 5] = s_U     [:,1]
    ret[:, 6] = s_L     [:,1]
    ret[:, 7] = s_U     [:,1]

    return ret


if step==4:
    
    sizex=1.3
    sizey=1.3
    
    plt.close('all')
    font_props = {"size":label_size}
    rc("font", **font_props)
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['mathtext.fontset']='stixsans'
    
    fig     = plt.figure(figsize=(print_size*sizex, print_size*sizey))
    plot    = plt.subplot2grid((5, 1), (0, 0), rowspan=3)
    plot_u  = plt.subplot2grid((5, 1), (3, 0), rowspan=2)
    plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
    
    plot.set_ylabel ( r'$\mathrm{q^{(\bar{p})} \, T_{\bar{p}}^{2.7}  \,  \, [GeV^{1.7}m^{-3}s^{-1}]}}$'     )
    plot.tick_params(labelbottom='off')
    
    plot.set_xscale ( 'log' )
    plot.set_yscale ( 'log' )
    
    plot.tick_params('both', length=20, width=2, which='major', top=True, right=True, direction='in', pad=10)
    plot.tick_params('both', length=10, width=1, which='minor', top=True, right=True, direction='in', pad=10)
    
    plot.grid(b=True, which='major', alpha=0.1, linestyle='-', linewidth=2)
    
    plot_u.set_xlabel ( r'$\mathrm{T_{\bar{p}} [GeV]}$'                                                       )
    plot_u.set_ylabel ( r'$k$'     )
    
    plot_u.set_xscale ( 'log' )
    plot_u.set_yscale ( 'linear' )
    
    plot_u.tick_params('both', length=20, width=2, which='major', top=True, right=True, direction='in', pad=10)
    plot_u.tick_params('both', length=10, width=1, which='minor', top=True, right=True, direction='in', pad=10)
    
    plot_u.grid(b=True, which='major', axis='x', alpha=0.1, linestyle='-', linewidth=2)
    
    f_01_01 = get_source_term( '1-1' )
    f_01_04 = get_source_term( '1-4' )
    f_04_01 = get_source_term( '4-1' )
    f_04_04 = get_source_term( '4-4' )
    
    f_12_01 = get_source_term( '12-1' )
    f_12_04 = get_source_term( '12-4' )
    f_14_01 = get_source_term( '14-1' )
    f_14_04 = get_source_term( '14-4' )
    f_16_01 = get_source_term( '16-1' )
    f_16_04 = get_source_term( '16-4' )
    
    
    f_07_01 = get_source_term(  '7-1' )
    f_07_04 = get_source_term(  '7-4' )
    f_09_01 = get_source_term(  '9-1' )
    f_09_04 = get_source_term(  '9-4' )
    f_11_01 = get_source_term( '11-1' )
    f_11_04 = get_source_term( '11-4' )
    
    f_20_01 = get_source_term( '20-1' )
    f_20_04 = get_source_term( '20-4' )
    f_24_01 = get_source_term( '24-1' )
    f_24_04 = get_source_term( '24-4' )
    f_28_01 = get_source_term( '28-1' )
    f_28_04 = get_source_term( '28-4' )
    f_56_01 = get_source_term( '56-1' )
    f_56_04 = get_source_term( '56-4' )
    
    f_01_12 = get_source_term( '1-12' )
    f_04_12 = get_source_term( '4-12' )
    f_01_14 = get_source_term( '1-14' )
    f_04_14 = get_source_term( '4-14' )
    f_01_16 = get_source_term( '1-16' )
    f_04_16 = get_source_term( '4-16' )
    f_01_20 = get_source_term( '1-20' )
    f_04_20 = get_source_term( '4-20' )
    f_01_24 = get_source_term( '1-24' )
    f_04_24 = get_source_term( '4-24' )
    f_01_28 = get_source_term( '1-28' )
    f_04_28 = get_source_term( '4-28' )
    f_01_56 = get_source_term( '1-56' )
    f_04_56 = get_source_term( '4-56' )

    
    f_onHeavy        = np.copy(f_01_12)
    f_onHeavy[:,1:] += f_04_12[:,1:]
    f_onHeavy[:,1:] += f_01_14[:,1:]
    f_onHeavy[:,1:] += f_04_14[:,1:]
    f_onHeavy[:,1:] += f_01_16[:,1:]
    f_onHeavy[:,1:] += f_04_16[:,1:]
    f_onHeavy[:,1:] += f_01_20[:,1:]
    f_onHeavy[:,1:] += f_04_20[:,1:]
    f_onHeavy[:,1:] += f_01_24[:,1:]
    f_onHeavy[:,1:] += f_04_24[:,1:]
    f_onHeavy[:,1:] += f_01_28[:,1:]
    f_onHeavy[:,1:] += f_04_28[:,1:]
    f_onHeavy[:,1:] += f_01_56[:,1:]
    f_onHeavy[:,1:] += f_04_56[:,1:]
    
    
    
    f_CNO          = np.copy(f_12_01)
    f_CNO  [:,1:] += f_12_04[:,1:]
    f_CNO  [:,1:] += f_14_01[:,1:]
    f_CNO  [:,1:] += f_14_04[:,1:]
    f_CNO  [:,1:] += f_16_01[:,1:]
    f_CNO  [:,1:] += f_16_04[:,1:]
    
    f_LiBeB          = np.copy(f_07_01)
    f_LiBeB[:,1:]   += f_07_04[:,1:]
    f_LiBeB[:,1:]   += f_09_01[:,1:]
    f_LiBeB[:,1:]   += f_09_04[:,1:]
    f_LiBeB[:,1:]   += f_11_01[:,1:]
    f_LiBeB[:,1:]   += f_11_04[:,1:]


    f_NeMgSi        = np.copy(f_20_01)
    f_NeMgSi[:,1:] += f_20_04[:,1:]
    f_NeMgSi[:,1:] += f_24_01[:,1:]
    f_NeMgSi[:,1:] += f_24_04[:,1:]
    f_NeMgSi[:,1:] += f_28_01[:,1:]
    f_NeMgSi[:,1:] += f_28_04[:,1:]

    f_Fe            = np.copy(f_56_01)
    f_Fe    [:,1:] += f_56_04[:,1:]

    f_all          = np.copy(f_CNO  )
    f_all  [:,1:] += f_LiBeB[:,1:]
    f_all  [:,1:] += f_NeMgSi[:,1:]
    f_all  [:,1:] += f_Fe   [:,1:]
    f_all  [:,1:] += f_01_01[:,1:]
    f_all  [:,1:] += f_04_01[:,1:]
    f_all  [:,1:] += f_01_04[:,1:]
    f_all  [:,1:] += f_04_04[:,1:]
    

    plot.   plot        ( f_01_01[:,0], f_01_01[:,1]     *np.power(f_01_01[:,0],2.7),                                            lw=3, zorder=1, dashes=( ), color=c1, label='$pp$')
    plot.   fill_between( f_01_01[:,0], f_01_01[:,2]     *np.power(f_01_01[:,0],2.7), f_01_01[:,3]*np.power(f_01_01[:,0],2.7),   lw=0, zorder=1,             color=c1, alpha=0.2 )
    plot.   plot        ( f_01_01[:,0], f_01_01[:,2]     *np.power(f_01_01[:,0],2.7),                                            lw=1, zorder=1, dashes=( ), color=c1            )
    plot.   plot        ( f_01_01[:,0],                                               f_01_01[:,3]*np.power(f_01_01[:,0],2.7),   lw=1, zorder=1, dashes=( ), color=c1            )

    plot.   plot        ( f_01_04[:,0], f_01_04[:,1]     *np.power(f_01_04[:,0],2.7),                                            lw=3, zorder=1, dashes=( ), color=c2, label='$pp$')
    plot.   fill_between( f_01_04[:,0], f_01_04[:,2]     *np.power(f_01_04[:,0],2.7), f_01_04[:,3]*np.power(f_01_04[:,0],2.7),   lw=0, zorder=1,             color=c2, alpha=0.2 )
    plot.   plot        ( f_01_04[:,0], f_01_04[:,2]     *np.power(f_01_04[:,0],2.7),                                            lw=1, zorder=1, dashes=( ), color=c2            )
    plot.   plot        ( f_01_04[:,0],                                               f_01_04[:,3]*np.power(f_01_04[:,0],2.7),   lw=1, zorder=1, dashes=( ), color=c2            )

    plot.   plot        ( f_04_01[:,0], f_04_01[:,1]     *np.power(f_04_01[:,0],2.7),                                            lw=3, zorder=1, dashes=( ), color=c3, label='$pp$')
    plot.   fill_between( f_04_01[:,0], f_04_01[:,2]     *np.power(f_04_01[:,0],2.7), f_04_01[:,3]*np.power(f_04_04[:,0],2.7),   lw=0, zorder=1,             color=c3, alpha=0.2 )
    plot.   plot        ( f_04_01[:,0], f_04_01[:,2]     *np.power(f_04_01[:,0],2.7),                                            lw=1, zorder=1, dashes=( ), color=c3            )
    plot.   plot        ( f_04_01[:,0],                                               f_04_01[:,3]*np.power(f_04_04[:,0],2.7),   lw=1, zorder=1, dashes=( ), color=c3            )

    plot.   plot        ( f_04_04[:,0], f_04_04[:,1]     *np.power(f_04_04[:,0],2.7),                                            lw=3, zorder=1, dashes=( ), color=c5, label='$pp$')
    plot.   fill_between( f_04_04[:,0], f_04_04[:,2]     *np.power(f_04_04[:,0],2.7), f_04_04[:,3]*np.power(f_04_04[:,0],2.7),   lw=0, zorder=1,             color=c5, alpha=0.2 )
    plot.   plot        ( f_04_04[:,0], f_04_04[:,2]     *np.power(f_04_04[:,0],2.7),                                            lw=1, zorder=1, dashes=( ), color=c5            )
    plot.   plot        ( f_04_04[:,0],                                               f_04_04[:,3]*np.power(f_04_04[:,0],2.7),   lw=1, zorder=1, dashes=( ), color=c5            )
    
    plot.   plot        ( f_CNO  [:,0], f_CNO  [:,1]     *np.power(f_CNO  [:,0],2.7),                                            lw=3, zorder=1, dashes=( ), color=c4, label='$pp$')
    plot.   fill_between( f_CNO  [:,0], f_CNO  [:,2]     *np.power(f_CNO  [:,0],2.7), f_CNO  [:,3]*np.power(f_CNO  [:,0],2.7),   lw=0, zorder=1,             color=c4, alpha=0.2 )
    plot.   plot        ( f_CNO  [:,0], f_CNO  [:,2]     *np.power(f_CNO  [:,0],2.7),                                            lw=1, zorder=1, dashes=( ), color=c4            )
    plot.   plot        ( f_CNO  [:,0],                                               f_CNO  [:,3]*np.power(f_CNO  [:,0],2.7),   lw=1, zorder=1, dashes=( ), color=c4            )
    
    plot.   plot        ( f_NeMgSi[:,0], f_NeMgSi[:,1]     *np.power(f_NeMgSi[:,0],2.7),                                              lw=3, zorder=1, dashes=( ), color=c8, label='$pp$')
    plot.   fill_between( f_NeMgSi[:,0], f_NeMgSi[:,2]     *np.power(f_NeMgSi[:,0],2.7), f_NeMgSi[:,3]*np.power(f_NeMgSi[:,0],2.7),   lw=0, zorder=1,             color=c8, alpha=0.2 )
    plot.   plot        ( f_NeMgSi[:,0], f_NeMgSi[:,2]     *np.power(f_NeMgSi[:,0],2.7),                                              lw=1, zorder=1, dashes=( ), color=c8            )
    plot.   plot        ( f_NeMgSi[:,0],                                                 f_NeMgSi[:,3]*np.power(f_NeMgSi[:,0],2.7),   lw=1, zorder=1, dashes=( ), color=c8            )
    
    plot.   plot        ( f_Fe   [:,0], f_Fe   [:,1]     *np.power(f_Fe   [:,0],2.7),                                            lw=3, zorder=1, dashes=( ), color=c7, label='$pp$')
    plot.   fill_between( f_Fe   [:,0], f_Fe   [:,2]     *np.power(f_Fe   [:,0],2.7), f_Fe   [:,3]*np.power(f_Fe   [:,0],2.7),   lw=0, zorder=1,             color=c7, alpha=0.2 )
    plot.   plot        ( f_Fe   [:,0], f_Fe   [:,2]     *np.power(f_Fe   [:,0],2.7),                                            lw=1, zorder=1, dashes=( ), color=c7            )
    plot.   plot        ( f_Fe   [:,0],                                               f_Fe   [:,3]*np.power(f_Fe   [:,0],2.7),   lw=1, zorder=1, dashes=( ), color=c7            )
    
    
    plot.   plot        ( f_LiBeB[:,0], f_LiBeB[:,1]     *np.power(f_LiBeB[:,0],2.7),                                            lw=3, zorder=1, dashes=( ), color=c6, label='$pp$')
    plot.   fill_between( f_LiBeB[:,0], f_LiBeB[:,2]     *np.power(f_LiBeB[:,0],2.7), f_LiBeB[:,3]*np.power(f_LiBeB[:,0],2.7),   lw=0, zorder=1,             color=c6, alpha=0.2 )
    plot.   plot        ( f_LiBeB[:,0], f_LiBeB[:,2]     *np.power(f_LiBeB[:,0],2.7),                                            lw=1, zorder=1, dashes=( ), color=c6            )
    plot.   plot        ( f_LiBeB[:,0],                                               f_LiBeB[:,3]*np.power(f_LiBeB[:,0],2.7),   lw=1, zorder=1, dashes=( ), color=c6            )
    
    plot.   plot        ( f_onHeavy[:,0], f_onHeavy[:,1]     *np.power(f_onHeavy[:,0],2.7),                                                lw=3, zorder=1, dashes=( ), color=c9, label='$pp$')
    plot.   fill_between( f_onHeavy[:,0], f_onHeavy[:,2]     *np.power(f_onHeavy[:,0],2.7), f_onHeavy[:,3]*np.power(f_onHeavy[:,0],2.7),   lw=0, zorder=1,             color=c9, alpha=0.2 )
    plot.   plot        ( f_onHeavy[:,0], f_onHeavy[:,2]     *np.power(f_onHeavy[:,0],2.7),                                                lw=1, zorder=1, dashes=( ), color=c9            )
    plot.   plot        ( f_onHeavy[:,0],                                                   f_onHeavy[:,3]*np.power(f_onHeavy[:,0],2.7),   lw=1, zorder=1, dashes=( ), color=c9            )
    
    
    plot.   plot        ( f_all  [:,0], f_all  [:,1]     *np.power(f_all  [:,0],2.7),                                            lw=3, zorder=1, dashes=( ), color=c0, label='$pp$')
    plot.   fill_between( f_all  [:,0], f_all  [:,2]     *np.power(f_all  [:,0],2.7), f_all  [:,3]*np.power(f_all  [:,0],2.7),   lw=0, zorder=1,             color=c0, alpha=0.2 )
    plot.   plot        ( f_all  [:,0], f_all  [:,2]     *np.power(f_all  [:,0],2.7),                                            lw=1, zorder=1, dashes=( ), color=c0            )
    plot.   plot        ( f_all  [:,0],                                               f_all  [:,3]*np.power(f_all  [:,0],2.7),   lw=1, zorder=1, dashes=( ), color=c0            )
    
    plot_u. plot        ( f_all  [:,0], f_all  [:,1]/f_all  [:,1]     ,                                            lw=3, zorder=1, dashes=( ), color=c0, label='$pp$')
    plot_u. fill_between( f_all  [:,0], f_all  [:,2]/f_all  [:,1]     , f_all  [:,3]/f_all  [:,1],                 lw=0, zorder=1,             color=c0, alpha=0.2 )
    plot_u. plot        ( f_all  [:,0], f_all  [:,2]/f_all  [:,1]     ,                                            lw=1, zorder=1, dashes=( ), color=c0            )
    plot_u. plot        ( f_all  [:,0],                                 f_all  [:,3]/f_all  [:,1],                 lw=1, zorder=1, dashes=( ), color=c0            )
    
    #plot_u. plot        ( f_all  [:,0], f_all  [:,2]/f_all  [:,1]+f_all  [:,4]/f_all  [:,1]-1,                     lw=1, zorder=1, dashes=( ), color=c0            )
    #plot_u. plot        ( f_all  [:,0], f_all  [:,3]/f_all  [:,1]+f_all  [:,5]/f_all  [:,1]-1,                     lw=1, zorder=1, dashes=( ), color=c0            )
    
    plot_u. plot        ( f_all  [:,0], f_all  [:,2]/f_all  [:,1]+f_all  [:,6]/f_all  [:,1]-1,                     lw=1, zorder=1, dashes=( ), color=c0            )
    plot_u. plot        ( f_all  [:,0], f_all  [:,3]/f_all  [:,1]+f_all  [:,7]/f_all  [:,1]-1,                     lw=1, zorder=1, dashes=( ), color=c0            )
    
    
    #f_err_all =  f_all[:,1]-f_all[:,2] #, f_all[:,3]-f_all[:,2] )
    #plot.fill_between( f_all  [:,0], f_err_all  *np.power(f_all  [:,0],2.7),  lw=0, zorder=1,             color=c0, alpha=0.2 )
    
    
    plot.text    (   8e2, 4.0e-21 , r'total',  color=c0,     rotation= +5, fontsize=label_size*0.6, ha='center' )
    plot.text    (   8e2, 2.0e-21 , r'$pp$',   color=c1,     rotation= +5, fontsize=label_size*0.6, ha='center' )
    plot.text    (   8e2, 8.2e-22 , r'He$p$',  color=c3,     rotation= +5, fontsize=label_size*0.6, ha='center' )
    plot.text    (   8e2, 3.3e-22 , r'$p$He',  color=c2,     rotation= +2, fontsize=label_size*0.6, ha='center' )
    plot.text    (   8e2, 1.3e-22 , r'HeHe',   color=c5,     rotation= +5, fontsize=label_size*0.6, ha='center' )
    plot.text    (   3e2, 2.6e-22 , r'CNO',    color=c4,     rotation= +7, fontsize=label_size*0.6, ha='center' )
    plot.text    (   3e2, 9.1e-23 , r'NeMgSi', color=c8,     rotation= +7, fontsize=label_size*0.6, ha='center' )
    plot.text    (   3e2, 4.5e-23 , r'Fe',     color=c7,     rotation= +7, fontsize=label_size*0.6, ha='center' )

    plot.text    (   3e2, 1.4e-23 , r'heavy ISM',  color=c9,     rotation= +4, fontsize=label_size*0.6, ha='center' )

    plot.text    (   3e2, 4.7e-24 , r'LiBeB',  color=c6,     rotation= -4, fontsize=label_size*0.6, ha='center' )

    plot  .set_xlim( (5e-1,  2e3)   )
    plot  .set_ylim( (1e-24, 1e-20) )

    plot_u.set_xlim( (5e-1,  2e3)   )
    plot_u.set_yticks([0.7,0.8,0.9, 1.0, 1.1, 1.2, 1.3] )
    plot_u.set_ylim( ( 0.7, 1.3 )   )

    plt.savefig('sourceterm_all.pdf')
    
    
    plot, fig = plot_1D (r'$\mathrm{T_{\bar{p}} [GeV]}$',      r'contribution', 'log', 'log' )
    
    
    
    plot.   plot        ( f_all  [:,0], f_all  [:,1] / f_all  [:,1]     ,                                            lw=3, zorder=1, dashes=( ), color=c0, label='$pp$')
    plot.   plot        ( f_01_01[:,0], f_01_01[:,1] / f_all  [:,1]     ,                                            lw=3, zorder=1, dashes=( ), color=c1, label='$pp$')
    plot.   plot        ( f_01_04[:,0], f_01_04[:,1] / f_all  [:,1]     ,                                            lw=3, zorder=1, dashes=( ), color=c2, label='$pp$')
    plot.   plot        ( f_04_01[:,0], f_04_01[:,1] / f_all  [:,1]     ,                                            lw=3, zorder=1, dashes=( ), color=c3, label='$pp$')
    plot.   plot        ( f_04_04[:,0], f_04_04[:,1] / f_all  [:,1]     ,                                            lw=3, zorder=1, dashes=( ), color=c5, label='$pp$')
    plot.   plot        ( f_CNO  [:,0], f_CNO  [:,1] / f_all  [:,1]     ,                                            lw=3, zorder=1, dashes=( ), color=c4, label='$pp$')
    plot.   plot        ( f_LiBeB[:,0], f_LiBeB[:,1] / f_all  [:,1]     ,                                            lw=3, zorder=1, dashes=( ), color=c6, label='$pp$')
    plot.   plot        (f_NeMgSi[:,0],f_NeMgSi[:,1] / f_all  [:,1]     ,                                            lw=3, zorder=1, dashes=( ), color=c8, label='$pp$')
    plot.   plot        ( f_Fe   [:,0], f_Fe   [:,1] / f_all  [:,1]     ,                                            lw=3, zorder=1, dashes=( ), color=c7, label='$pp$')

    plot.   plot        ( f_onHeavy[:,0], f_onHeavy[:,1] / f_all  [:,1]     ,                                        lw=3, zorder=1, dashes=( ), color=c9, label='$pp$')

    #plot.fill_between   ( f_all  [:,0], f_err_all  /f_all[:,1],  lw=0, zorder=1,             color=c0, alpha=0.2 )
    

    #plot.text    (   8e2, 4.0e-21 , r'total',  color=c0,     rotation= +5, fontsize=label_size*0.6, ha='center' )
    plot.text    (   8e2, 5.5e-1 , r'$pp$',   color=c1,     rotation= -2, fontsize=label_size*0.7, ha='center' )
    plot.text    (   8e2, 2.2e-1 , r'He$p$',  color=c3,     rotation= +3, fontsize=label_size*0.6, ha='center' )
    plot.text    (   8e2, 1.2e-1 , r'$p$He',  color=c2,     rotation= -3, fontsize=label_size*0.6, ha='center' )
    plot.text    (   8e2, 5.0e-2 , r'HeHe',   color=c5,     rotation= +3, fontsize=label_size*0.6, ha='center' )
    plot.text    (   8e2, 8.3e-2 , r'CNO',    color=c4,     rotation= +2, fontsize=label_size*0.6, ha='center' )
    plot.text    (   8e2, 2.8e-2 , r'NeMgSi', color=c8,     rotation= +2, fontsize=label_size*0.6, ha='center' )
    plot.text    (   8e2, 1.5e-2 , r'Fe',     color=c7,     rotation= +2, fontsize=label_size*0.6, ha='center' )
    plot.text    (   8e2, 4.0e-3 , r'heavy ISM',color=c9,   rotation= -2, fontsize=label_size*0.6, ha='center' )

    plot.text    (   4e0, 3.3e-3 , r'LiBeB',  color=c6,     rotation=  0, fontsize=label_size*0.6, ha='center' )


    plot  .set_xlim( (5e-1,  2e3)   )
    plot  .set_ylim( (2e-3, 2) )

    
    plt.savefig('contribution_all.pdf')


if step==5:
    #
    # This fuction assumes that all the results are obtained already, once for diMauro (cs: 0) and Winkler (cs: 1)
    # The results are stored in a directory "diMauro" and "Winkler" respectively.
    # In this case you can safely execute this step.
    #
    
    print 'XS_pp'
    plot_XS_comparison('standard', 'Winkler',    'diMauro',    'Param. II', 'Param. I')
    print 'XS_pHe'
    plot_XS_comparison('standard', 'Winkler',    'diMauro',    'Param. II', 'Param. I', '1-4')
    print 'XS_Hep'
    plot_XS_comparison('standard', 'Winkler',    'diMauro',    'Param. II', 'Param. I', '4-1')

    print 'source_term_pp'
    plot_sourceterm_comparison('standard', 'Winkler',    'diMauro',    'Param. II', 'Param. I')
    print 'source_term_pHe'
    plot_sourceterm_comparison('standard', 'Winkler',    'diMauro',    'Param. II', 'Param. I', '1-4')
    print 'source_term_Hep'
    plot_sourceterm_comparison('standard', 'Winkler',    'diMauro',    'Param. II', 'Param. I', '4-1')
    print 'source_term_HeHe'
    plot_sourceterm_comparison('standard', 'Winkler',    'diMauro',    'Param. II', 'Param. I', '4-4')

    #plot_sourceterm_comparison('standard', 'Winkler',    'diMauro',    'Param. II', 'Param. I', '12-1')
    #plot_sourceterm_comparison('standard', 'Winkler',    'diMauro',    'Param. II', 'Param. I', '12-4')
    #plot_sourceterm_comparison('standard', 'Winkler',    'diMauro',    'Param. II', 'Param. I', '16-1')
    #plot_sourceterm_comparison('standard', 'Winkler',    'diMauro',    'Param. II', 'Param. I', '16-4')


if step==6:
    data_01_00_01_00 = np.genfromtxt('sourceTerms_publish/dT_1-0_1-0_pbar_tot.txt')
    data_01_00_04_02 = np.genfromtxt('sourceTerms_publish/dT_1-0_4-2_pbar_tot.txt')
    
    data_02_01_01_00 = np.genfromtxt('sourceTerms_publish/dT_2-1_1-0_pbar_tot.txt')
    data_02_01_04_02 = np.genfromtxt('sourceTerms_publish/dT_2-1_4-2_pbar_tot.txt')
    
    data_03_02_01_00 = np.genfromtxt('sourceTerms_publish/dT_3-1_1-0_pbar_tot.txt')
    data_03_02_04_02 = np.genfromtxt('sourceTerms_publish/dT_3-1_4-2_pbar_tot.txt')
    
    data_04_02_01_00 = np.genfromtxt('sourceTerms_publish/dT_4-2_1-0_pbar_tot.txt')
    data_04_02_04_02 = np.genfromtxt('sourceTerms_publish/dT_4-2_4-2_pbar_tot.txt')
    
    data_12_06_01_00 = np.genfromtxt('sourceTerms_publish/dT_12-6_1-0_pbar_tot.txt')
    data_12_06_04_02 = np.genfromtxt('sourceTerms_publish/dT_12-6_4-2_pbar_tot.txt')
    
    data_13_07_01_00 = np.genfromtxt('sourceTerms_publish/dT_13-7_1-0_pbar_tot.txt')
    data_13_07_04_02 = np.genfromtxt('sourceTerms_publish/dT_13-7_4-2_pbar_tot.txt')
    
    data_14_07_01_00 = np.genfromtxt('sourceTerms_publish/dT_14-7_1-0_pbar_tot.txt')
    data_14_07_04_02 = np.genfromtxt('sourceTerms_publish/dT_14-7_4-2_pbar_tot.txt')
    
    data_15_08_01_00 = np.genfromtxt('sourceTerms_publish/dT_15-8_1-0_pbar_tot.txt')
    data_15_08_04_02 = np.genfromtxt('sourceTerms_publish/dT_15-8_4-2_pbar_tot.txt')
    
    data_16_08_01_00 = np.genfromtxt('sourceTerms_publish/dT_16-8_1-0_pbar_tot.txt')
    data_16_08_04_02 = np.genfromtxt('sourceTerms_publish/dT_16-8_4-2_pbar_tot.txt')
    
    data_17_09_01_00 = np.genfromtxt('sourceTerms_publish/dT_17-9_1-0_pbar_tot.txt')
    data_17_09_04_02 = np.genfromtxt('sourceTerms_publish/dT_17-9_4-2_pbar_tot.txt')
    
    data_18_10_01_00 = np.genfromtxt('sourceTerms_publish/dT_18-10_1-0_pbar_tot.txt')
    data_18_10_04_02 = np.genfromtxt('sourceTerms_publish/dT_18-10_4-2_pbar_tot.txt')
    


    toWrite = ''
    
    p=' I'
    if cs=='1':
        p = 'II'
    
    
    toWrite += '#\n'
    toWrite += '#        ************************* \n'
    toWrite += '#        * Sublementary material * \n'
    toWrite += '#        ************************* \n'
    toWrite += '#\n'
    toWrite += '#        * -------------------------------------------------------------------- * \n'
    toWrite += '#        * Production cross sections of cosmic antiprotons                      * \n'
    toWrite += '#        *            in the light of new data from NA49 and NA61 measurements  * \n'
    toWrite += '#        * -------------------------------------------------------------------- * \n'
    toWrite += '#\n'
    toWrite += '#          Michael Korsmeier, Fiorenza Donato, and Mattia Di Mauro      \n'
    toWrite += '#         *-------------------------------------------------------*     \n'
    toWrite += '#\n'
    toWrite += '#                    Published in: -                                    \n'
    toWrite += '#                    arXiv:        1802.????? [astro-ph.HE]             \n'
    toWrite += '#\n'
    toWrite += '#\n'
    toWrite += '#          IF YOU USE THIS TABLE, PLEASE CITE THe PAPER.                \n'
    toWrite += '#\n'
    toWrite += '#\n'
    toWrite += '#          We provide the energy differential cross section $d\sigma_{ij}/dT$        \n'
    toWrite += '#          for cosmic-ray (CR) component i and intersellar medium (ISM)              \n'
    toWrite += '#          component j. Here j are protons (Z=1,A=1) and Helium (Z=2,A=4).           \n'
    toWrite += '#          The first column contains the kinetic energy per nucleon of the           \n'
    toWrite += '#          incident CR $T_{projectile}$ and the second column the kinetic            \n'
    toWrite += '#          energy of the antiproton $T_{\bar{p}}, both in GeV.                       \n'
    toWrite += '#          The following columns contain the energy differential cross               \n'
    toWrite += '#          secions of various possible CR isotopes in units of m^2/GeV.              \n'
    toWrite += '#          These are total cross section, i.e. including antineutrons and            \n'
    toWrite += '#          antihyperons.                                                             \n'
    toWrite += '#\n'
    toWrite += '#\n'
    toWrite += '#               *-----------------------------------------------*            \n'
    toWrite += '#               | This corresponds to Param. '+str(p)+'-B in the paper. |    \n'
    toWrite += '#               *-----------------------------------------------*            \n'
    toWrite += '#\n'
    toWrite += '#          We recommend to use Param. II-B because of the better fit         \n'
    toWrite += '#          behaviour at high energies.                                       \n'
    toWrite += '#\n'
    toWrite += '#\n'
    toWrite += '#\n'

    toWrite += str('#').ljust(25*24-1, '*') + '\n'
    
    toWrite += '#  '+ str('T_{proj}/n [GeV]').ljust(21) + ' '
    toWrite +=        str('T_{pbar} [GeV]').ljust(24) + ' '
    
    toWrite += str('Z= 1, A= 1 (ISM p )').ljust(24) + ' '
    toWrite += str('Z= 1, A= 1 (ISM He)').ljust(24) + ' '
    
    toWrite += str('Z= 1, A= 2 (ISM p )').ljust(24) + ' '
    toWrite += str('Z= 1, A= 2 (ISM He)').ljust(24) + ' '
    
    toWrite += str('Z= 2, A= 3 (ISM p )').ljust(24) + ' '
    toWrite += str('Z= 2, A= 3 (ISM He)').ljust(24) + ' '
    
    toWrite += str('Z= 2, A= 4 (ISM p )').ljust(24) + ' '
    toWrite += str('Z= 2, A= 4 (ISM He)').ljust(24) + ' '
    
    toWrite += str('Z= 6, A=12 (ISM p )').ljust(24) + ' '
    toWrite += str('Z= 6, A=12 (ISM He)').ljust(24) + ' '
    
    toWrite += str('Z= 6, A=13 (ISM p )').ljust(24) + ' '
    toWrite += str('Z= 6, A=13 (ISM He)').ljust(24) + ' '
    
    toWrite += str('Z= 7, A=14 (ISM p )').ljust(24) + ' '
    toWrite += str('Z= 7, A=14 (ISM He)').ljust(24) + ' '
    
    toWrite += str('Z= 7, A=15 (ISM p )').ljust(24) + ' '
    toWrite += str('Z= 7, A=15 (ISM He)').ljust(24) + ' '
    
    toWrite += str('Z= 8, A=16 (ISM p )').ljust(24) + ' '
    toWrite += str('Z= 8, A=16 (ISM He)').ljust(24) + ' '
    
    toWrite += str('Z= 8, A=17 (ISM p )').ljust(24) + ' '
    toWrite += str('Z= 8, A=17 (ISM He)').ljust(24) + ' '
    
    toWrite += str('Z= 8, A=18 (ISM p )').ljust(24) + ' '
    toWrite += str('Z= 8, A=18 (ISM He)').ljust(24) + '\n'

    toWrite += str('#').ljust(25*24-1, '*') + '\n'
    
    
    
    for i in range(1,len(data_01_00_01_00[:,0])):
        for j in range(1,30*5+2): # len(data_01_00_01_00[0,:])
            toWrite += str(data_01_00_01_00[i,0]).ljust(24) + ' ' + str(data_01_00_01_00[0,j]).ljust(24) + ' '
            
            toWrite += str(data_01_00_01_00[i,j]).ljust(24) + ' '
            toWrite += str(data_01_00_04_02[i,j]).ljust(24) + ' '
            
            toWrite += str(data_02_01_01_00[i,j]).ljust(24) + ' '
            toWrite += str(data_02_01_04_02[i,j]).ljust(24) + ' '
            
            toWrite += str(data_03_02_01_00[i,j]).ljust(24) + ' '
            toWrite += str(data_03_02_04_02[i,j]).ljust(24) + ' '

            toWrite += str(data_04_02_01_00[i,j]).ljust(24) + ' '
            toWrite += str(data_04_02_04_02[i,j]).ljust(24) + ' '

            toWrite += str(data_12_06_01_00[i,j]).ljust(24) + ' '
            toWrite += str(data_12_06_04_02[i,j]).ljust(24) + ' '

            toWrite += str(data_13_07_01_00[i,j]).ljust(24) + ' '
            toWrite += str(data_13_07_04_02[i,j]).ljust(24) + ' '

            toWrite += str(data_14_07_01_00[i,j]).ljust(24) + ' '
            toWrite += str(data_14_07_04_02[i,j]).ljust(24) + ' '

            toWrite += str(data_15_08_01_00[i,j]).ljust(24) + ' '
            toWrite += str(data_15_08_04_02[i,j]).ljust(24) + ' '

            toWrite += str(data_16_08_01_00[i,j]).ljust(24) + ' '
            toWrite += str(data_16_08_04_02[i,j]).ljust(24) + ' '
            
            toWrite += str(data_17_09_01_00[i,j]).ljust(24) + ' '
            toWrite += str(data_17_09_04_02[i,j]).ljust(24) + ' '
            
            toWrite += str(data_18_10_01_00[i,j]).ljust(24) + ' '
            toWrite += str(data_18_10_04_02[i,j]).ljust(24) + '\n'

    f = open('XS_table.dat','w')
    f.write(toWrite)
    f.close()


if step==100:

#    plot2D(exp+'/dT_pbar_LAB.txt',         resfile=exp+'/XS.pdf'        , x_label=r'$\mathrm{T_{p}\quad [GeV]}$'  )
#    plot2D(exp+'/dT_pbar_LAB__exp.txt',    resfile=exp+'/XS__exp.pdf'   , x_label=r'$\mathrm{T_{p}\quad [GeV]}$'  )

    print exp

    source          = np.genfromtxt(exp+'/source.txt',         skip_header=1)
    source__exp     = np.genfromtxt(exp+'/source__exp.txt',    skip_header=1)

    plot, fig = plot_1D (r'$\mathrm{T_{\bar{p}} [GeV]}$',      r'$\mathrm{q^{(\bar{p})} \, T_{\bar{p}}^{2.7}  \,  \, [GeV^{1.7}m^{-3}s^{-1}]}}$', 'log', 'log' )

    plot.plot( source      [:,0], source      [:, 1]*np.power(source      [:,0],2.7) , color=c0, lw=3)
    plot.plot( source__exp [:,0], source__exp [:, 1]*np.power(source__exp [:,0],2.7) , color=c1, lw=3)

    plot.set_xlim( (1e-1, 1e4)         )
    plot.set_ylim( (2e-26, 5e-21)   )

    plt.savefig(exp+'/sourceterm_comparison.png')

    plot, fig = plot_1D (r'$\mathrm{T_{\bar{p}} [GeV]}$',      r'source term fraction', 'log', 'log' )
    
    frac = source__exp[:, 1]/source[:, 1] + 1e-10;
    plot.plot( source[:,0], frac, color=c0, lw=3)
    
    plot.set_xlim( (1e-1, 1e4)   )
    plot.set_ylim( (1e-5, 1e1)   )
    
    plt.savefig(exp+'/sourceterm_fraction.png')

    plot2D



if step==101:
    
    exp     = [ 'lhcb_pHe', 'lhcb_pHe_87', 'lhcb_pHe_43', 'lhcb_pHe_combined']
    label   = [ r'$\sqrt{s}= 100-120 \quad \mathrm{GeV}$', r'$\sqrt{s}= 80-100 \quad \mathrm{GeV}$', r'$\sqrt{s}= 40-50 \quad \mathrm{GeV}$', r'$\mathrm{combined}$', ]
    cc      = [ c3, c3, c3, c3 ]
    dashes  = [ (), (5,3), (7,3,2,3), (1.5,1.5) ]
    #dashes  = [ (), (12,4), (20,3,3,3), (3,3) ]
    
    plot, fig = plot_1D (r'$\mathrm{T_{\bar{p}} [GeV]}$',      r'source term contribution', 'log', 'log' )
    for i in range(len(exp)):
        source          = np.genfromtxt(exp[i]+'/source.txt',         skip_header=1)
        source__exp     = np.genfromtxt(exp[i]+'/source__exp.txt',    skip_header=1)
        frac = source__exp[:, 1]/source[:, 1] + 1e-10;
        plot.plot( source[:,0], frac , color=cc[i], dashes=dashes[i], lw=3, label=label[i])
    
    plot.set_xlim( (1e-1, 1e4)   )
    plot.set_ylim( (2e-5, 3e1)   )
    plot.axhline(y=1, lw=2, color='black')
    plot.legend( loc='upper center', bbox_to_anchor=(0.5, 0.97), frameon=False, ncol=2, fontsize=0.65*label_size )
    plot.text( 2e-1, 2e-1, 'LHCb', fontsize=0.65*label_size )
    plt.savefig('sourceterm_comparison_LHCb_pHe.pdf')


    exp     = [ 'lhcb_Hep', 'lhcb_Hep_87', 'lhcb_Hep_43', 'lhcb_Hep_combined']
    label   = [ r'$\sqrt{s}= 100-120 \quad \mathrm{GeV}$', r'$\sqrt{s}= 80-100 \quad \mathrm{GeV}$', r'$\sqrt{s}= 40-50 \quad \mathrm{GeV}$', r'$\mathrm{combined}$', ]
    cc      = [ c3, c3, c3, c3 ]
    dashes  = [ (), (5,3), (7,3,2,3), (1.5,1.5) ]

    plot, fig = plot_1D (r'$\mathrm{T_{\bar{p}} [GeV]}$',      r'source term contribution', 'log', 'log' )
    for i in range(len(exp)):
        source          = np.genfromtxt(exp[i]+'/source.txt',         skip_header=1)
        source__exp     = np.genfromtxt(exp[i]+'/source__exp.txt',    skip_header=1)
        frac = source__exp[:, 1]/source[:, 1] + 1e-10;
        plot.plot( source[:,0], frac , color=cc[i], dashes=dashes[i], lw=3, label=label[i])
    
    plot.set_xlim( (1e-1, 1e4)   )
    plot.set_ylim( (2e-5, 3e1)   )
    plot.axhline(y=1, lw=2, color='black')
    plot.legend( loc='upper center', bbox_to_anchor=(0.5, 0.97), frameon=False, ncol=2, fontsize=0.65*label_size )
    plot.text( 2e-1, 3e-1, 'LHCb', fontsize=0.65*label_size )
    plt.savefig('sourceterm_comparison_LHCb_Hep.pdf')



    exp     = ['na49', 'na49_large', 'na61', 'na61_large']
    label   = [r'NA49 $(\sqrt{s}= 15-20 \quad \mathrm{GeV})$', r'NA49 $(\sqrt{s}= 10-50 \quad \mathrm{GeV})$', 'NA61', r'NA61 $($scaled to $\sqrt{s}= 50 \quad \mathrm{GeV})$', ]
    cc      = [ c2, c2, c1, c1 ]
    dashes  = [ (), (5,3), (7,3,2,3), (1.5,1.5) ]

    plot, fig = plot_1D (r'$\mathrm{T_{\bar{p}} [GeV]}$',      r'source term contribution', 'log', 'log' )
    for i in range(len(exp)):
        source          = np.genfromtxt(exp[i]+'/source.txt',         skip_header=1)
        source__exp     = np.genfromtxt(exp[i]+'/source__exp.txt',    skip_header=1)
        frac = source__exp[:, 1]/source[:, 1] + 1e-10;
        plot.plot( source[:,0], frac , color=cc[i], dashes=dashes[i], lw=3, label=label[i])

    plot.set_xlim( (1e-1, 1e4)   )
    plot.set_ylim( (2e-2, 3e0)   )
    plot.axhline(y=1, lw=2, color='black')
    plot.legend( loc='upper center', bbox_to_anchor=(0.5, 0.97), frameon=False, ncol=2, fontsize=0.65*label_size )
    plt.savefig('sourceterm_comparison_NA49_NA61.pdf')



    exp     = [ 'na61', 'na61_large', 'na61_low', 'na61_all']
    label   = [ 'data', r'$(\sqrt{s}$ scaled to $50\; \mathrm{GeV})$', r'$(\sqrt{s}$ scaled to $6.3\; \mathrm{GeV} )$', r'$(\sqrt{s}$ scaled to $6.3$ and $50\; \mathrm{GeV})$', ]
    cc      = [ c1, c1, c4, c4 ]
    dashes  = [ (), (5,3), (7,3,2,3), (1.5,1.5) ]

    plot, fig = plot_1D (r'$\mathrm{T_{\bar{p}} [GeV]}$',      r'source term contribution', 'log', 'log' )
    for i in range(len(exp)):
        source          = np.genfromtxt(exp[i]+'/source.txt',         skip_header=1)
        source__exp     = np.genfromtxt(exp[i]+'/source__exp.txt',    skip_header=1)
        frac = source__exp[:, 1]/source[:, 1] + 1e-10;
        plot.plot( source[:,0], frac , color=cc[i], dashes=dashes[i], lw=3, label=label[i])

    plot.set_xlim( (1e-1, 1e4)   )
    plot.set_ylim( (2e-2, 3e0)   )
    plot.axhline(y=1, lw=2, color='black')
    plot.legend( loc='upper center', bbox_to_anchor=(0.5, 0.97), frameon=False, ncol=2, fontsize=0.65*label_size )
    plot.text( 2e-1, 6e-1, 'NA61', fontsize=0.65*label_size )
    plt.savefig('sourceterm_comparison_NA61.pdf')



    exp     = [ 'lhcb_pHe', 'na49_pC', 'lhcb_Hep', 'na49_Cp']
    label   = [ r'LHCb $p$He $(\sqrt{s}= 100-120 \quad \mathrm{GeV})$', r'NA49 $p$C $(\sqrt{s}= 15-20 \quad \mathrm{GeV})$', r'He$p$', r'C$p$' ]
    cc      = [ c3, c2, c3, c2 ]
    dashes  = [ (), (5,3), (7,3,2,3), (1.5,1.5) ]

    plot, fig = plot_1D (r'$\mathrm{T_{\bar{p}} [GeV]}$',      r'source term contribution', 'log', 'log' )
    for i in range(len(exp)):
        source          = np.genfromtxt(exp[i]+'/source.txt',         skip_header=1)
        source__exp     = np.genfromtxt(exp[i]+'/source__exp.txt',    skip_header=1)
        frac = source__exp[:, 1]/source[:, 1] + 1e-10;
        plot.plot( source[:,0], frac , color=cc[i], dashes=dashes[i], lw=3, label=label[i])

    plot.set_xlim( (1e-1, 1e4)   )
    plot.set_ylim( (2e-5, 2e1)   )
    plot.axhline(y=1, lw=2, color='black')
    plot.legend( loc='upper center', bbox_to_anchor=(0.5, 0.97), frameon=False, ncol=2, fontsize=0.65*label_size )
    plt.savefig('sourceterm_comparison_NA49_LHCb.pdf')


############################################################################################################################################################
############################################################################################################################################################
